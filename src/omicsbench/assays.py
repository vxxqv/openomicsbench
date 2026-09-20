"""Practical sequence-input profiles and local preprocessing plans."""
from __future__ import annotations

import shlex
from pathlib import Path


PROFILES = {
    "whole-genome": {
        "name": "Whole-genome DNA sequencing",
        "molecule": "dna",
        "inputs": ["paired or single-end FASTQ"],
        "checks": ["read structure", "per-base quality", "adapter presence", "duplication", "GC distribution", "pair integrity"],
        "local_scope": "Read validation, QC, trimming, filtering, sampling and k-mer comparison.",
        "downstream": "A reference-aware aligner and germline or somatic variant workflow are required for variant calls.",
    },
    "exome": {
        "name": "Whole-exome sequencing",
        "molecule": "dna",
        "inputs": ["paired FASTQ", "capture target BED for downstream analysis"],
        "checks": ["pair integrity", "per-base quality", "adapter presence", "duplication", "read length"],
        "local_scope": "Read validation, QC, trimming, filtering and deterministic subsets.",
        "downstream": "Target-aware alignment, coverage assessment and a variant caller are required for biological results.",
    },
    "targeted-panel": {
        "name": "Targeted DNA panel sequencing",
        "molecule": "dna",
        "inputs": ["paired FASTQ", "panel BED", "UMI layout when present"],
        "checks": ["pair integrity", "adapter presence", "quality", "read length", "UMI preservation"],
        "local_scope": "Read validation and preprocessing when UMIs are not altered.",
        "downstream": "UMI consensus, target coverage and low-frequency variant calling need panel-aware tools.",
    },
    "bulk-rna-seq": {
        "name": "Bulk RNA sequencing",
        "molecule": "dna",
        "inputs": ["paired or single-end FASTQ", "sample design", "reference genome or transcriptome"],
        "checks": ["read quality", "adapter presence", "pair integrity", "duplication", "read length"],
        "local_scope": "Read QC and preprocessing plus the existing count-matrix benchmarks.",
        "downstream": "Splice-aware alignment or transcript quantification and design-aware differential expression remain external.",
    },
    "single-cell-rna-seq": {
        "name": "Single-cell or single-nucleus RNA sequencing",
        "molecule": "dna",
        "inputs": ["FASTQ", "library chemistry", "barcode and UMI layout", "reference transcriptome"],
        "checks": ["pair integrity", "read length", "quality by position", "barcode-read preservation"],
        "local_scope": "Structural FASTQ validation, pair checks, QC and non-destructive sampling.",
        "downstream": "Chemistry-aware barcode correction, UMI counting and cell calling require a dedicated workflow.",
    },
    "atac-seq": {
        "name": "ATAC-seq",
        "molecule": "dna",
        "inputs": ["paired FASTQ", "reference genome", "blacklist regions"],
        "checks": ["pair integrity", "adapter presence", "quality", "read length", "duplication"],
        "local_scope": "Read validation, QC and adapter trimming.",
        "downstream": "Alignment, TSS enrichment, fragment-size QC and peak calling require genome-aware tools.",
    },
    "chip-seq": {
        "name": "ChIP-seq, CUT&RUN or CUT&Tag",
        "molecule": "dna",
        "inputs": ["FASTQ", "control relationship", "reference genome"],
        "checks": ["quality", "adapter presence", "pair integrity when paired", "duplication"],
        "local_scope": "Read validation, QC and preprocessing.",
        "downstream": "Control-aware alignment, signal tracks and broad or narrow peak calling remain external.",
    },
    "amplicon": {
        "name": "Marker-gene amplicon sequencing",
        "molecule": "dna",
        "inputs": ["paired or single-end FASTQ", "primer sequences", "marker and expected amplicon length"],
        "checks": ["primer presence", "pair integrity", "quality", "length", "ambiguity"],
        "local_scope": "Primer scans, trimming, filtering, sampling and k-mer inspection.",
        "downstream": "Denoising, chimera removal and taxonomic assignment need marker-specific reference data.",
    },
    "shotgun-metagenomics": {
        "name": "Shotgun metagenomic sequencing",
        "molecule": "dna",
        "inputs": ["paired or single-end FASTQ", "host reference when depletion is required"],
        "checks": ["quality", "adapter presence", "pair integrity", "complexity", "duplication"],
        "local_scope": "Read QC, filtering, sampling and k-mer similarity.",
        "downstream": "Host depletion, taxonomic classification and functional profiling require large external databases.",
    },
    "long-read-dna": {
        "name": "Long-read DNA sequencing",
        "molecule": "dna",
        "inputs": ["FASTA or FASTQ from Oxford Nanopore or PacBio"],
        "checks": ["read length distribution", "N50", "quality when available", "ambiguity"],
        "local_scope": "Long-read parsing, QC, filtering, sampling and sketches.",
        "downstream": "Platform-aware basecalling, assembly, polishing or long-read alignment remain external.",
    },
    "long-read-rna": {
        "name": "Long-read RNA or cDNA sequencing",
        "molecule": "dna",
        "inputs": ["FASTA or FASTQ", "library orientation", "reference genome and annotation"],
        "checks": ["read length distribution", "N50", "quality", "adapter presence"],
        "local_scope": "Sequence validation, QC, filtering, sampling and motif scans.",
        "downstream": "Splice-aware long-read alignment and isoform discovery remain external.",
    },
    "reference": {
        "name": "Genome or transcriptome reference",
        "molecule": "dna",
        "inputs": ["FASTA"],
        "checks": ["unique identifiers", "IUPAC symbols", "lengths", "N50", "ambiguity", "duplicate sequences"],
        "local_scope": "Validation, statistics, extraction, motifs, ORFs, translation and sketches.",
        "downstream": "Index construction is specific to the aligner or classifier being used.",
    },
    "rna-fasta": {
        "name": "RNA sequence collection",
        "molecule": "rna",
        "inputs": ["FASTA"],
        "checks": ["unique identifiers", "RNA IUPAC symbols", "lengths", "ambiguity", "duplicate sequences"],
        "local_scope": "Validation, statistics, extraction, motifs, reverse complements, back-transcription and translation.",
        "downstream": "Structure prediction and homology search are outside the built-in scope.",
    },
    "protein-fasta": {
        "name": "Protein sequence collection",
        "molecule": "protein",
        "inputs": ["FASTA"],
        "checks": ["unique identifiers", "amino-acid symbols", "lengths", "ambiguity", "duplicate sequences"],
        "local_scope": "Validation, statistics, filtering, extraction, sampling and deduplication.",
        "downstream": "Alignment, domain search, structure prediction and functional annotation remain external.",
    },
}


def list_profiles() -> list[dict]:
    return [{"id": key, "name": value["name"], "molecule": value["molecule"], "inputs": value["inputs"]} for key, value in PROFILES.items()]


def get_profile(profile_id: str) -> dict:
    if profile_id not in PROFILES:
        near = [key for key in PROFILES if profile_id.casefold() in key.casefold() or key.casefold() in profile_id.casefold()]
        raise ValueError(f"unknown assay profile {profile_id}. Try: {', '.join(near[:3]) or 'omicsbench assay list'}")
    return {"id": profile_id, **PROFILES[profile_id]}


def _command(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def build_plan(profile_id: str, inputs: list[Path], output_directory: Path, adapter: list[str] | None = None) -> dict:
    profile = get_profile(profile_id)
    if not inputs:
        raise ValueError("at least one input file is required")
    if len(inputs) > 2:
        raise ValueError("local plans accept one sequence file or one paired FASTQ set")
    output_directory = Path(output_directory)
    steps = []
    for index, path in enumerate(inputs, start=1):
        label = f"r{index}" if len(inputs) == 2 else "sequences"
        steps.append({"id": f"validate-{label}", "purpose": "Strict format and alphabet validation", "command": _command(["omicsbench", "seq", "validate", str(path), "--molecule", profile["molecule"]])})
        qc = ["omicsbench", "seq", "qc", str(path), "--molecule", profile["molecule"]]
        for sequence in adapter or []:
            qc.extend(["--adapter", sequence])
        steps.append({"id": f"qc-{label}", "purpose": "Composition, length, complexity and quality report", "command": _command(qc)})
    if len(inputs) == 2:
        steps.append({"id": "pair-check", "purpose": "Verify paired record counts and identifiers", "command": _command(["omicsbench", "seq", "pair-check", str(inputs[0]), str(inputs[1])])})
    if profile_id in {"whole-genome", "exome", "targeted-panel", "bulk-rna-seq", "single-cell-rna-seq", "atac-seq", "chip-seq", "amplicon", "shotgun-metagenomics"}:
        for index, path in enumerate(inputs, start=1):
            label = f"r{index}" if len(inputs) == 2 else "reads"
            target = output_directory / f"{label}.filtered.fastq.gz"
            command = ["omicsbench", "seq", "filter", str(path), str(target), "--molecule", "dna", "--min-length", "20", "--max-ambiguity", "0.05", "--min-mean-quality", "20"]
            steps.append({"id": f"filter-{label}", "purpose": "Remove very short, ambiguous or low-quality records", "command": _command(command)})
    return {"profile": profile, "inputs": [str(path) for path in inputs], "output_directory": str(output_directory), "steps": steps, "boundary": profile["downstream"]}
