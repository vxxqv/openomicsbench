"""Build the deterministic DNA, RNA, protein and paired-read v2 fixtures."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omicsbench.hashing import digest
from omicsbench.sequences import summarize


FIXTURES = {
    "sequence-001": {
        "title": "Synthetic nucleotide sequence validation set",
        "assay": "sequence_dna",
        "archetype": "DNA FASTA with IUPAC ambiguity and an exact duplicate",
        "molecule": "dna",
        "format": "fasta",
        "files": {
            "nano/sequences.fasta": (
                ">dna_alpha complete synthetic coding sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
                ">dna_beta ambiguous nucleotide sequence\nACGTRYSWKMBDHVNACGTACGT\n"
                ">dna_gamma short sequence\nGGGGCCCCAAAATTTT\n"
                ">dna_gamma_copy exact sequence duplicate\nGGGGCCCCAAAATTTT\n"
            )
        },
        "roles": ["sequence"],
        "limitations": ["Synthetic sequences exercise file handling and deterministic metrics; they do not represent a biological cohort."],
    },
    "sequence-002": {
        "title": "Synthetic RNA sequence validation set",
        "assay": "sequence_rna",
        "archetype": "RNA FASTA with coding, ambiguous and short records",
        "molecule": "rna",
        "format": "fasta",
        "files": {
            "nano/sequences.fasta": (
                ">rna_alpha synthetic transcript\nAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG\n"
                ">rna_beta ambiguous RNA\nACGURYSWKMBDHVNACGU\n"
                ">rna_gamma short noncoding sequence\nGGGGCCCCAAAAUUUU\n"
            )
        },
        "roles": ["sequence"],
        "limitations": ["Synthetic RNA records validate alphabet-aware operations and are not evidence for expression or structure."],
    },
    "sequence-003": {
        "title": "Synthetic protein sequence validation set",
        "assay": "sequence_protein",
        "archetype": "protein FASTA with standard, ambiguous and termination symbols",
        "molecule": "protein",
        "format": "fasta",
        "files": {
            "nano/sequences.fasta": (
                ">protein_alpha translated reference\nMAIVMGRXKGAR\n"
                ">protein_beta mixed amino acid sequence\nMKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\n"
                ">protein_gamma ambiguous residues\nACDEFGHIKLMNPQRSTVWYBXZJUO*\n"
            )
        },
        "roles": ["sequence"],
        "limitations": ["Synthetic proteins validate parsing and filtering only; no homology, structure or function is implied."],
    },
    "sequence-004": {
        "title": "Synthetic paired FASTQ validation set",
        "assay": "short_read_dna",
        "archetype": "paired-end reads with variable quality and exact read duplication",
        "molecule": "dna",
        "format": "fastq",
        "files": {
            "nano/reads_R1.fastq": (
                "@pair001/1 high quality\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
                "@pair002/1 mixed quality\nGGGGCCCCAAAA\n+\nIIIII5555!!!\n"
                "@pair003/1 ambiguous read\nACGTNNNNACGT\n+\n????????????\n"
                "@pair004/1 duplicate sequence\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
            ),
            "nano/reads_R2.fastq": (
                "@pair001/2 high quality\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n"
                "@pair002/2 mixed quality\nTTTTGGGGCCCC\n+\nIIIII5555!!!\n"
                "@pair003/2 ambiguous read\nTGCAAAAATGCA\n+\n????????????\n"
                "@pair004/2 duplicate sequence\nTGCATGCATGCA\n+\nIIIIIIIIIIII\n"
            ),
        },
        "roles": ["fastq_r1", "fastq_r2"],
        "limitations": ["Synthetic reads test structural QC and pair handling; they are not derived from an instrument or organism."],
    },
}


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")


def file_record(folder: Path, relative: str, tier: str, role: str, media_type: str) -> dict:
    path = folder / relative
    return {"path": relative, "tier": tier, "role": role, "media_type": media_type, "bytes": path.stat().st_size, "sha256": digest(path)}


def build(destination: Path) -> list[Path]:
    outputs = []
    for fixture_id, config in FIXTURES.items():
        folder = destination / "datasets" / "sequences" / fixture_id
        folder.mkdir(parents=True, exist_ok=True)
        for relative, content in config["files"].items():
            path = folder / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8", newline="\n")
        summaries = []
        for relative in config["files"]:
            report = summarize(folder / relative, config["molecule"])
            report.pop("path")
            summaries.append({"path": relative, "expected": report})
        profile = {"profile": "sequence-summary-v1", "files": summaries}
        paired = len(config["files"]) == 2
        if paired:
            profile["pairs"] = 4
        write_json(folder / "expected/validation.json", profile)
        inventory = []
        for (relative, _), role in zip(config["files"].items(), config["roles"]):
            inventory.append(file_record(folder, relative, "nano", role, "application/fastq" if config["format"] == "fastq" else "text/x-fasta"))
        inventory.append(file_record(folder, "expected/validation.json", "expected", "metrics", "application/json"))
        manifest = {
            "schema_version": "2.0",
            "id": fixture_id,
            "title": config["title"],
            "release": "2.0.0",
            "assay": config["assay"],
            "kind": "synthetic_fixture",
            "status": "validated",
            "archetype": config["archetype"],
            "organism": "Synthetic sequence model",
            "taxon_id": None,
            "source": {
                "repository": "Project generator",
                "accession": fixture_id,
                "url": "https://github.com/vxxqv/openomicsbench",
                "citation": f"OpenOmicsBench {fixture_id}, project-authored deterministic sequence fixture.",
                "retrieved": "2026-09-20",
                "sha256": None,
            },
            "rights": {
                "status": "GREEN",
                "license": "Apache-2.0",
                "evidence": "Generated from literals in the project workflow without third-party sequence observations.",
                "checked": "2026-09-20",
            },
            "reference": {"genome": "none", "annotation": "none", "namespace": "project synthetic identifiers", "compatibility": "not_applicable"},
            "samples": [{"sample_id": f"{fixture_id}-sample", "condition": "synthetic validation", "replicate": 1, "batch": "not_applicable", "strandedness": "not_applicable"}],
            "derivation": {
                "workflow": "workflows/build_sequence_fixtures.py",
                "commit": None,
                "seed": 0,
                "algorithm": "literal deterministic sequence fixtures",
                "parameters": {"records": sum(item["expected"]["records"] for item in summaries), "paired": paired},
                "versions": {"python": "3.12", "openomicsbench_contract": "2.0"},
            },
            "files": inventory,
            "validation": {"profile": "expected/validation.json", "baseline_version": "sequence-summary-v1", "metrics": []},
            "sequence": {"format": config["format"], "molecule": config["molecule"], "paired": paired, "quality_encoding": "phred33" if config["format"] == "fastq" else "not_applicable"},
            "limitations": config["limitations"],
        }
        write_json(folder / "manifest.json", manifest)
        outputs.append(folder)
    return outputs


if __name__ == "__main__":
    target = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT
    built = build(target)
    print(f"Built {len(built)} sequence fixtures in {target / 'datasets' / 'sequences'}")
