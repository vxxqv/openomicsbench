"""Build a small, fully synthetic paired-read classification test case."""

from __future__ import annotations

import json
import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omicsbench.hashing import digest
from omicsbench.sequences import reverse_complement, summarize
from omicsbench.taxonomy import exact_assignments


SEED = 20261101
READ_LENGTH = 50
INSERT_LENGTH = 130
POSITIVE_STARTS = (10, 30, 50, 70, 90, 110, 130, 420, 440, 460, 480, 500)
SHARED_STARTS = (270, 275, 280, 285, 290, 295)
BACKGROUND_STARTS = (10, 35, 60, 85, 110, 135)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")


def file_record(folder: Path, relative: str, tier: str, role: str, media_type: str) -> dict:
    path = folder / relative
    return {"path": relative, "tier": tier, "role": role, "media_type": media_type, "bytes": path.stat().st_size, "sha256": digest(path)}


def build(destination: Path) -> Path:
    folder = Path(destination) / "datasets/sequences/sequence-007"
    (folder / "nano").mkdir(parents=True, exist_ok=True)
    (folder / "expected").mkdir(parents=True, exist_ok=True)
    generator = random.Random(SEED)
    shared = "".join(generator.choices("ACGT", k=170))
    references = {}
    for code in (1, 2, 3):
        unique = "".join(generator.choices("ACGT", k=650))
        references[code] = unique[:260] + shared + unique[430:]
    background = "".join(generator.choices("ACGT", k=650))
    (folder / "nano/reference.fasta").write_text(
        "".join(f">mock_taxon_{code}\n{sequence}\n" for code, sequence in references.items()),
        encoding="utf-8", newline="\n",
    )
    read1 = []
    read2 = []
    truth = ["read_id\ttaxon_id"]

    def add_pair(name: str, sequence: str, start: int, code: int) -> None:
        left = sequence[start:start + READ_LENGTH]
        right = reverse_complement(sequence[start + INSERT_LENGTH - READ_LENGTH:start + INSERT_LENGTH], "dna")
        if len(left) != READ_LENGTH or len(right) != READ_LENGTH:
            raise ValueError("paired read lies outside its reference")
        read1.append(f"@{name}/1\n{left}\n+\n{'I' * READ_LENGTH}\n")
        read2.append(f"@{name}/2\n{right}\n+\n{'I' * READ_LENGTH}\n")
        truth.append(f"{name}\t{code}")

    for code, sequence in references.items():
        for index, start in enumerate(POSITIVE_STARTS, 1):
            add_pair(f"taxon{code}_{index:02d}", sequence, start, code)
    for index, start in enumerate(SHARED_STARTS, 1):
        add_pair(f"shared_{index:02d}", shared, start - 260, 0)
    for index, start in enumerate(BACKGROUND_STARTS, 1):
        add_pair(f"background_{index:02d}", background, start, 0)
    (folder / "nano/reads_R1.fastq").write_text("".join(read1), encoding="utf-8", newline="\n")
    (folder / "nano/reads_R2.fastq").write_text("".join(read2), encoding="utf-8", newline="\n")
    (folder / "expected/read-taxonomy.tsv").write_text("\n".join(truth) + "\n", encoding="utf-8", newline="\n")
    exact = exact_assignments(folder / "nano/reference.fasta", folder / "nano/reads_R1.fastq", folder / "nano/reads_R2.fastq")
    if list(exact.items()) != [(line.split("\t")[0], int(line.split("\t")[1])) for line in truth[1:]]:
        raise ValueError("synthetic assignments do not match the declared truth")
    summaries = []
    for relative in ("nano/reads_R1.fastq", "nano/reads_R2.fastq"):
        report = summarize(folder / relative, "dna")
        report.pop("path")
        summaries.append({"path": relative, "expected": report})
    profile = {
        "profile": "paired-read-taxonomy-v1",
        "files": summaries,
        "pairs": len(truth) - 1,
        "reference": "nano/reference.fasta",
        "read_truth": "expected/read-taxonomy.tsv",
        "taxonomy_namespace": "project-local mock codes, not NCBI taxids",
        "truth_counts": {str(code): count for code, count in sorted(Counter(exact.values()).items())},
        "minimum_scores": {"coverage": 1.0, "macro_f1": 0.95, "classified_precision": 0.95, "classified_recall": 0.95, "unclassified_recall": 1.0},
    }
    write_json(folder / "expected/validation.json", profile)
    inventory = [
        file_record(folder, "nano/reference.fasta", "nano", "reference", "text/x-fasta"),
        file_record(folder, "nano/reads_R1.fastq", "nano", "fastq_r1", "application/fastq"),
        file_record(folder, "nano/reads_R2.fastq", "nano", "fastq_r2", "application/fastq"),
        file_record(folder, "expected/read-taxonomy.tsv", "expected", "effects", "text/tab-separated-values"),
        file_record(folder, "expected/validation.json", "expected", "metrics", "application/json"),
    ]
    manifest = {
        "schema_version": "2.0", "id": "sequence-007", "title": "Synthetic paired-read taxonomy assignment benchmark",
        "release": "8.0.0", "assay": "shotgun_metagenomics", "kind": "synthetic_fixture", "status": "validated",
        "archetype": "three mock references, shared reads and absent-reference background pairs",
        "organism": "Synthetic community model", "taxon_id": None,
        "source": {"repository": "Project generator", "accession": "sequence-007", "url": "https://github.com/vxxqv/openomicsbench", "citation": "OpenOmicsBench sequence-007, project-authored deterministic read-classification fixture.", "retrieved": "2026-11-01", "sha256": None},
        "rights": {"status": "GREEN", "license": "Apache-2.0", "evidence": "All bases and reads are generated locally from a fixed seed without third-party sequence observations.", "checked": "2026-11-01"},
        "reference": {"genome": "mock_reference_v1", "annotation": "not_applicable", "namespace": "project-local mock taxon codes", "compatibility": "verified"},
        "samples": [{"sample_id": "mock-community-1", "condition": "three mock taxa with shared and background pairs", "replicate": 1, "batch": "not_applicable", "strandedness": "unstranded"}],
        "derivation": {"workflow": "workflows/build_taxonomy_fixture.py", "commit": None, "seed": SEED, "algorithm": "fixed-seed mock references and exact paired-read extraction", "parameters": {"reference_length": 650, "read_length": READ_LENGTH, "insert_length": INSERT_LENGTH, "positive_pairs_per_taxon": len(POSITIVE_STARTS), "shared_pairs": len(SHARED_STARTS), "background_pairs": len(BACKGROUND_STARTS)}, "versions": {"python": "3.12", "openomicsbench_contract": "2.0"}},
        "files": inventory,
        "validation": {"profile": "expected/validation.json", "baseline_version": "paired-read-taxonomy-v1", "metrics": []},
        "sequence": {"format": "fastq", "molecule": "dna", "paired": True, "quality_encoding": "phred33"},
        "limitations": ["The three reference labels are project-local mock codes, not NCBI taxonomy identifiers or biological species.", "Reads are error-free and short. The fixture tests result parsing, pair handling, exact assignment and unclassified cases, not real metagenomic profiling sensitivity.", "Shared-reference and absent-reference pairs have an unclassified truth label of 0."],
    }
    write_json(folder / "manifest.json", manifest)
    return folder


if __name__ == "__main__":
    target = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT
    built = build(target)
    print(f"Built {built}")
