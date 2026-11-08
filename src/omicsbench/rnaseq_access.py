"""Checked access to the biological bulk RNA-seq benchmarks."""

from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from itertools import islice
from pathlib import Path

import numpy as np

from .hashing import digest
from .rnaseq import read_counts


@dataclass(frozen=True)
class RNASeqData:
    """Feature-by-sample integer counts with ordered sample records."""

    id: str
    tier: str
    genes: tuple[str, ...]
    samples: tuple[str, ...]
    counts: np.ndarray
    metadata: tuple[dict, ...]
    source: dict
    limitations: tuple[str, ...]


def load_rnaseq(benchmark, tier: str = "pocket") -> RNASeqData:
    if benchmark.assay != "bulk_rna_seq" or benchmark.manifest["status"] != "validated" or benchmark.manifest["kind"] != "real":
        raise ValueError(f"{benchmark.id}: only validated biological bulk RNA-seq objects are supported")
    if tier not in {"pocket", "full"}:
        raise ValueError("tier must be 'pocket' or 'full'")
    role, file_tier = ("raw_counts", "pocket") if tier == "pocket" else ("baseline", "expected")
    genes, samples, counts = read_counts(benchmark.file(role=role, tier=file_tier))
    metadata = benchmark.manifest["samples"]
    if samples != [item["sample_id"] for item in metadata]:
        raise ValueError(f"{benchmark.id}: count columns and declared sample order differ")
    counts.setflags(write=False)
    return RNASeqData(benchmark.id, tier, tuple(genes), tuple(samples), counts,
                      tuple(dict(item) for item in metadata), dict(benchmark.manifest["source"]),
                      tuple(benchmark.manifest["limitations"]))


def export_rnaseq(benchmark, destination: str | Path, tier: str = "pocket") -> dict:
    """Write aligned counts, sample table and source details into a new directory."""
    data = load_rnaseq(benchmark, tier)
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f"{destination}: export directory already exists")
    if not destination.parent.is_dir():
        raise FileNotFoundError(f"{destination.parent}: parent directory does not exist")
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}-", dir=destination.parent))
    try:
        with (temporary / "counts.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", *data.samples])
            writer.writerows((gene, *row) for gene, row in zip(data.genes, data.counts.tolist()))
        fields = list(data.metadata[0])
        with (temporary / "samples.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(data.metadata)
        files = {name: {"bytes": (temporary / name).stat().st_size, "sha256": digest(temporary / name)}
                 for name in ("counts.tsv", "samples.tsv")}
        details = {"benchmark_id": data.id, "tier": data.tier, "orientation": "genes_by_samples",
                   "features": len(data.genes), "samples": len(data.samples), "source": data.source,
                   "limitations": list(data.limitations), "files": files}
        (temporary / "export.json").write_text(json.dumps(details, indent=2) + "\n", encoding="utf-8", newline="\n")
        os.replace(temporary, destination)
        return details
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def compare_counts(benchmark, path: str | Path, tier: str = "pocket", detail_limit: int = 20) -> dict:
    """Check an integer matrix against a declared count tier by gene and sample ID."""
    if detail_limit < 0:
        raise ValueError("detail_limit must be non-negative")
    reference = load_rnaseq(benchmark, tier)
    path = Path(path)
    genes, samples, submitted = read_counts(path)
    reference_genes = set(reference.genes)
    reference_samples = set(reference.samples)
    input_genes = set(genes)
    input_samples = set(samples)
    shared_genes = [gene for gene in reference.genes if gene in input_genes]
    shared_samples = [sample for sample in reference.samples if sample in input_samples]
    missing_genes = sorted(reference_genes - input_genes)
    extra_genes = sorted(input_genes - reference_genes)
    missing_samples = sorted(reference_samples - input_samples)
    extra_samples = sorted(input_samples - reference_samples)
    reference_gene_index = {gene: index for index, gene in enumerate(reference.genes)}
    reference_sample_index = {sample: index for index, sample in enumerate(reference.samples)}
    input_gene_index = {gene: index for index, gene in enumerate(genes)}
    input_sample_index = {sample: index for index, sample in enumerate(samples)}
    left = reference.counts[np.ix_([reference_gene_index[g] for g in shared_genes],
                                   [reference_sample_index[s] for s in shared_samples])]
    right = submitted[np.ix_([input_gene_index[g] for g in shared_genes],
                              [input_sample_index[s] for s in shared_samples])]
    different = left != right
    examples = [
        {"gene_id": shared_genes[row], "sample_id": shared_samples[column],
         "reference": int(left[row, column]), "submitted": int(right[row, column])}
        for row, column in islice(zip(*np.where(different)), detail_limit)
    ]
    sample_totals = [
        {"sample_id": sample, "reference_shared_gene_total": int(left[:, column].sum()),
         "submitted_shared_gene_total": int(right[:, column].sum())}
        for column, sample in enumerate(shared_samples)
    ]
    mismatch_count = int(np.count_nonzero(different))
    complete = not (missing_genes or extra_genes or missing_samples or extra_samples)
    return {
        "id": reference.id, "tier": reference.tier,
        "status": "pass" if complete and mismatch_count == 0 else "fail",
        "scope": "Exact integer count agreement after aligning gene and sample IDs; no biological method performance is assessed.",
        "submitted_sha256": digest(path),
        "reference_features": len(reference.genes), "submitted_features": len(genes),
        "reference_samples": len(reference.samples), "submitted_samples": len(samples),
        "shared_features": len(shared_genes), "shared_samples": len(shared_samples),
        "missing_genes": missing_genes[:detail_limit], "extra_genes": extra_genes[:detail_limit],
        "missing_samples": missing_samples[:detail_limit], "extra_samples": extra_samples[:detail_limit],
        "missing_gene_count": len(missing_genes), "extra_gene_count": len(extra_genes),
        "missing_sample_count": len(missing_samples), "extra_sample_count": len(extra_samples),
        "mismatched_cells": mismatch_count,
        "mismatched_genes": int(np.count_nonzero(np.any(different, axis=1))),
        "mismatched_samples": int(np.count_nonzero(np.any(different, axis=0))),
        "examples": examples, "sample_totals": sample_totals,
    }


def verify_export(benchmark, directory: str | Path) -> dict:
    """Check an exported RNA-seq directory against its manifest and source object."""
    directory = Path(directory)
    details = json.loads((directory / "export.json").read_text(encoding="utf-8"))
    if not isinstance(details, dict) or details.get("tier") not in {"pocket", "full"}:
        raise ValueError(f"{directory}: invalid export.json")
    data = load_rnaseq(benchmark, details["tier"])
    issues = []
    expected_names = {"counts.tsv", "samples.tsv", "export.json"}
    actual_names = {path.name for path in directory.iterdir()}
    if actual_names != expected_names:
        issues.append(f"Export file set differs: missing={sorted(expected_names - actual_names)}, extra={sorted(actual_names - expected_names)}")
    expected_fields = {"benchmark_id": data.id, "orientation": "genes_by_samples",
                       "features": len(data.genes), "samples": len(data.samples),
                       "source": data.source, "limitations": list(data.limitations)}
    for field, expected in expected_fields.items():
        if details.get(field) != expected:
            issues.append(f"export.json field differs: {field}")
    declared_files = details.get("files")
    if not isinstance(declared_files, dict) or set(declared_files) != {"counts.tsv", "samples.tsv"}:
        issues.append("export.json has an invalid file inventory")
        declared_files = {}
    for name in ("counts.tsv", "samples.tsv"):
        path = directory / name
        if not path.is_file():
            issues.append(f"Missing {name}")
            continue
        record = declared_files.get(name)
        if not isinstance(record, dict) or record.get("bytes") != path.stat().st_size or record.get("sha256") != digest(path):
            issues.append(f"{name} differs from its recorded size or SHA-256")
    samples_path = directory / "samples.tsv"
    if samples_path.is_file():
        try:
            with samples_path.open("r", encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                rows = list(reader)
                fields = list(data.metadata[0])
                if reader.fieldnames != fields or len(rows) != len(data.metadata) or any(
                    row != {field: str(expected[field]) for field in fields}
                    for row, expected in zip(rows, data.metadata)
                ):
                    issues.append("samples.tsv differs from the declared sample records or order")
        except (UnicodeError, csv.Error):
            issues.append("samples.tsv cannot be read as a UTF-8 tab-separated table")
    comparison = None
    counts_path = directory / "counts.tsv"
    if counts_path.is_file():
        try:
            comparison = compare_counts(benchmark, counts_path, details["tier"], detail_limit=5)
            if comparison["status"] != "pass":
                issues.append("counts.tsv does not exactly match the declared count tier")
        except (ValueError, UnicodeError):
            issues.append("counts.tsv is not a valid integer count matrix")
    return {"id": data.id, "tier": data.tier, "status": "pass" if not issues else "fail",
            "issues": issues, "count_comparison": comparison}
