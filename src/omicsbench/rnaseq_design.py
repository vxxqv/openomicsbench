"""Check count matrices, sample sheets and differential-expression contrasts together."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np

from .hashing import digest


def _table(path: Path):
    suffixes = path.suffixes
    compressed = suffixes[-1:] == [".gz"]
    suffix = suffixes[-2] if compressed and len(suffixes) > 1 else path.suffix
    if suffix not in {".tsv", ".csv"}:
        raise ValueError(f"{path}: expected .tsv or .csv, optionally gzip-compressed")
    opener = gzip.open if compressed else open
    return opener(path, "rt", encoding="utf-8-sig", newline=""), "\t" if suffix == ".tsv" else ","


def _matrix(path: Path) -> dict:
    stream, delimiter = _table(path)
    with stream:
        reader = csv.reader(stream, delimiter=delimiter)
        header = next(reader, None)
        if not header or len(header) < 3 or header[0] not in {"gene_id", "Geneid", "GeneID"}:
            raise ValueError(f"{path}: expected gene_id followed by at least two sample columns")
        samples = header[1:]
        if any(not item or item != item.strip() for item in samples) or len(samples) != len(set(samples)):
            raise ValueError(f"{path}: sample IDs are empty, padded or duplicated")
        totals = [0] * len(samples)
        features = set()
        for line, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"{path}:{line}: row has {len(row)} columns, expected {len(header)}")
            gene = row[0]
            if not gene or gene != gene.strip() or gene in features:
                raise ValueError(f"{path}:{line}: empty, padded or duplicate gene ID")
            features.add(gene)
            for index, cell in enumerate(row[1:]):
                if not cell.isascii() or not cell.isdigit():
                    raise ValueError(f"{path}:{line}: counts must be non-negative integers")
                totals[index] += int(cell)
        if not features:
            raise ValueError(f"{path}: count matrix has no features")
    return {"samples": samples, "features": len(features), "totals": totals,
            "sha256": digest(path), "bytes": path.stat().st_size}


def _rows(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    stream, delimiter = _table(path)
    with stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)) or not required.issubset(fields):
            raise ValueError(f"{path}: required columns are {', '.join(sorted(required))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: table has no rows")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"{path}: inconsistent row width")
    return fields, rows


def _rank(samples: list[dict[str, str]], variable: str, target: str,
          blocks: list[str]) -> tuple[int, int]:
    matrix = [[1.0, float(row[variable] == target)] for row in samples]
    for block in blocks:
        levels = sorted({row[block] for row in samples})
        for level in levels[1:]:
            for vector, row in zip(matrix, samples):
                vector.append(float(row[block] == level))
    design = np.asarray(matrix, dtype=np.float64)
    return int(np.linalg.matrix_rank(design)), int(design.shape[1])


def check_rnaseq_design(counts: str | Path, samples: str | Path, contrasts: str | Path,
                        benchmark=None) -> dict:
    """Report structural mismatches and identifiability before a DE workflow runs."""
    counts, samples, contrasts = map(Path, (counts, samples, contrasts))
    matrix = _matrix(counts)
    sample_fields, sample_rows = _rows(samples, {"sample"})
    contrast_fields, contrast_rows = _rows(contrasts, {"id", "variable", "reference", "target"})
    sample_ids = [row["sample"] for row in sample_rows]
    if any(not value or value != value.strip() for value in sample_ids) or len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"{samples}: sample IDs are empty, padded or duplicated")
    contrast_ids = [row["id"] for row in contrast_rows]
    if any(not value or value != value.strip() for value in contrast_ids) or len(contrast_ids) != len(set(contrast_ids)):
        raise ValueError(f"{contrasts}: contrast IDs are empty, padded or duplicated")
    missing = sorted(set(matrix["samples"]) - set(sample_ids))
    extra = sorted(set(sample_ids) - set(matrix["samples"]))
    problems = []
    if missing or extra:
        problems.append("Matrix and sample sheet IDs differ")
    empty_libraries = [name for name, total in zip(matrix["samples"], matrix["totals"]) if total == 0]
    if empty_libraries:
        problems.append("One or more samples have zero total counts")
    case_reports = []
    for row in contrast_rows:
        issues = []
        target_count = reference_count = 0
        rank = columns = None
        variable, target, reference = (row[key] for key in ("variable", "target", "reference"))
        blocks = [value.strip() for value in row.get("blocking", "").split(";") if value.strip()]
        if not variable or variable not in sample_fields:
            issues.append(f"Unknown contrast variable: {variable}")
        if not target or not reference or target == reference:
            issues.append("Reference and target must be distinct nonempty levels")
        if len(blocks) != len(set(blocks)) or variable in blocks:
            issues.append("Blocking columns must be unique and distinct from the contrast variable")
        for block in blocks:
            if block not in sample_fields:
                issues.append(f"Unknown blocking column: {block}")
        if not issues:
            levels = {item[variable] for item in sample_rows}
            if {target, reference} - levels:
                issues.append("Reference or target is absent from the sample sheet")
            if levels - {target, reference}:
                issues.append("This checker requires the contrast variable to have exactly the two declared levels")
            subset = [item for item in sample_rows if item[variable] in {target, reference}]
            target_count = sum(item[variable] == target for item in subset)
            reference_count = sum(item[variable] == reference for item in subset)
            if min(target_count, reference_count) < 2:
                issues.append("Each contrast level needs at least two samples")
            if any(not item[block] for item in subset for block in blocks):
                issues.append("Blocking values must be present for every sample")
            if not issues:
                rank, columns = _rank(subset, variable, target, blocks)
                if rank < columns:
                    issues.append("Contrast and blocking variables are confounded")
                if len(subset) <= columns:
                    issues.append("Design has no residual degrees of freedom")
        case_reports.append({"id": row["id"], "status": "fail" if issues else "pass",
                             "target": target, "reference": reference, "variable": variable,
                             "blocking": blocks, "target_samples": target_count,
                             "reference_samples": reference_count, "rank": rank, "columns": columns,
                             "issues": issues})
    if benchmark is not None:
        from .rnaseq_pack import _design
        source = benchmark.rnaseq()
        expected = _design(benchmark)
        expected_samples = {row["sample_id"]: row for row in source.metadata}
        if set(matrix["samples"]) != set(source.samples):
            problems.append(f"{benchmark.id}: matrix sample IDs differ from the benchmark")
        if set(sample_ids) != set(source.samples):
            problems.append(f"{benchmark.id}: sample sheet IDs differ from the benchmark")
        for row in sample_rows:
            original = expected_samples.get(row["sample"])
            if original and any(row.get(key) != str(original[key]) for key in
                                ("condition", "batch", "replicate", "strandedness") if key in row):
                problems.append(f"{benchmark.id}: sample metadata differs for {row['sample']}")
        expected_contrast = {key: expected[key] for key in ("id", "variable", "reference", "target", "blocking")}
        if len(contrast_rows) != 1 or any(contrast_rows[0].get(key, "") != value for key, value in expected_contrast.items()):
            problems.append(f"{benchmark.id}: contrast differs from the declared reference")
    failed = int(bool(problems)) + sum(item["status"] == "fail" for item in case_reports)
    total = 1 + len(case_reports)
    return {
        "operation": "check_rnaseq_design", "benchmark_id": benchmark.id if benchmark else None,
        "inputs": {"counts": {"path": str(counts), "sha256": matrix["sha256"], "bytes": matrix["bytes"]},
                   "samples": {"path": str(samples), "sha256": digest(samples)},
                   "contrasts": {"path": str(contrasts), "sha256": digest(contrasts)}},
        "features": matrix["features"], "matrix_samples": len(matrix["samples"]),
        "sample_sheet_rows": len(sample_rows), "sample_order_matches": matrix["samples"] == sample_ids,
        "missing_sample_ids": missing, "extra_sample_ids": extra,
        "empty_libraries": empty_libraries, "sample_totals": dict(zip(matrix["samples"], matrix["totals"])),
        "issues": problems, "contrasts": case_reports,
        "summary": {"status": "fail" if failed else "pass", "total": total,
                    "passed": total - failed,
                    "failed": failed, "skipped": 0},
        "scope": "Checks IDs, integer counts and a two-level fixed-effect design; it does not run differential expression.",
    }
