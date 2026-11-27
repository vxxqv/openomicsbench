"""Descriptive QC for integer bulk RNA-seq count matrices."""

from __future__ import annotations

import csv
from pathlib import Path

from .rnaseq_design import _matrix, _table, check_rnaseq_design


def summarize_counts(counts: str | Path, samples: str | Path | None = None,
                     contrasts: str | Path | None = None, benchmark=None,
                     mitochondrial_prefix: str = "MT-") -> dict:
    """Summarize counts without imposing biological QC thresholds."""
    if (samples is None) != (contrasts is None):
        raise ValueError("samples and contrasts must be supplied together")
    if benchmark is not None and samples is None:
        raise ValueError("benchmark comparison requires samples and contrasts")
    if not mitochondrial_prefix:
        raise ValueError("mitochondrial prefix must not be empty")
    counts = Path(counts)
    matrix = _matrix(counts)
    names = matrix["samples"]
    detected = [0] * len(names)
    maxima = [0] * len(names)
    mitochondrial = [0] * len(names)
    mitochondrial_features = 0
    stream, delimiter = _table(counts)
    with stream:
        reader = csv.reader(stream, delimiter=delimiter)
        next(reader)
        for row in reader:
            is_mitochondrial = row[0].startswith(mitochondrial_prefix)
            mitochondrial_features += is_mitochondrial
            for i, cell in enumerate(row[1:]):
                value = int(cell)
                detected[i] += value > 0
                maxima[i] = max(maxima[i], value)
                if is_mitochondrial:
                    mitochondrial[i] += value
    observations = []
    for i, name in enumerate(names):
        total = matrix["totals"][i]
        observations.append({
            "sample": name,
            "library_size": total,
            "detected_features": detected[i],
            "zero_fraction": (matrix["features"] - detected[i]) / matrix["features"],
            "largest_feature_fraction": maxima[i] / total if total else None,
            "mitochondrial_fraction": mitochondrial[i] / total if mitochondrial_features and total else None,
        })
    positive = [item for item in matrix["totals"] if item > 0]
    design = (check_rnaseq_design(counts, samples, contrasts, benchmark)
              if samples is not None and contrasts is not None else None)
    failed = bool(len(positive) != len(names) or
                  (design is not None and design["summary"]["status"] == "fail"))
    return {
        "operation": "rnaseq_qc", "benchmark_id": benchmark.id if benchmark else None,
        "counts": str(counts), "counts_sha256": matrix["sha256"],
        "features": matrix["features"], "samples": len(names),
        "mitochondrial_prefix": mitochondrial_prefix,
        "mitochondrial_features_matched": mitochondrial_features,
        "library_size_ratio": max(positive) / min(positive) if positive else None,
        "sample_qc": observations, "design_check": design,
        "summary": {"status": "fail" if failed else "pass", "total": 1,
                    "passed": 0 if failed else 1, "failed": 1 if failed else 0,
                    "skipped": 0},
        "scope": "Descriptive count-level QC. No universal biological thresholds, normalization, differential expression, or read-level checks are applied.",
    }
