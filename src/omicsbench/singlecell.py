"""Strict single-cell matrix checks and label-invariant clustering scores."""
from __future__ import annotations

import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

from .hashing import contained


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).lower().endswith(".gz") else path.open("r", encoding="utf-8", newline="")


def read_lines(path: Path) -> list[str]:
    with _open_text(Path(path)) as stream:
        values = [line.rstrip("\r\n") for line in stream]
    if not values or any(not value for value in values):
        raise ValueError(f"{path}: expected non-empty lines")
    return values


def read_features(path: Path, expected_columns: int | None = None) -> list[tuple[str, str]]:
    rows = []
    for line_number, line in enumerate(read_lines(Path(path)), start=1):
        fields = line.split("\t")
        if expected_columns is not None and len(fields) != expected_columns:
            raise ValueError(f"{path}: line {line_number} has {len(fields)} fields; expected {expected_columns}")
        if len(fields) not in {2, 3} or not fields[0] or not fields[1]:
            raise ValueError(f"{path}: invalid feature row at line {line_number}")
        rows.append((fields[0], fields[1]))
    identifiers = [row[0] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"{path}: feature identifiers must be unique")
    return rows


def read_matrix_market(path: Path) -> dict:
    with _open_text(Path(path)) as stream:
        try:
            header = next(stream).rstrip("\r\n")
        except StopIteration:
            raise ValueError(f"{path}: empty Matrix Market file") from None
        if header != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"{path}: expected an integer coordinate Matrix Market matrix")
        line_number = 1
        dimensions = None
        entries = []
        for raw in stream:
            line_number += 1
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            fields = line.split()
            if dimensions is None:
                if len(fields) != 3:
                    raise ValueError(f"{path}: invalid dimension line")
                try:
                    dimensions = tuple(int(value) for value in fields)
                except ValueError:
                    raise ValueError(f"{path}: matrix dimensions must be integers") from None
                if any(value < 0 for value in dimensions):
                    raise ValueError(f"{path}: matrix dimensions cannot be negative")
                continue
            if len(fields) != 3:
                raise ValueError(f"{path}: invalid coordinate at line {line_number}")
            try:
                row, column, count = (int(value) for value in fields)
            except ValueError:
                raise ValueError(f"{path}: coordinates and counts must be integers at line {line_number}") from None
            entries.append((row, column, count))
    if dimensions is None:
        raise ValueError(f"{path}: matrix dimensions are missing")
    rows, columns, declared_nnz = dimensions
    if len(entries) != declared_nnz:
        raise ValueError(f"{path}: declared {declared_nnz} entries but found {len(entries)}")
    coordinates = set()
    for row, column, count in entries:
        if not (1 <= row <= rows and 1 <= column <= columns):
            raise ValueError(f"{path}: coordinate outside the declared matrix shape")
        if count <= 0:
            raise ValueError(f"{path}: sparse counts must be positive integers")
        if (row, column) in coordinates:
            raise ValueError(f"{path}: duplicate matrix coordinate {row},{column}")
        coordinates.add((row, column))
    return {"features": rows, "cells": columns, "nonzero": declared_nnz, "entries": entries}


def summarize_matrix(matrix_path: Path, barcodes_path: Path, features_path: Path, mitochondrial_prefix: str = "MT-", feature_columns: int | None = None) -> dict:
    matrix = read_matrix_market(Path(matrix_path))
    barcodes = read_lines(Path(barcodes_path))
    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"{barcodes_path}: barcodes must be unique")
    features = read_features(Path(features_path), feature_columns)
    if matrix["cells"] != len(barcodes) or matrix["features"] != len(features):
        raise ValueError("matrix dimensions do not match the barcode and feature tables")
    totals = [0] * len(barcodes)
    detected = [0] * len(barcodes)
    mitochondrial = [0] * len(barcodes)
    for row, column, count in matrix["entries"]:
        totals[column - 1] += count
        detected[column - 1] += 1
        if features[row - 1][1].upper().startswith(mitochondrial_prefix.upper()):
            mitochondrial[column - 1] += count
    ordered_totals = sorted(totals)
    ordered_detected = sorted(detected)
    midpoint = len(barcodes) // 2
    median = lambda values: float(values[midpoint]) if len(values) % 2 else (values[midpoint - 1] + values[midpoint]) / 2
    total_counts = sum(totals)
    return {
        "cells": len(barcodes),
        "features": len(features),
        "nonzero": matrix["nonzero"],
        "total_counts": total_counts,
        "counts_per_cell": {"minimum": min(totals), "median": median(ordered_totals), "maximum": max(totals)},
        "features_per_cell": {"minimum": min(detected), "median": median(ordered_detected), "maximum": max(detected)},
        "mitochondrial_fraction": sum(mitochondrial) / total_counts if total_counts else 0.0,
        "zero_count_cells": sum(value == 0 for value in totals),
    }


def read_labels(path: Path, columns: dict[str, str] | None = None) -> dict[str, str]:
    mapping = {"cell_id": "cell_id", "label": "label", **(columns or {})}
    delimiter = "," if ".csv" in str(path).lower() else "\t"
    with _open_text(Path(path)) as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        if reader.fieldnames is None or any(source not in reader.fieldnames for source in mapping.values()):
            raise ValueError(f"{path}: expected columns {mapping['cell_id']} and {mapping['label']}")
        labels = {}
        for line_number, row in enumerate(reader, start=2):
            cell = (row.get(mapping["cell_id"]) or "").strip()
            label = (row.get(mapping["label"]) or "").strip()
            if not cell or not label:
                raise ValueError(f"{path}: blank cell or label at line {line_number}")
            if cell in labels:
                raise ValueError(f"{path}: duplicate cell identifier {cell}")
            labels[cell] = label
    if not labels:
        raise ValueError(f"{path}: no label rows")
    return labels


def _comb2(value: int) -> float:
    return value * (value - 1) / 2


def clustering_metrics(truth: list[str], predicted: list[str]) -> dict[str, float | int]:
    if len(truth) != len(predicted) or not truth:
        raise ValueError("truth and predicted labels require the same non-zero length")
    table = Counter(zip(truth, predicted))
    true_counts = Counter(truth)
    predicted_counts = Counter(predicted)
    pairs = _comb2(len(truth))
    index = sum(_comb2(value) for value in table.values())
    true_sum = sum(_comb2(value) for value in true_counts.values())
    predicted_sum = sum(_comb2(value) for value in predicted_counts.values())
    if pairs == 0:
        ari = 1.0
    else:
        expected = true_sum * predicted_sum / pairs
        denominator = (true_sum + predicted_sum) / 2 - expected
        ari = (index - expected) / denominator if denominator else (1.0 if truth == predicted else 0.0)
    total = len(truth)
    mutual_information = 0.0
    for (true_label, predicted_label), count in table.items():
        mutual_information += count / total * math.log(count * total / (true_counts[true_label] * predicted_counts[predicted_label]))
    true_entropy = -sum(count / total * math.log(count / total) for count in true_counts.values())
    predicted_entropy = -sum(count / total * math.log(count / total) for count in predicted_counts.values())
    entropy_mean = (true_entropy + predicted_entropy) / 2
    nmi = mutual_information / entropy_mean if entropy_mean else 1.0
    by_cluster = defaultdict(Counter)
    for true_label, predicted_label in zip(truth, predicted):
        by_cluster[predicted_label][true_label] += 1
    purity = sum(max(values.values()) for values in by_cluster.values()) / total
    return {"adjusted_rand_index": ari, "normalized_mutual_information": nmi, "purity": purity, "truth_clusters": len(true_counts), "predicted_clusters": len(predicted_counts)}


def compare_labels(model, folder: Path, results: Path, detail_limit: int = 20, columns: dict[str, str] | None = None) -> dict:
    if model.assay != "single_cell_rna_seq" or model.validation is None:
        raise ValueError(f"{model.id}: no single-cell clustering truth")
    import json
    profile = json.loads(contained(folder, model.validation.profile).read_text(encoding="utf-8"))
    truth = read_labels(contained(folder, profile["label_truth"]))
    observed = read_labels(Path(results), columns)
    missing = sorted(set(truth) - set(observed))
    unexpected = sorted(set(observed) - set(truth))
    shared = [cell for cell in truth if cell in observed]
    metrics = clustering_metrics([truth[cell] for cell in shared], [observed[cell] for cell in shared]) if shared else {"adjusted_rand_index": 0.0, "normalized_mutual_information": 0.0, "purity": 0.0, "truth_clusters": len(set(truth.values())), "predicted_clusters": 0}
    metrics["coverage"] = len(shared) / len(truth)
    thresholds = profile["clustering_comparison"]
    reasons = []
    for name in ("adjusted_rand_index", "normalized_mutual_information", "purity", "coverage"):
        minimum = thresholds[f"minimum_{name}"]
        if metrics[name] < minimum:
            reasons.append(f"{name}={metrics[name]:.6g} is below {minimum}")
    if unexpected:
        reasons.append(f"{len(unexpected)} unexpected cells")
    return {
        "id": model.id,
        "status": "fail" if reasons else "pass",
        "input": str(results),
        "input_columns": {"cell_id": "cell_id", "label": "label", **(columns or {})},
        "cells": {"expected": len(truth), "observed": len(observed), "shared": len(shared), "missing": len(missing), "unexpected": len(unexpected)},
        "metrics": metrics,
        "thresholds": thresholds,
        "missing_examples": missing[:detail_limit],
        "unexpected_examples": unexpected[:detail_limit],
        "reasons": reasons,
        "scope": "label-invariant recovery of synthetic cell states; this does not validate biological cell-type annotation",
    }
