"""Separate clustering, exact annotation and embedding diagnostics for small cell sets."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from .hashing import contained
from .singlecell import compare_labels, read_labels


def read_embedding(path: Path) -> dict[str, tuple[float, ...]]:
    path = Path(path)
    delimiter = "," if ".csv" in path.name.lower() else "\t"
    points: dict[str, tuple[float, ...]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        fields = reader.fieldnames or []
        if len(fields) < 3 or fields[0] != "cell_id" or len(set(fields)) != len(fields):
            raise ValueError(f"{path}: embedding needs cell_id and at least two distinct coordinate columns")
        for line, row in enumerate(reader, start=2):
            cell = (row.get("cell_id") or "").strip()
            if not cell or cell in points or None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{line}: blank, duplicate or malformed cell row")
            try:
                values = tuple(float(row[field]) for field in fields[1:])
            except ValueError:
                raise ValueError(f"{path}:{line}: embedding coordinates must be numeric") from None
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path}:{line}: embedding coordinates must be finite")
            points[cell] = values
    if not points:
        raise ValueError(f"{path}: no embedding rows")
    return points


def _batch_labels(path: Path) -> dict[str, str]:
    batches = {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if "batch" not in (reader.fieldnames or []):
            return {}
        for row in reader:
            batches[row["cell_id"]] = row["batch"]
    return batches


def _embedding_report(points: dict[str, tuple[float, ...]], truth: dict[str, str], batches: dict[str, str], k: int) -> dict:
    if set(points) != set(truth):
        raise ValueError("embedding cell IDs must exactly match benchmark truth")
    if len(points) < 3 or not 1 <= k < len(points):
        raise ValueError("neighbor count must be positive and smaller than the cell count")
    if len(points) > 5000:
        raise ValueError("embedding diagnostics are limited to 5,000 cells")
    ordered = list(truth)
    n = len(ordered)
    matrix = np.asarray([points[cell] for cell in ordered], dtype=float)
    labels = [truth[cell] for cell in ordered]
    batch_labels = [batches.get(cell) for cell in ordered]
    if batches and any(not value for value in batch_labels):
        raise ValueError("batch labels are missing for some cells")
    same_state = []
    same_batch = []
    for index in range(len(ordered)):
        distances = np.sum((matrix - matrix[index]) ** 2, axis=1)
        distances[index] = np.inf
        neighbors = np.lexsort((np.arange(n), distances))[:k]
        same_state.append(sum(labels[j] == labels[index] for j in neighbors) / k)
        if batches:
            same_batch.append(sum(batch_labels[j] == batch_labels[index] for j in neighbors) / k)
    state_counts = Counter(labels)
    report = {
        "cells": n,
        "dimensions": matrix.shape[1],
        "neighbors": k,
        "mean_same_state_neighbor_fraction": float(np.mean(same_state)),
        "random_state_neighbor_fraction": sum(count * (count - 1) for count in state_counts.values()) / (n * (n - 1)),
        "scope": "Descriptive local-neighborhood conservation, not proof of biological cell-type accuracy.",
    }
    if batches:
        batch_counts = Counter(batch_labels)
        report["mean_same_batch_neighbor_fraction"] = float(np.mean(same_batch))
        report["random_batch_neighbor_fraction"] = sum(count * (count - 1) for count in batch_counts.values()) / (n * (n - 1))
        report["batch_scope"] = "Compare the observed and random batch fractions alongside state conservation; neither is a standalone quality score."
    return report


def assess_cells(model, folder: Path, clusters: Path, *, annotation: Path | None = None,
                 embedding: Path | None = None, neighbors: int = 10) -> dict:
    if model.assay != "single_cell_rna_seq" or model.validation is None:
        raise ValueError(f"{model.id}: no single-cell assessment reference")
    profile = json.loads(contained(folder, model.validation.profile).read_text(encoding="utf-8"))
    truth_path = contained(folder, profile["label_truth"])
    truth = read_labels(truth_path)
    clustering = compare_labels(model, folder, clusters)
    report = {"operation": "cell_assess", "id": model.id, "clustering": clustering}
    reasons = list(clustering["reasons"])
    if annotation is not None:
        names = read_labels(annotation)
        missing = sorted(set(truth) - set(names))
        extra = sorted(set(names) - set(truth))
        label_counts = Counter(truth.values())
        correct = Counter(truth[cell] for cell in truth if names.get(cell) == truth[cell])
        report["annotation"] = {
            "input": str(annotation),
            "exact_label_accuracy": sum(correct.values()) / len(truth),
            "per_label_recall": {label: correct[label] / count for label, count in sorted(label_counts.items())},
            "missing_cells": len(missing), "unexpected_cells": len(extra),
            "scope": "Exact declared-label agreement; unlike clustering, label names are meaningful.",
        }
        if missing or extra:
            reasons.append("annotation cell IDs differ from benchmark truth")
        if sum(correct.values()) != len(truth):
            reasons.append("annotation labels differ from the declared reference")
    if embedding is not None:
        points = read_embedding(embedding)
        batches = _batch_labels(truth_path)
        report["embedding"] = _embedding_report(points, truth, batches, neighbors)
    report["status"] = "fail" if reasons else "pass"
    report["reasons"] = reasons
    return report
