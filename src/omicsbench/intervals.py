"""Strict BED-family parsing and overlap scores for genomic peak sets."""
from __future__ import annotations

import gzip
import math
from dataclasses import dataclass
from pathlib import Path

from .hashing import contained


@dataclass(frozen=True, order=True)
class Interval:
    chrom: str
    start: int
    end: int
    name: str = "."
    summit: int | None = None


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).lower().endswith(".gz") else Path(path).open("r", encoding="utf-8", newline="")


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with _open_text(Path(path)) as stream:
        for number, raw in enumerate(stream, 1):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0] or fields[0] in sizes:
                raise ValueError(f"{path}: invalid chromosome-size row at line {number}")
            try:
                size = int(fields[1])
            except ValueError:
                raise ValueError(f"{path}: chromosome size is not an integer at line {number}") from None
            if size <= 0:
                raise ValueError(f"{path}: chromosome size must be positive at line {number}")
            sizes[fields[0]] = size
    if not sizes:
        raise ValueError(f"{path}: no chromosome sizes")
    return sizes


def read_intervals(path: Path, chrom_sizes: dict[str, int] | None = None) -> tuple[list[Interval], int]:
    intervals: list[Interval] = []
    field_count = None
    seen = set()
    with _open_text(Path(path)) as stream:
        for number, raw in enumerate(stream, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            if line.startswith(("track", "browser")):
                raise ValueError(f"{path}: track and browser headers are not accepted")
            fields = line.split("\t")
            if len(fields) not in {3, 6, 10}:
                raise ValueError(f"{path}: line {number} must contain 3, 6 or 10 tab-separated fields")
            if field_count is None:
                field_count = len(fields)
            elif len(fields) != field_count:
                raise ValueError(f"{path}: BED field count changes at line {number}")
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise ValueError(f"{path}: coordinates are not integers at line {number}") from None
            if not chrom or start < 0 or end <= start:
                raise ValueError(f"{path}: invalid positive half-open interval at line {number}")
            if chrom_sizes is not None and (chrom not in chrom_sizes or end > chrom_sizes[chrom]):
                raise ValueError(f"{path}: interval is outside the declared reference at line {number}")
            name, summit = (fields[3], None) if len(fields) >= 6 else (".", None)
            if len(fields) >= 6:
                try:
                    score = int(fields[4])
                except ValueError:
                    raise ValueError(f"{path}: BED score is not an integer at line {number}") from None
                if score < 0 or score > 1000 or fields[5] not in {"+", "-", "."}:
                    raise ValueError(f"{path}: invalid BED score or strand at line {number}")
            if len(fields) == 10:
                try:
                    signal, pvalue, qvalue = (float(value) for value in fields[6:9])
                    offset = int(fields[9])
                except ValueError:
                    raise ValueError(f"{path}: invalid narrowPeak fields at line {number}") from None
                if any(not math.isfinite(value) or value < -1 for value in (signal, pvalue, qvalue)) or offset < -1 or offset >= end - start:
                    raise ValueError(f"{path}: narrowPeak values are outside their allowed range at line {number}")
                summit = start + offset if offset >= 0 else None
            key = (chrom, start, end, name)
            if key in seen:
                raise ValueError(f"{path}: duplicate interval at line {number}")
            seen.add(key)
            intervals.append(Interval(chrom, start, end, name, summit))
    if not intervals:
        raise ValueError(f"{path}: no intervals")
    return intervals, int(field_count)


def _merged(intervals: list[Interval]) -> list[Interval]:
    output: list[Interval] = []
    for current in sorted(intervals):
        if output and current.chrom == output[-1].chrom and current.start <= output[-1].end:
            previous = output[-1]
            output[-1] = Interval(previous.chrom, previous.start, max(previous.end, current.end))
        else:
            output.append(Interval(current.chrom, current.start, current.end))
    return output


def _intersection_bp(first: list[Interval], second: list[Interval]) -> int:
    left, right, total = _merged(first), _merged(second), 0
    i = j = 0
    while i < len(left) and j < len(right):
        if left[i].chrom < right[j].chrom or (left[i].chrom == right[j].chrom and left[i].end <= right[j].start):
            i += 1
        elif right[j].chrom < left[i].chrom or (left[i].chrom == right[j].chrom and right[j].end <= left[i].start):
            j += 1
        else:
            total += max(0, min(left[i].end, right[j].end) - max(left[i].start, right[j].start))
            if left[i].end <= right[j].end:
                i += 1
            else:
                j += 1
    return total


def summarize_intervals(path: Path, chrom_sizes_path: Path | None = None) -> dict:
    sizes = read_chrom_sizes(chrom_sizes_path) if chrom_sizes_path else None
    values, fields = read_intervals(Path(path), sizes)
    merged = _merged(values)
    widths = sorted(item.end - item.start for item in values)
    middle = len(widths) // 2
    median = float(widths[middle]) if len(widths) % 2 else (widths[middle - 1] + widths[middle]) / 2
    return {
        "format": {3: "BED3", 6: "BED6", 10: "narrowPeak"}[fields],
        "intervals": len(values),
        "contigs": len({item.chrom for item in values}),
        "span_bp": sum(widths),
        "union_bp": sum(item.end - item.start for item in merged),
        "overlapping_intervals": len(values) - len(merged),
        "width_bp": {"minimum": min(widths), "median": median, "maximum": max(widths)},
        "sorted": values == sorted(values),
    }


def interval_metrics(truth: list[Interval], calls: list[Interval], reciprocal_overlap: float) -> dict:
    truth_hit, call_hit = set(), set()
    summit_distances = []
    for ti, expected in enumerate(truth):
        for ci, observed in enumerate(calls):
            if expected.chrom != observed.chrom:
                continue
            overlap = max(0, min(expected.end, observed.end) - max(expected.start, observed.start))
            if overlap / (expected.end - expected.start) >= reciprocal_overlap and overlap / (observed.end - observed.start) >= reciprocal_overlap:
                truth_hit.add(ti)
                call_hit.add(ci)
                if expected.summit is not None and observed.summit is not None:
                    summit_distances.append(abs(expected.summit - observed.summit))
    precision = len(call_hit) / len(calls)
    recall = len(truth_hit) / len(truth)
    intersection = _intersection_bp(truth, calls)
    truth_bp = sum(item.end - item.start for item in _merged(truth))
    call_bp = sum(item.end - item.start for item in _merged(calls))
    base_precision = intersection / call_bp
    base_recall = intersection / truth_bp
    f1 = lambda p, r: 2 * p * r / (p + r) if p + r else 0.0
    return {
        "interval_precision": precision,
        "interval_recall": recall,
        "interval_f1": f1(precision, recall),
        "base_precision": base_precision,
        "base_recall": base_recall,
        "base_f1": f1(base_precision, base_recall),
        "truth_intervals": len(truth),
        "called_intervals": len(calls),
        "matched_truth": len(truth_hit),
        "matched_calls": len(call_hit),
        "intersection_bp": intersection,
        "median_summit_distance_bp": (sorted(summit_distances)[len(summit_distances) // 2] if summit_distances else None),
    }


def compare_intervals(model, folder: Path, results: Path) -> dict:
    if model.intervals is None or model.validation is None:
        raise ValueError(f"{model.id}: no genomic interval truth")
    import json
    profile = json.loads(contained(folder, model.validation.profile).read_text(encoding="utf-8"))
    sizes = read_chrom_sizes(contained(folder, profile["chrom_sizes"]))
    truth, _ = read_intervals(contained(folder, profile["interval_truth"]), sizes)
    calls, _ = read_intervals(Path(results), sizes)
    metrics = interval_metrics(truth, calls, model.intervals.reciprocal_overlap)
    thresholds = profile["peak_comparison"]
    reasons = []
    for name in ("interval_precision", "interval_recall", "interval_f1", "base_precision", "base_recall", "base_f1"):
        minimum = thresholds[f"minimum_{name}"]
        if metrics[name] < minimum:
            reasons.append(f"{name}={metrics[name]:.6g} is below {minimum}")
    return {"id": model.id, "status": "fail" if reasons else "pass", "input": str(results), "metrics": metrics, "thresholds": thresholds, "reciprocal_overlap": model.intervals.reciprocal_overlap, "reasons": reasons, "scope": "coordinate and overlap recovery against deterministic synthetic peak truth; biological peak validity is outside scope"}
