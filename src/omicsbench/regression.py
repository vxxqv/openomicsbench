"""Detect status and metric regressions between benchmark reports."""
from __future__ import annotations

import json
import math
from pathlib import Path

from .suite import _summary


TRACKED_METRICS = ("spearman_logfc", "top_k_jaccard", "sign_concordance")


def _reject_constant(value: str):
    raise ValueError(f"non-finite JSON value {value}")


def read_report(path: Path) -> dict:
    try:
        report = json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=_reject_constant)
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ValueError(f"{path}: invalid JSON report: {error}") from None
    if not isinstance(report, dict) or report.get("operation") not in {"validate", "compare", "matrix"}:
        raise ValueError(f"{path}: unsupported OpenOmicsBench report")
    if not isinstance(report.get("results"), list):
        raise ValueError(f"{path}: report results must be a list")
    return report


def _case_key(operation: str, item: dict) -> str:
    dataset_id = item.get("id")
    if not isinstance(dataset_id, str) or not dataset_id:
        raise ValueError("report result is missing a benchmark ID")
    if operation == "matrix":
        method = item.get("method")
        if not isinstance(method, str) or not method:
            raise ValueError(f"{dataset_id}: matrix result is missing its method")
        return f"{method}/{dataset_id}"
    return dataset_id


def _cases(report: dict) -> dict[str, dict]:
    cases = {}
    for item in report["results"]:
        if not isinstance(item, dict):
            raise ValueError("report results must contain objects")
        key = _case_key(report["operation"], item)
        if key in cases:
            raise ValueError(f"duplicate report case: {key}")
        cases[key] = item
    return cases


def _measurements(item: dict) -> dict[str, float]:
    details = item.get("details")
    if not isinstance(details, dict):
        return {}
    metrics = details.get("metrics", {})
    genes = details.get("genes", {})
    values = {name: metrics.get(name) for name in TRACKED_METRICS}
    values["coverage"] = genes.get("coverage")
    result = {}
    for name, value in values.items():
        if value is None:
            continue
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
            raise ValueError(f"{item.get('id', 'case')}: {name} must be finite")
        result[name] = float(value)
    return result


def compare_reports(
    baseline_path: Path,
    candidate_path: Path,
    absolute_tolerance: float = 0.0,
    allow_missing: bool = False,
) -> dict:
    if not math.isfinite(absolute_tolerance) or absolute_tolerance < 0:
        raise ValueError("absolute tolerance must be finite and zero or greater")
    baseline = read_report(baseline_path)
    candidate = read_report(candidate_path)
    if baseline["operation"] != candidate["operation"]:
        raise ValueError("baseline and candidate report operations differ")
    before = _cases(baseline)
    after = _cases(candidate)
    results = []
    for key in sorted(before):
        prior = before[key]
        current = after.get(key)
        reasons = []
        deltas = {}
        if current is None:
            if not allow_missing:
                reasons.append("case is missing from the candidate report")
            candidate_status = "missing"
        else:
            candidate_status = current.get("status")
            if prior.get("status") in {"pass", "skipped"} and candidate_status == "fail":
                reasons.append(f"status changed from {prior.get('status')} to fail")
            old_values = _measurements(prior)
            new_values = _measurements(current)
            for name in sorted(set(old_values) & set(new_values)):
                delta = new_values[name] - old_values[name]
                deltas[name] = delta
                if delta < -absolute_tolerance:
                    reasons.append(f"{name} decreased by {abs(delta):.6g}")
        results.append({
            "id": key,
            "status": "fail" if reasons else "pass",
            "baseline_status": prior.get("status"),
            "candidate_status": candidate_status,
            "metric_deltas": deltas,
            "reason": "; ".join(reasons),
        })
    added = sorted(set(after) - set(before))
    return {
        "operation": "regress",
        "baseline": str(Path(baseline_path)),
        "candidate": str(Path(candidate_path)),
        "source_operation": baseline["operation"],
        "absolute_tolerance": absolute_tolerance,
        "allow_missing": allow_missing,
        "summary": {**_summary(results), "added": len(added), "missing": len(set(before) - set(after))},
        "added_cases": added,
        "results": results,
    }
