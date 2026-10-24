"""Policy gates for scientific results, benchmark coverage and resources."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from .models import StrictModel
from .regression import read_report
from .suite import _summary


class MetricFloor(StrictModel):
    metric: Literal[
        "spearman_logfc", "top_k_jaccard", "sign_concordance", "coverage",
        "precision", "recall", "f1", "protein_recall", "residue_coverage",
        "adjusted_rand_index", "normalized_mutual_information", "purity",
        "interval_precision", "interval_recall", "interval_f1",
        "base_precision", "base_recall", "base_f1",
    ]
    minimum: float = Field(ge=0, le=1)
    assay: str | None = Field(default=None, min_length=1)
    comparator: str | None = Field(default=None, min_length=1)


class GatePolicy(StrictModel):
    schema_version: str = Field(pattern=r"^1\.0$")
    name: str = Field(min_length=1, max_length=128)
    required_benchmarks: list[str] = Field(default_factory=list)
    minimum_pass_rate: float = Field(default=1.0, ge=0, le=1)
    minimum_assay_pass_rate: dict[str, float] = Field(default_factory=dict)
    require_method_receipt: bool = False
    require_run_fingerprint: bool = True
    maximum_runtime_seconds: float | None = Field(default=None, ge=0)
    maximum_peak_memory_mb: float | None = Field(default=None, ge=0)
    metric_floors: list[MetricFloor] = Field(default_factory=list)

    @field_validator("required_benchmarks")
    @classmethod
    def unique_benchmarks(cls, value):
        if any(not item or item.strip() != item or any(character.isspace() for character in item) for item in value):
            raise ValueError("required benchmark IDs must be non-empty tokens")
        if len(value) != len(set(value)):
            raise ValueError("required benchmark IDs must be unique")
        return value

    @field_validator("minimum_assay_pass_rate")
    @classmethod
    def assay_rates(cls, value):
        for assay, rate in value.items():
            if not assay or assay.strip() != assay or not math.isfinite(rate) or not 0 <= rate <= 1:
                raise ValueError("assay pass rates require non-empty assay names and values between zero and one")
        return value

    @model_validator(mode="after")
    def unique_metric_scopes(self):
        scopes = [(item.metric, item.assay, item.comparator) for item in self.metric_floors]
        if len(scopes) != len(set(scopes)):
            raise ValueError("metric floor scopes must be unique")
        return self


def read_policy(path: Path) -> GatePolicy:
    try:
        return GatePolicy.model_validate_json(Path(path).read_text(encoding="utf-8"))
    except (ValueError, UnicodeDecodeError) as error:
        raise ValueError(f"{path}: invalid gate policy: {error}") from None


def _check(identifier: str, passed: bool, reason: str, details: dict | None = None) -> dict:
    item = {"id": identifier, "status": "pass" if passed else "fail", "reason": "" if passed else reason}
    if details:
        item["details"] = details
    return item


def _metric(item: dict, name: str):
    details = item.get("details", {})
    if name == "coverage":
        return details.get("metrics", {}).get("coverage", details.get("genes", {}).get("coverage"))
    return details.get("metrics", {}).get(name)


def evaluate_policy(report_path: Path, policy_path: Path) -> dict:
    report = read_report(Path(report_path))
    if report["operation"] != "evaluate":
        raise ValueError("gate policies require a suite evaluate report")
    policy = read_policy(Path(policy_path))
    cases = {item.get("id"): item for item in report["results"] if isinstance(item, dict)}
    if None in cases or len(cases) != len(report["results"]):
        raise ValueError("evaluation report contains missing or duplicate benchmark IDs")
    checks = []
    fingerprint = report.get("run_fingerprint")
    valid_fingerprint = isinstance(fingerprint, str) and len(fingerprint) == 71 and fingerprint.startswith("sha256:")
    if policy.require_run_fingerprint:
        checks.append(_check("run-fingerprint", valid_fingerprint, "evaluation report has no valid SHA-256 run fingerprint"))

    summary = report.get("summary", {})
    total, passed = summary.get("total"), summary.get("passed")
    if not isinstance(total, int) or not isinstance(passed, int) or total < 1 or not 0 <= passed <= total:
        raise ValueError("evaluation report has an invalid summary")
    pass_rate = passed / total
    checks.append(_check(
        "overall-pass-rate",
        pass_rate >= policy.minimum_pass_rate,
        f"pass rate {pass_rate:.6g} is below {policy.minimum_pass_rate:.6g}",
        {"observed": pass_rate, "minimum": policy.minimum_pass_rate},
    ))

    for benchmark in policy.required_benchmarks:
        case = cases.get(benchmark)
        ok = case is not None and case.get("status") == "pass"
        observed = "missing" if case is None else case.get("status")
        checks.append(_check(f"benchmark:{benchmark}", ok, f"required benchmark is {observed}", {"observed": observed}))

    assays = summary.get("assays", {})
    for assay, minimum in sorted(policy.minimum_assay_pass_rate.items()):
        values = assays.get(assay)
        if not isinstance(values, dict) or not values.get("total"):
            checks.append(_check(f"assay:{assay}", False, "assay is absent from the evaluation"))
            continue
        observed = values.get("passed", 0) / values["total"]
        checks.append(_check(
            f"assay:{assay}", observed >= minimum, f"assay pass rate {observed:.6g} is below {minimum:.6g}",
            {"observed": observed, "minimum": minimum},
        ))

    receipt = report.get("method")
    if policy.require_method_receipt:
        checks.append(_check("method-receipt", isinstance(receipt, dict), "method receipt is required"))
    for key, maximum, label in (
        ("runtime_seconds", policy.maximum_runtime_seconds, "runtime-seconds"),
        ("peak_memory_mb", policy.maximum_peak_memory_mb, "peak-memory-mb"),
    ):
        if maximum is None:
            continue
        observed = receipt.get(key) if isinstance(receipt, dict) else None
        ok = isinstance(observed, (int, float)) and not isinstance(observed, bool) and math.isfinite(observed) and observed <= maximum
        reason = f"{key} is missing" if observed is None else f"{key} {observed:.6g} exceeds {maximum:.6g}"
        checks.append(_check(label, ok, reason, {"observed": observed, "maximum": maximum}))

    for requirement in policy.metric_floors:
        matches = [
            item for item in report["results"]
            if (requirement.assay is None or item.get("assay") == requirement.assay)
            and (requirement.comparator is None or item.get("comparator") == requirement.comparator)
        ]
        scope = ":".join(value or "all" for value in (requirement.metric, requirement.assay, requirement.comparator))
        measured = [(item["id"], _metric(item, requirement.metric)) for item in matches]
        measured = [(identifier, value) for identifier, value in measured if isinstance(value, (int, float)) and not isinstance(value, bool)]
        if not measured:
            checks.append(_check(f"metric:{scope}", False, "no matching metric values were found"))
            continue
        failing = [(identifier, float(value)) for identifier, value in measured if not math.isfinite(value) or value < requirement.minimum]
        checks.append(_check(
            f"metric:{scope}", not failing,
            "; ".join(f"{identifier}={value:.6g}" for identifier, value in failing) + f" below {requirement.minimum:.6g}",
            {"minimum": requirement.minimum, "cases": len(measured), "failed": len(failing)},
        ))

    gate_summary = _summary(checks)
    gate_summary.update({"policy": policy.name, "source_cases": total, "source_pass_rate": pass_rate})
    return {
        "operation": "gate",
        "policy_path": str(Path(policy_path)),
        "policy": policy.model_dump(mode="json"),
        "source_report": str(Path(report_path)),
        "source_run_fingerprint": fingerprint,
        "summary": gate_summary,
        "results": checks,
    }
