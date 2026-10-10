"""Cross-run histories for versioned evaluation reports."""
from __future__ import annotations

from pathlib import Path

from .evaluate import _content_fingerprint
from .matrix import METHOD_NAME
from .regression import _measurements, read_report
from .suite import _summary


def parse_run(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError("run must use NAME=PATH")
    name, raw_path = value.split("=", 1)
    if not METHOD_NAME.fullmatch(name):
        raise ValueError(f"{name or '(empty)'}: invalid run name")
    if not raw_path:
        raise ValueError(f"{name}: report path is empty")
    return name, Path(raw_path)


def _case_map(report: dict) -> dict[str, dict]:
    cases = {}
    for item in report["results"]:
        identifier = item.get("id") if isinstance(item, dict) else None
        if not isinstance(identifier, str) or not identifier:
            raise ValueError("evaluation result is missing a benchmark ID")
        if identifier in cases:
            raise ValueError(f"duplicate evaluation result: {identifier}")
        cases[identifier] = item
    return cases


def trend_reports(runs: list[tuple[str, Path]]) -> dict:
    if len(runs) < 2:
        raise ValueError("a trend report requires at least two --run entries")
    names = [name for name, _ in runs]
    if len(names) != len(set(names)):
        raise ValueError("run names must be unique")
    loaded = []
    run_summaries = []
    for name, path in runs:
        report = read_report(path)
        if report["operation"] != "evaluate":
            raise ValueError(f"{path}: trend inputs must be suite evaluate reports")
        fingerprint = report.get("run_fingerprint")
        if not isinstance(fingerprint, str) or not fingerprint.startswith("sha256:") or len(fingerprint) != 71:
            raise ValueError(f"{path}: evaluation report has no valid run fingerprint")
        cases = _case_map(report)
        summary = report["summary"]
        total = summary.get("total", 0)
        run_summaries.append({
            "name": name,
            "path": str(Path(path)),
            "run_fingerprint": fingerprint,
            "status": summary.get("status"),
            "passed": summary.get("passed"),
            "failed": summary.get("failed"),
            "total": total,
            "pass_rate": summary.get("passed", 0) / total if total else 0.0,
            "method": report.get("method"),
        })
        loaded.append((name, cases))

    identifiers = sorted(set().union(*(set(cases) for _, cases in loaded)))
    results = []
    for identifier in identifiers:
        history = []
        for name, cases in loaded:
            item = cases.get(identifier)
            history.append({
                "run": name,
                "status": item.get("status") if item else "missing",
                "assay": item.get("assay") if item else None,
                "comparator": item.get("comparator") if item else None,
                "metrics": _measurements(item) if item else {},
                "input_sha256": item.get("input_sha256") if item else None,
                "benchmark_manifest_sha256": item.get("benchmark_manifest_sha256") if item else None,
            })
        first, latest = history[0], history[-1]
        common = sorted(set(first["metrics"]) & set(latest["metrics"]))
        deltas = {name: latest["metrics"][name] - first["metrics"][name] for name in common}
        latest_case = loaded[-1][1].get(identifier)
        status = latest_case.get("status") if latest_case else "fail"
        reason = "benchmark is missing from the latest run" if latest_case is None else latest_case.get("reason", "")
        results.append({
            "id": identifier,
            "assay": latest.get("assay") or first.get("assay") or "",
            "comparator": latest.get("comparator") or first.get("comparator") or "",
            "status": status,
            "reason": reason,
            "details": {"metrics": latest["metrics"]},
            "metric_deltas": deltas,
            "history": history,
        })

    summary = _summary(results)
    summary.update({"runs": len(runs), "benchmarks": len(results), "first_run": names[0], "latest_run": names[-1]})
    fingerprint = _content_fingerprint({
        "contract": "openomicsbench-trend-v1",
        "runs": [{"name": item["name"], "run_fingerprint": item["run_fingerprint"]} for item in run_summaries],
    })
    return {
        "operation": "trend",
        "run_fingerprint": fingerprint,
        "fingerprint_contract": "openomicsbench-trend-v1",
        "summary": summary,
        "runs": run_summaries,
        "results": results,
    }
