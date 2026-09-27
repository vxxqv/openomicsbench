"""Collection-wide validation and comparison reports for CI systems."""
from __future__ import annotations

import json
import os
import csv
import html
import io
import xml.etree.ElementTree as ET
from pathlib import Path

from .compare import compare
from .registry import registry
from .validate import validate


RESULT_SUFFIXES = (".csv", ".tsv", ".csv.gz", ".tsv.gz")


def _selected(root: Path, assay: str | None, dataset_ids: list[str] | None):
    records = registry(root)
    requested = list(dict.fromkeys(dataset_ids or []))
    unknown = [dataset_id for dataset_id in requested if dataset_id not in records]
    if unknown:
        raise ValueError(f"unknown dataset IDs: {', '.join(unknown)}")
    selected = [(model, folder) for model, folder in records.values() if not assay or model.assay == assay]
    if requested:
        wanted = set(requested)
        selected = [(model, folder) for model, folder in selected if model.id in wanted]
        excluded = [dataset_id for dataset_id in requested if dataset_id not in {model.id for model, _ in selected}]
        if excluded:
            raise ValueError(f"dataset IDs do not match assay {assay}: {', '.join(excluded)}")
    if not selected:
        raise ValueError("no datasets matched the suite selection")
    return selected


def _summary(results: list[dict]) -> dict:
    counts = {
        "passed": sum(item["status"] == "pass" for item in results),
        "failed": sum(item["status"] == "fail" for item in results),
        "skipped": sum(item["status"] == "skipped" for item in results),
    }
    return {"status": "fail" if counts["failed"] else "pass", "total": len(results), **counts}


def validate_suite(root: Path, assay: str | None = None, dataset_ids: list[str] | None = None) -> dict:
    results = []
    for model, folder in _selected(Path(root), assay, dataset_ids):
        try:
            details = validate(model, folder)
            if details["status"] == "link_only":
                results.append({"id": model.id, "status": "skipped", "reason": "link-only object", "details": details})
            else:
                results.append({"id": model.id, "status": "pass", "details": details})
        except (ValueError, OSError, KeyError) as error:
            results.append({"id": model.id, "status": "fail", "reason": str(error)})
    return {"operation": "validate", "selection": {"assay": assay, "ids": dataset_ids or []}, "summary": _summary(results), "results": results}


def _result_file(directory: Path, dataset_id: str) -> Path | None:
    matches = [directory / f"{dataset_id}{suffix}" for suffix in RESULT_SUFFIXES if (directory / f"{dataset_id}{suffix}").is_file()]
    if len(matches) > 1:
        raise ValueError(f"{dataset_id}: multiple result files found: {', '.join(path.name for path in matches)}")
    return matches[0] if matches else None


def compare_suite(root: Path, results_directory: Path, dataset_ids: list[str] | None = None, detail_limit: int = 20) -> dict:
    directory = Path(results_directory)
    if not directory.is_dir():
        raise ValueError(f"{directory}: results directory does not exist")
    selected = _selected(Path(root), "bulk_rna_seq", dataset_ids)
    comparable = [(model, folder) for model, folder in selected if model.kind == "real" and model.validation is not None]
    if not comparable:
        raise ValueError("no selected datasets have a differential-expression comparison reference")
    results = []
    for model, folder in comparable:
        try:
            path = _result_file(directory, model.id)
            if path is None:
                results.append({"id": model.id, "status": "fail", "reason": f"missing {model.id}.csv or {model.id}.tsv result file"})
                continue
            details = compare(model, folder, path, detail_limit)
            result = {"id": model.id, "status": details["status"], "input": str(path), "details": details}
            if details["status"] == "fail":
                result["reason"] = "; ".join(details["reasons"])
            results.append(result)
        except (ValueError, OSError, KeyError) as error:
            results.append({"id": model.id, "status": "fail", "reason": str(error)})
    return {
        "operation": "compare",
        "selection": {"assay": "bulk_rna_seq", "ids": dataset_ids or []},
        "results_directory": str(directory),
        "summary": _summary(results),
        "results": results,
    }


def markdown_report(report: dict) -> str:
    summary = report["summary"]
    lines = [
        "# OpenOmicsBench suite report",
        "",
        f"Operation: `{report['operation']}`",
        "",
        f"Status: **{summary['status'].upper()}**",
        "",
        f"Passed: {summary['passed']}  ",
        f"Failed: {summary['failed']}  ",
        f"Skipped: {summary['skipped']}",
        "",
        "| Benchmark | Status | Detail |",
        "|---|---|---|",
    ]
    for item in report["results"]:
        detail = item.get("reason", "")
        if not detail and report["operation"] == "compare":
            metrics = item["details"]["metrics"]
            detail = f"Spearman {metrics['spearman_logfc']:.3f}; top-k Jaccard {metrics['top_k_jaccard']:.3f}; sign agreement {metrics['sign_concordance']:.3f}"
        detail = detail.replace("|", "\\|").replace("\n", " ")
        lines.append(f"| `{item['id']}` | {item['status']} | {detail} |")
    return "\n".join(lines) + "\n"


def junit_report(report: dict) -> str:
    summary = report["summary"]
    suite = ET.Element(
        "testsuite",
        name=f"openomicsbench.{report['operation']}",
        tests=str(summary["total"]),
        failures=str(summary["failed"]),
        skipped=str(summary["skipped"]),
    )
    for item in report["results"]:
        case = ET.SubElement(suite, "testcase", classname="openomicsbench", name=item["id"])
        if item["status"] == "fail":
            failure = ET.SubElement(case, "failure", message=item.get("reason", "benchmark failed"))
            failure.text = item.get("reason", "benchmark failed")
        elif item["status"] == "skipped":
            ET.SubElement(case, "skipped", message=item.get("reason", "skipped"))
        output = ET.SubElement(case, "system-out")
        output.text = json.dumps(item, sort_keys=True, allow_nan=False)
    ET.indent(suite, space="  ")
    return ET.tostring(suite, encoding="unicode", xml_declaration=True) + "\n"


def csv_report(report: dict) -> str:
    fields = (
        "method", "benchmark", "status", "reason", "baseline_status", "candidate_status",
        "spearman_logfc", "top_k_jaccard", "sign_concordance", "coverage",
        "delta_spearman_logfc", "delta_top_k_jaccard", "delta_sign_concordance", "delta_coverage",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for item in report["results"]:
        details = item.get("details", {})
        metrics = details.get("metrics", {})
        genes = details.get("genes", {})
        deltas = item.get("metric_deltas", {})
        writer.writerow({
            "method": item.get("method", ""),
            "benchmark": item["id"],
            "status": item["status"],
            "reason": item.get("reason", ""),
            "baseline_status": item.get("baseline_status", ""),
            "candidate_status": item.get("candidate_status", ""),
            "spearman_logfc": metrics.get("spearman_logfc", ""),
            "top_k_jaccard": metrics.get("top_k_jaccard", ""),
            "sign_concordance": metrics.get("sign_concordance", ""),
            "coverage": genes.get("coverage", ""),
            "delta_spearman_logfc": deltas.get("spearman_logfc", ""),
            "delta_top_k_jaccard": deltas.get("top_k_jaccard", ""),
            "delta_sign_concordance": deltas.get("sign_concordance", ""),
            "delta_coverage": deltas.get("coverage", ""),
        })
    return stream.getvalue()


def html_report(report: dict) -> str:
    summary = report["summary"]
    rows = []
    for item in report["results"]:
        details = item.get("details", {})
        metrics = details.get("metrics", {})
        values = []
        for name in ("spearman_logfc", "top_k_jaccard", "sign_concordance"):
            value = metrics.get(name)
            if value is not None:
                values.append(f"{name.replace('_', ' ')} {value:.3f}")
        detail = item.get("reason") or "; ".join(values)
        benchmark = f"{item.get('method')}/{item['id']}" if item.get("method") else item["id"]
        rows.append(
            f"<tr><td><code>{html.escape(benchmark)}</code></td>"
            f"<td class=\"{html.escape(item['status'])}\">{html.escape(item['status'])}</td>"
            f"<td>{html.escape(detail or '')}</td></tr>"
        )
    leaderboard = ""
    if report.get("leaderboard"):
        leaders = "".join(
            f"<tr><td>{index}</td><td><code>{html.escape(row['name'])}</code></td>"
            f"<td>{row['passed']}/{row['total']}</td><td>{row['pass_rate']:.1%}</td></tr>"
            for index, row in enumerate(report["leaderboard"], start=1)
        )
        leaderboard = f"<h2>Methods</h2><table><thead><tr><th>Rank</th><th>Method</th><th>Passed</th><th>Pass rate</th></tr></thead><tbody>{leaders}</tbody></table>"
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>OpenOmicsBench report</title><style>
body{font:16px system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#111;background:#fff}
table{border-collapse:collapse;width:100%;margin:1rem 0 2rem}th,td{border:1px solid #bbb;padding:.55rem;text-align:left;vertical-align:top}
th{background:#eee}.pass{color:#126b2e;font-weight:700}.fail{color:#a01414;font-weight:700}.skipped{color:#555;font-weight:700}
code{font-family:ui-monospace,monospace}dl{display:grid;grid-template-columns:max-content auto;gap:.35rem 1rem}dt{font-weight:700}dd{margin:0}
</style></head><body>""" + (
        f"<h1>OpenOmicsBench suite report</h1><dl><dt>Operation</dt><dd><code>{html.escape(report['operation'])}</code></dd>"
        f"<dt>Status</dt><dd class=\"{html.escape(summary['status'])}\">{html.escape(summary['status'].upper())}</dd>"
        f"<dt>Passed</dt><dd>{summary['passed']}</dd><dt>Failed</dt><dd>{summary['failed']}</dd>"
        f"<dt>Skipped</dt><dd>{summary['skipped']}</dd></dl>{leaderboard}"
        f"<h2>Cases</h2><table><thead><tr><th>Benchmark</th><th>Status</th><th>Detail</th></tr></thead><tbody>{''.join(rows)}</tbody></table>"
        "</body></html>\n"
    )


def _atomic_text(path: Path, text: str, force: bool) -> None:
    path = Path(path)
    if path.exists() and not force:
        raise ValueError(f"{path}: output already exists; use --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_reports(
    report: dict,
    json_path: Path | None = None,
    markdown_path: Path | None = None,
    junit_path: Path | None = None,
    force: bool = False,
    csv_path: Path | None = None,
    html_path: Path | None = None,
) -> dict:
    destinations = [Path(path) for path in (json_path, markdown_path, junit_path, csv_path, html_path) if path is not None]
    duplicates = {path for path in destinations if destinations.count(path) > 1}
    if duplicates:
        raise ValueError(f"report paths must be distinct: {', '.join(str(path) for path in sorted(duplicates))}")
    existing = [path for path in destinations if path.exists()]
    if existing and not force:
        raise ValueError(f"{existing[0]}: output already exists; use --force to replace it")
    if json_path is not None:
        _atomic_text(Path(json_path), json.dumps(report, indent=2, allow_nan=False) + "\n", force)
    if markdown_path is not None:
        _atomic_text(Path(markdown_path), markdown_report(report), force)
    if junit_path is not None:
        _atomic_text(Path(junit_path), junit_report(report), force)
    if csv_path is not None:
        _atomic_text(Path(csv_path), csv_report(report), force)
    if html_path is not None:
        _atomic_text(Path(html_path), html_report(report), force)
    return {
        "json": str(Path(json_path)) if json_path is not None else None,
        "markdown": str(Path(markdown_path)) if markdown_path is not None else None,
        "junit": str(Path(junit_path)) if junit_path is not None else None,
        "csv": str(Path(csv_path)) if csv_path is not None else None,
        "html": str(Path(html_path)) if html_path is not None else None,
    }
