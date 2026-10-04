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
        "| Benchmark | Assay | Status | Detail |",
        "|---|---|---|---|",
    ]
    for item in report["results"]:
        detail = _case_detail(item)
        detail = detail.replace("|", "\\|").replace("\n", " ")
        benchmark = f"{item.get('method')}/{item['id']}" if item.get("method") else item["id"]
        lines.append(f"| `{benchmark}` | {item.get('assay', '')} | {item['status']} | {detail} |")
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
        "method", "benchmark", "assay", "comparator", "status", "reason", "baseline_status", "candidate_status",
        "spearman_logfc", "top_k_jaccard", "sign_concordance", "coverage", "precision", "recall", "f1",
        "true_positive", "false_positive", "false_negative",
        "delta_spearman_logfc", "delta_top_k_jaccard", "delta_sign_concordance", "delta_coverage",
        "delta_precision", "delta_recall", "delta_f1",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for item in report["results"]:
        details = item.get("details", {})
        metrics = details.get("metrics", {})
        genes = details.get("genes", {})
        counts = details.get("counts", {})
        deltas = item.get("metric_deltas", {})
        writer.writerow({
            "method": item.get("method", ""),
            "benchmark": item["id"],
            "assay": item.get("assay", ""),
            "comparator": item.get("comparator", ""),
            "status": item["status"],
            "reason": item.get("reason", ""),
            "baseline_status": item.get("baseline_status", ""),
            "candidate_status": item.get("candidate_status", ""),
            "spearman_logfc": metrics.get("spearman_logfc", ""),
            "top_k_jaccard": metrics.get("top_k_jaccard", ""),
            "sign_concordance": metrics.get("sign_concordance", ""),
            "coverage": genes.get("coverage", ""),
            "precision": metrics.get("precision", ""),
            "recall": metrics.get("recall", ""),
            "f1": metrics.get("f1", ""),
            "true_positive": counts.get("true_positive", ""),
            "false_positive": counts.get("false_positive", ""),
            "false_negative": counts.get("false_negative", ""),
            "delta_spearman_logfc": deltas.get("spearman_logfc", ""),
            "delta_top_k_jaccard": deltas.get("top_k_jaccard", ""),
            "delta_sign_concordance": deltas.get("sign_concordance", ""),
            "delta_coverage": deltas.get("coverage", ""),
            "delta_precision": deltas.get("precision", ""),
            "delta_recall": deltas.get("recall", ""),
            "delta_f1": deltas.get("f1", ""),
        })
    return stream.getvalue()


METRIC_LABELS = {
    "spearman_logfc": "Spearman",
    "top_k_jaccard": "top-k Jaccard",
    "sign_concordance": "sign agreement",
    "precision": "precision",
    "recall": "recall",
    "f1": "F1",
}


def _case_metrics(item: dict) -> list[tuple[str, float]]:
    metrics = item.get("details", {}).get("metrics", {})
    return [(METRIC_LABELS[name], metrics[name]) for name in METRIC_LABELS if metrics.get(name) is not None]


def _case_detail(item: dict) -> str:
    if item.get("reason"):
        return item["reason"]
    values = [f"{label} {value:.3f}" for label, value in _case_metrics(item)]
    coverage = item.get("details", {}).get("genes", {}).get("coverage")
    if coverage is not None:
        values.append(f"coverage {coverage:.3f}")
    return "; ".join(values)


def _status_chart(summary: dict) -> str:
    total = summary.get("total", 0) or 1
    passed = 100 * summary.get("passed", 0) / total
    failed = 100 * summary.get("failed", 0) / total
    return (
        '<svg class="status-chart" viewBox="0 0 100 12" role="img" aria-label="case status distribution">'
        f'<rect x="0" y="0" width="{passed:.6f}" height="12" class="bar-pass"/>'
        f'<rect x="{passed:.6f}" y="0" width="{failed:.6f}" height="12" class="bar-fail"/>'
        f'<rect x="{passed + failed:.6f}" y="0" width="{100 - passed - failed:.6f}" height="12" class="bar-skip"/>'
        "</svg>"
    )


def _leaderboard(report: dict) -> str:
    if not report.get("leaderboard"):
        return ""
    rows = []
    for index, row in enumerate(report["leaderboard"], start=1):
        receipt = row.get("receipt") or {}
        runtime = receipt.get("runtime_seconds")
        memory = receipt.get("peak_memory_mb")
        assays = row.get("assays", {})
        assay_text = ", ".join(
            f"{name}: {values['passed']}/{values['total']}" for name, values in assays.items()
        )
        rows.append(
            f"<tr><td>{index}</td><td><code>{html.escape(row['name'])}</code></td>"
            f"<td>{row['passed']}/{row['total']}</td><td>{row['pass_rate']:.1%}</td>"
            f"<td>{html.escape(assay_text)}</td>"
            f"<td>{html.escape(f'{runtime:.3f} s' if runtime is not None else '')}</td>"
            f"<td>{html.escape(f'{memory:.3f} MB' if memory is not None else '')}</td></tr>"
        )
    ranking = html.escape(report.get("ranking", ""))
    note = f"<p class=\"note\">{ranking}</p>" if ranking else ""
    return (
        "<section><h2>Method leaderboard</h2>" + note
        + "<div class=\"table-wrap\"><table><thead><tr><th>Rank</th><th>Method</th><th>Passed</th><th>Pass rate</th>"
        + "<th>Assays</th><th>Runtime</th><th>Peak memory</th></tr></thead><tbody>"
        + "".join(rows) + "</tbody></table></div></section>"
    )


def _matrix(report: dict) -> str:
    if not report.get("leaderboard") or not any(item.get("method") for item in report["results"]):
        return ""
    methods = [row["name"] for row in report["leaderboard"]]
    benchmarks = list(dict.fromkeys(item["id"] for item in report["results"]))
    cases = {(item["method"], item["id"]): item for item in report["results"]}
    rows = []
    for method in methods:
        cells = []
        for benchmark in benchmarks:
            item = cases.get((method, benchmark))
            status = item.get("status", "missing") if item else "missing"
            detail = _case_detail(item) if item else "case absent"
            cells.append(f'<td class="cell-{html.escape(status)}" title="{html.escape(detail)}">{html.escape(status)}</td>')
        rows.append(f"<tr><th><code>{html.escape(method)}</code></th>{''.join(cells)}</tr>")
    headers = "".join(f"<th><code>{html.escape(benchmark)}</code></th>" for benchmark in benchmarks)
    return (
        "<section><h2>Coverage matrix</h2><div class=\"table-wrap\"><table class=\"matrix\"><thead><tr><th>Method</th>"
        + headers + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div></section>"
    )


def html_report(report: dict) -> str:
    summary = report["summary"]
    rows = []
    for item in report["results"]:
        detail = _case_detail(item)
        benchmark = f"{item.get('method')}/{item['id']}" if item.get("method") else item["id"]
        metrics = "".join(f'<span class="metric">{html.escape(label)} <strong>{value:.3f}</strong></span>' for label, value in _case_metrics(item))
        rows.append(
            f"<tr><td><code>{html.escape(benchmark)}</code></td>"
            f"<td>{html.escape(item.get('assay', ''))}</td>"
            f"<td><span class=\"status {html.escape(item['status'])}\">{html.escape(item['status'])}</span></td>"
            f"<td>{metrics or html.escape(detail or '')}</td></tr>"
        )
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>OpenOmicsBench report</title><style>
*{box-sizing:border-box}body{font:15px/1.5 system-ui,sans-serif;max-width:1240px;margin:0 auto;padding:2rem;color:#171717;background:#f5f5f3}
h1{font-size:2rem;margin:.2rem 0}h2{margin:2rem 0 .5rem}header,section{background:#fff;border:1px solid #d3d3cf;border-radius:10px;padding:1.25rem;margin-bottom:1rem}
.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-size:.75rem;font-weight:700;color:#555}.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:.75rem;margin:1.25rem 0}
.card{border:1px solid #d8d8d3;border-radius:8px;padding:.8rem}.card strong{display:block;font-size:1.65rem}.status-chart{width:100%;height:12px;border-radius:6px;overflow:hidden;background:#ddd}
.bar-pass{fill:#267a45}.bar-fail{fill:#b33a3a}.bar-skip{fill:#888}table{border-collapse:collapse;width:100%}th,td{border-bottom:1px solid #ddd;padding:.65rem;text-align:left;vertical-align:top}th{background:#f0f0ed;font-size:.8rem;text-transform:uppercase;letter-spacing:.04em}
.table-wrap{overflow:auto}.status{display:inline-block;border-radius:999px;padding:.15rem .55rem;font-weight:700}.pass,.cell-pass{color:#145b30;background:#e5f4e9}.fail,.cell-fail{color:#8b1f1f;background:#fae6e6}.skipped,.cell-skipped,.cell-missing{color:#444;background:#e9e9e6}
.metric{display:inline-block;margin:0 .35rem .25rem 0;padding:.2rem .5rem;border:1px solid #ccc;border-radius:4px;white-space:nowrap}.matrix td{text-align:center;font-weight:700}.note{color:#555;max-width:80ch}code{font-family:ui-monospace,monospace}
@media(max-width:720px){body{padding:.75rem}.cards{grid-template-columns:repeat(2,1fr)}header,section{padding:.85rem}}
@media(prefers-color-scheme:dark){body{color:#eee;background:#151515}header,section{background:#1e1e1e;border-color:#444}.card{border-color:#555}th{background:#2b2b2b}th,td{border-color:#444}.note,.eyebrow{color:#bbb}.metric{border-color:#666}.pass,.cell-pass{color:#b9efc9;background:#173d26}.fail,.cell-fail{color:#ffc4c4;background:#4a1d1d}.skipped,.cell-skipped,.cell-missing{color:#ddd;background:#333}}
</style></head><body>""" + (
        f"<header><div class=\"eyebrow\">OpenOmicsBench</div><h1>Benchmark report</h1><p>Operation <code>{html.escape(report['operation'])}</code></p>"
        f"<div class=\"cards\"><div class=\"card\"><span>Status</span><strong class=\"{html.escape(summary['status'])}\">{html.escape(summary['status'].upper())}</strong></div>"
        f"<div class=\"card\"><span>Passed</span><strong>{summary['passed']}</strong></div><div class=\"card\"><span>Failed</span><strong>{summary['failed']}</strong></div>"
        f"<div class=\"card\"><span>Total</span><strong>{summary['total']}</strong></div></div>{_status_chart(summary)}</header>"
        f"{_leaderboard(report)}{_matrix(report)}"
        f"<section><h2>Benchmark cases</h2><div class=\"table-wrap\"><table><thead><tr><th>Benchmark</th><th>Assay</th><th>Status</th><th>Metrics or reason</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></section>"
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
