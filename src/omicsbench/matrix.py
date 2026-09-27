"""Compare several analysis methods across the RNA-seq benchmark collection."""
from __future__ import annotations

import re
from pathlib import Path

from .suite import _summary, compare_suite


METHOD_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
SCORED_METRICS = ("spearman_logfc", "top_k_jaccard", "sign_concordance")


def parse_method(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"{value}: method must use NAME=RESULTS_DIRECTORY")
    name, directory = value.split("=", 1)
    if not METHOD_NAME.fullmatch(name):
        raise ValueError(f"{name}: method name must use letters, numbers, dots, underscores or hyphens")
    if not directory:
        raise ValueError(f"{value}: results directory is empty")
    return name, Path(directory)


def _method_summary(name: str, directory: Path, results: list[dict]) -> dict:
    summary = _summary(results)
    completed = [item for item in results if "details" in item]
    mean_metrics = {
        metric: (sum(item["details"]["metrics"][metric] for item in completed) / len(completed) if completed else None)
        for metric in SCORED_METRICS
    }
    mean_coverage = (
        sum(item["details"]["genes"]["coverage"] for item in completed) / len(completed)
        if completed else None
    )
    return {
        "name": name,
        "results_directory": str(directory),
        **summary,
        "pass_rate": summary["passed"] / summary["total"] if summary["total"] else 0.0,
        "mean_metrics": mean_metrics,
        "mean_coverage": mean_coverage,
    }


def compare_matrix(
    root: Path,
    methods: list[tuple[str, Path]],
    dataset_ids: list[str] | None = None,
    detail_limit: int = 20,
) -> dict:
    if len(methods) < 2:
        raise ValueError("a benchmark matrix requires at least two --method entries")
    names = [name for name, _ in methods]
    if len(names) != len(set(names)):
        raise ValueError("method names must be unique")
    cells = []
    method_summaries = []
    benchmark_ids = None
    for name, directory in methods:
        if not METHOD_NAME.fullmatch(name):
            raise ValueError(f"{name}: invalid method name")
        report = compare_suite(root, directory, dataset_ids, detail_limit)
        ids = [item["id"] for item in report["results"]]
        if benchmark_ids is None:
            benchmark_ids = ids
        elif ids != benchmark_ids:
            raise ValueError(f"{name}: benchmark selection differs from the first method")
        method_results = []
        for item in report["results"]:
            cell = {"method": name, **item}
            cells.append(cell)
            method_results.append(item)
        method_summaries.append(_method_summary(name, directory, method_results))

    def rank_key(item: dict):
        spearman = item["mean_metrics"]["spearman_logfc"]
        return (-item["passed"], item["failed"], -(spearman if spearman is not None else -1.0), item["name"])

    leaderboard = sorted(method_summaries, key=rank_key)
    return {
        "operation": "matrix",
        "selection": {"assay": "bulk_rna_seq", "ids": dataset_ids or []},
        "summary": {**_summary(cells), "methods": len(methods), "benchmarks": len(benchmark_ids or [])},
        "leaderboard": leaderboard,
        "results": cells,
    }
