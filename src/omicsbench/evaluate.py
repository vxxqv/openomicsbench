"""Unified evaluation of RNA effects and exact synthetic SNV calls."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

from pydantic import Field, field_validator

from .compare import compare
from .matrix import METHOD_NAME
from .models import StrictModel
from .registry import registry
from .suite import RESULT_SUFFIXES, _summary
from .variants import compare_variants


VCF_SUFFIXES = (".vcf", ".vcf.gz")


class MethodReceipt(StrictModel):
    schema_version: str = Field(pattern=r"^1\.0$")
    name: str = Field(min_length=1, max_length=128)
    version: str = Field(min_length=1, max_length=128)
    description: str | None = Field(default=None, max_length=500)
    command: str | None = Field(default=None, max_length=4000)
    container: str | None = Field(default=None, max_length=500)
    source_revision: str | None = Field(default=None, max_length=128)
    runtime_seconds: float | None = Field(default=None, ge=0)
    peak_memory_mb: float | None = Field(default=None, ge=0)
    threads: int | None = Field(default=None, gt=0, strict=True)
    parameters: dict[str, str | int | float | bool | None] = Field(default_factory=dict)

    @field_validator("runtime_seconds", "peak_memory_mb")
    @classmethod
    def finite_resource_value(cls, value):
        if value is not None and not math.isfinite(value):
            raise ValueError("resource measurements must be finite")
        return value


def read_method_receipt(directory: Path) -> dict | None:
    path = Path(directory) / "method.json"
    if not path.exists():
        return None
    if not path.is_file():
        raise ValueError(f"{path}: method receipt is not a file")
    try:
        receipt = MethodReceipt.model_validate_json(path.read_text(encoding="utf-8"))
    except (ValueError, UnicodeDecodeError) as error:
        raise ValueError(f"{path}: invalid method receipt: {error}") from None
    return {"path": str(path), **receipt.model_dump(mode="json")}


def _candidate_files(directory: Path, dataset_id: str, suffixes: tuple[str, ...]) -> list[Path]:
    return [directory / f"{dataset_id}{suffix}" for suffix in suffixes if (directory / f"{dataset_id}{suffix}").is_file()]


def _result_file(directory: Path, dataset_id: str, suffixes: tuple[str, ...]) -> Path | None:
    matches = _candidate_files(directory, dataset_id, suffixes)
    if len(matches) > 1:
        raise ValueError(f"{dataset_id}: multiple result files found: {', '.join(path.name for path in matches)}")
    return matches[0] if matches else None


def _evaluator(model) -> tuple[str, tuple[str, ...]] | None:
    if model.assay == "bulk_rna_seq" and model.kind == "real" and model.validation is not None:
        return "differential_expression", RESULT_SUFFIXES
    if model.assay == "whole_genome_dna_seq" and model.validation is not None:
        return "exact_snv", VCF_SUFFIXES
    return None


def _selection(root: Path, dataset_ids: list[str] | None) -> list[tuple]:
    records = registry(Path(root))
    requested = list(dict.fromkeys(dataset_ids or []))
    unknown = [dataset_id for dataset_id in requested if dataset_id not in records]
    if unknown:
        raise ValueError(f"unknown dataset IDs: {', '.join(unknown)}")
    if requested:
        unsupported = [dataset_id for dataset_id in requested if _evaluator(records[dataset_id][0]) is None]
        if unsupported:
            raise ValueError(f"no result comparator for: {', '.join(unsupported)}")
        return [records[dataset_id] for dataset_id in requested]
    selected = [record for record in records.values() if _evaluator(record[0]) is not None]
    if not selected:
        raise ValueError("the collection has no comparable benchmarks")
    return selected


def _assay_summary(results: list[dict]) -> dict:
    assays = {}
    for assay in sorted({item["assay"] for item in results}):
        assays[assay] = _summary([item for item in results if item["assay"] == assay])
    return assays


def evaluate_suite(
    root: Path,
    results_directory: Path,
    dataset_ids: list[str] | None = None,
    detail_limit: int = 20,
) -> dict:
    directory = Path(results_directory)
    if not directory.is_dir():
        raise ValueError(f"{directory}: results directory does not exist")
    if detail_limit < 0:
        raise ValueError("detail limit must be zero or greater")
    receipt = read_method_receipt(directory)
    results = []
    for model, folder in _selection(Path(root), dataset_ids):
        comparator, suffixes = _evaluator(model)  # type: ignore[misc]
        try:
            path = _result_file(directory, model.id, suffixes)
            if path is None:
                expected = " or ".join(f"{model.id}{suffix}" for suffix in suffixes)
                results.append({
                    "id": model.id,
                    "assay": model.assay,
                    "comparator": comparator,
                    "status": "fail",
                    "reason": f"missing result file; expected {expected}",
                })
                continue
            details = compare(model, folder, path, detail_limit) if comparator == "differential_expression" else compare_variants(model, folder, path, detail_limit)
            item = {
                "id": model.id,
                "assay": model.assay,
                "comparator": comparator,
                "status": details["status"],
                "input": str(path),
                "details": details,
            }
            if details["status"] == "fail":
                item["reason"] = "; ".join(details["reasons"])
            results.append(item)
        except (ValueError, OSError, KeyError) as error:
            results.append({
                "id": model.id,
                "assay": model.assay,
                "comparator": comparator,
                "status": "fail",
                "reason": str(error),
            })
    summary = _summary(results)
    summary["assays"] = _assay_summary(results)
    return {
        "operation": "evaluate",
        "selection": {"ids": dataset_ids or [], "default": "all comparable benchmarks"},
        "results_directory": str(directory),
        "method": receipt,
        "summary": summary,
        "results": results,
    }


def _leader(name: str, directory: Path, report: dict) -> dict:
    summary = report["summary"]
    receipt = report.get("method")
    return {
        "name": name,
        "results_directory": str(directory),
        "receipt": receipt,
        "status": summary["status"],
        "total": summary["total"],
        "passed": summary["passed"],
        "failed": summary["failed"],
        "skipped": summary["skipped"],
        "pass_rate": summary["passed"] / summary["total"] if summary["total"] else 0.0,
        "assays": summary["assays"],
    }


def evaluate_matrix(
    root: Path,
    methods: list[tuple[str, Path]],
    dataset_ids: list[str] | None = None,
    detail_limit: int = 20,
) -> dict:
    if len(methods) < 2:
        raise ValueError("a multi-assay matrix requires at least two --method entries")
    names = [name for name, _ in methods]
    if len(names) != len(set(names)):
        raise ValueError("method names must be unique")
    for name in names:
        if not METHOD_NAME.fullmatch(name):
            raise ValueError(f"{name}: invalid method name")
    cells = []
    leaders = []
    benchmark_ids = None
    for name, directory in methods:
        report = evaluate_suite(root, directory, dataset_ids, detail_limit)
        ids = [item["id"] for item in report["results"]]
        if benchmark_ids is None:
            benchmark_ids = ids
        elif ids != benchmark_ids:
            raise ValueError(f"{name}: benchmark selection differs from the first method")
        cells.extend({"method": name, **item} for item in report["results"])
        leaders.append(_leader(name, directory, report))
    leaderboard = sorted(leaders, key=lambda item: (-item["passed"], item["failed"], item["name"]))
    summary = _summary(cells)
    summary.update({"methods": len(methods), "benchmarks": len(benchmark_ids or []), "assays": _assay_summary(cells)})
    return {
        "operation": "evaluate_matrix",
        "selection": {"ids": dataset_ids or [], "default": "all comparable benchmarks"},
        "ranking": "passed cases descending, failed cases ascending, then method name; assay metrics are not pooled",
        "summary": summary,
        "leaderboard": leaderboard,
        "results": cells,
    }
