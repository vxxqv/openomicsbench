"""Find external result files in a prepared workspace without guessing their meaning."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

from .hashing import contained, digest
from .results import (
    REQUIRED_COLUMNS,
    RESULT_SUFFIXES,
    comparator_for,
    create_result_manifest,
    parse_column,
    parse_result,
    result_suffix,
)
from .registry import registry
from .workspace import verify_workspace


HEADER_NAMES = {
    "differential_expression": {
        "gene_id": ("gene_id", "Geneid", "GeneID", "gene", "Gene", "feature_id"),
        "log2_fold_change": ("log2_fold_change", "log2FoldChange", "log2FC", "logFC"),
    },
    "peptide_identification": {
        "protein_id": ("protein_id", "Protein", "protein", "ProteinId"),
        "peptide": ("peptide", "Peptide", "Sequence", "sequence"),
        "q_value": ("q_value", "QValue", "qvalue", "q-value"),
    },
    "single_cell_clustering": {
        "cell_id": ("cell_id", "Barcode", "barcode", "cell", "Cell"),
        "label": ("label", "Cluster", "cluster", "leiden", "Louvain"),
    },
}


def _header(path: Path) -> list[str]:
    suffix = result_suffix(path)
    delimiter = "," if suffix.startswith(".csv") else "\t"
    opener = gzip.open if suffix.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        fields = next(csv.reader(stream, delimiter=delimiter), None)
    if not fields or len(fields) != len(set(fields)):
        raise ValueError("the table header is empty or contains duplicate columns")
    return fields


def _columns(path: Path, comparator: str, overrides: dict[str, str]) -> dict[str, str]:
    if comparator not in HEADER_NAMES:
        if overrides:
            raise ValueError("this result format has no remappable columns")
        return {}
    fields = _header(path)
    names = HEADER_NAMES[comparator]
    unknown = sorted(set(overrides) - set(names))
    if unknown:
        raise ValueError(f"unsupported column mappings: {', '.join(unknown)}")
    selected = {}
    for logical, aliases in names.items():
        if logical in overrides:
            source = overrides[logical]
            if source not in fields:
                raise ValueError(f"mapped {logical} column is absent: {source}")
        else:
            matches = [alias for alias in aliases if alias in fields]
            if len(matches) > 1:
                raise ValueError(f"ambiguous {logical} columns: {', '.join(matches)}; use --column")
            source = matches[0] if matches else None
        if source is None:
            if logical in REQUIRED_COLUMNS[comparator]:
                raise ValueError(f"no {logical} column; use --column ID:{logical}=HEADER")
            continue
        if source in selected.values():
            raise ValueError("two result fields use the same source column")
        selected[logical] = source
    return selected


def _belongs(path: Path, identifier: str) -> bool:
    return path.name.startswith((identifier + ".", identifier + "_", identifier + "-")) or identifier in path.parts[:-1]


def discover_workspace_results(
    root: Path,
    workspace: Path,
    result_values: list[str] | None = None,
    column_values: list[str] | None = None,
    write: bool = False,
    force: bool = False,
) -> dict:
    workspace = Path(workspace).resolve()
    verified = verify_workspace(root, workspace)
    directory = workspace / "results"
    if not directory.is_dir():
        raise ValueError(f"{directory}: results directory is missing")
    if force and not write:
        raise ValueError("--force requires --write")
    identifiers = verified["benchmarks"]
    explicit = dict(parse_result(value) for value in result_values or [])
    if len(explicit) != len(result_values or []):
        raise ValueError("result mappings must use unique benchmark IDs")
    if set(explicit) - set(identifiers):
        raise ValueError(f"result mapping is outside the workspace: {', '.join(sorted(set(explicit) - set(identifiers)))}")
    mapped: dict[str, dict[str, str]] = {}
    for identifier, logical, source in map(parse_column, column_values or []):
        if identifier not in identifiers:
            raise ValueError(f"{identifier}: column mapping is outside the workspace")
        if logical in mapped.setdefault(identifier, {}):
            raise ValueError(f"{identifier}: duplicate {logical} column mapping")
        mapped[identifier][logical] = source
    records = registry(Path(root))
    candidates = sorted((path for path in directory.rglob("*") if path.is_file() and path.name not in {"results.json", "method.json"}), key=lambda path: path.as_posix())
    cases = []
    result_arguments = []
    column_arguments = []
    for identifier in identifiers:
        comparator = comparator_for(records[identifier][0])
        assert comparator is not None
        allowed = RESULT_SUFFIXES[comparator]
        if identifier in explicit:
            relative = explicit[identifier]
            path = contained(directory, relative)
            matches = [path] if path.is_file() else []
        else:
            matches = [path for path in candidates if _belongs(path.relative_to(directory), identifier) and result_suffix(path) in allowed]
        item = {"id": identifier, "comparator": comparator}
        try:
            if not matches:
                raise ValueError("no matching result file; use --result ID=RELATIVE_PATH")
            if len(matches) > 1:
                raise ValueError("multiple matching files: " + ", ".join(path.relative_to(directory).as_posix() for path in matches))
            path = matches[0]
            relative = path.relative_to(directory).as_posix()
            if result_suffix(path) not in allowed:
                raise ValueError(f"unsupported result suffix: {relative}")
            path = contained(directory, relative)
            columns = _columns(path, comparator, mapped.get(identifier, {}))
            item.update(status="pass", path=relative, sha256=digest(path), columns=columns)
            result_arguments.append(f"{identifier}={relative}")
            column_arguments.extend(f"{identifier}:{logical}={source}" for logical, source in columns.items())
        except (OSError, UnicodeError, ValueError) as error:
            item.update(status="fail", reason=str(error))
        cases.append(item)
    failed = sum(item["status"] == "fail" for item in cases)
    report = {
        "operation": "discover_workspace_results",
        "workspace": str(workspace),
        "summary": {"status": "fail" if failed else "pass", "total": len(cases), "passed": len(cases) - failed, "failed": failed, "skipped": 0},
        "results": cases,
    }
    if write:
        if failed:
            raise ValueError("result discovery failed; no manifest was written: " + "; ".join(f"{item['id']}: {item['reason']}" for item in cases if item["status"] == "fail"))
        indexed = create_result_manifest(root, directory, result_arguments, column_arguments, force)
        report["result_manifest"] = indexed["result_manifest"]
    return report
