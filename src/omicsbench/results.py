"""Hash-locked mappings from benchmark IDs to pipeline result files."""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator

from .hashing import contained, digest
from .models import StrictModel, safe_relative
from .registry import default_root, registry


MANIFEST_NAME = "results.json"
RESULT_ID = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
COLUMN_KEY = re.compile(r"^[a-z][a-z0-9_]{0,63}$")
RESULT_SUFFIXES = {
    "differential_expression": (".csv", ".csv.gz", ".tsv", ".tsv.gz"),
    "exact_snv": (".vcf", ".vcf.gz"),
    "peptide_identification": (".tsv", ".tsv.gz", ".peptides.tsv", ".peptides.tsv.gz"),
}
DEFAULT_COLUMNS = {
    "differential_expression": {"gene_id": "gene_id", "log2_fold_change": "log2_fold_change"},
    "exact_snv": {},
    "peptide_identification": {"protein_id": "protein_id", "peptide": "peptide", "q_value": "q_value"},
}
REQUIRED_COLUMNS = {
    "differential_expression": {"gene_id", "log2_fold_change"},
    "exact_snv": set(),
    "peptide_identification": {"protein_id", "peptide"},
}


def comparator_for(model) -> str | None:
    if model.assay == "bulk_rna_seq" and model.kind == "real" and model.validation is not None:
        return "differential_expression"
    if model.assay == "whole_genome_dna_seq" and model.validation is not None:
        return "exact_snv"
    if model.assay == "sequence_protein" and model.validation is not None and model.validation.baseline_version == "peptide-identification-v1":
        return "peptide_identification"
    return None


def result_suffix(path: Path) -> str:
    suffixes = [value.lower() for value in Path(path).suffixes]
    return "".join(suffixes[-2:]) if suffixes[-1:] == [".gz"] else (suffixes[-1] if suffixes else "")


def effective_columns(comparator: str, overrides: dict[str, str] | None = None) -> dict[str, str]:
    if comparator not in DEFAULT_COLUMNS:
        raise ValueError(f"unsupported result comparator: {comparator}")
    overrides = dict(overrides or {})
    allowed = set(DEFAULT_COLUMNS[comparator])
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise ValueError(f"{comparator}: unsupported column mappings: {', '.join(unknown)}")
    if comparator == "exact_snv" and overrides:
        raise ValueError("exact_snv: VCF columns cannot be remapped")
    columns = {**DEFAULT_COLUMNS[comparator], **overrides}
    required = REQUIRED_COLUMNS[comparator]
    if not required.issubset(columns):
        raise ValueError(f"{comparator}: missing required column mappings")
    required_sources = [columns[key] for key in sorted(required)]
    if len(required_sources) != len(set(required_sources)):
        raise ValueError(f"{comparator}: required fields must use different source columns")
    return columns


class ResultEntry(StrictModel):
    benchmark_id: str = Field(pattern=r"^[a-z0-9][a-z0-9-]{0,63}$")
    path: str
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    columns: dict[str, str] = Field(default_factory=dict)

    _path = field_validator("path")(safe_relative)

    @field_validator("columns")
    @classmethod
    def valid_columns(cls, value):
        for key, source in value.items():
            if not COLUMN_KEY.fullmatch(key) or not source or source.strip() != source or any(character in source for character in "\r\n\t"):
                raise ValueError("column mappings require simple logical keys and non-empty source headers")
        return value


class ResultManifest(StrictModel):
    schema_version: Literal["1.0"]
    results: list[ResultEntry] = Field(min_length=1)

    @model_validator(mode="after")
    def unique_benchmarks(self):
        identifiers = [entry.benchmark_id for entry in self.results]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("result manifest benchmark IDs must be unique")
        return self


def parse_result(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise ValueError("result mapping must use BENCHMARK=RELATIVE_PATH")
    benchmark_id, relative = value.split("=", 1)
    if not RESULT_ID.fullmatch(benchmark_id):
        raise ValueError(f"{benchmark_id or '(empty)'}: invalid benchmark ID")
    safe_relative(relative)
    return benchmark_id, relative


def parse_column(value: str) -> tuple[str, str, str]:
    if "=" not in value or ":" not in value.split("=", 1)[0]:
        raise ValueError("column mapping must use BENCHMARK:FIELD=SOURCE_COLUMN")
    left, source = value.split("=", 1)
    benchmark_id, logical = left.split(":", 1)
    if not RESULT_ID.fullmatch(benchmark_id) or not COLUMN_KEY.fullmatch(logical):
        raise ValueError("column mapping has an invalid benchmark ID or logical field")
    if not source or source.strip() != source or any(character in source for character in "\r\n\t"):
        raise ValueError("source column must be a non-empty header")
    return benchmark_id, logical, source


def _record_contract(records: dict, benchmark_id: str) -> tuple[object, str]:
    if benchmark_id not in records:
        raise ValueError(f"unknown benchmark ID: {benchmark_id}")
    model = records[benchmark_id][0]
    comparator = comparator_for(model)
    if comparator is None:
        raise ValueError(f"{benchmark_id}: no result comparator")
    return model, comparator


def _validate_entry(directory: Path, entry: ResultEntry, records: dict) -> dict:
    model, comparator = _record_contract(records, entry.benchmark_id)
    path = contained(directory, entry.path)
    if not path.is_file():
        raise ValueError(f"{entry.benchmark_id}: result file is missing: {entry.path}")
    suffix = result_suffix(path)
    if suffix not in RESULT_SUFFIXES[comparator]:
        allowed = ", ".join(RESULT_SUFFIXES[comparator])
        raise ValueError(f"{entry.benchmark_id}: expected one of {allowed}, found {entry.path}")
    observed = digest(path)
    if observed != entry.sha256:
        raise ValueError(f"{entry.benchmark_id}: SHA-256 mismatch for {entry.path}")
    columns = effective_columns(comparator, entry.columns)
    return {
        "id": entry.benchmark_id,
        "assay": model.assay,
        "comparator": comparator,
        "path": path,
        "relative_path": entry.path,
        "bytes": path.stat().st_size,
        "sha256": observed,
        "columns": columns,
    }


def read_result_manifest(directory: Path, records: dict | None = None) -> tuple[dict, dict[str, dict]] | None:
    directory = Path(directory)
    path = directory / MANIFEST_NAME
    if not path.exists():
        return None
    if not path.is_file():
        raise ValueError(f"{path}: result manifest is not a file")
    try:
        manifest = ResultManifest.model_validate_json(path.read_text(encoding="utf-8"))
    except (ValueError, UnicodeDecodeError) as error:
        raise ValueError(f"{path}: invalid result manifest: {error}") from None
    records = records or registry(default_root())
    resolved = [_validate_entry(directory, entry, records) for entry in manifest.results]
    return (
        {
            "path": str(path),
            "sha256": digest(path),
            "schema_version": manifest.schema_version,
            "entries": len(resolved),
        },
        {item["id"]: item for item in resolved},
    )


def create_result_manifest(
    root: Path,
    directory: Path,
    result_values: list[str],
    column_values: list[str],
    force: bool = False,
) -> dict:
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"{directory}: results directory does not exist")
    pairs = [parse_result(value) for value in result_values]
    if not pairs:
        raise ValueError("at least one --result mapping is required")
    identifiers = [identifier for identifier, _ in pairs]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("result mappings must use unique benchmark IDs")
    columns: dict[str, dict[str, str]] = {}
    for benchmark_id, logical, source in map(parse_column, column_values):
        if benchmark_id not in identifiers:
            raise ValueError(f"{benchmark_id}: column mapping has no result mapping")
        if logical in columns.setdefault(benchmark_id, {}):
            raise ValueError(f"{benchmark_id}: duplicate column mapping for {logical}")
        columns[benchmark_id][logical] = source

    records = registry(Path(root))
    entries = []
    for benchmark_id, relative in pairs:
        _, comparator = _record_contract(records, benchmark_id)
        path = contained(directory, relative)
        if not path.is_file():
            raise ValueError(f"{benchmark_id}: result file is missing: {relative}")
        suffix = result_suffix(path)
        if suffix not in RESULT_SUFFIXES[comparator]:
            allowed = ", ".join(RESULT_SUFFIXES[comparator])
            raise ValueError(f"{benchmark_id}: expected one of {allowed}, found {relative}")
        overrides = columns.get(benchmark_id, {})
        effective_columns(comparator, overrides)
        entries.append(ResultEntry(benchmark_id=benchmark_id, path=relative, sha256=digest(path), columns=overrides))
    manifest = ResultManifest(schema_version="1.0", results=entries)
    target = directory / MANIFEST_NAME
    if target.exists() and not force:
        raise ValueError(f"{target}: output already exists; use --force to replace it")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    try:
        temporary.write_text(json.dumps(manifest.model_dump(mode="json"), indent=2) + "\n", encoding="utf-8", newline="\n")
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    loaded = read_result_manifest(directory, records)
    if loaded is None:
        raise ValueError("result manifest could not be reloaded")
    info, resolved = loaded
    items = [
        {
            "id": item["id"],
            "assay": item["assay"],
            "comparator": item["comparator"],
            "status": "pass",
            "path": item["relative_path"],
            "bytes": item["bytes"],
            "sha256": item["sha256"],
            "columns": item["columns"],
        }
        for item in resolved.values()
    ]
    return {
        "operation": "index_results",
        "result_manifest": info,
        "summary": {"status": "pass", "total": len(items), "passed": len(items), "failed": 0, "skipped": 0},
        "results": items,
    }
