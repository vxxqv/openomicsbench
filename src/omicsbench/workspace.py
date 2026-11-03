"""Prepare and grade a local workspace for independently produced results."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path

from . import __version__
from .evaluate import evaluate_suite
from .evidence import create_evidence_crate
from .hashing import contained, digest
from .models import safe_relative
from .registry import registry
from .results import comparator_for, effective_columns
from .suite import write_reports


RESULT_NAMES = {
    "differential_expression": ".tsv",
    "exact_snv": ".vcf",
    "peptide_identification": ".peptides.tsv",
    "single_cell_clustering": ".cells.tsv",
    "genomic_intervals": ".narrowPeak",
    "read_taxonomy": ".taxonomy.tsv",
}


def _selected(root: Path, identifiers: list[str]):
    if not identifiers:
        raise ValueError("select at least one benchmark with --id")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("benchmark IDs must be unique")
    records = registry(root)
    selected = []
    for identifier in identifiers:
        if identifier not in records:
            raise ValueError(f"unknown benchmark ID: {identifier}")
        model, folder = records[identifier]
        comparator = comparator_for(model)
        if comparator is None or model.status != "validated":
            raise ValueError(f"{identifier}: no validated result comparator")
        inputs = [item for item in model.files if item.tier in {"nano", "pocket"} or (item.tier == "metadata" and item.role == "samples")]
        if not inputs:
            raise ValueError(f"{identifier}: no distributable input files")
        selected.append((model, folder, comparator, inputs))
    return selected


def prepare_workspace(root: Path, destination: Path, identifiers: list[str]) -> dict:
    root = Path(root).resolve()
    destination = Path(destination).absolute()
    if destination.exists():
        raise ValueError(f"{destination}: output already exists; choose a new directory")
    if not destination.parent.is_dir():
        raise ValueError(f"{destination.parent}: parent directory does not exist")
    selected = _selected(root, identifiers)
    cases = []
    for model, folder, comparator, inputs in selected:
        files = []
        for item in inputs:
            source = contained(folder, item.path)
            if not source.is_file() or source.stat().st_size != item.bytes or digest(source) != item.sha256:
                raise ValueError(f"{model.id}: input differs from its declared inventory: {item.path}")
            files.append({
                "path": f"inputs/{model.id}/{item.path}",
                "role": item.role,
                "bytes": item.bytes,
                "sha256": item.sha256,
            })
        cases.append({
            "id": model.id,
            "title": model.title,
            "assay": model.assay,
            "comparator": comparator,
            "benchmark_manifest_sha256": digest(folder / "manifest.json"),
            "suggested_result": f"results/{model.id}{RESULT_NAMES[comparator]}",
            "expected_columns": effective_columns(comparator),
            "source_accession": model.source.accession,
            "source_citation": model.source.citation,
            "rights": model.rights.model_dump(mode="json"),
            "samples": [sample.model_dump(mode="json") for sample in model.samples],
            "limitations": model.limitations,
            "inputs": files,
        })
    payload = {"schema_version": "1.0", "package_version": __version__, "cases": cases}
    temporary = Path(tempfile.mkdtemp(prefix=destination.name + ".tmp-", dir=destination.parent))
    if temporary.resolve().parent != destination.parent.resolve():
        raise ValueError("temporary workspace escaped its destination directory")
    try:
        for (model, folder, _, inputs), case in zip(selected, cases):
            for item, entry in zip(inputs, case["inputs"]):
                target = contained(temporary, entry["path"])
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(contained(folder, item.path), target)
                if digest(target) != item.sha256:
                    raise ValueError(f"{model.id}: copied input failed its checksum")
        (temporary / "results").mkdir()
        (temporary / "workspace.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
        os.replace(temporary, destination)
    except Exception:
        if temporary.exists() and temporary.resolve().parent == destination.parent.resolve():
            shutil.rmtree(temporary)
        raise
    return {"operation": "prepare_workspace", "directory": str(destination), "benchmarks": identifiers, "files": sum(len(case["inputs"]) for case in cases), "summary": {"status": "pass", "total": len(cases), "passed": len(cases), "failed": 0, "skipped": 0}}


def verify_workspace(root: Path, directory: Path) -> dict:
    root = Path(root).resolve()
    directory = Path(directory).resolve()
    manifest_path = directory / "workspace.json"
    if not manifest_path.is_file():
        raise ValueError(f"{manifest_path}: workspace manifest is missing")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (ValueError, UnicodeDecodeError) as error:
        raise ValueError(f"{manifest_path}: invalid JSON: {error}") from None
    if manifest.get("schema_version") != "1.0" or not isinstance(manifest.get("cases"), list) or not manifest["cases"]:
        raise ValueError(f"{manifest_path}: unsupported or empty workspace")
    records = registry(root)
    seen_ids = set()
    expected_files = set()
    for case in manifest["cases"]:
        if not isinstance(case, dict) or not isinstance(case.get("id"), str):
            raise ValueError("workspace case has no benchmark ID")
        identifier = case["id"]
        if identifier in seen_ids or identifier not in records:
            raise ValueError(f"{identifier}: duplicate or unknown workspace benchmark")
        seen_ids.add(identifier)
        model, folder = records[identifier]
        if case.get("benchmark_manifest_sha256") != digest(folder / "manifest.json") or case.get("comparator") != comparator_for(model):
            raise ValueError(f"{identifier}: benchmark contract differs from the installed collection")
        expected_fields = {
            "title": model.title,
            "assay": model.assay,
            "expected_columns": effective_columns(case["comparator"]),
            "source_accession": model.source.accession,
            "source_citation": model.source.citation,
            "rights": model.rights.model_dump(mode="json"),
            "samples": [sample.model_dump(mode="json") for sample in model.samples],
            "limitations": model.limitations,
        }
        if any(case.get(key) != value for key, value in expected_fields.items()):
            raise ValueError(f"{identifier}: workspace provenance or result contract differs from the collection")
        if case.get("suggested_result") != f"results/{identifier}{RESULT_NAMES[case['comparator']]}":
            raise ValueError(f"{identifier}: suggested result path differs from its comparator")
        if not isinstance(case.get("inputs"), list) or not case["inputs"]:
            raise ValueError(f"{identifier}: workspace inputs are missing")
        declared = {f"inputs/{identifier}/{item.path}": item for item in model.files if item.tier in {"nano", "pocket"} or (item.tier == "metadata" and item.role == "samples")}
        observed = set()
        for entry in case["inputs"]:
            if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
                raise ValueError(f"{identifier}: malformed input entry")
            relative = safe_relative(entry["path"])
            if relative in observed or relative not in declared:
                raise ValueError(f"{identifier}: duplicate or undeclared workspace input: {relative}")
            observed.add(relative)
            expected_files.add(relative)
            item = declared[relative]
            path = contained(directory, relative)
            if not path.is_file() or path.stat().st_size != item.bytes or digest(path) != item.sha256:
                raise ValueError(f"{identifier}: workspace input failed its checksum: {relative}")
            if entry.get("bytes") != item.bytes or entry.get("sha256") != item.sha256 or entry.get("role") != item.role:
                raise ValueError(f"{identifier}: workspace inventory differs: {relative}")
        if observed != set(declared):
            raise ValueError(f"{identifier}: workspace input inventory is incomplete")
    input_root = directory / "inputs"
    actual_files = {path.relative_to(directory).as_posix() for path in input_root.rglob("*") if path.is_file()} if input_root.is_dir() else set()
    if actual_files != expected_files:
        raise ValueError(f"workspace input file set differs: missing={sorted(expected_files - actual_files)}, extra={sorted(actual_files - expected_files)}")
    count = len(seen_ids)
    return {"operation": "verify_workspace", "directory": str(directory), "manifest_sha256": digest(manifest_path), "benchmarks": [case["id"] for case in manifest["cases"]], "files": len(expected_files), "summary": {"status": "pass", "total": count, "passed": count, "failed": 0, "skipped": 0}}


def grade_workspace(root: Path, directory: Path, force: bool = False) -> dict:
    directory = Path(directory).resolve()
    verified = verify_workspace(root, directory)
    results = directory / "results"
    if not results.is_dir():
        raise ValueError(f"{results}: results directory is missing")
    outputs = {name: directory / name for name in ("report.json", "report.html", "report.csv", "report.xml", "evidence.zip")}
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing and not force:
        raise ValueError(f"{existing[0]}: output already exists; use --force to replace the previous report")
    report = evaluate_suite(root, results, verified["benchmarks"])
    report["workspace"] = {"path": str(directory / "workspace.json"), "manifest_sha256": verified["manifest_sha256"], "input_files": verified["files"]}
    write_reports(report, json_path=outputs["report.json"], html_path=outputs["report.html"], csv_path=outputs["report.csv"], junit_path=outputs["report.xml"], force=force)
    create_evidence_crate(report, outputs["evidence.zip"], force=force)
    return {"operation": "grade_workspace", "directory": str(directory), "summary": report["summary"], "reports": {name: str(path) for name, path in outputs.items()}}
