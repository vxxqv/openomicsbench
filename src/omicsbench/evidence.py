"""Deterministic, verifiable evidence crates for multi-assay evaluations."""
from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import zipfile
from pathlib import Path, PurePosixPath

from . import __version__
from .suite import csv_report, html_report


ZIP_TIME = (1980, 1, 1, 0, 0, 0)
SAFE_METHOD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
REQUIRED = {"report.json", "report.csv", "report.html", "ro-crate-metadata.json", "checksums.sha256"}


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _safe_name(name: str) -> str:
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name or str(path) != name:
        raise ValueError(f"unsafe evidence path: {name}")
    return name


def _add_source(files: dict[str, bytes], archive_name: str, source: Path) -> None:
    name = _safe_name(archive_name)
    if name in files:
        raise ValueError(f"duplicate evidence path: {name}")
    path = Path(source)
    if not path.is_file():
        raise ValueError(f"{path}: evidence input is missing")
    files[name] = path.read_bytes()


def _portable_report(report: dict) -> tuple[dict, dict[str, bytes]]:
    if report.get("operation") not in {"evaluate", "evaluate_matrix"}:
        raise ValueError("evidence crates require an evaluate or evaluate_matrix report")
    portable = copy.deepcopy(report)
    portable.pop("reports", None)
    files: dict[str, bytes] = {}
    for item in portable.get("results", []):
        source_value = item.get("input")
        if not source_value:
            continue
        source = Path(source_value)
        method = item.get("method")
        prefix = f"inputs/{method}" if method else "inputs"
        if method and not SAFE_METHOD.fullmatch(method):
            raise ValueError(f"invalid method name in report: {method}")
        archive_name = f"{prefix}/{source.name}"
        if archive_name not in files:
            _add_source(files, archive_name, source)
        item["input"] = archive_name
        details = item.get("details")
        if isinstance(details, dict) and details.get("input"):
            details["input"] = archive_name

    if portable.get("operation") == "evaluate":
        portable["results_directory"] = "inputs"
        receipt = portable.get("method")
        if receipt:
            source = Path(receipt["path"])
            archive_name = "method.json"
            _add_source(files, archive_name, source)
            receipt["path"] = archive_name
    else:
        for leader in portable.get("leaderboard", []):
            name = leader["name"]
            leader["results_directory"] = f"inputs/{name}"
            receipt = leader.get("receipt")
            if receipt:
                archive_name = f"methods/{name}/method.json"
                _add_source(files, archive_name, Path(receipt["path"]))
                receipt["path"] = archive_name
    return portable, files


def _media_type(name: str) -> str:
    lower = name.lower()
    if lower.endswith(".json"):
        return "application/json"
    if lower.endswith(".csv"):
        return "text/csv"
    if lower.endswith(".html"):
        return "text/html"
    if lower.endswith(".vcf") or lower.endswith(".vcf.gz"):
        return "text/vcf"
    if lower.endswith(".tsv") or lower.endswith(".tsv.gz"):
        return "text/tab-separated-values"
    return "application/octet-stream"


def _crate_metadata(files: dict[str, bytes], report: dict) -> bytes:
    parts = [
        {"@id": "ro-crate-metadata.json", "@type": "CreativeWork", "about": {"@id": "./"}, "conformsTo": {"@id": "https://w3id.org/ro/crate/1.3"}},
        {
            "@id": "./",
            "@type": "Dataset",
            "name": "OpenOmicsBench evaluation evidence",
            "description": "Portable reports and exact submitted result files for an OpenOmicsBench evaluation.",
            "license": {"@id": "https://spdx.org/licenses/Apache-2.0"},
            "hasPart": [{"@id": name} for name in sorted(files)],
            "mentions": {"@id": "#openomicsbench"},
        },
        {
            "@id": "#openomicsbench",
            "@type": "SoftwareApplication",
            "name": "OpenOmicsBench",
            "softwareVersion": __version__,
            "url": "https://github.com/vxxqv/openomicsbench",
        },
    ]
    for name, data in sorted(files.items()):
        parts.append({
            "@id": name,
            "@type": "File",
            "name": Path(name).name,
            "encodingFormat": _media_type(name),
            "contentSize": str(len(data)),
            "identifier": f"sha256:{_sha256(data)}",
        })
    metadata = {"@context": "https://w3id.org/ro/crate/1.3/context", "@graph": parts}
    return (json.dumps(metadata, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(_safe_name(name), ZIP_TIME)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o100644 << 16
    info.create_system = 3
    return info


def create_evidence_crate(report: dict, output: Path, force: bool = False) -> dict:
    output = Path(output)
    if output.suffix.lower() != ".zip":
        raise ValueError("evidence output must end in .zip")
    if output.exists() and not force:
        raise ValueError(f"{output}: output already exists; use --force to replace it")
    portable, files = _portable_report(report)
    files["report.json"] = (json.dumps(portable, indent=2, allow_nan=False) + "\n").encode("utf-8")
    files["report.csv"] = csv_report(portable).encode("utf-8")
    files["report.html"] = html_report(portable).encode("utf-8")
    files["ro-crate-metadata.json"] = _crate_metadata(files, portable)
    checksum_lines = [f"{_sha256(data)}  {name}" for name, data in sorted(files.items())]
    files["checksums.sha256"] = ("\n".join(checksum_lines) + "\n").encode("utf-8")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for name, data in sorted(files.items()):
                archive.writestr(_zip_info(name), data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    data = output.read_bytes()
    return {
        "path": str(output),
        "files": len(files),
        "inputs": sum(name.startswith("inputs/") for name in files),
        "bytes": len(data),
        "sha256": _sha256(data),
        "format": "RO-Crate 1.3 compatible ZIP",
    }


def _read_json(archive: zipfile.ZipFile, name: str) -> dict:
    try:
        value = json.loads(archive.read(name).decode("utf-8"))
    except (KeyError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{name}: invalid JSON: {error}") from None
    if not isinstance(value, dict):
        raise ValueError(f"{name}: expected a JSON object")
    return value


def verify_evidence_crate(path: Path) -> dict:
    path = Path(path)
    try:
        archive = zipfile.ZipFile(path, "r")
    except (OSError, zipfile.BadZipFile) as error:
        raise ValueError(f"{path}: invalid evidence ZIP: {error}") from None
    with archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError(f"{path}: duplicate ZIP members")
        for info in archive.infolist():
            _safe_name(info.filename)
            mode = (info.external_attr >> 16) & 0o170000
            if info.is_dir() or mode == 0o120000:
                raise ValueError(f"{path}: directories and symbolic links are not allowed: {info.filename}")
        missing = sorted(REQUIRED - set(names))
        if missing:
            raise ValueError(f"{path}: missing required evidence files: {', '.join(missing)}")
        try:
            lines = archive.read("checksums.sha256").decode("utf-8").splitlines()
        except UnicodeDecodeError as error:
            raise ValueError(f"{path}: checksums.sha256 is not UTF-8") from error
        expected = {}
        for line_number, line in enumerate(lines, start=1):
            parts = line.split("  ", 1)
            if len(parts) != 2 or not re.fullmatch(r"[0-9a-f]{64}", parts[0]):
                raise ValueError(f"{path}: invalid checksum line {line_number}")
            name = _safe_name(parts[1])
            if name in expected:
                raise ValueError(f"{path}: duplicate checksum entry {name}")
            expected[name] = parts[0]
        members = set(names) - {"checksums.sha256"}
        if set(expected) != members:
            raise ValueError(f"{path}: checksum inventory does not match ZIP members")
        for name, checksum in expected.items():
            if _sha256(archive.read(name)) != checksum:
                raise ValueError(f"{path}: checksum mismatch for {name}")
        report = _read_json(archive, "report.json")
        if report.get("operation") not in {"evaluate", "evaluate_matrix"} or not isinstance(report.get("results"), list):
            raise ValueError(f"{path}: unsupported evaluation report")
        metadata = _read_json(archive, "ro-crate-metadata.json")
        if metadata.get("@context") != "https://w3id.org/ro/crate/1.3/context" or not isinstance(metadata.get("@graph"), list):
            raise ValueError(f"{path}: invalid RO-Crate 1.3 metadata")
        identifiers = {entry.get("@id") for entry in metadata["@graph"] if isinstance(entry, dict)}
        if not {"./", "ro-crate-metadata.json"}.issubset(identifiers):
            raise ValueError(f"{path}: RO-Crate root entities are missing")
    data = path.read_bytes()
    return {
        "status": "pass",
        "path": str(path),
        "operation": report["operation"],
        "files": len(names),
        "inputs": sum(name.startswith("inputs/") for name in names),
        "bytes": len(data),
        "sha256": _sha256(data),
        "ro_crate": "1.3",
    }
