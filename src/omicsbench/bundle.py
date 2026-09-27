"""Create and verify portable, deterministic benchmark bundles."""
from __future__ import annotations

import hashlib
import json
import os
import zipfile
from pathlib import Path, PurePosixPath

from . import __version__
from .models import Dataset
from .registry import registry


BUNDLE_FORMAT = "openomicsbench-bundle-v1"
SUPPORT_FILES = ("CITATION.cff", "LICENSE", "LICENSE-METADATA", "NOTICE")
FIXED_TIME = (1980, 1, 1, 0, 0, 0)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _safe_name(value: str) -> str:
    if "\\" in value or "\x00" in value:
        raise ValueError(f"unsafe bundle path: {value}")
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"unsafe bundle path: {value}")
    return path.as_posix()


def _selected(root: Path, dataset_ids: list[str] | None, assay: str | None):
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
        raise ValueError("no datasets matched the bundle selection")
    return selected


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(_safe_name(name), FIXED_TIME)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = 0o100644 << 16
    return info


def create_bundle(
    root: Path,
    destination: Path,
    dataset_ids: list[str] | None = None,
    assay: str | None = None,
    force: bool = False,
) -> dict:
    root = Path(root).resolve()
    destination = Path(destination)
    if destination.exists() and not force:
        raise ValueError(f"{destination}: output already exists; use --force to replace it")
    selected = _selected(root, dataset_ids, assay)
    files: dict[str, bytes] = {}
    for relative in SUPPORT_FILES:
        path = root / relative
        if not path.is_file():
            raise ValueError(f"{path}: required bundle support file is missing")
        files[relative] = path.read_bytes()
    for model, folder in selected:
        prefix = folder.resolve().relative_to(root).as_posix()
        manifest_path = (folder / "manifest.json").resolve()
        included = [manifest_path, *(folder / record.path for record in model.files)]
        for path in included:
            resolved = path.resolve()
            try:
                relative = resolved.relative_to(root).as_posix()
            except ValueError:
                raise ValueError(f"{path}: bundle input escapes the collection root") from None
            if not resolved.is_file():
                raise ValueError(f"{relative}: declared bundle file is missing")
            data = resolved.read_bytes()
            if resolved != manifest_path:
                record = next(item for item in model.files if (folder / item.path).resolve() == resolved)
                if len(data) != record.bytes or _sha256(data) != record.sha256:
                    raise ValueError(f"{relative}: declared size or SHA-256 does not match")
            files[_safe_name(relative)] = data
        if f"{prefix}/manifest.json" not in files:
            raise ValueError(f"{model.id}: manifest was not added to the bundle")

    file_records = [
        {"path": name, "bytes": len(files[name]), "sha256": _sha256(files[name])}
        for name in sorted(files)
    ]
    manifest = {
        "format": BUNDLE_FORMAT,
        "software_version": __version__,
        "selection": {"assay": assay, "ids": [model.id for model, _ in selected] if dataset_ids else []},
        "datasets": [model.id for model, _ in selected],
        "files": file_records,
    }
    manifest_data = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            archive.writestr(_zip_info("bundle.json"), manifest_data, compresslevel=9)
            for name in sorted(files):
                archive.writestr(_zip_info(name), files[name], compresslevel=9)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    result = verify_bundle(destination)
    result["path"] = str(destination)
    return result


def verify_bundle(path: Path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"{path}: bundle does not exist")
    try:
        with zipfile.ZipFile(path) as archive:
            infos = [info for info in archive.infolist() if not info.is_dir()]
            names = [_safe_name(info.filename) for info in infos]
            if len(names) != len(set(names)):
                raise ValueError(f"{path}: bundle contains duplicate paths")
            if names.count("bundle.json") != 1:
                raise ValueError(f"{path}: bundle.json is missing or duplicated")
            for info in infos:
                if ((info.external_attr >> 16) & 0o170000) == 0o120000:
                    raise ValueError(f"{path}: symbolic links are not permitted")
            try:
                manifest = json.loads(archive.read("bundle.json").decode("utf-8"))
            except (json.JSONDecodeError, UnicodeDecodeError) as error:
                raise ValueError(f"{path}: invalid bundle.json: {error}") from None
            if manifest.get("format") != BUNDLE_FORMAT:
                raise ValueError(f"{path}: unsupported bundle format")
            datasets = manifest.get("datasets")
            records = manifest.get("files")
            if not isinstance(datasets, list) or not datasets or len(datasets) != len(set(datasets)):
                raise ValueError(f"{path}: dataset list is empty or duplicated")
            if not isinstance(records, list):
                raise ValueError(f"{path}: file inventory is missing")
            declared = {}
            for record in records:
                if not isinstance(record, dict):
                    raise ValueError(f"{path}: invalid file inventory record")
                name = _safe_name(record.get("path", ""))
                if name == "bundle.json" or name in declared:
                    raise ValueError(f"{path}: duplicate or reserved inventory path {name}")
                if not isinstance(record.get("bytes"), int) or record["bytes"] < 0:
                    raise ValueError(f"{path}: invalid byte count for {name}")
                sha256 = record.get("sha256", "")
                if not isinstance(sha256, str) or len(sha256) != 64 or any(char not in "0123456789abcdef" for char in sha256):
                    raise ValueError(f"{path}: invalid SHA-256 for {name}")
                declared[name] = record
            expected = set(declared) | {"bundle.json"}
            if set(names) != expected:
                missing = sorted(expected - set(names))
                extra = sorted(set(names) - expected)
                raise ValueError(f"{path}: archive inventory differs: missing={missing}, extra={extra}")
            for name, record in declared.items():
                data = archive.read(name)
                if len(data) != record["bytes"] or _sha256(data) != record["sha256"]:
                    raise ValueError(f"{path}: size or SHA-256 differs for {name}")
            manifest_ids = []
            for name in sorted(declared):
                if name.startswith("datasets/") and name.endswith("/manifest.json"):
                    model = Dataset.model_validate_json(archive.read(name))
                    manifest_ids.append(model.id)
            if sorted(manifest_ids) != sorted(datasets):
                raise ValueError(f"{path}: dataset manifests differ from the dataset list")
    except zipfile.BadZipFile as error:
        raise ValueError(f"{path}: invalid ZIP archive: {error}") from None
    data = path.read_bytes()
    return {
        "status": "pass",
        "format": BUNDLE_FORMAT,
        "software_version": manifest.get("software_version"),
        "datasets": datasets,
        "dataset_count": len(datasets),
        "file_count": len(declared),
        "archive_bytes": len(data),
        "sha256": _sha256(data),
    }
