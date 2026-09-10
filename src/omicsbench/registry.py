import difflib
import json
import os
from pathlib import Path
from .models import Dataset

def _contains_collection(root: Path) -> bool:
    return next((root / "datasets").glob("**/manifest.json"), None) is not None

def default_root() -> Path:
    configured = os.environ.get("OMICSBENCH_ROOT")
    if configured:
        return Path(configured)
    current = Path.cwd()
    if _contains_collection(current):
        return current
    packaged = Path(__file__).resolve().parent / "_collection"
    if _contains_collection(packaged):
        return packaged
    checkout = Path(__file__).resolve().parents[2]
    if _contains_collection(checkout):
        return checkout
    return packaged

def registry(root: Path) -> dict[str, tuple[Dataset, Path]]:
    result = {}
    for path in sorted((root / "datasets").glob("**/manifest.json")):
        model = Dataset.model_validate_json(path.read_text(encoding="utf-8"))
        if model.id in result:
            raise ValueError(f"Duplicate dataset ID: {model.id}")
        if path.parent.name != model.id:
            raise ValueError(f"{path}: folder name must match {model.id}")
        result[model.id] = (model, path.parent)
    return result

def lookup(root: Path, dataset_id: str):
    records = registry(root)
    if dataset_id not in records:
        near = difflib.get_close_matches(dataset_id, records, n=3)
        raise ValueError(f"Unknown dataset {dataset_id}. Try: {', '.join(near) or 'omicsbench list'}")
    return records[dataset_id]

def write_registry(root: Path, destination: Path):
    destination.write_text(json.dumps([m.model_dump(mode="json") for m, _ in registry(root).values()], indent=2) + "\n", encoding="utf-8", newline="\n")
