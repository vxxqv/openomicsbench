import json
from pathlib import Path
from .hashing import digest, contained
from .download import transfer
from .registry import lookup, registry

def get(root: Path, cache: Path, dataset_id: str, tier: str) -> Path:
    model, folder = lookup(root, dataset_id)
    if model.status == "original_link" or model.rights.status != "GREEN":
        raise ValueError(f"{dataset_id} is link-only ({model.rights.status}). Source: {model.source.url}")
    files = [f for f in model.files if f.tier == tier]
    if not files:
        raise ValueError(f"{dataset_id} has no {tier} tier. Run info to see available files.")
    manifest_hash = digest(folder / "manifest.json")
    destination = cache / model.id / model.release / manifest_hash[:16]
    marker = destination / f"{tier}.complete.json"
    marker.unlink(missing_ok=True)
    for f in files:
        local = contained(folder, f.path)
        source = local if local.exists() else f.url
        if source is None:
            raise ValueError(f"{dataset_id}/{f.path}: neither local file nor remote URL is available")
        transfer(source, contained(destination, f.path), f.sha256, f.bytes)
    marker.write_text(json.dumps({"id":model.id,"release":model.release,"manifest_sha256":manifest_hash,"tier":tier,"files":[f.path for f in files]}, indent=2) + "\n", encoding="utf-8", newline="\n")
    return destination

def verify_cache(root: Path, cache: Path) -> dict:
    records = registry(root)
    checked = 0
    for marker in cache.glob("*/*/*/*.complete.json"):
        receipt = json.loads(marker.read_text(encoding="utf-8"))
        model, folder = records[receipt["id"]]
        if receipt["release"] != model.release or receipt["manifest_sha256"] != digest(folder / "manifest.json"):
            raise ValueError(f"{marker}: cache receipt does not match the current manifest")
        expected = {f.path:f for f in model.files if f.tier == receipt["tier"]}
        if set(receipt["files"]) != set(expected) or len(receipt["files"]) != len(expected):
            raise ValueError(f"{marker}: incomplete cache receipt")
        for name, f in expected.items():
            p = contained(marker.parent, name)
            if not p.is_file() or p.stat().st_size != f.bytes or digest(p) != f.sha256:
                raise ValueError(f"{p}: cached file is missing or corrupt; retrieve the tier again")
            checked += 1
    return {"verified_files":checked,"partial_files":len(list(cache.glob("**/.partial-*")))}
