"""Write the deterministic file inventory for a release candidate."""
import json
from pathlib import Path
from omicsbench.hashing import digest

root = Path(__file__).resolve().parents[1]
metadata = json.loads((root / "release/metadata-input.json").read_text(encoding="utf-8"))
destination = root / "release/file-manifest.json"
excluded = {destination.resolve()}
excluded_directories = {".git", ".snakemake", "__pycache__", "build", "staging"}
files = []
for path in sorted(root.rglob("*")):
    if not path.is_file() or path.resolve() in excluded or any(part in excluded_directories for part in path.parts) or path.suffix == ".pyc":
        continue
    files.append({
        "path": path.relative_to(root).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": digest(path),
    })
payload = {
    "version": metadata["software_version"],
    "created": metadata["candidate_date"],
    "scope": "Files in the reviewed version 1 release candidate; generated caches and staged source downloads are excluded.",
    "self_excluded": "release/file-manifest.json",
    "file_count": len(files),
    "files": files,
}
destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
print(f"Inventoried {len(files)} files.")
