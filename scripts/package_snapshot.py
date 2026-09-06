"""Write a deterministic inventory for the development handoff."""
import json
from pathlib import Path
from omicsbench.hashing import digest

root = Path(__file__).resolve().parents[1]
destination = root / "release/development-file-manifest.json"
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
    "snapshot": "1.0.0.dev1",
    "created": "2026-09-06",
    "scope": "Project-authored development handoff; staged third-party source data excluded.",
    "self_excluded": "release/development-file-manifest.json",
    "file_count": len(files),
    "files": files,
}
destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
print(f"Inventoried {len(files)} files.")
