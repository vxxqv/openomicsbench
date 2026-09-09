"""Write the deterministic file inventory for a release candidate."""
import json
import hashlib
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
metadata = json.loads((root / "release/metadata-input.json").read_text(encoding="utf-8"))
destination = root / "release/file-manifest.json"
excluded = {destination.resolve()}
files = []
git = os.environ.get("OPENOMICSBENCH_GIT", "git")
listed = subprocess.run(
    [git, "-c", f"safe.directory={root.as_posix()}", "ls-files", "-z"],
    cwd=root, check=True, capture_output=True,
).stdout.decode("utf-8").split("\0")
for relative in sorted(value for value in listed if value):
    path = root / relative
    if not path.is_file() or path.resolve() in excluded:
        continue
    blob = subprocess.run(
        [git, "-c", f"safe.directory={root.as_posix()}", "show", f":{relative}"],
        cwd=root, check=True, capture_output=True,
    ).stdout
    files.append({
        "path": path.relative_to(root).as_posix(),
        "bytes": len(blob),
        "sha256": hashlib.sha256(blob).hexdigest(),
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
