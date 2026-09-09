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
release_commit = metadata.get("release_commit")
command = [git, "-c", f"safe.directory={root.as_posix()}"]
command += ["ls-tree", "-r", "--name-only", "-z", release_commit] if release_commit else ["ls-files", "-z"]
listed = subprocess.run(command, cwd=root, check=True, capture_output=True).stdout.decode("utf-8").split("\0")
for relative in sorted(value for value in listed if value):
    path = root / relative
    if path.resolve() in excluded:
        continue
    reference = release_commit or ""
    object_name = f"{reference}:{relative}" if reference else f":{relative}"
    blob = subprocess.run(
        [git, "-c", f"safe.directory={root.as_posix()}", "show", object_name],
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
