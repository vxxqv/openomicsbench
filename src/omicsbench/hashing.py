from hashlib import sha256
from pathlib import Path

def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return sha256_stream(stream)

def sha256_stream(stream) -> str:
    h = sha256()
    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
        h.update(chunk)
    return h.hexdigest()

def contained(root: Path, relative: str) -> Path:
    from .models import safe_relative
    safe_relative(relative)
    root = root.resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"Path escapes dataset directory: {relative}")
    return path
