"""Verified transfers; incomplete files never become completed cache entries."""
import os
import tempfile
import urllib.request
from pathlib import Path
from .hashing import digest

def transfer(source: str | Path, target: Path, expected_hash: str, expected_bytes: int):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size == expected_bytes and digest(target) == expected_hash:
        return target
    fd, temp_name = tempfile.mkstemp(prefix=".partial-", dir=target.parent)
    temp = Path(temp_name)
    try:
        if isinstance(source, Path):
            stream = source.open("rb")
        else:
            if not source.startswith("https://"):
                raise ValueError("Remote downloads require HTTPS")
            stream = urllib.request.urlopen(source, timeout=60)
            if not stream.geturl().startswith("https://"):
                stream.close()
                raise ValueError("Download redirected away from HTTPS")
        with os.fdopen(fd, "wb") as out, stream:
            fd = None
            total = 0
            while chunk := stream.read(1024 * 1024):
                total += len(chunk)
                if total > expected_bytes:
                    raise ValueError(f"{target.name}: download exceeds declared size")
                out.write(chunk)
            out.flush()
            os.fsync(out.fileno())
        if temp.stat().st_size != expected_bytes or digest(temp) != expected_hash:
            raise ValueError(f"{target.name}: checksum or size mismatch; source was not cached")
        os.replace(temp, target)
        return target
    finally:
        if fd is not None:
            os.close(fd)
        temp.unlink(missing_ok=True)
