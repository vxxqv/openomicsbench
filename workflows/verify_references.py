"""Pin Ensembl references and verify pocket gene identifiers against release GTFs."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
import urllib.request
from pathlib import Path
from urllib.parse import urlsplit

from omicsbench.hashing import digest

ROOT = Path(__file__).resolve().parents[1]
GENE_ID = re.compile(r'(?:^|;\s*)gene_id "([^"]+)";')


def fetch(url: str, limit: int = 2_000_000) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as response:
        if response.geturl() != url:
            raise ValueError(f"Unexpected redirect from {url}")
        body = response.read(limit + 1)
    if len(body) > limit:
        raise ValueError(f"Response exceeds {limit} bytes: {url}")
    return body


def head(url: str) -> int:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.geturl() != url:
            raise ValueError(f"Unexpected redirect from {url}")
        return int(response.headers["Content-Length"])


def bsd_sum(path: Path) -> tuple[int, int]:
    checksum = 0
    size = 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) | ((checksum & 1) << 15))
                checksum = (checksum + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def download(url: str, target: Path, max_bytes: int = 200_000_000) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=".partial-", dir=target.parent)
    temporary = Path(temporary_name)
    total = 0
    sha256 = hashlib.sha256()
    try:
        with os.fdopen(fd, "wb") as output, urllib.request.urlopen(url, timeout=120) as response:
            fd = None
            if response.geturl() != url:
                raise ValueError(f"Unexpected redirect from {url}")
            while chunk := response.read(1024 * 1024):
                total += len(chunk)
                if total > max_bytes:
                    raise ValueError(f"Download exceeds {max_bytes} bytes: {url}")
                sha256.update(chunk)
                output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, target)
    finally:
        if fd is not None:
            os.close(fd)
        temporary.unlink(missing_ok=True)
    return {"bytes": total, "sha256": sha256.hexdigest()}


def verify_checksum_record(resource: dict) -> dict:
    body = fetch(resource["checksums_url"])
    filename = Path(urlsplit(resource["url"]).path).name
    line = f'{resource["bsd_sum"]:05d} {resource["blocks_1024"]:5d} {filename}'
    if line not in body.decode("ascii").splitlines():
        raise ValueError(f"Official CHECKSUMS does not contain the pinned record for {filename}")
    return {
        "url": resource["checksums_url"],
        "bytes": len(body),
        "sha256": hashlib.sha256(body).hexdigest(),
        "record": line,
    }


def gtf_genes(path: Path) -> set[str]:
    genes = set()
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = GENE_ID.search(fields[8])
            if not match:
                raise ValueError("A GTF gene row lacks gene_id")
            genes.add(match.group(1))
    if not genes:
        raise ValueError("The downloaded GTF contains no gene identifiers")
    return genes


def pocket_genes(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "gene_id":
            raise ValueError(f"Unexpected pocket header: {path}")
        genes = [row[0] for row in reader]
    if len(genes) != len(set(genes)):
        raise ValueError(f"Pocket gene identifiers are duplicated: {path}")
    return genes


def run(config_path: Path, staging: Path) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    profiles = []
    verified_accessions = []
    for profile in config["profiles"]:
        fasta_record = verify_checksum_record(profile["fasta"])
        fasta_bytes = head(profile["fasta"]["url"])
        gtf_record = verify_checksum_record(profile["gtf"])
        filename = Path(urlsplit(profile["gtf"]["url"]).path).name
        target = staging / profile["id"] / filename
        checksum, blocks = bsd_sum(target) if target.is_file() else (None, None)
        if (checksum, blocks) == (profile["gtf"]["bsd_sum"], profile["gtf"]["blocks_1024"]):
            transfer = {"bytes": target.stat().st_size, "sha256": digest(target)}
        else:
            transfer = download(profile["gtf"]["url"], target)
            checksum, blocks = bsd_sum(target)
        if (checksum, blocks) != (profile["gtf"]["bsd_sum"], profile["gtf"]["blocks_1024"]):
            raise ValueError(f"Downloaded GTF fails its official BSD checksum: {profile['id']}")
        annotation_genes = gtf_genes(target)
        checks = []
        for object_key in profile["objects"]:
            slug = object_key.lower()
            diagnostic_path = ROOT / f"staging/{slug}-diagnostic/diagnostic-report.json"
            diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
            size = diagnostic["smallest_passing_features"]
            genes = pocket_genes(ROOT / f"staging/{slug}-diagnostic/candidates/{size}/counts.tsv")
            missing = sorted(set(genes) - annotation_genes)
            checks.append({
                "object_key": object_key,
                "pocket_features": len(genes),
                "matched_features": len(genes) - len(missing),
                "unmatched_features": missing,
                "status": "pass" if not missing else "fail",
            })
            if not missing:
                verified_accessions.append(object_key)
        profiles.append({
            "id": profile["id"],
            "provider": profile["provider"],
            "release": profile["release"],
            "organism": profile["organism"],
            "assembly": profile["assembly"],
            "fasta": {
                "url": profile["fasta"]["url"],
                "bytes": fasta_bytes,
                "official_bsd_sum": profile["fasta"]["bsd_sum"],
                "official_blocks_1024": profile["fasta"]["blocks_1024"],
                "checksums_evidence": fasta_record,
                "downloaded": False,
            },
            "gtf": {
                "url": profile["gtf"]["url"],
                **transfer,
                "official_bsd_sum": checksum,
                "official_blocks_1024": blocks,
                "checksums_evidence": gtf_record,
                "gene_identifiers": len(annotation_genes),
                "gzip_test": "pass",
            },
            "pocket_checks": checks,
        })
    expected = sorted(object_key for profile in config["profiles"] for object_key in profile["objects"])
    report = {
        "schema_version": "1.0",
        "status": "pass" if sorted(verified_accessions) == expected else "fail",
        "checked": config["checked"],
        "verified_objects": sorted(verified_accessions),
        "excluded_sources": config.get("excluded_sources", []),
        "profiles": profiles,
        "config_sha256": digest(config_path),
        "workflow_sha256": digest(Path(__file__)),
    }
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "curation/reference-sources.json")
    parser.add_argument("--staging", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.config.resolve(), args.staging.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "pass" else 1)
