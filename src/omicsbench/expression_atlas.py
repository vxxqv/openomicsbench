"""Discovery and staging for Expression Atlas bulk RNA-seq resources."""
import csv
import hashlib
import json
import os
import re
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urljoin, urlsplit

ORIGIN = "https://www.ebi.ac.uk"
ROOT = f"{ORIGIN}/gxa/"
ACCESSION = re.compile(r"^E-[A-Z0-9]+-[0-9]+$")
RESOURCE_TYPES = {
    "icon-raw-counts": "raw_counts",
    "icon-experiment-design": "experiment_design",
}


@dataclass(frozen=True)
class AtlasResource:
    role: str
    url: str
    description: str


def normalize_accession(value: str) -> str:
    accession = value.strip().upper()
    if not ACCESSION.fullmatch(accession):
        raise ValueError("Expression Atlas accession must look like E-MTAB-8572")
    return accession


def _official_url(accession: str, value: str) -> str:
    url = urljoin(ROOT, value)
    parsed = urlsplit(url)
    prefix = f"/gxa/experiments-content/{accession}/"
    if parsed.scheme != "https" or parsed.netloc != "www.ebi.ac.uk" or not parsed.path.startswith(prefix):
        raise ValueError("Expression Atlas returned a resource outside the requested experiment")
    if parsed.query or parsed.fragment or "/../" in parsed.path:
        raise ValueError("Expression Atlas returned an unsafe resource URL")
    return url


def parse_catalogue(accession: str, payload: bytes) -> dict[str, AtlasResource]:
    accession = normalize_accession(accession)
    try:
        rows = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("Expression Atlas resource catalogue is not valid JSON") from exc
    if not isinstance(rows, list):
        raise ValueError("Expression Atlas resource catalogue must be a list")
    found: dict[str, AtlasResource] = {}
    for row in rows:
        if not isinstance(row, dict) or row.get("type") not in RESOURCE_TYPES:
            continue
        role = RESOURCE_TYPES[row["type"]]
        if role in found:
            raise ValueError(f"Expression Atlas returned more than one {role} resource")
        if not isinstance(row.get("url"), str) or not isinstance(row.get("description"), str):
            raise ValueError(f"Expression Atlas {role} metadata is incomplete")
        found[role] = AtlasResource(role, _official_url(accession, row["url"]), row["description"].strip())
    missing = sorted(set(RESOURCE_TYPES.values()) - set(found))
    if missing:
        raise ValueError(f"Expression Atlas resource catalogue is missing: {', '.join(missing)}")
    return found


def discover(accession: str, opener=urllib.request.urlopen) -> dict[str, AtlasResource]:
    accession = normalize_accession(accession)
    url = f"{ROOT}json/experiments/{accession}/resources/DATA"
    with opener(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60) as response:
        if response.geturl() != url:
            raise ValueError("Expression Atlas catalogue redirected unexpectedly")
        payload = response.read(1_000_001)
    if len(payload) > 1_000_000:
        raise ValueError("Expression Atlas resource catalogue exceeds 1 MB")
    return parse_catalogue(accession, payload)


def _stage_one(resource: AtlasResource, target: Path, opener, max_bytes: int) -> dict:
    fd, temp_name = tempfile.mkstemp(prefix=".partial-", dir=target.parent)
    temporary = Path(temp_name)
    total = 0
    checksum = hashlib.sha256()
    try:
        with os.fdopen(fd, "wb") as output, opener(resource.url, timeout=120) as response:
            fd = None
            if response.geturl() != resource.url:
                raise ValueError(f"{resource.role} redirected unexpectedly")
            while chunk := response.read(1024 * 1024):
                total += len(chunk)
                if total > max_bytes:
                    raise ValueError(f"{resource.role} exceeds the staging limit")
                checksum.update(chunk)
                output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, target)
        return {"role": resource.role, "path": target.name, "url": resource.url, "bytes": total, "sha256": checksum.hexdigest()}
    finally:
        if fd is not None:
            os.close(fd)
        temporary.unlink(missing_ok=True)


def _check_counts(path: Path) -> tuple[int, list[str]]:
    genes = 0
    with path.open(encoding="utf-8", newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        header = next(rows, None)
        if not header or header[:2] != ["Gene ID", "Gene Name"] or len(header) < 4:
            raise ValueError("Raw counts need Gene ID, Gene Name and at least two samples")
        samples = header[2:]
        if len(samples) != len(set(samples)) or any(not sample for sample in samples):
            raise ValueError("Raw-count sample identifiers must be present and unique")
        for number, row in enumerate(rows, start=2):
            if len(row) != len(header) or not row[0]:
                raise ValueError(f"Raw-count row {number} has the wrong shape")
            try:
                counts = [int(value) for value in row[2:]]
            except ValueError as exc:
                raise ValueError(f"Raw-count row {number} contains a non-integer value") from exc
            if any(value < 0 for value in counts):
                raise ValueError(f"Raw-count row {number} contains a negative value")
            genes += 1
    if genes == 0:
        raise ValueError("Raw-count table has no genes")
    return genes, samples


def _check_design(path: Path) -> tuple[int, list[str], list[str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "Run" not in reader.fieldnames or "Analysed" not in reader.fieldnames:
            raise ValueError("Experiment design needs Run and Analysed columns")
        runs = [row["Run"].strip() for row in reader if row["Analysed"].strip().lower() == "yes"]
    if len(runs) != len(set(runs)) or any(not run for run in runs):
        raise ValueError("Analysed run identifiers must be present and unique")
    return len(runs), runs, reader.fieldnames


def stage(accession: str, destination: Path, opener=urllib.request.urlopen, max_bytes: int = 1_000_000_000) -> dict:
    accession = normalize_accession(accession)
    resources = discover(accession, opener)
    destination.mkdir(parents=True, exist_ok=True)
    paths = {
        "raw_counts": destination / f"{accession}-raw-counts.tsv",
        "experiment_design": destination / f"{accession}-experiment-design.tsv",
    }
    records = [_stage_one(resources[role], paths[role], opener, max_bytes) for role in ("raw_counts", "experiment_design")]
    genes, count_samples = _check_counts(paths["raw_counts"])
    analysed, design_samples, columns = _check_design(paths["experiment_design"])
    if set(count_samples) != set(design_samples):
        raise ValueError("Raw-count samples do not match analysed experiment-design runs")
    return {
        "accession": accession,
        "genes": genes,
        "samples": len(count_samples),
        "design_columns": columns,
        "resources": records,
    }
