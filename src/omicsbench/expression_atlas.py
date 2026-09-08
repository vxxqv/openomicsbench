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
import numpy as np

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


def load_counts(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    genes: list[str] = []
    matrix: list[list[int]] = []
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
            genes.append(row[0])
            matrix.append(counts)
    if not genes:
        raise ValueError("Raw-count table has no genes")
    if len(genes) != len(set(genes)):
        raise ValueError("Raw-count gene identifiers must be unique")
    values = np.asarray(matrix, dtype=np.int64)
    if np.any(values.sum(axis=0, dtype=np.float64) <= 0):
        raise ValueError("Raw-count table contains an empty sample library")
    return genes, samples, values


def _check_design(path: Path) -> tuple[int, list[str], list[str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "Run" not in reader.fieldnames or "Analysed" not in reader.fieldnames:
            raise ValueError("Experiment design needs Run and Analysed columns")
        runs = [row["Run"].strip() for row in reader if row["Analysed"].strip().lower() == "yes"]
    if len(runs) != len(set(runs)) or any(not run for run in runs):
        raise ValueError("Analysed run identifiers must be present and unique")
    return len(runs), runs, reader.fieldnames


def load_design_fields(path: Path, fields: list[str]) -> dict[str, dict[str, str]]:
    if not fields or len(fields) != len(set(fields)) or any(not field for field in fields):
        raise ValueError("Design fields must be present and unique")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = ["Run", "Analysed", *fields]
        missing = [field for field in required if not reader.fieldnames or field not in reader.fieldnames]
        if missing:
            raise ValueError(f"Experiment design is missing columns: {', '.join(missing)}")
        design: dict[str, dict[str, str]] = {}
        for row in reader:
            if row["Analysed"].strip().lower() != "yes":
                continue
            run = row["Run"].strip()
            values = {field: row[field].strip() for field in fields}
            if not run or any(not value for value in values.values()):
                raise ValueError("Analysed run identifiers and requested design values must be present")
            if run in design:
                raise ValueError("Analysed run identifiers must be unique")
            design[run] = values
    if len(design) < 4:
        raise ValueError("Experiment design needs at least four analysed runs")
    return design


def load_design(path: Path, factor: str) -> dict[str, str]:
    fields = load_design_fields(path, [factor])
    design = {run: values[factor] for run, values in fields.items()}
    if len(set(design.values())) != 2:
        raise ValueError("The diagnostic requires exactly two factor values")
    return design


def stage(accession: str, destination: Path, opener=urllib.request.urlopen, max_bytes: int = 1_000_000_000) -> dict:
    accession = normalize_accession(accession)
    resources = discover(accession, opener)
    destination.mkdir(parents=True, exist_ok=True)
    paths = {
        "raw_counts": destination / f"{accession}-raw-counts.tsv",
        "experiment_design": destination / f"{accession}-experiment-design.tsv",
    }
    records = [_stage_one(resources[role], paths[role], opener, max_bytes) for role in ("raw_counts", "experiment_design")]
    genes, count_samples, _ = load_counts(paths["raw_counts"])
    analysed, design_samples, columns = _check_design(paths["experiment_design"])
    missing = sorted(set(design_samples) - set(count_samples))
    if missing:
        raise ValueError("Raw-count samples are missing analysed experiment-design runs")
    excluded = [sample for sample in count_samples if sample not in set(design_samples)]
    return {
        "accession": accession,
        "genes": len(genes),
        "samples": len(design_samples),
        "raw_count_columns": len(count_samples),
        "excluded_count_columns": excluded,
        "design_columns": columns,
        "resources": records,
    }
