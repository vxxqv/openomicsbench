"""Screen Expression Atlas candidates without downloading count matrices."""
import argparse
import csv
import hashlib
import io
import json
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

from omicsbench.expression_atlas import discover, normalize_accession


def fetch(url: str, limit: int) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as response:
        if response.geturl() != url:
            raise ValueError(f"Unexpected redirect from {url}")
        body = response.read(limit + 1)
    if len(body) > limit:
        raise ValueError(f"Response exceeds {limit} bytes: {url}")
    return body


def attributes(rows: list[dict]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for row in rows:
        if isinstance(row, dict) and isinstance(row.get("name"), str) and isinstance(row.get("value"), str):
            result.setdefault(row["name"], []).append(row["value"])
    return result


def screen(accession: str) -> dict:
    accession = normalize_accession(accession)
    resources = discover(accession)
    raw = resources["raw_counts"]
    request = urllib.request.Request(raw.url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.geturl() != raw.url:
            raise ValueError(f"Unexpected redirect from {raw.url}")
        raw_bytes = int(response.headers["Content-Length"])
    design_resource = resources["experiment_design"]
    design_body = fetch(design_resource.url, 5_000_000)
    reader = csv.DictReader(io.StringIO(design_body.decode("utf-8")), delimiter="\t")
    if not reader.fieldnames or "Run" not in reader.fieldnames or "Analysed" not in reader.fieldnames:
        raise ValueError(f"{accession}: design lacks Run or Analysed")
    rows = [row for row in reader if row["Analysed"].strip().lower() == "yes"]
    runs = [row["Run"].strip() for row in rows]
    if len(runs) != len(set(runs)) or any(not run for run in runs):
        raise ValueError(f"{accession}: analysed run identifiers are missing or duplicated")
    factors = {}
    for column in reader.fieldnames:
        if column.startswith("Factor Value["):
            counts = Counter(row[column].strip() for row in rows)
            factors[column] = dict(sorted(counts.items()))
    study_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}"
    study = json.loads(fetch(study_url, 5_000_000))
    top = attributes(study.get("attributes", []))
    section = attributes(study.get("section", {}).get("attributes", []))
    return {
        "accession": accession,
        "title": (top.get("Title") or section.get("Title") or [None])[0],
        "release_date": (top.get("ReleaseDate") or [None])[0],
        "study_types": section.get("Study type", []),
        "organisms": section.get("Organism", []),
        "analysed_runs": len(runs),
        "factors": factors,
        "raw_counts": {"url": raw.url, "bytes": raw_bytes},
        "experiment_design": {
            "url": design_resource.url,
            "bytes": len(design_body),
            "sha256": hashlib.sha256(design_body).hexdigest(),
        },
        "biostudies_url": f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions", nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = []
    failures = []
    for value in args.accessions:
        try:
            results.append(screen(value))
        except (KeyError, OSError, ValueError, json.JSONDecodeError) as exc:
            failures.append({"accession": value, "error": str(exc)})
    payload = {"checked": date.today().isoformat(), "screened": results, "failures": failures}
    text = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")
