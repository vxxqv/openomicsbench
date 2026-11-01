"""Small, strictly scoped read-level taxonomy fixtures and comparisons."""

from __future__ import annotations

import gzip
import json
import os
import re
from collections import Counter
from itertools import chain, zip_longest
from pathlib import Path

from .hashing import contained
from .sequences import read_records, reverse_complement


REFERENCE_ID = re.compile(r"^mock_taxon_([1-9][0-9]*)$")


def _lines(path: Path):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        for number, raw in enumerate(stream, 1):
            line = raw.rstrip("\r\n")
            if not line:
                raise ValueError(f"{path}:{number}: empty assignment row")
            yield number, line


def read_assignments(path: Path) -> dict[str, int]:
    """Read canonical TSV or Kraken 2's five-column per-sequence format."""
    path = Path(path)
    rows = iter(_lines(path))
    try:
        number, first = next(rows)
    except StopIteration:
        raise ValueError(f"{path}: assignment file is empty") from None
    canonical = first == "read_id\ttaxon_id"
    if not canonical:
        rows = chain([(number, first)], rows)
    assignments = {}
    for number, line in rows:
        fields = line.split("\t")
        if canonical:
            if len(fields) != 2:
                raise ValueError(f"{path}:{number}: expected read_id and taxon_id")
            read_id, taxon_text = fields
        else:
            if len(fields) != 5 or fields[0] not in {"C", "U"}:
                raise ValueError(f"{path}:{number}: expected canonical TSV or five-column Kraken output")
            status, read_id, taxon_text, length, kmer_list = fields
            if not length or not kmer_list:
                raise ValueError(f"{path}:{number}: incomplete Kraken output")
            if (status == "U") != (taxon_text == "0"):
                raise ValueError(f"{path}:{number}: Kraken status and taxon code disagree")
        if not read_id or read_id.strip() != read_id or any(char.isspace() for char in read_id):
            raise ValueError(f"{path}:{number}: invalid read ID")
        if not taxon_text.isascii() or not taxon_text.isdecimal():
            raise ValueError(f"{path}:{number}: taxon code must be a non-negative integer")
        taxon_id = int(taxon_text)
        if read_id in assignments:
            raise ValueError(f"{path}:{number}: duplicate read ID {read_id}")
        assignments[read_id] = taxon_id
    if not assignments:
        raise ValueError(f"{path}: no assignments")
    return assignments


def _references(path: Path) -> dict[int, str]:
    references = {}
    for record in read_records(path, "fasta"):
        match = REFERENCE_ID.fullmatch(record.identifier)
        if not match:
            raise ValueError(f"{path}: reference ID must be mock_taxon_<positive integer>")
        code = int(match.group(1))
        if code in references:
            raise ValueError(f"{path}: duplicate taxon code {code}")
        if set(record.sequence) - set("ACGT"):
            raise ValueError(f"{path}: reference contains ambiguous bases")
        references[code] = record.sequence
    if not references:
        raise ValueError(f"{path}: no reference sequences")
    return references


def _pairs(read1: Path, read2: Path):
    seen = set()
    for left, right in zip_longest(read_records(read1, "fastq"), read_records(read2, "fastq")):
        if left is None or right is None:
            raise ValueError("paired FASTQ files have different record counts")
        if not left.identifier.endswith("/1") or not right.identifier.endswith("/2"):
            raise ValueError("paired FASTQ IDs must end in /1 and /2")
        identifier = left.identifier[:-2]
        if right.identifier[:-2] != identifier or identifier in seen:
            raise ValueError("paired FASTQ IDs are unordered or repeated")
        seen.add(identifier)
        yield identifier, left.sequence, right.sequence


def exact_assignments(reference: Path, read1: Path, read2: Path) -> dict[str, int]:
    """Assign only when both error-free mates match exactly one mock reference."""
    references = _references(Path(reference))
    assignments = {}
    for identifier, left, right in _pairs(Path(read1), Path(read2)):
        reverse = reverse_complement(right, "dna")
        matches = [code for code, sequence in references.items() if left in sequence and reverse in sequence]
        assignments[identifier] = matches[0] if len(matches) == 1 else 0
    return assignments


def classify_exact(reference: Path, read1: Path, read2: Path, output: Path, force: bool = False) -> dict:
    """Write a canonical TSV for a small mock-reference teaching case."""
    output = Path(output)
    if output.exists() and not force:
        raise ValueError(f"{output}: output already exists; use --force to replace it")
    assignments = exact_assignments(reference, read1, read2)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        temporary.write_text(
            "read_id\ttaxon_id\n" + "".join(f"{name}\t{code}\n" for name, code in assignments.items()),
            encoding="utf-8", newline="\n",
        )
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    counts = Counter(assignments.values())
    return {"operation": "exact_mock_taxonomy", "scope": "Exact, error-free matches to project-local mock references only", "output": str(output), "reads": len(assignments), "assigned": len(assignments) - counts[0], "unclassified": counts[0]}


def compare_taxonomy(model, folder: Path, submitted: Path, detail_limit: int = 20) -> dict:
    if model.assay != "shotgun_metagenomics" or model.validation is None:
        raise ValueError(f"{model.id}: no read-level taxonomy truth set")
    if detail_limit < 0:
        raise ValueError("detail limit must be zero or greater")
    profile = json.loads(contained(folder, model.validation.profile).read_text(encoding="utf-8"))
    truth_name = profile.get("read_truth")
    if truth_name not in {item.path for item in model.files}:
        raise ValueError(f"{model.id}: read truth is undeclared")
    truth = read_assignments(contained(folder, truth_name))
    calls = read_assignments(Path(submitted))
    missing = sorted(truth.keys() - calls.keys())
    unexpected = sorted(calls.keys() - truth.keys())
    shared = truth.keys() & calls.keys()
    codes = sorted(set(truth.values()) - {0})
    per_taxon = {}
    for code in codes:
        tp = sum(truth[name] == code and calls[name] == code for name in shared)
        fp = sum(calls[name] == code and truth[name] != code for name in shared)
        fn = sum(truth[name] == code and calls[name] != code for name in shared) + sum(truth[name] == code for name in missing)
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_taxon[str(code)] = {"support": sum(value == code for value in truth.values()), "true_positive": tp, "false_positive": fp, "false_negative": fn, "precision": precision, "recall": recall, "f1": f1}
    correct_classified = sum(truth[name] == calls[name] and truth[name] != 0 for name in shared)
    called_classified = sum(calls[name] != 0 for name in shared)
    truth_classified = sum(value != 0 for value in truth.values())
    truth_unclassified = sum(value == 0 for value in truth.values())
    correctly_unclassified = sum(truth[name] == 0 and calls[name] == 0 for name in shared)
    metrics = {
        "coverage": len(shared) / len(truth),
        "accuracy": sum(truth[name] == calls[name] for name in shared) / len(truth),
        "macro_f1": sum(item["f1"] for item in per_taxon.values()) / len(per_taxon) if per_taxon else 0.0,
        "classified_precision": correct_classified / called_classified if called_classified else 0.0,
        "classified_recall": correct_classified / truth_classified if truth_classified else 0.0,
        "unclassified_recall": correctly_unclassified / truth_unclassified if truth_unclassified else 1.0,
    }
    limits = profile.get("minimum_scores", {})
    reasons = []
    if missing:
        reasons.append(f"missing {len(missing)} read assignments")
    if unexpected:
        reasons.append(f"unexpected {len(unexpected)} read assignments")
    for name, minimum in limits.items():
        if name not in metrics or metrics[name] < minimum:
            reasons.append(f"{name}={metrics.get(name, 0):.5f} is below {minimum}")
    return {
        "id": model.id,
        "status": "fail" if reasons else "pass",
        "scope": "Read-level assignments to project-local mock taxon codes; not biological classification accuracy",
        "metrics": metrics,
        "per_taxon": per_taxon,
        "truth_counts": dict(sorted(Counter(truth.values()).items())),
        "missing": missing[:detail_limit],
        "unexpected": unexpected[:detail_limit],
        "reasons": reasons,
    }
