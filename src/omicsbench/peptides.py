"""Exact peptide-identification comparison for synthetic protein fixtures."""
from __future__ import annotations

import csv
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path

from .hashing import contained, digest
from .sequences import read_records


CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True, order=True, slots=True)
class Peptide:
    protein_id: str
    sequence: str

    def label(self) -> str:
        return f"{self.protein_id}:{self.sequence}"


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix.lower() == ".gz" else path.open("r", encoding="utf-8", newline="")


def _suffix(path: Path) -> str:
    suffixes = [value.lower() for value in path.suffixes]
    if suffixes[-2:] == [".tsv", ".gz"]:
        return ".tsv.gz"
    if suffixes[-1:] == [".tsv"]:
        return ".tsv"
    raise ValueError(f"{path}: expected a .tsv or .tsv.gz peptide table")


def read_peptides(path: Path, maximum_q_value: float, columns: dict[str, str] | None = None) -> tuple[set[Peptide], dict]:
    """Read unique protein-peptide pairs, optionally filtering by q-value."""
    path = Path(path)
    _suffix(path)
    if not math.isfinite(maximum_q_value) or not 0 <= maximum_q_value <= 1:
        raise ValueError("maximum q-value must be between zero and one")
    names = {"protein_id": "protein_id", "peptide": "peptide", "q_value": "q_value", **(columns or {})}
    if set(names) != {"protein_id", "peptide", "q_value"} or names["protein_id"] == names["peptide"]:
        raise ValueError("peptide columns must map protein_id, peptide and q_value")
    peptides: set[Peptide] = set()
    rows = filtered = 0
    with _open_text(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if names["protein_id"] not in fields or names["peptide"] not in fields:
            raise ValueError(f"{path}: peptide table requires {names['protein_id']} and {names['peptide']} columns")
        q_values = names["q_value"] in fields
        for line_number, row in enumerate(reader, start=2):
            rows += 1
            protein_id = (row.get(names["protein_id"]) or "").strip()
            sequence = (row.get(names["peptide"]) or "").strip().upper()
            if not protein_id or any(character.isspace() for character in protein_id):
                raise ValueError(f"{path}:{line_number}: protein_id must be a non-empty token")
            if not sequence or set(sequence) - CANONICAL_AMINO_ACIDS:
                raise ValueError(f"{path}:{line_number}: peptide must contain canonical amino-acid letters only")
            if q_values:
                try:
                    q_value = float(row.get(names["q_value"]) or "")
                except ValueError:
                    raise ValueError(f"{path}:{line_number}: q_value must be numeric") from None
                if not math.isfinite(q_value) or not 0 <= q_value <= 1:
                    raise ValueError(f"{path}:{line_number}: q_value must be finite and between zero and one")
                if q_value > maximum_q_value:
                    filtered += 1
                    continue
            peptide = Peptide(protein_id, sequence)
            if peptide in peptides:
                raise ValueError(f"{path}:{line_number}: duplicate identification {peptide.label()}")
            peptides.add(peptide)
    return peptides, {
        "rows": rows,
        "assessed_identifications": len(peptides),
        "filtered_identifications": filtered,
        "q_value_column": q_values,
        "columns": names,
        "maximum_q_value": maximum_q_value,
    }


def _declared_file(model, folder: Path, role: str, path: str | None = None) -> Path:
    records = [record for record in model.files if record.role == role and (path is None or record.path == path)]
    if len(records) != 1:
        raise ValueError(f"{model.id}: expected one declared {role} file")
    record = records[0]
    target = contained(folder, record.path)
    if not target.is_file() or target.stat().st_size != record.bytes or digest(target) != record.sha256:
        raise ValueError(f"{model.id}: {record.path} is missing or corrupt")
    return target


def _truth(model, folder: Path) -> tuple[set[Peptide], dict[Peptide, tuple[int, int]], dict[str, str], dict]:
    if model.assay != "sequence_protein" or model.validation is None:
        raise ValueError(f"{model.id}: no peptide-identification comparison reference")
    profile_path = contained(folder, model.validation.profile)
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    truth_relative = profile.get("peptide_truth")
    reference_relative = profile.get("protein_reference")
    comparison = profile.get("peptide_comparison")
    if not isinstance(truth_relative, str) or not isinstance(reference_relative, str) or not isinstance(comparison, dict):
        raise ValueError(f"{model.id}: peptide comparison profile is incomplete")
    truth_path = _declared_file(model, folder, "effects", truth_relative)
    reference_path = _declared_file(model, folder, "sequence", reference_relative)
    proteins = {record.identifier: record.sequence.upper().replace("*", "") for record in read_records(reference_path, "fasta")}
    if not proteins:
        raise ValueError(f"{model.id}: protein reference contains no records")
    truth: set[Peptide] = set()
    coordinates: dict[Peptide, tuple[int, int]] = {}
    with truth_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        expected = ["protein_id", "peptide", "start", "end", "missed_cleavages"]
        if reader.fieldnames != expected:
            raise ValueError(f"{model.id}: peptide truth columns must be {', '.join(expected)}")
        for line_number, row in enumerate(reader, start=2):
            protein_id = row["protein_id"]
            sequence = row["peptide"].upper()
            try:
                start, end, missed = int(row["start"]), int(row["end"]), int(row["missed_cleavages"])
            except ValueError:
                raise ValueError(f"{model.id}: peptide truth line {line_number} has invalid coordinates") from None
            if protein_id not in proteins or start < 0 or end <= start or end > len(proteins[protein_id]):
                raise ValueError(f"{model.id}: peptide truth line {line_number} has invalid protein coordinates")
            if proteins[protein_id][start:end] != sequence:
                raise ValueError(f"{model.id}: peptide truth line {line_number} disagrees with the protein reference")
            if missed < 0 or set(sequence) - CANONICAL_AMINO_ACIDS:
                raise ValueError(f"{model.id}: peptide truth line {line_number} is invalid")
            peptide = Peptide(protein_id, sequence)
            if peptide in truth:
                raise ValueError(f"{model.id}: duplicate truth identification {peptide.label()}")
            truth.add(peptide)
            coordinates[peptide] = (start, end)
    if len(truth) != profile.get("peptide_count"):
        raise ValueError(f"{model.id}: peptide truth count disagrees with the validation profile")
    return truth, coordinates, proteins, comparison


def _covered_residues(peptides: set[Peptide], coordinates: dict[Peptide, tuple[int, int]]) -> set[tuple[str, int]]:
    covered = set()
    for peptide in peptides:
        start, end = coordinates[peptide]
        covered.update((peptide.protein_id, position) for position in range(start, end))
    return covered


def compare_peptides(model, folder: Path, result_path: Path, detail_limit: int = 20, columns: dict[str, str] | None = None) -> dict:
    if detail_limit < 0:
        raise ValueError("detail limit must be zero or greater")
    truth, coordinates, proteins, comparison = _truth(model, folder)
    maximum_q_value = float(comparison.get("maximum_q_value", 0.01))
    submitted, input_summary = read_peptides(Path(result_path), maximum_q_value, columns)
    true_positive = truth & submitted
    false_positive = submitted - truth
    false_negative = truth - submitted
    precision = len(true_positive) / len(submitted) if submitted else 0.0
    recall = len(true_positive) / len(truth) if truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    truth_proteins = {item.protein_id for item in truth}
    matched_proteins = {item.protein_id for item in true_positive}
    protein_recall = len(matched_proteins) / len(truth_proteins) if truth_proteins else 0.0
    truth_residues = _covered_residues(truth, coordinates)
    matched_residues = _covered_residues(true_positive, coordinates)
    residue_coverage = len(matched_residues) / len(truth_residues) if truth_residues else 0.0
    metrics = {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "protein_recall": protein_recall,
        "residue_coverage": residue_coverage,
    }
    thresholds = {
        name: float(comparison[f"minimum_{name}"])
        for name in ("precision", "recall", "protein_recall", "residue_coverage")
    }
    checks = {name: metrics[name] >= minimum for name, minimum in thresholds.items()}
    reasons = [f"{name} is below {minimum}" for name, minimum in thresholds.items() if not checks[name]]
    return {
        "id": model.id,
        "status": "pass" if not reasons else "fail",
        "comparison": {
            "mode": "exact protein-peptide pair",
            "key": "protein_id and unmodified peptide sequence",
            "scope": "project-authored in-silico tryptic digest",
            "q_value_policy": f"q_value <= {maximum_q_value} when the column is present",
            "outside_scope": ["spectra", "modifications", "charge states", "protein inference", "FDR calibration"],
        },
        "reference": {"proteins": len(proteins), "truth_identifications": len(truth)},
        "input": str(Path(result_path)),
        "input_summary": input_summary,
        "counts": {
            "truth": len(truth),
            "submitted": len(submitted),
            "true_positive": len(true_positive),
            "false_positive": len(false_positive),
            "false_negative": len(false_negative),
            "truth_proteins": len(truth_proteins),
            "matched_proteins": len(matched_proteins),
        },
        "metrics": metrics,
        "thresholds": thresholds,
        "checks": checks,
        "false_positive_examples": [item.label() for item in sorted(false_positive)[:detail_limit]],
        "false_negative_examples": [item.label() for item in sorted(false_negative)[:detail_limit]],
        "examples_truncated": len(false_positive) > detail_limit or len(false_negative) > detail_limit,
        "reasons": reasons,
    }
