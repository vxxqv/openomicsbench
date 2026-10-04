"""Exact small-variant comparison for the bundled synthetic DNA benchmark."""
from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

from .hashing import contained, digest
from .sequences import read_records


@dataclass(frozen=True, order=True, slots=True)
class Allele:
    contig: str
    position: int
    ref: str
    alt: str

    def label(self) -> str:
        return f"{self.contig}:{self.position}:{self.ref}>{self.alt}"


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix.lower() == ".gz" else path.open("r", encoding="utf-8", newline="")


def read_vcf(path: Path) -> tuple[set[Allele], dict]:
    """Read PASS and unfiltered SNV alleles from a VCF 4.x file.

    Multi-allelic records are expanded into individual alleles. Filtered records
    are counted but not assessed. Genotypes, phase and haplotype equivalence are
    deliberately outside this exact-match comparator.
    """
    path = Path(path)
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-2:] == [".vcf", ".gz"]:
        pass
    elif suffixes[-1:] != [".vcf"]:
        raise ValueError(f"{path}: expected a .vcf or .vcf.gz file")
    alleles: set[Allele] = set()
    fileformat = None
    columns = None
    records = assessed = filtered = 0
    with _open_text(path) as stream:
        for line_number, raw in enumerate(stream, start=1):
            line = raw.rstrip("\r\n")
            if line.startswith("##fileformat="):
                fileformat = line.split("=", 1)[1]
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                columns = line.split("\t")
                if columns[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                    raise ValueError(f"{path}:{line_number}: VCF columns are not in the required order")
                continue
            if not line:
                continue
            if line.startswith("#"):
                raise ValueError(f"{path}:{line_number}: unexpected VCF header line")
            if columns is None:
                raise ValueError(f"{path}:{line_number}: variant record appears before #CHROM header")
            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: VCF record has fewer than eight columns")
            records += 1
            contig, position_text, _, ref, alt_text, _, filter_value, _ = fields[:8]
            if not contig or any(character.isspace() for character in contig):
                raise ValueError(f"{path}:{line_number}: invalid contig")
            try:
                position = int(position_text)
            except ValueError:
                raise ValueError(f"{path}:{line_number}: POS must be a positive integer") from None
            if position < 1:
                raise ValueError(f"{path}:{line_number}: POS must be a positive integer")
            ref = ref.upper()
            alts = [value.upper() for value in alt_text.split(",")]
            if not ref or ref == "." or any(not alt or alt == "." for alt in alts):
                raise ValueError(f"{path}:{line_number}: REF and ALT must be present")
            if filter_value not in {"PASS", "."}:
                filtered += len(alts)
                continue
            for alt in alts:
                if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                    raise ValueError(f"{path}:{line_number}: this benchmark accepts single-nucleotide A/C/G/T substitutions only")
                allele = Allele(contig, position, ref, alt)
                if allele in alleles:
                    raise ValueError(f"{path}:{line_number}: duplicate allele {allele.label()}")
                alleles.add(allele)
                assessed += 1
    if fileformat is None or not fileformat.startswith("VCFv4."):
        raise ValueError(f"{path}: missing a supported ##fileformat=VCFv4.x declaration")
    if columns is None:
        raise ValueError(f"{path}: missing #CHROM header")
    return alleles, {
        "fileformat": fileformat,
        "records": records,
        "assessed_alleles": assessed,
        "filtered_alleles": filtered,
        "filter_policy": "FILTER is PASS or .",
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


def _truth(model, folder: Path) -> tuple[set[Allele], str, str, dict]:
    if model.assay != "whole_genome_dna_seq" or model.validation is None:
        raise ValueError(f"{model.id}: no small-variant comparison reference")
    profile_path = contained(folder, model.validation.profile)
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    truth_relative = profile.get("variant_truth")
    reference_relative = profile.get("reference")
    comparison = profile.get("variant_comparison")
    if not isinstance(truth_relative, str) or not isinstance(reference_relative, str) or not isinstance(comparison, dict):
        raise ValueError(f"{model.id}: variant comparison profile is incomplete")
    truth_path = _declared_file(model, folder, "effects", truth_relative)
    reference_path = _declared_file(model, folder, "reference", reference_relative)
    references = list(read_records(reference_path, "fasta"))
    if len(references) != 1:
        raise ValueError(f"{model.id}: variant benchmark requires one reference sequence")
    reference = references[0]
    truth: set[Allele] = set()
    with truth_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["position", "ref", "alt"]:
            raise ValueError(f"{model.id}: truth table columns must be position, ref and alt")
        for line_number, row in enumerate(reader, start=2):
            try:
                position = int(row["position"])
            except (TypeError, ValueError):
                raise ValueError(f"{model.id}: truth line {line_number} has an invalid position") from None
            ref = row["ref"].upper()
            alt = row["alt"].upper()
            allele = Allele(reference.identifier, position, ref, alt)
            if position < 1 or position > len(reference.sequence) or reference.sequence[position - 1] != ref:
                raise ValueError(f"{model.id}: truth allele {allele.label()} disagrees with the reference")
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                raise ValueError(f"{model.id}: truth contains a non-SNV allele")
            if allele in truth:
                raise ValueError(f"{model.id}: duplicate truth allele {allele.label()}")
            truth.add(allele)
    if len(truth) != profile.get("variant_count"):
        raise ValueError(f"{model.id}: truth allele count disagrees with the validation profile")
    return truth, reference.identifier, reference.sequence, comparison


def compare_variants(model, folder: Path, result_path: Path, detail_limit: int = 20) -> dict:
    if detail_limit < 0:
        raise ValueError("detail limit must be zero or greater")
    truth, reference_name, reference_sequence, comparison = _truth(model, folder)
    submitted, input_summary = read_vcf(Path(result_path))
    wrong_contigs = sorted({allele.contig for allele in submitted if allele.contig != reference_name})
    if wrong_contigs:
        raise ValueError(f"{model.id}: VCF contig must be {reference_name}; found {', '.join(wrong_contigs)}")
    for allele in submitted:
        if allele.position > len(reference_sequence) or reference_sequence[allele.position - 1] != allele.ref:
            raise ValueError(f"{model.id}: submitted allele {allele.label()} disagrees with the reference")

    true_positive = truth & submitted
    false_positive = submitted - truth
    false_negative = truth - submitted
    precision = len(true_positive) / len(submitted) if submitted else 0.0
    recall = len(true_positive) / len(truth) if truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    metrics = {"precision": precision, "recall": recall, "f1": f1}
    thresholds = {
        "precision": float(comparison["minimum_precision"]),
        "recall": float(comparison["minimum_recall"]),
    }
    checks = {name: metrics[name] >= minimum for name, minimum in thresholds.items()}
    reasons = [f"{name} is below {minimum}" for name, minimum in thresholds.items() if not checks[name]]
    return {
        "id": model.id,
        "status": "pass" if not reasons else "fail",
        "comparison": {
            "mode": "exact allele",
            "key": "contig, one-based position, REF and ALT",
            "scope": "single-nucleotide substitutions",
            "genotypes": "not assessed",
            "representation": "no haplotype reconciliation",
        },
        "reference": {"contig": reference_name, "truth_alleles": len(truth)},
        "input": str(Path(result_path)),
        "input_summary": input_summary,
        "counts": {
            "truth": len(truth),
            "submitted": len(submitted),
            "true_positive": len(true_positive),
            "false_positive": len(false_positive),
            "false_negative": len(false_negative),
        },
        "metrics": metrics,
        "thresholds": thresholds,
        "checks": checks,
        "false_positive_examples": [allele.label() for allele in sorted(false_positive)[:detail_limit]],
        "false_negative_examples": [allele.label() for allele in sorted(false_negative)[:detail_limit]],
        "examples_truncated": len(false_positive) > detail_limit or len(false_negative) > detail_limit,
        "reasons": reasons,
    }
