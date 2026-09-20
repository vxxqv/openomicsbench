"""Small deterministic amino-acid analyses without external services."""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from .sequences import read_records, validate_letters


RESIDUE_MASS = {
    "A":71.0788,"R":156.1875,"N":114.1038,"D":115.0886,"C":103.1388,"E":129.1155,"Q":128.1307,
    "G":57.0519,"H":137.1411,"I":113.1594,"L":113.1594,"K":128.1741,"M":131.1926,"F":147.1766,
    "P":97.1167,"S":87.0782,"T":101.1051,"W":186.2132,"Y":163.1760,"V":99.1326,"U":150.0388,"O":255.3134,
}
HYDROPATHY = {
    "I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"G":-0.4,"T":-0.7,"S":-0.8,
    "W":-0.9,"Y":-1.3,"P":-1.6,"H":-3.2,"E":-3.5,"Q":-3.5,"D":-3.5,"N":-3.5,"K":-3.9,"R":-4.5,
}
POSITIVE_PKA = {"N_term": 9.69, "H": 6.00, "K": 10.50, "R": 12.50}
NEGATIVE_PKA = {"C_term": 2.34, "C": 8.33, "D": 3.86, "E": 4.25, "Y": 10.07}


def net_charge(sequence: str, ph: float) -> float:
    if not 0 <= ph <= 14:
        raise ValueError("pH must be between 0 and 14")
    value = sequence.upper().replace("*", "")
    counts = Counter(value)
    positive = 1 / (1 + 10 ** (ph - POSITIVE_PKA["N_term"]))
    positive += sum(counts[residue] / (1 + 10 ** (ph - pka)) for residue, pka in POSITIVE_PKA.items() if residue != "N_term")
    negative = 1 / (1 + 10 ** (NEGATIVE_PKA["C_term"] - ph))
    negative += sum(counts[residue] / (1 + 10 ** (pka - ph)) for residue, pka in NEGATIVE_PKA.items() if residue != "C_term")
    return positive - negative


def isoelectric_point(sequence: str) -> float:
    low, high = 0.0, 14.0
    for _ in range(80):
        middle = (low + high) / 2
        if net_charge(sequence, middle) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def protein_record_stats(identifier: str, sequence: str, ph: float = 7.0) -> dict:
    validate_letters(sequence, "protein")
    value = sequence.upper().replace("*", "")
    counts = Counter(value)
    unknown = sorted(set(value) - set(RESIDUE_MASS))
    known_hydropathy = [HYDROPATHY[residue] for residue in value if residue in HYDROPATHY]
    mass = sum(RESIDUE_MASS[residue] for residue in value) + 18.01528 if value and not unknown else None
    return {
        "id": identifier,
        "length": len(value),
        "composition": dict(sorted(counts.items())),
        "unknown_mass_residues": unknown,
        "molecular_weight_da": mass,
        "net_charge": net_charge(value, ph),
        "charge_ph": ph,
        "estimated_isoelectric_point": isoelectric_point(value),
        "mean_hydropathy": sum(known_hydropathy) / len(known_hydropathy) if known_hydropathy else None,
        "aromatic_fraction": sum(counts[residue] for residue in "FWY") / len(value) if value else 0.0,
        "extinction_coefficient_reduced_280nm": counts["W"] * 5500 + counts["Y"] * 1490,
    }


def protein_report(path: Path, ph: float = 7.0) -> dict:
    records = list(read_records(path))
    if not records:
        raise ValueError("protein input contains no records")
    results = [protein_record_stats(record.identifier, record.sequence, ph) for record in records]
    return {"input": str(path), "records": len(results), "charge_ph": ph, "proteins": results}


def cleavage_sites(sequence: str, enzyme: str) -> list[int]:
    value = sequence.upper().replace("*", "")
    rules = {
        "trypsin": set("KR"),
        "lys-c": {"K"},
        "arg-c": {"R"},
        "chymotrypsin": set("FWYLM"),
    }
    if enzyme not in rules:
        raise ValueError(f"unknown enzyme: {enzyme}")
    sites = [0]
    for index, residue in enumerate(value[:-1], start=1):
        if residue in rules[enzyme] and value[index] != "P":
            sites.append(index)
    if value and value[-1] in rules[enzyme]:
        sites.append(len(value))
    elif sites[-1] != len(value):
        sites.append(len(value))
    return sorted(set(sites))


def digest_sequence(identifier: str, sequence: str, enzyme: str, missed_cleavages: int = 0, minimum_length: int = 1, maximum_length: int | None = None) -> list[dict]:
    validate_letters(sequence, "protein")
    if missed_cleavages < 0 or minimum_length < 1 or maximum_length is not None and maximum_length < minimum_length:
        raise ValueError("digestion length or missed-cleavage settings are invalid")
    value = sequence.upper().replace("*", "")
    sites = cleavage_sites(value, enzyme)
    peptides = []
    for left_index in range(len(sites) - 1):
        for missed in range(missed_cleavages + 1):
            right_index = left_index + missed + 1
            if right_index >= len(sites):
                break
            start, end = sites[left_index], sites[right_index]
            length = end - start
            if length < minimum_length or maximum_length is not None and length > maximum_length:
                continue
            peptides.append({"id": identifier, "start": start, "end": end, "sequence": value[start:end], "length": length, "missed_cleavages": missed})
    return peptides


def digest_report(path: Path, enzyme: str, missed_cleavages: int = 0, minimum_length: int = 1, maximum_length: int | None = None) -> dict:
    peptides = []
    records = 0
    for record in read_records(path):
        records += 1
        peptides.extend(digest_sequence(record.identifier, record.sequence, enzyme, missed_cleavages, minimum_length, maximum_length))
    return {"input": str(path), "enzyme": enzyme, "coordinate_system": "zero-based half-open", "records": records, "peptide_count": len(peptides), "peptides": peptides}
