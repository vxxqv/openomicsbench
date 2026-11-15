"""Describe native differential-expression results without treating agreement as truth."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from .compare import read_effects
from .hashing import digest
from .rnaseq_design import _table
from .rnaseq_pack import EFFECT_HEADERS, GENE_HEADERS, _design

ADJUSTED_HEADERS = ("padj", "FDR", "adj.P.Val", "qvalue", "q_value", "adjusted_p_value")


def _column(headers: list[str], explicit: str | None, aliases: tuple[str, ...], label: str) -> str:
    if explicit is not None:
        if explicit not in headers:
            raise ValueError(f"{label} column is absent: {explicit}")
        return explicit
    found = [name for name in aliases if name in headers]
    if len(found) != 1:
        raise ValueError(f"{label} column is {'ambiguous' if found else 'missing'}; specify it explicitly")
    return found[0]


def _number(raw: str, path: Path, line: int, label: str) -> float | None:
    if raw.strip() in {"", "NA", "NaN", "nan", "NULL"}:
        return None
    try:
        value = float(raw)
    except ValueError:
        raise ValueError(f"{path}:{line}: {label} is not numeric") from None
    if not math.isfinite(value):
        raise ValueError(f"{path}:{line}: {label} must be finite or missing")
    return value


def assess_de(benchmark, source: str | Path, *, alpha: float = 0.05,
              gene_column: str | None = None, effect_column: str | None = None,
              adjusted_column: str | None = None, reverse_sign: bool = False) -> dict:
    """Report significance calls and effect agreement with the archived reference."""
    if not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be strictly between zero and one")
    source = Path(source)
    stream, delimiter = _table(source)
    with stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        headers = reader.fieldnames
        if not headers or len(headers) != len(set(headers)):
            raise ValueError("result table has an empty or duplicate header")
        gene_name = _column(headers, gene_column, GENE_HEADERS, "gene")
        effect_name = _column(headers, effect_column, EFFECT_HEADERS, "effect")
        adjusted_name = _column(headers, adjusted_column, ADJUSTED_HEADERS, "adjusted p-value")
        if len({gene_name, effect_name, adjusted_name}) != 3:
            raise ValueError("gene, effect and adjusted p-value columns must differ")
        seen = set()
        effects = {}
        tested = significant = missing_adjusted = missing_effect = 0
        significant_positive = significant_negative = 0
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{source}:{line}: inconsistent row width")
            gene = row[gene_name].strip()
            if not gene or gene in seen:
                raise ValueError(f"{source}:{line}: empty or duplicate gene ID")
            seen.add(gene)
            effect = _number(row[effect_name], source, line, "effect")
            adjusted = _number(row[adjusted_name], source, line, "adjusted p-value")
            if adjusted is not None and not 0 <= adjusted <= 1:
                raise ValueError(f"{source}:{line}: adjusted p-value must be between zero and one")
            if effect is None:
                missing_effect += 1
            else:
                effects[gene] = -effect if reverse_sign else effect
            if adjusted is None:
                missing_adjusted += 1
            else:
                tested += 1
                if adjusted <= alpha:
                    significant += 1
                    if effect is not None:
                        significant_positive += int(effects[gene] > 0)
                        significant_negative += int(effects[gene] < 0)
    if not seen:
        raise ValueError("result table has no genes")
    reference = read_effects(benchmark.file(role="effects", tier="expected"))
    shared = set(effects) & set(reference)
    nonzero = [gene for gene in shared if effects[gene] != 0 and reference[gene] != 0]
    concordant = sum((effects[gene] > 0) == (reference[gene] > 0) for gene in nonzero)
    design = _design(benchmark)
    return {
        "operation": "assess_de", "benchmark_id": benchmark.id,
        "input": str(source), "input_sha256": digest(source),
        "columns": {"gene_id": gene_name, "log2_fold_change": effect_name,
                    "adjusted_p_value": adjusted_name},
        "alpha": alpha, "rows": len(seen), "tested_rows": tested,
        "missing_adjusted_p_values": missing_adjusted,
        "missing_effects": missing_effect,
        "significant_rows": significant,
        "significant_fraction_of_tested": significant / tested if tested else None,
        "significant_positive_effects": significant_positive,
        "significant_negative_effects": significant_negative,
        "significant_missing_or_zero_effects": significant - significant_positive - significant_negative,
        "reference_effects": len(reference), "shared_finite_effects": len(shared),
        "shared_nonzero_effects": len(nonzero),
        "effect_direction_agreement": concordant / len(nonzero) if nonzero else None,
        "contrast": {"target": design["target"], "reference": design["reference"]},
        "sign_reversed": reverse_sign,
        "scope": "Counts calls at the supplied adjusted-p threshold and compares effect direction with one archived DESeq2 analysis. No independent differential-expression truth set or false-discovery calibration is implied.",
    }
