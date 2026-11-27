"""Evaluate processed protein-peptide tables without implying spectral validation."""

from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from .peptides import Peptide, read_peptides


def _open(path: Path):
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if path.suffix == '.gz' else path.open(encoding='utf-8', newline='')


def _mapping(path: Path | None) -> dict[str, set[str]]:
    if path is None:
        return {}
    mapping = defaultdict(set)
    with Path(path).open(encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames != ['source_id', 'reference_id']:
            raise ValueError('mapping must have source_id and reference_id columns')
        for number, row in enumerate(reader, 2):
            source, target = row['source_id'].strip(), row['reference_id'].strip()
            if not source or not target:
                raise ValueError(f'{path}:{number}: blank accession mapping')
            mapping[source].add(target)
    return dict(mapping)


def _quantities(path: Path) -> dict[Peptide, float] | None:
    with _open(Path(path)) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if 'quantity' not in (reader.fieldnames or []):
            return None
        if not {'protein_id', 'peptide'} <= set(reader.fieldnames or []):
            raise ValueError(f'{path}: quantity table needs protein_id and peptide')
        amounts = {}
        for number, row in enumerate(reader, 2):
            key = Peptide(row['protein_id'].strip(), row['peptide'].strip().upper())
            try:
                value = float(row['quantity'])
            except (ValueError, TypeError):
                raise ValueError(f'{path}:{number}: invalid quantity') from None
            if not math.isfinite(value) or value <= 0 or key in amounts:
                raise ValueError(f'{path}:{number}: quantity must be positive, finite and unique')
            amounts[key] = value
        return amounts


def _counts(truth: set[Peptide], calls: set[Peptide]) -> dict:
    matched = truth & calls
    precision = len(matched) / len(calls) if calls else 0.0
    recall = len(matched) / len(truth) if truth else 0.0
    truth_proteins = {item.protein_id for item in truth}
    matched_proteins = {item.protein_id for item in matched}
    return {'truth_pairs': len(truth), 'called_pairs': len(calls),
            'matched_pairs': len(matched), 'pair_precision': precision,
            'pair_recall': recall,
            'pair_f1': 2*precision*recall/(precision+recall) if precision+recall else 0.0,
            'truth_proteins': len(truth_proteins), 'matched_proteins': len(matched_proteins),
            'protein_recovery': len(matched_proteins)/len(truth_proteins) if truth_proteins else 0.0,
            'sequence_recovery': len({item.sequence for item in matched}) /
                                 len({item.sequence for item in truth}) if truth else 0.0}


def assess_peptides(truth_table: Path, calls_table: Path, *, mapping_table: Path | None = None,
                    q_values: tuple[float, ...] = (0.001, 0.01, 0.05)) -> dict:
    if not q_values or any(not math.isfinite(q) or not 0 <= q <= 1 for q in q_values):
        raise ValueError('q thresholds must be finite and between zero and one')
    aliases = _mapping(mapping_table)
    results = {}
    ambiguous = set()
    for threshold in q_values:
        truth, truth_input = read_peptides(truth_table, threshold)
        submitted, call_input = read_peptides(calls_table, threshold)
        mapped = set()
        ambiguous_at_threshold = set()
        for pair in submitted:
            destinations = aliases.get(pair.protein_id, {pair.protein_id})
            if len(destinations) != 1:
                ambiguous_at_threshold.add(pair.protein_id)
                continue
            mapped.add(Peptide(next(iter(destinations)), pair.sequence))
        ambiguous.update(ambiguous_at_threshold)
        if len(mapped) != len(submitted) - sum(pair.protein_id in ambiguous_at_threshold for pair in submitted):
            raise ValueError('accession mapping collapses distinct submitted identifications')
        results[str(threshold)] = {'identification': _counts(truth, mapped),
                                   'truth_input': truth_input, 'calls_input': call_input,
                                   'ambiguous_accessions_excluded': len(ambiguous_at_threshold)}
    truth_amounts, call_amounts = _quantities(truth_table), _quantities(calls_table)
    quantitative = None
    if truth_amounts is not None and call_amounts is not None:
        mapped_amounts = {}
        for pair, value in call_amounts.items():
            destinations = aliases.get(pair.protein_id, {pair.protein_id})
            if len(destinations) == 1:
                key = Peptide(next(iter(destinations)), pair.sequence)
                if key in mapped_amounts:
                    raise ValueError('accession mapping collapses quantities')
                mapped_amounts[key] = value
        shared = truth_amounts.keys() & mapped_amounts.keys()
        errors = sorted(abs(math.log2(mapped_amounts[key]/truth_amounts[key])) for key in shared)
        middle = len(errors)//2
        median = (errors[middle] if len(errors)%2 else (errors[middle-1]+errors[middle])/2) if errors else None
        quantitative = {'matched_pairs': len(shared),
                        'median_absolute_log2_ratio_error': median,
                        'scope': 'Positive processed quantities on exact mapped pairs; no normalization, missing-value imputation or protein-level quantification.'}
    return {'operation': 'peptide_assess', 'thresholds': results,
            'ambiguous_accessions': sorted(ambiguous), 'quantification': quantitative,
            'scope': 'Processed identification agreement and protein/sequence recovery. No raw-spectrum scoring, modification handling, FDR recalibration or protein inference.'}
