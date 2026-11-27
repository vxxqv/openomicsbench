"""Descriptive peak-set checks, separate from ENCODE-style biological validation."""

from __future__ import annotations

import json
import math
from pathlib import Path

from .intervals import interval_metrics, read_intervals


def _supplied_qc(path: Path) -> dict:
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(data, dict) or not data.get('source'):
        raise ValueError('QC JSON must name the source of the supplied measurements')
    for name in ('frip', 'tss_enrichment'):
        if name in data:
            value = data[name]
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
                raise ValueError(f'{name} must be a finite number')
            if name == 'frip' and not 0 <= value <= 1 or name == 'tss_enrichment' and value < 0:
                raise ValueError(f'{name} is outside its valid range')
    return data


def assess_intervals(peaks: Path, *, replicate: Path | None = None,
                     reference: Path | None = None, qc_json: Path | None = None,
                     idr_json: Path | None = None, reciprocal_overlap: float = 0.5) -> dict:
    if not 0 < reciprocal_overlap <= 1:
        raise ValueError('reciprocal overlap must be within (0, 1]')
    calls, _ = read_intervals(Path(peaks))
    widths = sorted(peak.end - peak.start for peak in calls)
    report = {'operation': 'interval_assess', 'peaks': str(peaks),
              'peak_count': len(calls),
              'width_bp': {'min': widths[0], 'median': widths[len(widths)//2], 'max': widths[-1]},
              'reciprocal_overlap': reciprocal_overlap}
    if replicate is not None:
        other, _ = read_intervals(Path(replicate))
        report['replicate_overlap'] = interval_metrics(other, calls, reciprocal_overlap)
        report['replicate_scope'] = 'Symmetric overlap of two supplied peak sets; not IDR or evidence of biological reproducibility by itself.'
    if reference is not None:
        truth, _ = read_intervals(Path(reference))
        report['reference_recovery'] = interval_metrics(truth, calls, reciprocal_overlap)
        bins = {}
        for name, predicate in [
            ('under_200bp', lambda width: width < 200),
            ('200_to_999bp', lambda width: 200 <= width < 1000),
            ('1000bp_or_more', lambda width: width >= 1000),
        ]:
            expected = [peak for peak in truth if predicate(peak.end - peak.start)]
            observed = [peak for peak in calls if predicate(peak.end - peak.start)]
            if expected:
                bins[name] = (interval_metrics(expected, observed, reciprocal_overlap) if observed else
                              {'truth_intervals': len(expected), 'called_intervals': 0,
                               'interval_precision': 0.0, 'interval_recall': 0.0, 'interval_f1': 0.0})
        report['width_sensitivity'] = bins
        report['reference_scope'] = 'Coordinate overlap with a user-supplied reference, not confirmation of biological peak validity.'
    if qc_json is not None:
        report['supplied_qc'] = _supplied_qc(qc_json)
    if idr_json is not None:
        idr = json.loads(Path(idr_json).read_text(encoding='utf-8'))
        if not isinstance(idr, dict) or not idr.get('source') or not isinstance(idr.get('idr_peaks'), int) or idr['idr_peaks'] < 0:
            raise ValueError('IDR JSON must contain source and nonnegative idr_peaks')
        report['supplied_idr'] = idr
        report['idr_scope'] = 'Imported external IDR summary; this command does not perform IDR analysis.'
    return report
