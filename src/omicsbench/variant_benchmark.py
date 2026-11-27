"""Region-limited exact representation comparison of small VCF alleles.

This is an input audit and baseline, not haplotype-aware variant benchmarking.
"""

from __future__ import annotations

import gzip
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_regions(path: Path) -> dict[str, list[tuple[int, int]]]:
    regions = defaultdict(list)
    with _open(Path(path)) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected BED3 or more columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise ValueError(f"{path}:{line_number}: invalid BED coordinates") from None
            if not fields[0] or start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            regions[fields[0]].append((start, end))
    if not regions:
        raise ValueError(f"{path}: no benchmark intervals")
    merged = {}
    for chrom, values in regions.items():
        output = []
        for start, end in sorted(values):
            if output and start <= output[-1][1]:
                output[-1] = (output[-1][0], max(end, output[-1][1]))
            else:
                output.append((start, end))
        merged[chrom] = output
    return merged


def read_alleles(path: Path) -> tuple[set[tuple[str, int, str, str]], dict]:
    alleles = set()
    format_seen = columns_seen = False
    filtered = 0
    with _open(Path(path)) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith('##fileformat=VCFv4.'):
                format_seen = True
            if line.startswith('##') or not line.strip():
                continue
            if line.startswith('#CHROM\t'):
                columns_seen = True
                continue
            if not columns_seen:
                raise ValueError(f"{path}:{line_number}: variant before VCF columns")
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: fewer than eight VCF columns")
            chrom, pos_text, _, ref, alternate, _, filter_value, _ = fields[:8]
            try:
                pos = int(pos_text)
            except ValueError:
                raise ValueError(f"{path}:{line_number}: invalid VCF position") from None
            if pos < 1 or not chrom or not ref or not alternate:
                raise ValueError(f"{path}:{line_number}: invalid VCF allele")
            alts = alternate.upper().split(',')
            if filter_value not in {'PASS', '.'}:
                filtered += len(alts)
                continue
            for alt in alts:
                ref_upper = ref.upper()
                if not alt or not set(ref_upper + alt) <= set('ACGT') or ref_upper == alt:
                    raise ValueError(f"{path}:{line_number}: only literal A/C/G/T alleles are supported")
                allele = (chrom, pos, ref_upper, alt)
                if allele in alleles:
                    raise ValueError(f"{path}:{line_number}: duplicate allele")
                alleles.add(allele)
    if not format_seen or not columns_seen:
        raise ValueError(f"{path}: missing VCF 4.x header")
    return alleles, {'passed_alleles': len(alleles), 'filtered_alleles': filtered}


def _inside(allele, regions):
    chrom, pos, ref, _ = allele
    start, end = pos - 1, pos - 1 + len(ref)
    values = regions.get(chrom, [])
    index = bisect_right(values, (start, float('inf'))) - 1
    return index >= 0 and end <= values[index][1]


def _metrics(truth: set, calls: set) -> dict:
    found = len(truth & calls)
    precision = found / len(calls) if calls else 0.0
    recall = found / len(truth) if truth else 0.0
    return {'truth': len(truth), 'calls': len(calls), 'true_positive': found,
            'false_positive': len(calls - truth), 'false_negative': len(truth - calls),
            'precision': precision, 'recall': recall,
            'f1': 2 * precision * recall / (precision + recall) if precision + recall else 0.0}


def benchmark_variants(truth_vcf: Path, calls_vcf: Path, regions_bed: Path,
                       *, assembly: str, truth_assembly: str) -> dict:
    if not assembly or not truth_assembly or assembly != truth_assembly:
        raise ValueError('declared call and truth reference assemblies must match')
    regions = read_regions(regions_bed)
    truth, truth_input = read_alleles(truth_vcf)
    calls, calls_input = read_alleles(calls_vcf)
    in_truth = {allele for allele in truth if _inside(allele, regions)}
    in_calls = {allele for allele in calls if _inside(allele, regions)}
    if not in_truth:
        raise ValueError('no truth alleles fall fully inside benchmark regions')
    by_type = {}
    for name, predicate in [('SNV', lambda a: len(a[2]) == len(a[3]) == 1),
                            ('indel', lambda a: len(a[2]) != len(a[3])),
                            ('MNV', lambda a: len(a[2]) == len(a[3]) > 1)]:
        by_type[name] = _metrics({a for a in in_truth if predicate(a)},
                                 {a for a in in_calls if predicate(a)})
    return {'operation': 'variant_benchmark', 'assembly': assembly,
            'truth_input': truth_input, 'calls_input': calls_input,
            'regions': sum(map(len, regions.values())),
            'outside_regions': {'truth': len(truth - in_truth), 'calls': len(calls - in_calls)},
            'overall': _metrics(in_truth, in_calls), 'by_type': by_type,
            'scope': 'Exact CHROM/POS/REF/ALT matching after full REF-span containment in supplied BED regions. No left normalization, haplotype equivalence, genotype comparison or confidence-region generation; use a representation-aware comparator for those tasks.'}
