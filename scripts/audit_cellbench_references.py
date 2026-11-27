"""Check CellBench feature identifiers against versioned source-era annotations."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

from audit_cellbench_sources import STEMS, sha256


GENE_ID = re.compile(r'gene_id "([^\"]+)"')
GENE_NAME = re.compile(r'gene_name "([^\"]+)"')


def annotation(path: Path) -> tuple[set[str], set[str]]:
    ids: set[str] = set()
    names: set[str] = set()
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.split('\t', 8)
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            identifier = GENE_ID.search(fields[8])
            name = GENE_NAME.search(fields[8])
            if identifier:
                ids.add(identifier.group(1).split('.')[0])
            if name:
                names.add(name.group(1))
    if not ids:
        raise ValueError(f'{path}: no gene records found')
    return ids, names


def features(path: Path) -> set[str]:
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream)
        next(reader)
        return {row[0] for row in reader}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sources', type=Path)
    parser.add_argument('--ensembl88', required=True, type=Path)
    parser.add_argument('--gencode26', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    ensembl_ids, ensembl_names = annotation(args.ensembl88)
    gencode_ids, gencode_names = annotation(args.gencode26)
    records = []
    for stem in STEMS:
        source = args.sources / f'{stem}.count.csv.gz'
        values = features(source)
        expected = ensembl_names if stem == 'sc_10x_5cl' else ensembl_ids
        secondary = gencode_names if stem == 'sc_10x_5cl' else gencode_ids
        controls = {value for value in values if value.startswith('ERCC-')}
        test = values - controls
        missing = sorted(test - expected)
        records.append({
            'stem': stem, 'source_sha256': sha256(source), 'feature_count': len(values),
            'ERCC_controls_excluded': len(controls), 'annotation_matched': len(test) - len(missing),
            'annotation_unmatched': len(missing), 'unmatched_examples': missing[:25],
            'secondary_gencode26_matched': len(test & secondary),
            'candidate_annotation': 'Ensembl GRCh38 release 88 (GRCh38.p10)',
        })
    report = {
        'scope': 'Identifier overlap check only; this does not prove the original processing pipeline used these exact GTF files.',
        'annotations': {
            'ensembl88': {'sha256': sha256(args.ensembl88), 'gene_ids': len(ensembl_ids), 'gene_names': len(ensembl_names)},
            'gencode26': {'sha256': sha256(args.gencode26), 'gene_ids': len(gencode_ids), 'gene_names': len(gencode_names)},
        },
        'cases': records,
    }
    payload = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding='utf-8', newline='\n')
    print(payload)


if __name__ == '__main__':
    main()
