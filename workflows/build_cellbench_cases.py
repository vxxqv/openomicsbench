"""Build compact, source-checked CellBench label-recovery cases.

Use only with the pinned source files and GTF recorded by the v15 intake audit.
The cases report agreement with source cell-line labels; they do not establish
general cell-type truth or make the seven CellBench files independent studies.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import random
import re
import shutil
from datetime import date
from pathlib import Path

from omicsbench.hashing import digest
from omicsbench.singlecell import summarize_matrix


ROOT = Path(__file__).resolve().parents[1]
GEO = {
    'sc_10x': ('singlecell-002', 'GSM3022245', 'Three-line 10x Chromium cell mixture'),
    'sc_celseq2': ('singlecell-003', 'GSM3336845', 'Three-line CEL-seq2 cell mixture'),
    'sc_dropseq': ('singlecell-004', 'GSM3336849', 'Three-line Drop-seq cell mixture'),
    'sc_celseq2_5cl_p1': ('singlecell-005', 'GSM3618022', 'Five-line CEL-seq2 mixture, plate 1'),
    'sc_celseq2_5cl_p2': ('singlecell-006', 'GSM3618023', 'Five-line CEL-seq2 mixture, plate 2'),
    'sc_celseq2_5cl_p3': ('singlecell-007', 'GSM3618024', 'Five-line CEL-seq2 mixture, plate 3'),
}
GENE_ID = re.compile(r'gene_id "([^\"]+)"')
GENE_NAME = re.compile(r'gene_name "([^\"]+)"')


def read_gtf(path: Path) -> dict[str, str]:
    names = {}
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
                names[identifier.group(1).split('.')[0]] = name.group(1) if name else identifier.group(1)
    if not names:
        raise ValueError('GTF contains no gene records')
    return names


def source_labels(path: Path) -> dict[str, str]:
    labels = {}
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream)
        header = next(reader)
        label_column = header.index('cell_line_demuxlet') + 1
        for number, row in enumerate(reader, 2):
            if len(row) != len(header) + 1 or not row[0] or not row[label_column] or row[0] in labels:
                raise ValueError(f'{path}:{number}: malformed cell-line label')
            labels[row[0]] = row[label_column]
    return labels


def record(path: Path, folder: Path, tier: str, role: str, media: str) -> dict:
    return {'path': path.relative_to(folder).as_posix(), 'tier': tier, 'role': role,
            'media_type': media, 'bytes': path.stat().st_size, 'sha256': digest(path)}


def build(stem: str, sources: Path, names: dict[str, str], source_hashes: dict,
          reference_hash: str, workflow_commit: str, license_path: Path) -> dict:
    dataset_id, geo, title = GEO[stem]
    count_path = sources / f'{stem}.count.csv.gz'
    metadata_path = sources / f'{stem}.metadata.csv.gz'
    for key, path in [('counts', count_path), ('metadata', metadata_path)]:
        if digest(path) != source_hashes[key]['sha256']:
            raise ValueError(f'{stem}: {key} source hash differs from intake audit')
    labels = source_labels(metadata_path)
    with gzip.open(count_path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream)
        all_cells = next(reader)
    if set(all_cells) != set(labels):
        raise ValueError(f'{stem}: source count/metadata cell IDs differ')
    if stem == 'sc_10x':
        chosen = set(random.Random(15).sample(all_cells, 300))
        cells = [cell for cell in all_cells if cell in chosen]
    else:
        cells = all_cells
    selected_cells = set(cells)
    columns = [index for index, cell in enumerate(all_cells, 1) if cell in selected_cells]
    target = ROOT / 'datasets' / 'singlecell' / dataset_id
    if target.exists():
        raise ValueError(f'{target}: refusing to overwrite existing case')
    (target / 'pocket').mkdir(parents=True)
    (target / 'expected').mkdir()
    (target / 'provenance').mkdir()
    matrix_path = target / 'pocket' / 'matrix.mtx.gz'
    barcode_path = target / 'pocket' / 'barcodes.tsv'
    feature_path = target / 'pocket' / 'features.tsv'
    truth_path = target / 'expected' / 'cell-labels.tsv'
    license_target = target / 'provenance' / 'source-license.txt'
    shutil.copyfile(license_path, license_target)
    barcode_path.write_text(''.join(f'{cell}\n' for cell in cells), encoding='utf-8', newline='\n')
    truth_path.write_text('cell_id\tlabel\n' + ''.join(f'{cell}\t{labels[cell]}\n' for cell in cells), encoding='utf-8', newline='\n')

    kept_features: list[tuple[str, str]] = []
    source_nonzero = source_total = 0
    removed = []
    with gzip.open(count_path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream)
        next(reader)
        for number, row in enumerate(reader, 2):
            gene = row[0]
            if gene not in names and not gene.startswith('ERCC-'):
                removed.append(gene)
                continue
            kept_features.append((gene, names.get(gene, gene)))
            for column in columns:
                value = int(row[column])
                if value < 0:
                    raise ValueError(f'{stem}:{number}: negative count')
                source_nonzero += value > 0
                source_total += value
    if not kept_features or len(set(gene for gene, _ in kept_features)) != len(kept_features):
        raise ValueError(f'{stem}: duplicate or empty feature set')
    if len(removed) > 25:
        raise ValueError(f'{stem}: more than 25 features absent from the declared annotation')
    feature_path.write_text(''.join(f'{gene}\t{name}\tGene Expression\n' for gene, name in kept_features),
                            encoding='utf-8', newline='\n')

    with matrix_path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='\n') as output:
                output.write('%%MatrixMarket matrix coordinate integer general\n')
                output.write(f'{len(kept_features)} {len(cells)} {source_nonzero}\n')
                row_index = 0
                with gzip.open(count_path, 'rt', encoding='utf-8', newline='') as stream:
                    reader = csv.reader(stream)
                    next(reader)
                    for row in reader:
                        if row[0] not in names and not row[0].startswith('ERCC-'):
                            continue
                        row_index += 1
                        for cell_index, source_column in enumerate(columns, 1):
                            value = int(row[source_column])
                            if value:
                                output.write(f'{row_index} {cell_index} {value}\n')
    summary = summarize_matrix(matrix_path, barcode_path, feature_path, 'MT-', 3)
    if summary['nonzero'] != source_nonzero or summary['total_counts'] != source_total:
        raise ValueError(f'{stem}: transformed matrix differs from selected source counts')

    profile_path = target / 'expected' / 'validation.json'
    profile = {
        'profile': 'cellbench-cell-line-reference-v1',
        'matrix': 'pocket/matrix.mtx.gz', 'barcodes': 'pocket/barcodes.tsv',
        'features': 'pocket/features.tsv', 'label_truth': 'expected/cell-labels.tsv',
        'expected_summary': summary,
        'clustering_comparison': {
            'minimum_adjusted_rand_index': -1.0, 'minimum_normalized_mutual_information': 0.0,
            'minimum_purity': 0.0, 'minimum_coverage': 1.0,
        },
        'threshold_rationale': 'Report label recovery without an arbitrary performance floor; require complete cell coverage.',
    }
    profile_path.write_text(json.dumps(profile, indent=2) + '\n', encoding='utf-8', newline='\n')
    provenance_path = target / 'provenance' / 'source-and-transform.json'
    provenance = {
        'dataset_id': dataset_id, 'source_stem': stem, 'source_commit': source_hashes['commit'],
        'source_count_sha256': source_hashes['counts']['sha256'],
        'source_metadata_sha256': source_hashes['metadata']['sha256'],
        'reference_gtf_sha256': reference_hash, 'workflow_commit': workflow_commit,
        'cell_selection': 'All source cells' if stem != 'sc_10x' else '300 cells sampled without labels using Python Random(15), retained in source order',
        'selected_cells': len(cells), 'selected_features': len(kept_features),
        'excluded_unmatched_feature_ids': removed,
        'label_origin': 'The upstream cell_line_demuxlet metadata column; SNP-based reference labels, not direct cell-type truth.',
        'source_nonzero_selected': source_nonzero, 'source_total_counts_selected': source_total,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8', newline='\n')

    files = [
        record(matrix_path, target, 'pocket', 'matrix', 'application/gzip'),
        record(barcode_path, target, 'pocket', 'barcodes', 'text/tab-separated-values'),
        record(feature_path, target, 'pocket', 'features', 'text/tab-separated-values'),
        record(truth_path, target, 'expected', 'cell_metadata', 'text/tab-separated-values'),
        record(profile_path, target, 'expected', 'metrics', 'application/json'),
        record(provenance_path, target, 'metadata', 'provenance', 'application/json'),
        record(license_target, target, 'metadata', 'license', 'text/plain'),
    ]
    today = date.today().isoformat()
    model = {
        'schema_version': '3.0', 'id': dataset_id, 'title': title, 'release': '15.0.0',
        'assay': 'single_cell_rna_seq', 'kind': 'real', 'status': 'validated',
        'archetype': 'Known cell-line mixture measured by ' + stem.replace('_', ' '),
        'organism': 'Homo sapiens', 'taxon_id': 9606,
        'source': {'repository': 'CellBench / GEO', 'accession': geo,
                   'url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo}',
                   'citation': 'Tian L et al. Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments. Nature Methods 2019. doi:10.1038/s41592-019-0425-8.',
                   'retrieved': today, 'sha256': source_hashes['counts']['sha256']},
        'rights': {'status': 'GREEN', 'license': 'MIT',
                   'evidence': f'Count and metadata CSVs are in LuyiTian/sc_mixology at commit {source_hashes["commit"]}; the repository LICENSE grants redistribution with notice, copied into this case. Original source authors generated the data; derived pocket is attributed and source hashes are pinned.',
                   'checked': today},
        'reference': {'genome': 'GRCh38.p10', 'annotation': 'Ensembl GRCh38 release 88 (feature-ID compatibility check)',
                      'namespace': 'Ensembl gene IDs plus ERCC controls', 'compatibility': 'verified'},
        'samples': [{'sample_id': stem, 'condition': 'cell-line mixture', 'replicate': 1,
                     'batch': stem, 'strandedness': 'unknown'}],
        'derivation': {'workflow': 'workflows/build_cellbench_cases.py', 'commit': workflow_commit,
                       'seed': 15 if stem == 'sc_10x' else 0, 'algorithm': 'source-count-preserving subset and Matrix Market conversion',
                       'parameters': {'cells': len(cells), 'features': len(kept_features), 'excluded_features': len(removed)},
                       'versions': {'python': '3.12', 'openomicsbench_contract': '3.0', 'ensembl_gtf': '88'}},
        'files': files,
        'validation': {'profile': 'expected/validation.json', 'baseline_version': 'cellbench-cell-line-reference-v1', 'metrics': []},
        'sequence': None, 'single_cell': {'format': 'matrix_market', 'orientation': 'features_by_cells',
                                          'feature_columns': 3, 'mitochondrial_prefix': 'MT-'},
        'intervals': None,
        'limitations': [
            'Cell-line labels come from upstream demuxlet metadata, not direct cell-type measurements.',
            'All six cases are from one source study; CEL-seq2 five-line plates are technical repeats, not independent cohorts.',
            'The pocket preserves selected cell counts and all features matching the declared Ensembl annotation; excluded unmatched IDs are recorded in provenance.',
            'Gene-ID overlap confirms feature compatibility with Ensembl 88 but does not prove the original aligner used this exact downloaded GTF file.',
        ],
    }
    (target / 'manifest.json').write_text(json.dumps(model, indent=2) + '\n', encoding='utf-8', newline='\n')
    return {'id': dataset_id, 'cells': len(cells), 'features': len(kept_features), 'excluded_features': len(removed), 'summary': summary}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_directory', type=Path)
    parser.add_argument('--gtf', required=True, type=Path)
    parser.add_argument('--workflow-commit', required=True)
    args = parser.parse_args()
    if not re.fullmatch(r'[0-9a-f]{40}', args.workflow_commit):
        raise ValueError('workflow commit must be a full SHA-1')
    audited = json.loads((ROOT / 'curation/v15-cellbench-source-audit.json').read_text(encoding='utf-8'))
    source_audit = {case['stem']: {**case['source_files'], 'commit': audited['upstream_commit']} for case in audited['cases']}
    reference_hash = digest(args.gtf)
    expected_reference = json.loads((ROOT / 'curation/v15-cellbench-reference-audit.json').read_text(encoding='utf-8'))
    if reference_hash != expected_reference['annotations']['ensembl88']['sha256']:
        raise ValueError('GTF checksum differs from reference audit')
    names = read_gtf(args.gtf)
    license_path = ROOT / 'workflows/cellbench-source-license.txt'
    if not license_path.is_file():
        raise ValueError('pinned source license is missing')
    for stem in GEO:
        print(json.dumps(build(stem, args.source_directory, names, source_audit[stem], reference_hash,
                               args.workflow_commit, license_path)))


if __name__ == '__main__':
    main()
