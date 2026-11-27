"""Optional AnnData adapter for explicit, reviewable cell-label tables."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def export_h5ad(source: Path, output: Path, *, label_column: str,
                batch_column: str | None = None, embedding_key: str | None = None) -> dict:
    try:
        import anndata
    except ImportError:
        raise ValueError('AnnData support is optional; install openomicsbench with anndata separately') from None
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError(f'{output}: output directory already exists')
    adata = anndata.read_h5ad(source, backed='r')
    try:
        ids = [str(item) for item in adata.obs_names]
        if not ids or len(ids) != len(set(ids)) or any(not item.strip() for item in ids):
            raise ValueError('AnnData observation IDs must be nonempty and unique')
        for column in (label_column, batch_column):
            if column is not None and column not in adata.obs:
                raise ValueError(f'AnnData obs column is missing: {column}')
        labels = [str(value) for value in adata.obs[label_column]]
        batches = [str(value) for value in adata.obs[batch_column]] if batch_column else None
        if any(not value.strip() or value.lower() in {'nan', 'none'} for value in labels):
            raise ValueError('AnnData labels must be present for every cell')
        if batches and any(not value.strip() or value.lower() in {'nan', 'none'} for value in batches):
            raise ValueError('AnnData batch values must be present for every cell')
        coordinates = None
        if embedding_key:
            if embedding_key not in adata.obsm:
                raise ValueError(f'AnnData embedding is missing: {embedding_key}')
            coordinates = np.asarray(adata.obsm[embedding_key], dtype=float)
            if coordinates.ndim != 2 or coordinates.shape[0] != len(ids) or coordinates.shape[1] < 2:
                raise ValueError('AnnData embedding must have one row per cell and at least two dimensions')
            if not np.isfinite(coordinates).all():
                raise ValueError('AnnData embedding coordinates must be finite')
        output.mkdir(parents=True)
        labels_path = output / 'labels.tsv'
        with labels_path.open('w', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(['cell_id', 'label'] + (['batch'] if batches else []))
            for index, cell in enumerate(ids):
                writer.writerow([cell, labels[index]] + ([batches[index]] if batches else []))
        embedding_path = None
        if coordinates is not None:
            embedding_path = output / 'embedding.tsv'
            with embedding_path.open('w', encoding='utf-8', newline='') as stream:
                writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
                writer.writerow(['cell_id'] + [f'dimension_{j+1}' for j in range(coordinates.shape[1])])
                for cell, row in zip(ids, coordinates):
                    writer.writerow([cell] + list(row))
        return {'operation': 'cell_export_h5ad', 'cells': len(ids),
                'labels': str(labels_path), 'embedding': str(embedding_path) if embedding_path else None,
                'scope': 'Exports declared AnnData obs labels and optional obsm coordinates; does not infer cell types or calculate a quality score.'}
    finally:
        if getattr(adata, 'file', None) is not None:
            adata.file.close()
