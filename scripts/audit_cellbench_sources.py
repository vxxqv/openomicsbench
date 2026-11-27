"""Audit pinned CellBench processed CSVs before proposing collection objects.

This verifies file structure and provenance inputs. It does not certify a
biological benchmark or establish permission to redistribute a derived file.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path


UPSTREAM_COMMIT = "7d00ffb0bd043467f365a9358d7312f294f13c9f"
STEMS = ("sc_10x", "sc_celseq2", "sc_dropseq", "sc_10x_5cl")
GENE_ID = re.compile(r"^ENSG\d+(?:\.\d+)?$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_pair(directory: Path, stem: str) -> dict:
    counts_path = directory / f"{stem}.count.csv.gz"
    metadata_path = directory / f"{stem}.metadata.csv.gz"
    if not counts_path.is_file() or not metadata_path.is_file():
        raise ValueError(f"{stem}: source pair is incomplete")

    with gzip.open(metadata_path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        if "cell_line_demuxlet" not in header:
            raise ValueError(f"{stem}: expected demuxlet reference labels")
        label_col = header.index("cell_line_demuxlet") + 1
        cell_ids: list[str] = []
        labels: Counter[str] = Counter()
        for row_number, row in enumerate(reader, 2):
            if len(row) != len(header) + 1 or not row[0] or not row[label_col]:
                raise ValueError(f"{stem}: malformed metadata row {row_number}")
            cell_ids.append(row[0])
            labels[row[label_col]] += 1
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError(f"{stem}: duplicate metadata cell identifiers")

    with gzip.open(counts_path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        count_cells = next(reader)
        if not count_cells or len(count_cells) != len(set(count_cells)):
            raise ValueError(f"{stem}: missing or duplicate count-matrix cells")
        if set(count_cells) != set(cell_ids):
            raise ValueError(f"{stem}: count-matrix and metadata cell IDs differ")
        gene_ids: set[str] = set()
        namespaces: Counter[str] = Counter()
        duplicated_genes: set[str] = set()
        feature_rows = 0
        nonzero = total = 0
        zero_genes = 0
        for row_number, row in enumerate(reader, 2):
            if len(row) != len(count_cells) + 1 or not row[0]:
                raise ValueError(f"{stem}: malformed gene row {row_number}")
            if row[0] in gene_ids:
                duplicated_genes.add(row[0])
            gene_ids.add(row[0])
            feature_rows += 1
            namespaces["Ensembl" if GENE_ID.fullmatch(row[0]) else
                       "ERCC" if row[0].startswith("ERCC-") else "other"] += 1
            row_total = 0
            for value in row[1:]:
                if not value.isascii() or not value.isdecimal():
                    raise ValueError(f"{stem}: invalid count at gene row {row_number}")
                number = int(value)
                row_total += number
                nonzero += number > 0
            total += row_total
            zero_genes += row_total == 0
    if not gene_ids or total == 0:
        raise ValueError(f"{stem}: empty count matrix")

    return {
        "stem": stem,
        "cells": len(count_cells),
        "features": feature_rows,
        "unique_feature_ids": len(gene_ids),
        "duplicated_feature_ids": len(duplicated_genes),
        "feature_namespaces": dict(sorted(namespaces.items())),
        "nonzero_entries": nonzero,
        "total_counts": total,
        "zero_count_genes": zero_genes,
        "cell_labels": dict(sorted(labels.items())),
        "matrix_metadata_cell_ids_match": True,
        "source_files": {
            name: {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "url": f"https://raw.githubusercontent.com/LuyiTian/sc_mixology/{UPSTREAM_COMMIT}/data/csv/{path.name}",
            }
            for name, path in (("counts", counts_path), ("metadata", metadata_path))
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "source_structure_checked",
        "checked": date.today().isoformat(),
        "upstream_repository": "https://github.com/LuyiTian/sc_mixology",
        "upstream_commit": UPSTREAM_COMMIT,
        "scope": "Count and metadata integrity only. Assembly, annotation, rights, and scoring remain separate release gates.",
        "cases": [audit_pair(args.directory, stem) for stem in STEMS],
    }
    payload = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload)


if __name__ == "__main__":
    main()
