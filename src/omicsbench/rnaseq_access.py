"""Checked access to the biological bulk RNA-seq benchmarks."""

from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .hashing import digest
from .rnaseq import read_counts


@dataclass(frozen=True)
class RNASeqData:
    """Feature-by-sample integer counts with ordered sample records."""

    id: str
    tier: str
    genes: tuple[str, ...]
    samples: tuple[str, ...]
    counts: np.ndarray
    metadata: tuple[dict, ...]
    source: dict
    limitations: tuple[str, ...]


def load_rnaseq(benchmark, tier: str = "pocket") -> RNASeqData:
    if benchmark.assay != "bulk_rna_seq" or benchmark.manifest["status"] != "validated" or benchmark.manifest["kind"] != "real":
        raise ValueError(f"{benchmark.id}: only validated biological bulk RNA-seq objects are supported")
    if tier not in {"pocket", "full"}:
        raise ValueError("tier must be 'pocket' or 'full'")
    role, file_tier = ("raw_counts", "pocket") if tier == "pocket" else ("baseline", "expected")
    genes, samples, counts = read_counts(benchmark.file(role=role, tier=file_tier))
    metadata = benchmark.manifest["samples"]
    if samples != [item["sample_id"] for item in metadata]:
        raise ValueError(f"{benchmark.id}: count columns and declared sample order differ")
    counts.setflags(write=False)
    return RNASeqData(benchmark.id, tier, tuple(genes), tuple(samples), counts,
                      tuple(dict(item) for item in metadata), dict(benchmark.manifest["source"]),
                      tuple(benchmark.manifest["limitations"]))


def export_rnaseq(benchmark, destination: str | Path, tier: str = "pocket") -> dict:
    """Write aligned counts, sample table and source details into a new directory."""
    data = load_rnaseq(benchmark, tier)
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f"{destination}: export directory already exists")
    if not destination.parent.is_dir():
        raise FileNotFoundError(f"{destination.parent}: parent directory does not exist")
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}-", dir=destination.parent))
    try:
        with (temporary / "counts.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", *data.samples])
            writer.writerows((gene, *row) for gene, row in zip(data.genes, data.counts.tolist()))
        fields = list(data.metadata[0])
        with (temporary / "samples.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(data.metadata)
        files = {name: {"bytes": (temporary / name).stat().st_size, "sha256": digest(temporary / name)}
                 for name in ("counts.tsv", "samples.tsv")}
        details = {"benchmark_id": data.id, "tier": data.tier, "orientation": "genes_by_samples",
                   "features": len(data.genes), "samples": len(data.samples), "source": data.source,
                   "limitations": list(data.limitations), "files": files}
        (temporary / "export.json").write_text(json.dumps(details, indent=2) + "\n", encoding="utf-8", newline="\n")
        os.replace(temporary, destination)
        return details
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
