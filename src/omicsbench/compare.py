"""Compare submitted differential-expression effects with a bundled reference."""
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np

from .hashing import contained, digest
from .rnaseq import correlation, ranked


METRICS = ("spearman_logfc", "top_k_jaccard", "sign_concordance")


def read_effects(path: Path) -> dict[str, float]:
    suffixes = path.suffixes
    data_suffix = suffixes[-2] if suffixes and suffixes[-1] == ".gz" and len(suffixes) > 1 else path.suffix
    if data_suffix not in {".csv", ".tsv"}:
        raise ValueError(f"{path}: expected a .csv or .tsv result file")
    delimiter = "," if data_suffix == ".csv" else "\t"
    opener = gzip.open if path.suffix == ".gz" else open
    effects = {}
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        required = {"gene_id", "log2_fold_change"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: required columns are gene_id and log2_fold_change")
        for line, row in enumerate(reader, start=2):
            gene = (row.get("gene_id") or "").strip()
            if not gene:
                raise ValueError(f"{path}:{line}: gene_id is empty")
            if gene in effects:
                raise ValueError(f"{path}:{line}: duplicate gene_id {gene}")
            try:
                value = float(row["log2_fold_change"])
            except (TypeError, ValueError):
                raise ValueError(f"{path}:{line}: log2_fold_change is not numeric") from None
            if not math.isfinite(value):
                raise ValueError(f"{path}:{line}: log2_fold_change must be finite")
            effects[gene] = value
    if not effects:
        raise ValueError(f"{path}: result file contains no effects")
    return effects


def top_genes(effects: dict[str, float], k: int) -> set[str]:
    return set(sorted(effects, key=lambda gene: (-abs(effects[gene]), gene))[:k])


def compare(model, folder: Path, result_path: Path, detail_limit: int = 20) -> dict:
    if detail_limit < 0:
        raise ValueError("detail limit must be zero or greater")
    if model.kind != "real" or model.validation is None:
        raise ValueError(f"{model.id}: no differential-expression comparison reference")
    profile = json.loads(contained(folder, model.validation.profile).read_text(encoding="utf-8"))
    relative = profile.get("reference_effects")
    if not relative:
        raise ValueError(f"{model.id}: no differential-expression comparison reference")
    record = next((item for item in model.files if item.path == relative), None)
    reference_path = contained(folder, relative)
    if record is None or not reference_path.is_file() or reference_path.stat().st_size != record.bytes or digest(reference_path) != record.sha256:
        raise ValueError(f"{model.id}: comparison reference is missing or corrupt")

    reference = read_effects(reference_path)
    submitted = read_effects(Path(result_path))
    reference_genes = set(reference)
    submitted_genes = set(submitted)
    shared = sorted(reference_genes & submitted_genes)
    missing = sorted(reference_genes - submitted_genes)
    unexpected = sorted(submitted_genes - reference_genes)
    top_k = int(profile["top_k"])
    if len(shared) < top_k:
        raise ValueError(f"{model.id}: at least {top_k} shared finite effects are required")

    reference_shared = np.asarray([reference[gene] for gene in shared], dtype=float)
    submitted_shared = np.asarray([submitted[gene] for gene in shared], dtype=float)
    eligible = np.abs(reference_shared) >= 0.1
    if not np.any(eligible):
        raise ValueError(f"{model.id}: no shared reference effects meet the sign filter")
    reference_top = top_genes(reference, top_k)
    submitted_top = top_genes({gene: submitted[gene] for gene in shared}, top_k)
    metrics = {
        "spearman_logfc": correlation(ranked(reference_shared), ranked(submitted_shared)),
        "top_k_jaccard": len(reference_top & submitted_top) / len(reference_top | submitted_top),
        "sign_concordance": float(np.mean(np.sign(reference_shared[eligible]) == np.sign(submitted_shared[eligible]))),
    }
    thresholds = {metric.name: metric.minimum for metric in model.validation.metrics if metric.name in METRICS}
    checks = {name: metrics[name] >= thresholds[name] for name in METRICS}
    reasons = [f"{name} is below {thresholds[name]}" for name in METRICS if not checks[name]]
    if unexpected:
        reasons.append("submitted gene IDs are absent from the reference")
    evidence_path = contained(folder, profile["deseq2_evidence"])
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    return {
        "id": model.id,
        "status": "pass" if not reasons else "fail",
        "reference": {
            "method": "DESeq2 full-source analysis",
            "version": evidence["runtime"]["deseq2"],
            "contrast": evidence["design"]["contrast"],
        },
        "input": str(Path(result_path)),
        "genes": {
            "reference": len(reference),
            "submitted": len(submitted),
            "shared": len(shared),
            "coverage": len(shared) / len(reference),
            "missing": len(missing),
            "unexpected": len(unexpected),
            "missing_examples": missing[:detail_limit],
            "unexpected_examples": unexpected[:detail_limit],
            "examples_truncated": len(missing) > detail_limit or len(unexpected) > detail_limit,
        },
        "top_k": top_k,
        "metrics": metrics,
        "thresholds": thresholds,
        "checks": checks,
        "reasons": reasons,
    }
