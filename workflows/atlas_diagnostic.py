"""Evaluate deterministic reductions of one staged Expression Atlas experiment."""
import argparse
import csv
import json
import platform
from collections import Counter
from pathlib import Path
import numpy as np

from omicsbench.expression_atlas import load_counts, load_design, stage
from omicsbench.hashing import digest
from omicsbench.models import Sample
from omicsbench.rnaseq import contrast_selection_order, diagnostic, preservation

ROOT = Path(__file__).resolve().parents[1]


def write_counts(path: Path, genes: list[str], samples: list[str], counts: np.ndarray):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])
        writer.writerows([[gene, *map(int, row)] for gene, row in zip(genes, counts, strict=True)])


def run(config_path: Path, staging: Path) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    source = staging / "source"
    transfer = stage(config["accession"], source)
    genes, sample_ids, counts = load_counts(source / f'{config["accession"]}-raw-counts.tsv')
    design = load_design(source / f'{config["accession"]}-experiment-design.tsv', config["factor_column"])
    if set(sample_ids) != set(design):
        raise ValueError("Count columns and analysed design runs differ")
    seen: Counter[str] = Counter()
    samples = []
    for sample_id in sample_ids:
        condition = design[sample_id]
        seen[condition] += 1
        samples.append(Sample(sample_id=sample_id, condition=condition, replicate=seen[condition], batch="not_reported", strandedness="unknown"))
    full = diagnostic(counts, samples, normalization=config["normalization"])
    order = contrast_selection_order(genes, counts, full["effects"], config["seed"], config["effect_anchors"])
    curve = []
    selected = None
    for size in config["candidate_sizes"]:
        indices = sorted(order[:size])
        names = [genes[index] for index in indices]
        reduced = diagnostic(counts[indices], samples, normalization=config["normalization"])
        metrics = preservation(full, reduced, genes, names, config["top_k"])
        candidate = staging / f"candidates/{size}/counts.tsv"
        write_counts(candidate, names, sample_ids, counts[indices])
        passed = all(metrics[name] >= threshold for name, threshold in config["minimum"].items())
        curve.append({
            "features": size,
            "bytes": candidate.stat().st_size,
            "size_ratio": candidate.stat().st_size / transfer["resources"][0]["bytes"],
            "metrics": metrics,
            "pass": passed,
        })
        if selected is None and passed:
            selected = size
    report = {
        "accession": config["accession"],
        "status": "diagnostic_only",
        "source_features": len(genes),
        "source_samples": len(sample_ids),
        "conditions": dict(seen),
        "sample_order": sample_ids,
        "contrast": full["contrast"],
        "baseline": config["baseline"],
        "normalization": config["normalization"],
        "full_size_factors": full["size_factors"].tolist(),
        "candidate_curve": curve,
        "smallest_passing_features": selected,
        "rights": config["rights"],
        "reference": config["reference"],
        "release_eligible": False,
        "source_resources": transfer["resources"],
        "versions": {"python": platform.python_version(), "numpy": np.__version__},
        "config_sha256": digest(config_path),
        "workflow_sha256": digest(Path(__file__)),
    }
    (staging / "diagnostic-report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "workflows/expression-atlas-source.json")
    parser.add_argument("--staging", type=Path, required=True)
    args = parser.parse_args()
    if args.staging.resolve().is_relative_to((ROOT / "datasets").resolve()):
        parser.error("Use a staging directory outside datasets until the candidate is certified.")
    print(json.dumps(run(args.config.resolve(), args.staging.resolve()), indent=2))
