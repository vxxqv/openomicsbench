"""Evaluate deterministic reductions of one staged Expression Atlas experiment."""
import argparse
import csv
import json
import platform
from collections import Counter
from pathlib import Path
import numpy as np

from omicsbench.expression_atlas import load_counts, load_design_fields, stage
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
    genes, all_sample_ids, counts = load_counts(source / f'{config["accession"]}-raw-counts.tsv')
    factor = config["factor_column"]
    batch_column = config.get("batch_column")
    fields = [factor] + ([batch_column] if batch_column else [])
    design = load_design_fields(source / f'{config["accession"]}-experiment-design.tsv', fields)
    if not set(design) <= set(all_sample_ids):
        raise ValueError("Count columns are missing analysed design runs")
    analysed_columns = [index for index, sample_id in enumerate(all_sample_ids) if sample_id in design]
    sample_ids = [all_sample_ids[index] for index in analysed_columns]
    counts = counts[:, analysed_columns]
    condition_values = config.get("condition_values")
    if condition_values is not None:
        if not isinstance(condition_values, list) or len(condition_values) != 2 or len(set(condition_values)) != 2 or any(not value for value in condition_values):
            raise ValueError("condition_values must contain two distinct non-empty values")
        allowed = set(condition_values)
        observed = {values[factor] for values in design.values()}
        if not allowed <= observed:
            raise ValueError("condition_values contains a value absent from the experiment design")
        selected_columns = [index for index, sample_id in enumerate(sample_ids) if design[sample_id][factor] in allowed]
        sample_ids = [sample_ids[index] for index in selected_columns]
        counts = counts[:, selected_columns]
    else:
        if len({values[factor] for values in design.values()}) != 2:
            raise ValueError("The diagnostic requires exactly two factor values or an explicit condition_values selection")
    seen: Counter[str] = Counter()
    batches: Counter[str] = Counter()
    samples = []
    for sample_id in sample_ids:
        condition = design[sample_id][factor]
        batch = design[sample_id][batch_column] if batch_column else "not_reported"
        seen[condition] += 1
        batches[batch] += 1
        samples.append(Sample(sample_id=sample_id, condition=condition, replicate=seen[condition], batch=batch, strandedness="unknown"))
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
        "source_samples": len(all_sample_ids),
        "analysed_samples": len(design),
        "selected_samples": len(sample_ids),
        "excluded_count_columns": transfer["excluded_count_columns"],
        "deduplicated_count_columns": transfer["deduplicated_count_columns"],
        "conditions": dict(seen),
        "factor_column": factor,
        "batch_column": batch_column,
        "batches": dict(batches),
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
