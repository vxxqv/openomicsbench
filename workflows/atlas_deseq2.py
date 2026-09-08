"""Confirm an Expression Atlas pocket against its full matrix with DESeq2."""
import argparse
import csv
import json
import math
import platform
import subprocess
from pathlib import Path

import numpy as np

from omicsbench.expression_atlas import load_counts, load_design_fields, select_design_samples
from omicsbench.hashing import digest
from omicsbench.rnaseq import correlation, ranked

ROOT = Path(__file__).resolve().parents[1]


def write_counts(path: Path, genes: list[str], samples: list[str], counts: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])
        writer.writerows([[gene, *map(int, row)] for gene, row in zip(genes, counts, strict=True)])


def design_inputs(config: dict, diagnostic_dir: Path):
    accession = config["accession"]
    source = diagnostic_dir / "source"
    genes, count_samples, counts = load_counts(source / f"{accession}-raw-counts.tsv")
    factor = config["factor_column"]
    batch_column = config.get("batch_column")
    subset = config.get("subset", {})
    fields = list(dict.fromkeys([factor] + ([batch_column] if batch_column else []) + list(subset)))
    design = load_design_fields(source / f"{accession}-experiment-design.tsv", fields)
    samples = select_design_samples(count_samples, design, factor, config.get("condition_values"), subset)
    lookup = {sample: index for index, sample in enumerate(count_samples)}
    counts = counts[:, [lookup[sample] for sample in samples]]
    metadata = [
        (sample, design[sample][factor], design[sample][batch_column] if batch_column else "not_reported")
        for sample in samples
    ]
    return genes, samples, counts, metadata


def read_effects(path: Path) -> dict[str, float]:
    effects = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            gene = row.get("") or row.get("gene_id")
            try:
                value = float(row["log2FoldChange"])
            except (KeyError, TypeError, ValueError):
                continue
            if gene and math.isfinite(value):
                effects[gene] = value
    return effects


def sample_distances(path: Path) -> np.ndarray:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if len(header) < 3:
            raise ValueError("DESeq2 normalized counts need at least two samples")
        values = np.asarray([[float(value) for value in row[1:]] for row in reader], dtype=float)
    transformed = np.log2(values + 1).T
    distances = np.sqrt(((transformed[:, None, :] - transformed[None, :, :]) ** 2).sum(axis=2))
    return distances[np.triu_indices(len(header) - 1, 1)]


def preservation(full_dir: Path, pocket_dir: Path, pocket_genes: list[str], top_k: int) -> tuple[dict, int]:
    full = read_effects(full_dir / "differential_expression.tsv")
    pocket = read_effects(pocket_dir / "differential_expression.tsv")
    shared = [gene for gene in pocket_genes if gene in full and gene in pocket]
    if len(shared) < top_k:
        raise ValueError("Too few finite shared DESeq2 effects for the declared top-k metric")
    a = np.asarray([full[gene] for gene in shared])
    b = np.asarray([pocket[gene] for gene in shared])
    top = lambda values: set(sorted(range(len(values)), key=lambda index: (-abs(values[index]), shared[index]))[:top_k])
    aa, bb = top(a), top(b)
    eligible = np.abs(a) >= 0.1
    if not np.any(eligible):
        raise ValueError("No DESeq2 effects satisfy the sign-concordance filter")
    metrics = {
        "spearman_logfc": correlation(ranked(a), ranked(b)),
        "top_k_jaccard": len(aa & bb) / len(aa | bb),
        "distance_correlation": correlation(sample_distances(full_dir / "normalized_counts.tsv"), sample_distances(pocket_dir / "normalized_counts.tsv")),
        "sign_concordance": float(np.mean(np.sign(a[eligible]) == np.sign(b[eligible]))),
    }
    return metrics, len(shared)


def artifacts(directory: Path) -> dict:
    records = {}
    for path in sorted(directory.iterdir()):
        if path.is_file():
            records[path.name] = {"bytes": path.stat().st_size, "sha256": digest(path)}
    return records


def run(config_path: Path, diagnostic_dir: Path, output: Path, conda: Path, environment: Path) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    diagnostic_path = diagnostic_dir / "diagnostic-report.json"
    diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
    pocket_size = diagnostic["smallest_passing_features"]
    if pocket_size is None:
        raise ValueError("The diagnostic has no passing pocket")
    genes, samples, counts, metadata = design_inputs(config, diagnostic_dir)
    input_dir = output / "inputs"
    write_counts(input_dir / "full-counts.tsv", genes, samples, counts)
    with (input_dir / "samples.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "condition", "batch"])
        writer.writerows(metadata)
    pocket_input = diagnostic_dir / f"candidates/{pocket_size}/counts.tsv"
    pocket_genes, pocket_samples, _ = load_counts_for_baseline(pocket_input)
    if pocket_samples != samples:
        raise ValueError("Pocket and full matrices have different sample order")
    baseline = ROOT / "workflows/deseq2_baseline.R"
    for label, counts_path in (("full", input_dir / "full-counts.tsv"), ("pocket", pocket_input)):
        destination = output / label
        command = [str(conda), "run", "-p", str(environment), "Rscript", str(baseline), str(counts_path), str(input_dir / "samples.tsv"), str(destination)]
        completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        if completed.returncode:
            raise RuntimeError(f"DESeq2 {label} baseline failed: {completed.stderr.strip() or completed.stdout.strip()}")
    metrics, finite = preservation(output / "full", output / "pocket", pocket_genes, config["top_k"])
    passed = all(metrics[name] >= minimum for name, minimum in config["minimum"].items())
    contrast = (output / "full/contrast.log").read_text(encoding="utf-8").strip()
    report = {
        "schema_version": "1.0",
        "accession": config["accession"],
        "object_key": config.get("object_key", config["accession"]),
        "status": "pass" if passed else "fail",
        "release_eligible": False,
        "release_gates_remaining": ["exact reference-file checksums", "dataset attribution package", "independent clean-environment reproduction"],
        "design": {
            "formula": "~ batch + condition" if len({row[2] for row in metadata}) > 1 else "~ condition",
            "factor_column": config["factor_column"],
            "batch_column": config.get("batch_column"),
            "subset": config.get("subset", {}),
            "contrast": contrast,
            "samples": len(samples),
            "full_features": len(genes),
            "pocket_features": len(pocket_genes),
            "finite_shared_effects": finite,
        },
        "thresholds": config["minimum"],
        "metrics": metrics,
        "runtime": {"r": "4.5.3", "deseq2": "1.50.2", "python": platform.python_version(), "numpy": np.__version__},
        "inputs": {
            "diagnostic_report_sha256": digest(diagnostic_path),
            "config_sha256": digest(config_path),
            "full_counts_sha256": digest(input_dir / "full-counts.tsv"),
            "pocket_counts_sha256": digest(pocket_input),
            "samples_sha256": digest(input_dir / "samples.tsv"),
        },
        "artifacts": {"full": artifacts(output / "full"), "pocket": artifacts(output / "pocket")},
        "workflow": {"python_sha256": digest(Path(__file__)), "r_sha256": digest(baseline)},
    }
    (output / "baseline-report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    return report


def load_counts_for_baseline(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "gene_id" or len(header) < 3:
            raise ValueError("Pocket counts need gene_id and at least two samples")
        rows = list(reader)
    genes = [row[0] for row in rows]
    counts = np.asarray([[int(value) for value in row[1:]] for row in rows], dtype=np.int64)
    return genes, header[1:], counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--diagnostic-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--conda", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.config.resolve(), args.diagnostic_dir.resolve(), args.output.resolve(), args.conda.resolve(), args.environment.resolve()), indent=2))
