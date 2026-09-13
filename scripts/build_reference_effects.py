"""Package finite full-source DESeq2 effects for result comparison."""
import csv
import gzip
import io
import json
import math
from pathlib import Path

from omicsbench.hashing import digest


ROOT = Path(__file__).resolve().parents[1]
RUNS = {
    "rnaseq-002": "e-mtab-8572-deseq2",
    "rnaseq-003": "e-mtab-6866-deseq2",
    "rnaseq-004": "e-geod-33979-deseq2",
    "rnaseq-005": "e-mtab-567-deseq2",
    "rnaseq-006": "e-mtab-8845-deseq2",
    "rnaseq-007": "e-mtab-9206-deseq2",
    "rnaseq-008": "e-mtab-10322-deseq2",
    "rnaseq-009": "e-mtab-8845-genotype-deseq2",
    "rnaseq-010": "e-mtab-8845-infection-wild-deseq2",
    "rnaseq-011": "e-mtab-8845-infection-gsnor1-deseq2",
    "rnaseq-012": "e-mtab-10322-mdx-deseq2",
    "rnaseq-013": "e-mtab-10322-beta-deseq2",
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")


def extract(source: Path, destination: Path) -> int:
    rows = []
    with source.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            gene = row.get("") or row.get("gene_id")
            try:
                effect = float(row["log2FoldChange"])
            except (KeyError, TypeError, ValueError):
                continue
            if gene and math.isfinite(effect):
                rows.append((gene, format(effect, ".17g")))
    if len(rows) < 50 or len({gene for gene, _ in rows}) != len(rows):
        raise ValueError(f"{source}: reference effects are missing or duplicated")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(["gene_id", "log2_fold_change"])
                writer.writerows(rows)
    return len(rows)


def update_manifest(folder: Path, output: Path) -> None:
    profile_path = folder / "expected/validation.json"
    profile = read_json(profile_path)
    profile["reference_effects"] = "expected/reference-effects.tsv.gz"
    write_json(profile_path, profile)

    manifest_path = folder / "manifest.json"
    manifest = read_json(manifest_path)
    manifest["release"] = "1.1.0"
    records = {record["path"]: record for record in manifest["files"]}
    records["expected/reference-effects.tsv.gz"] = {
        "path": "expected/reference-effects.tsv.gz",
        "tier": "expected",
        "role": "effects",
        "media_type": "application/gzip",
        "bytes": output.stat().st_size,
        "sha256": digest(output),
    }
    validation = records["expected/validation.json"]
    validation["bytes"] = profile_path.stat().st_size
    validation["sha256"] = digest(profile_path)
    manifest["files"] = sorted(records.values(), key=lambda record: record["path"])
    write_json(manifest_path, manifest)


def main() -> None:
    total = 0
    for dataset_id, run in RUNS.items():
        folder = ROOT / "datasets/rnaseq" / dataset_id
        evidence = read_json(folder / "expected/deseq2.json")
        source = ROOT / "staging" / run / "full/differential_expression.tsv"
        expected = evidence["artifacts"]["full"]["differential_expression.tsv"]
        if not source.is_file() or source.stat().st_size != expected["bytes"] or digest(source) != expected["sha256"]:
            raise ValueError(f"{dataset_id}: full DESeq2 result does not match its evidence record")
        output = folder / "expected/reference-effects.tsv.gz"
        count = extract(source, output)
        update_manifest(folder, output)
        total += count
        print(f"{dataset_id}: {count} finite effects")
    print(f"Packaged {total} reference effects for {len(RUNS)} biological objects.")


if __name__ == "__main__":
    main()
