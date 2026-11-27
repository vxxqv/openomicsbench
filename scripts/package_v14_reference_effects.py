"""Add checked full-source DESeq2 reference effects to the v14 RNA-seq objects."""

import argparse
from pathlib import Path

from build_reference_effects import extract, read_json, write_json
from omicsbench.hashing import digest


ROOT = Path(__file__).resolve().parents[1]
OBJECTS = {
    "rnaseq-014": "e-mtab-7925",
    "rnaseq-015": "e-mtab-9639",
    "rnaseq-016": "e-mtab-10374",
    "rnaseq-017": "e-mtab-3101",
}


def package(dataset_id: str, slug: str, staging_root: Path) -> int:
    folder = ROOT / "datasets/rnaseq" / dataset_id
    source = staging_root / f"openomicsbench-v14-{slug}-deseq2/full/differential_expression.tsv"
    evidence = read_json(folder / "expected/deseq2.json")
    expected = evidence["artifacts"]["full"]["differential_expression.tsv"]
    if not source.is_file() or source.stat().st_size != expected["bytes"] or digest(source) != expected["sha256"]:
        raise ValueError(f"{dataset_id}: full DESeq2 result differs from its evidence record")

    output = folder / "expected/reference-effects.tsv.gz"
    count = extract(source, output)
    profile_path = folder / "expected/validation.json"
    profile = read_json(profile_path)
    profile["reference_effects"] = "expected/reference-effects.tsv.gz"
    write_json(profile_path, profile)

    card_path = folder / "README.md"
    card = card_path.read_text(encoding="utf-8")
    effects_line = "- `expected/reference-effects.tsv.gz` contains finite full-source DESeq2 effects for result comparison.\n"
    if effects_line not in card:
        anchor = "- `expected/source-counts.tsv` contains the full source matrix for the declared samples.\n"
        if anchor not in card:
            raise ValueError(f"{dataset_id}: README file list has changed")
        card_path.write_text(card.replace(anchor, anchor + effects_line), encoding="utf-8", newline="\n")

    manifest_path = folder / "manifest.json"
    manifest = read_json(manifest_path)
    files = {record["path"]: record for record in manifest["files"]}
    files["expected/reference-effects.tsv.gz"] = {
        "path": "expected/reference-effects.tsv.gz",
        "tier": "expected",
        "role": "effects",
        "media_type": "application/gzip",
        "bytes": output.stat().st_size,
        "sha256": digest(output),
    }
    files["expected/validation.json"]["bytes"] = profile_path.stat().st_size
    files["expected/validation.json"]["sha256"] = digest(profile_path)
    files["README.md"]["bytes"] = card_path.stat().st_size
    files["README.md"]["sha256"] = digest(card_path)
    manifest["files"] = sorted(files.values(), key=lambda record: record["path"])
    write_json(manifest_path, manifest)
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-root", type=Path, required=True)
    args = parser.parse_args()
    for dataset_id, slug in OBJECTS.items():
        count = package(dataset_id, slug, args.staging_root)
        print(f"{dataset_id}: {count} finite effects")


if __name__ == "__main__":
    main()
