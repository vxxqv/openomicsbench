"""Audit collection identity, sample reuse, file reuse and repeated prose."""
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from omicsbench.registry import registry


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "release/collection-audit.json"


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def contrast(folder: Path) -> dict:
    profile = json.loads((folder / "expected/validation.json").read_text(encoding="utf-8"))
    return profile["contrast"]


def markdown_paragraphs() -> dict[str, list[str]]:
    locations = defaultdict(list)
    for path in sorted(ROOT.rglob("*.md")):
        if any(part in {".git", ".snakemake", "build", "staging"} for part in path.parts):
            continue
        in_code = False
        current = []
        for line in path.read_text(encoding="utf-8").splitlines() + [""]:
            if line.startswith("```"):
                in_code = not in_code
                continue
            if in_code or line.startswith("#") or line.startswith("-"):
                continue
            if line.strip():
                current.append(line.strip())
                continue
            if current:
                paragraph = re.sub(r"\s+", " ", " ".join(current)).strip()
                if len(paragraph) >= 100:
                    locations[paragraph].append(path.relative_to(ROOT).as_posix())
                current = []
    return locations


def main() -> None:
    records = [
        (model, folder)
        for model, folder in registry(ROOT).values()
        if model.kind == "real" and model.status == "validated"
    ]
    blockers = []
    ids = [model.id for model, _ in records]
    titles = [model.title.casefold() for model, _ in records]
    if len(ids) != len(set(ids)):
        blockers.append("Dataset identifiers are repeated.")
    if len(titles) != len(set(titles)):
        blockers.append("Dataset titles are repeated.")

    object_keys = {}
    contrasts = {}
    samples = {}
    files_by_hash = defaultdict(list)
    for model, folder in records:
        attribution = json.loads((folder / "attribution.json").read_text(encoding="utf-8"))
        key = attribution["object_key"]
        if key in object_keys:
            blockers.append(f"Object key is repeated: {key}.")
        object_keys[key] = model.id
        contrasts[model.id] = contrast(folder)
        samples[model.id] = {sample.sample_id for sample in model.samples}
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                files_by_hash[sha256(path)].append(
                    {"dataset_id": model.id, "path": path.relative_to(folder).as_posix()}
                )

    overlap = []
    for index, (left, _) in enumerate(records):
        for right, _ in records[index + 1 :]:
            shared = sorted(samples[left.id] & samples[right.id])
            if not shared:
                continue
            same_source = left.source.accession == right.source.accession
            if not same_source:
                blockers.append(
                    f"{left.id} and {right.id} share sample identifiers across different accessions."
                )
            overlap.append(
                {
                    "left": left.id,
                    "right": right.id,
                    "source_accession": left.source.accession if same_source else None,
                    "shared_samples": len(shared),
                    "left_samples": len(samples[left.id]),
                    "right_samples": len(samples[right.id]),
                    "classification": "same-study contrast or declared subset" if same_source else "conflict",
                }
            )

    exact_reuse = []
    by_id = {model.id: (model, folder) for model, folder in records}
    for digest_value, group in sorted(files_by_hash.items()):
        dataset_ids = {item["dataset_id"] for item in group}
        if len(dataset_ids) < 2:
            continue
        paths = {item["path"] for item in group}
        classification = None
        if paths == {"rights.json"}:
            classification = "shared provider licence record"
        elif paths == {"expected/source-counts.tsv"}:
            models = [by_id[dataset_id][0] for dataset_id in sorted(dataset_ids)]
            same_accession = len({model.source.accession for model in models}) == 1
            same_samples = len({tuple(sample.sample_id for sample in model.samples) for model in models}) == 1
            distinct_contrasts = len(
                {json.dumps(contrasts[model.id], sort_keys=True) for model in models}
            ) == len(models)
            if same_accession and same_samples and distinct_contrasts:
                classification = "same-study source matrix used for distinct declared contrasts"
        if classification is None:
            blockers.append(
                "Unclassified exact file reuse: "
                + ", ".join(f"{item['dataset_id']}/{item['path']}" for item in group)
            )
            classification = "conflict"
        exact_reuse.append(
            {"sha256": digest_value, "files": group, "classification": classification}
        )

    repeated_paragraphs = []
    for paragraph, paths in markdown_paragraphs().items():
        unique_paths = sorted(set(paths))
        if len(unique_paths) < 2:
            continue
        classification = "conflict"
        dataset_ids = []
        for relative in unique_paths:
            parts = Path(relative).parts
            if len(parts) == 4 and parts[:2] == ("datasets", "rnaseq") and parts[3] == "README.md":
                dataset_ids.append(parts[2])
        if len(dataset_ids) == len(unique_paths) and paragraph.startswith("Source citation:"):
            accessions = {by_id[dataset_id][0].source.accession for dataset_id in dataset_ids}
            if len(accessions) == 1:
                classification = "required citation shared by objects from one source study"
        if classification == "conflict":
            blockers.append("Unclassified repeated long prose: " + ", ".join(unique_paths))
        repeated_paragraphs.append(
            {"files": unique_paths, "text": paragraph, "classification": classification}
        )

    report = {
        "schema_version": "1.0",
        "status": "pass" if not blockers else "fail",
        "biological_objects": len(records),
        "source_studies": len({model.source.accession for model, _ in records}),
        "unique_object_keys": len(object_keys),
        "unique_titles": len(set(titles)),
        "shared_sample_sets": overlap,
        "exact_file_reuse": exact_reuse,
        "repeated_long_prose": repeated_paragraphs,
        "blockers": blockers,
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(1 if blockers else 0)


if __name__ == "__main__":
    main()
