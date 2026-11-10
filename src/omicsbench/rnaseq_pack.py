"""Portable RNA-seq test inputs and explicit result-table conversion."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import shutil
import tempfile
from pathlib import Path

from .compare import compare, read_effects
from .hashing import digest


GENE_HEADERS = ("gene_id", "Geneid", "GeneID", "gene", "Gene", "feature_id", "Unnamed: 0", "")
EFFECT_HEADERS = ("log2_fold_change", "log2FoldChange", "log2FC", "logFC")
PACK_FILES = ("counts.tsv", "samples.tsv", "contrasts.tsv")
SAMPLE_FIELDS = ("sample", "condition", "batch", "replicate", "strandedness")


def _design(benchmark) -> dict:
    data = benchmark.rnaseq()
    evidence_path = benchmark.file(role="provenance", tier="expected")
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    design = evidence.get("design", {})
    contrast = design.get("contrast", "")
    if not isinstance(contrast, str) or contrast.count(" versus ") != 1:
        raise ValueError(f"{benchmark.id}: the reference contrast cannot be converted unambiguously")
    target, reference = contrast.split(" versus ")
    conditions = {row["condition"] for row in data.metadata}
    if not target or not reference or {target, reference} != conditions:
        raise ValueError(f"{benchmark.id}: reference contrast and sample conditions differ")
    formula = design.get("formula")
    if formula not in {"~ condition", "~ batch + condition"}:
        raise ValueError(f"{benchmark.id}: unsupported reference design: {formula}")
    blocking = "batch" if formula == "~ batch + condition" else ""
    if blocking and any(not row.get("batch") or row["batch"] == "not_reported" for row in data.metadata):
        raise ValueError(f"{benchmark.id}: reference design requires a declared batch for every sample")
    return {
        "id": f"{benchmark.id}_contrast", "variable": "condition", "reference": reference,
        "target": target, "blocking": blocking, "formula": formula,
        "evidence_sha256": digest(evidence_path),
    }


def _write_one(benchmark, destination: Path, tier: str) -> dict:
    data = benchmark.rnaseq(tier)
    design = _design(benchmark)
    destination.mkdir()
    with (destination / "counts.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *data.samples])
        writer.writerows((gene, *(int(value) for value in row)) for gene, row in zip(data.genes, data.counts))
    with (destination / "samples.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(SAMPLE_FIELDS)
        writer.writerows((row["sample_id"], row["condition"], row["batch"],
                          row["replicate"], row["strandedness"]) for row in data.metadata)
    with (destination / "contrasts.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("id", "variable", "reference", "target", "blocking"))
        writer.writerow(tuple(design[key] for key in ("id", "variable", "reference", "target", "blocking")))
    files = {name: {"bytes": (destination / name).stat().st_size, "sha256": digest(destination / name)}
             for name in PACK_FILES}
    record = {
        "schema_version": "1.0", "benchmark_id": benchmark.id, "tier": tier,
        "features": len(data.genes), "samples": len(data.samples),
        "design": design, "source": data.source, "rights": benchmark.manifest["rights"],
        "limitations": list(data.limitations), "files": files,
        "scope": "Checked input data and reference contrast; this pack does not establish pipeline accuracy.",
    }
    (destination / "pack.json").write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    return {"id": benchmark.id, "pack_sha256": digest(destination / "pack.json"), "files": 4,
            "features": len(data.genes), "samples": len(data.samples), "contrast": design["id"]}


def create_rnaseq_pack(collection, destination: str | Path, identifiers: list[str] | None = None,
                       tier: str = "pocket") -> dict:
    """Write a checked, multi-study input pack without exposing expected results."""
    if tier not in {"pocket", "full"}:
        raise ValueError("tier must be 'pocket' or 'full'")
    if identifiers is None:
        identifiers = [row["id"] for row in collection.list(assay="bulk_rna_seq", status="validated", kind="real")]
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ValueError("select at least one unique benchmark ID")
    identifiers = sorted(identifiers)
    benchmarks = [collection.get(identifier) for identifier in identifiers]
    for benchmark in benchmarks:
        _design(benchmark)
    destination = Path(destination).absolute()
    if destination.exists():
        raise FileExistsError(f"{destination}: output already exists")
    if not destination.parent.is_dir():
        raise FileNotFoundError(f"{destination.parent}: parent directory does not exist")
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}-", dir=destination.parent))
    try:
        cases = [_write_one(benchmark, temporary / benchmark.id, tier) for benchmark in benchmarks]
        index = {"schema_version": "1.0", "tier": tier, "cases": cases,
                 "scope": "Input pack for workflow testing; expected results remain in the installed collection."}
        (temporary / "pack-index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8", newline="\n")
        os.replace(temporary, destination)
        return {"operation": "create_rnaseq_pack", "directory": str(destination), "tier": tier,
                "cases": cases, "files": 1 + sum(case["files"] for case in cases),
                "summary": {"status": "pass", "total": len(cases), "passed": len(cases), "failed": 0, "skipped": 0}}
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def verify_rnaseq_pack(collection, directory: str | Path) -> dict:
    """Recheck every pack file, its IDs, contrast and original count tier."""
    directory = Path(directory).resolve()
    index = json.loads((directory / "pack-index.json").read_text(encoding="utf-8"))
    if not isinstance(index, dict) or index.get("schema_version") != "1.0" or index.get("tier") not in {"pocket", "full"}:
        raise ValueError(f"{directory}: invalid pack-index.json")
    cases = index.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError(f"{directory}: pack index has no cases")
    identifiers = [case.get("id") for case in cases if isinstance(case, dict)]
    if (len(identifiers) != len(cases) or any(not isinstance(identifier, str) for identifier in identifiers)
            or len(set(identifiers)) != len(cases)):
        raise ValueError(f"{directory}: pack index has duplicate or malformed cases")
    if set(path.name for path in directory.iterdir()) != {"pack-index.json", *identifiers}:
        raise ValueError(f"{directory}: pack directory inventory differs")
    results = []
    for case in cases:
        identifier = case["id"]
        benchmark = collection.get(identifier)
        folder = directory / identifier
        if not folder.is_dir() or set(path.name for path in folder.iterdir()) != {*PACK_FILES, "pack.json"}:
            raise ValueError(f"{identifier}: pack file inventory differs")
        record_path = folder / "pack.json"
        record = json.loads(record_path.read_text(encoding="utf-8"))
        data = benchmark.rnaseq(index["tier"])
        design = _design(benchmark)
        expected = {
            "id": identifier, "pack_sha256": digest(record_path), "files": 4,
            "features": len(data.genes), "samples": len(data.samples), "contrast": design["id"],
        }
        if case != expected:
            raise ValueError(f"{identifier}: pack index entry differs")
        fields = {"schema_version": "1.0", "benchmark_id": identifier, "tier": index["tier"],
                  "features": len(data.genes), "samples": len(data.samples), "design": design,
                  "source": data.source, "rights": benchmark.manifest["rights"],
                  "limitations": list(data.limitations)}
        if any(record.get(key) != value for key, value in fields.items()):
            raise ValueError(f"{identifier}: pack metadata differs from the installed collection")
        if set(record) != {*fields, "files", "scope"} or record["scope"] != "Checked input data and reference contrast; this pack does not establish pipeline accuracy.":
            raise ValueError(f"{identifier}: pack metadata has unsupported fields or scope")
        if record.get("files") is None or set(record["files"]) != set(PACK_FILES):
            raise ValueError(f"{identifier}: pack file inventory is invalid")
        for name in PACK_FILES:
            path = folder / name
            declared = record["files"].get(name)
            if declared != {"bytes": path.stat().st_size, "sha256": digest(path)}:
                raise ValueError(f"{identifier}: {name} failed its declared checksum")
        comparison = benchmark.compare_counts(folder / "counts.tsv", index["tier"], detail_limit=0)
        if comparison["status"] != "pass":
            raise ValueError(f"{identifier}: counts differ from the installed collection")
        with (folder / "samples.tsv").open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            rows = list(reader)
        expected_rows = [{"sample": row["sample_id"], "condition": row["condition"],
                          "batch": str(row["batch"]), "replicate": str(row["replicate"]),
                          "strandedness": row["strandedness"]}
                         for row in data.metadata]
        if reader.fieldnames != list(SAMPLE_FIELDS) or rows != expected_rows:
            raise ValueError(f"{identifier}: sample sheet differs from the declared order or design")
        with (folder / "contrasts.tsv").open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            contrasts = list(reader)
        contrast_keys = ["id", "variable", "reference", "target", "blocking"]
        if reader.fieldnames != contrast_keys or contrasts != [{key: design[key] for key in contrast_keys}]:
            raise ValueError(f"{identifier}: contrast differs from the reference design")
        results.append({"id": identifier, "status": "pass", "features": len(data.genes), "samples": len(data.samples)})
    return {"operation": "verify_rnaseq_pack", "directory": str(directory), "tier": index["tier"],
            "cases": results, "files": 1 + 4 * len(results),
            "summary": {"status": "pass", "total": len(results), "passed": len(results), "failed": 0, "skipped": 0}}


def normalize_effects(benchmark, source: str | Path, destination: str | Path, *,
                      gene_column: str | None = None, effect_column: str | None = None,
                      reverse_sign: bool = False) -> dict:
    """Convert a DESeq2, edgeR or limma style table to the canonical result columns."""
    design = _design(benchmark)
    source = Path(source)
    suffixes = source.suffixes
    compressed = suffixes[-1:] == [".gz"]
    suffix = suffixes[-2] if compressed and len(suffixes) > 1 else source.suffix
    if suffix not in {".csv", ".tsv"}:
        raise ValueError("effect input must be CSV, TSV or a gzip-compressed equivalent")
    opener = gzip.open if compressed else open
    with opener(source, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="," if suffix == ".csv" else "\t")
        headers = reader.fieldnames
        if not headers or len(headers) != len(set(headers)):
            raise ValueError("effect table has an empty or duplicate header")
        def choose(explicit, aliases, label):
            if explicit is not None:
                if explicit not in headers:
                    raise ValueError(f"{label} column is absent: {explicit}")
                return explicit
            found = [name for name in aliases if name in headers]
            if len(found) != 1:
                raise ValueError(f"{label} column is {'ambiguous' if found else 'missing'}; specify it explicitly")
            return found[0]
        gene_name = choose(gene_column, GENE_HEADERS, "gene")
        effect_name = choose(effect_column, EFFECT_HEADERS, "effect")
        if gene_name == effect_name:
            raise ValueError("gene and effect columns must differ")
        output = []
        seen = set()
        dropped = 0
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{source}:{line}: inconsistent row width")
            gene = row[gene_name].strip()
            if not gene or gene in seen:
                raise ValueError(f"{source}:{line}: empty or duplicate gene ID")
            seen.add(gene)
            raw = row[effect_name].strip()
            if raw in {"", "NA", "NaN", "nan", "NULL"}:
                dropped += 1
                continue
            try:
                value = float(raw)
            except ValueError:
                raise ValueError(f"{source}:{line}: effect is not numeric") from None
            if not math.isfinite(value):
                raise ValueError(f"{source}:{line}: effect must be finite or missing")
            output.append((gene, -value if reverse_sign else value))
    if not output:
        raise ValueError("effect table has no finite effects")
    reference = read_effects(benchmark.file(role="effects", tier="expected"))
    overlap = len(set(gene for gene, _ in output) & set(reference))
    destination = Path(destination).absolute()
    if destination.exists():
        raise FileExistsError(f"{destination}: output already exists")
    if not destination.parent.is_dir():
        raise FileNotFoundError(f"{destination.parent}: parent directory does not exist")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=destination.parent,
                                         prefix=f".{destination.name}-", delete=False) as stream:
            temporary = Path(stream.name)
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(("gene_id", "log2_fold_change"))
            writer.writerows((gene, format(value, ".17g")) for gene, value in output)
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return {"operation": "normalize_effects", "benchmark_id": benchmark.id,
            "source_sha256": digest(source), "output": str(destination), "output_sha256": digest(destination),
            "columns": {"gene_id": gene_name, "log2_fold_change": effect_name},
            "finite_effects": len(output), "missing_effects_dropped": dropped,
            "shared_reference_genes": overlap, "reference_genes": len(reference),
            "contrast": {"target": design["target"], "reference": design["reference"]},
            "sign_reversed": reverse_sign,
            "scope": "Format conversion only; compare or grade the output to assess agreement."}


def score_native_effects(benchmark, source: str | Path, *, gene_column: str | None = None,
                         effect_column: str | None = None, reverse_sign: bool = False,
                         detail_limit: int = 20) -> dict:
    """Compare a native result table without asking the user to rename or copy it."""
    with tempfile.TemporaryDirectory(prefix="omicsbench-effects-") as temporary:
        normalized = Path(temporary) / "effects.tsv"
        conversion = normalize_effects(benchmark, source, normalized, gene_column=gene_column,
                                       effect_column=effect_column, reverse_sign=reverse_sign)
        comparison = compare(benchmark._model, benchmark._folder, normalized, detail_limit)
    conversion.pop("output")
    conversion.pop("output_sha256")
    comparison["input"] = str(Path(source))
    return {"operation": "score_native_effects", "benchmark_id": benchmark.id,
            "status": comparison["status"], "conversion": conversion, "comparison": comparison}
