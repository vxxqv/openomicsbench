"""Build the certified biological count-matrix collection from staged Atlas runs."""
import csv
import json
import os
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from omicsbench.expression_atlas import load_counts, load_design_fields, select_design_samples
from omicsbench.hashing import digest
from omicsbench.models import Dataset

ROOT = Path(__file__).resolve().parents[1]
COLLECTION_ROOT = (ROOT / "datasets/rnaseq").resolve()
GIT = os.environ.get("OPENOMICSBENCH_GIT", "git")
METRICS = [
    {
        "name": "spearman_logfc",
        "minimum": 0.9,
        "definition": "Average-tie rank correlation of median-ratio-normalized contrast effects on the shared selected feature universe; higher is better.",
    },
    {
        "name": "top_k_jaccard",
        "minimum": 0.6,
        "definition": "Jaccard overlap of the top 50 absolute contrast effects on the shared selected feature universe, with ties settled by gene ID; higher is better.",
    },
    {
        "name": "distance_correlation",
        "minimum": 0.9,
        "definition": "Pearson correlation of sample-pair distances on log2 median-ratio-normalized counts for the full and pocket matrices; higher is better.",
    },
    {
        "name": "sign_concordance",
        "minimum": 0.9,
        "definition": "Fraction of matching contrast-effect signs where the absolute full-source effect is at least 0.1; higher is better.",
    },
]


def git_commit(path: str) -> str:
    result = subprocess.run(
        [GIT, "-c", f"safe.directory={ROOT.as_posix()}", "log", "-1", "--format=%H", "--", path],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    commit = result.stdout.strip()
    if len(commit) != 40:
        raise ValueError(f"No committed provenance found for {path}")
    return commit


def write_counts(path: Path, genes: list[str], samples: list[str], counts) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])
        writer.writerows([[gene, *map(int, row)] for gene, row in zip(genes, counts, strict=True)])


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")


def file_record(folder: Path, relative: str, tier: str, role: str, media_type: str) -> dict:
    path = folder / relative
    return {"path": relative, "tier": tier, "role": role, "media_type": media_type, "bytes": path.stat().st_size, "sha256": digest(path)}


def selected_inputs(config: dict, diagnostic_dir: Path):
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
    return genes, samples, counts[:, [lookup[sample] for sample in samples]], design


def reference_profiles() -> tuple[dict, str]:
    path = ROOT / "evidence/reference-verification.json"
    evidence = json.loads(path.read_text(encoding="utf-8"))
    if evidence["status"] != "pass":
        raise ValueError("Reference verification has not passed")
    profiles = {}
    for profile in evidence["profiles"]:
        for check in profile["pocket_checks"]:
            if check["status"] != "pass":
                raise ValueError(f"Reference check failed for {check['object_key']}")
            profiles[check["object_key"]] = profile
    return profiles, digest(path)


def build(entry: dict, release: str, checked: str, references: dict, reference_hash: str) -> Path:
    dataset_id = entry["dataset_id"]
    object_key = entry["object_key"]
    slug = object_key.lower()
    target = ROOT / f"datasets/rnaseq/{dataset_id}"
    partial = target.with_name(f".{dataset_id}.partial")
    if target.resolve().parent != COLLECTION_ROOT or partial.resolve().parent != COLLECTION_ROOT:
        raise ValueError(f"Collection path escapes datasets/rnaseq: {dataset_id}")
    if target.exists() or partial.exists():
        raise FileExistsError(f"Refusing to replace an existing collection object: {target}")
    partial.mkdir(parents=True)
    try:
        config_path = ROOT / entry["config"]
        config = json.loads(config_path.read_text(encoding="utf-8"))
        diagnostic_dir = ROOT / f"staging/{slug}-diagnostic"
        diagnostic_path = ROOT / f"evidence/{slug}-diagnostic.json"
        diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
        baseline_path = ROOT / f"evidence/{slug}-deseq2.json"
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        intake_path = ROOT / f"curation/intake/{slug}.json"
        intake = json.loads(intake_path.read_text(encoding="utf-8"))
        if baseline["status"] != "pass" or diagnostic["smallest_passing_features"] is None:
            raise ValueError(f"Unqualified object: {object_key}")
        genes, samples, counts, design = selected_inputs(config, diagnostic_dir)
        factor = config["factor_column"]
        batch_column = config.get("batch_column")
        write_counts(partial / "expected/source-counts.tsv", genes, samples, counts)
        pocket_source = diagnostic_dir / f"candidates/{diagnostic['smallest_passing_features']}/counts.tsv"
        (partial / "pocket").mkdir()
        shutil.copyfile(pocket_source, partial / "pocket/counts.tsv")
        seen = Counter()
        sample_records = []
        for sample in samples:
            condition = design[sample][factor]
            seen[condition] += 1
            sample_records.append({
                "sample_id": sample,
                "condition": condition,
                "replicate": seen[condition],
                "batch": design[sample][batch_column] if batch_column else "not_reported",
                "strandedness": "unknown",
            })
        write_json(partial / "samples.json", sample_records)
        validation = {
            "baseline_counts": "expected/source-counts.tsv",
            "sample_ids": samples,
            "shapes": {"pocket": [diagnostic["smallest_passing_features"], len(samples)]},
            "top_k": config["top_k"],
            "normalization": config["normalization"],
            "contrast": {
                "factor_column": factor,
                "condition_values": config.get("condition_values"),
                "batch_column": batch_column,
                "subset": config.get("subset", {}),
            },
            "deseq2_evidence": "expected/deseq2.json",
        }
        write_json(partial / "expected/validation.json", validation)
        shutil.copyfile(baseline_path, partial / "expected/deseq2.json")
        profile = references[object_key]
        reference = {
            "provider": profile["provider"],
            "release": profile["release"],
            "assembly": profile["assembly"],
            "genome": profile["fasta"],
            "annotation": profile["gtf"],
            "pocket_check": next(check for check in profile["pocket_checks"] if check["object_key"] == object_key),
            "collection_evidence_sha256": reference_hash,
        }
        write_json(partial / "reference.json", reference)
        attribution = {
            "dataset_id": dataset_id,
            "object_key": object_key,
            "source_study": entry["citation"],
            "source_accession": config["accession"],
            "source_url": f"https://www.ebi.ac.uk/gxa/experiments/{config['accession']}",
            "provider": "Expression Atlas, EMBL-EBI",
            "provider_citation": "Expression Atlas in 2026: enabling FAIR and open expression data through community collaboration and integration. Nucleic Acids Research. DOI 10.1093/nar/gkaf1238.",
            "curator": "Vivaan Patni",
            "repository_account": "vxxqv",
            "changes": "Selected a deterministic feature subset and retained the source integer counts for the declared samples and contrast.",
        }
        write_json(partial / "attribution.json", attribution)
        rights = {
            "status": "GREEN",
            "license": "Creative Commons Attribution 4.0 International",
            "license_id": "CC-BY-4.0",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "evidence_url": "https://www.ebi.ac.uk/gxa/licence.html",
            "credit": "Expression Atlas, EMBL-EBI",
            "attribution_file": "attribution.json",
            "checked": checked,
        }
        write_json(partial / "rights.json", rights)
        workflow_path = "workflows/atlas_diagnostic.py"
        provenance = {
            "dataset_id": dataset_id,
            "object_key": object_key,
            "algorithm": "contrast-anchored deterministic feature selection",
            "seed": config["seed"],
            "source_resources": diagnostic["source_resources"],
            "selection": {
                "config": entry["config"],
                "config_sha256": digest(config_path),
                "config_commit": git_commit(entry["config"]),
                "workflow": workflow_path,
                "workflow_sha256": digest(ROOT / workflow_path),
                "workflow_commit": git_commit(workflow_path),
                "candidate_sizes": config["candidate_sizes"],
                "selected_features": diagnostic["smallest_passing_features"],
                "factor_column": factor,
                "batch_column": batch_column,
                "condition_values": config.get("condition_values"),
                "subset": config.get("subset", {}),
            },
            "diagnostic": {"path": f"evidence/{slug}-diagnostic.json", "sha256": digest(diagnostic_path)},
            "deseq2": {"path": f"evidence/{slug}-deseq2.json", "sha256": digest(baseline_path)},
            "reference": {"path": "evidence/reference-verification.json", "sha256": reference_hash},
            "intake": {"path": f"curation/intake/{slug}.json", "sha256": digest(intake_path)},
        }
        write_json(partial / "provenance/transform.json", provenance)
        files = [
            file_record(partial, "attribution.json", "metadata", "documentation", "application/json"),
            file_record(partial, "expected/deseq2.json", "expected", "provenance", "application/json"),
            file_record(partial, "expected/source-counts.tsv", "expected", "baseline", "text/tab-separated-values"),
            file_record(partial, "expected/validation.json", "expected", "metrics", "application/json"),
            file_record(partial, "pocket/counts.tsv", "pocket", "raw_counts", "text/tab-separated-values"),
            file_record(partial, "provenance/transform.json", "metadata", "provenance", "application/json"),
            file_record(partial, "reference.json", "metadata", "reference", "application/json"),
            file_record(partial, "rights.json", "metadata", "license", "application/json"),
            file_record(partial, "samples.json", "metadata", "samples", "application/json"),
        ]
        manifest = {
            "schema_version": "1.0",
            "id": dataset_id,
            "title": entry["title"],
            "release": release,
            "assay": "bulk_rna_seq",
            "kind": "real",
            "status": "validated",
            "archetype": entry["archetype"],
            "organism": intake["organism"],
            "taxon_id": intake["taxon_id"],
            "source": {
                "repository": "Expression Atlas and BioStudies",
                "accession": config["accession"],
                "url": f"https://www.ebi.ac.uk/gxa/experiments/{config['accession']}",
                "citation": entry["citation"],
                "retrieved": checked,
                "sha256": diagnostic["source_resources"][0]["sha256"],
            },
            "rights": {
                "status": "GREEN",
                "license": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
                "evidence": "Expression Atlas applies CC BY 4.0 to copyrightable website material and requires provider citation; rights.json records the licence and attribution.",
                "checked": checked,
            },
            "reference": {
                "genome": f"{profile['provider']} release {profile['release']} {profile['assembly']}; {profile['fasta']['url']}; official checksum {profile['fasta']['checksums_evidence']['record']}",
                "annotation": f"{profile['gtf']['url']}; SHA-256 {profile['gtf']['sha256']}",
                "namespace": "Ensembl gene identifiers",
                "compatibility": "verified",
            },
            "samples": sample_records,
            "derivation": {
                "workflow": workflow_path,
                "commit": git_commit(workflow_path),
                "seed": config["seed"],
                "algorithm": "contrast-anchored deterministic feature selection",
                "parameters": {
                    "candidate_sizes": config["candidate_sizes"],
                    "selected_features": diagnostic["smallest_passing_features"],
                    "top_k": config["top_k"],
                    "factor_column": factor,
                    "batch_column": batch_column or "none",
                    "condition_values": json.dumps(config.get("condition_values")),
                    "subset": json.dumps(config.get("subset", {}), sort_keys=True),
                },
                "versions": {"python": diagnostic["versions"]["python"], "numpy": diagnostic["versions"]["numpy"], "r": "4.5.3", "deseq2": "1.50.2"},
            },
            "files": files,
            "validation": {
                "profile": "expected/validation.json",
                "baseline_version": "deseq2-1.50.2-and-median-ratio-envelope-v1",
                "metrics": METRICS,
            },
            "limitations": [
                "The object starts from an archive count matrix and does not revalidate read alignment or quantification.",
                "The pocket is intended for workflow testing and method validation, not as a substitute for the full study.",
                "Strandedness is recorded as unknown because it is not required for count-matrix validation.",
            ],
        }
        Dataset.model_validate(manifest)
        write_json(partial / "manifest.json", manifest)
        os.replace(partial, target)
        return target
    except Exception:
        if partial.exists():
            shutil.rmtree(partial)
        raise


if __name__ == "__main__":
    plan = json.loads((ROOT / "curation/collection-plan.json").read_text(encoding="utf-8"))
    references, reference_hash = reference_profiles()
    built = [str(build(entry, plan["release"], plan["checked"], references, reference_hash).relative_to(ROOT)) for entry in plan["objects"]]
    print(json.dumps({"built": built}, indent=2))
