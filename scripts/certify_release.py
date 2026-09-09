"""Verify the v1 collection and write a fail-closed release decision."""
import argparse
import json
import re
from pathlib import Path

from omicsbench.hashing import digest
from omicsbench.registry import registry
from omicsbench.validate import validate


ROOT = Path(__file__).resolve().parents[1]
SHA256 = re.compile(r"^[0-9a-f]{64}$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
PUBLICATION_FIELDS = ("version_doi", "external_user_trial", "post_upload_verification")


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def declared_path(model, folder: Path, relative: str, blockers: list[str]) -> Path | None:
    if relative not in {record.path for record in model.files}:
        blockers.append(f"{model.id}: {relative} is not declared in the manifest.")
        return None
    path = folder / relative
    if not path.is_file():
        blockers.append(f"{model.id}: {relative} is missing.")
        return None
    return path


def check_recorded_file(record: dict, label: str, blockers: list[str]) -> None:
    relative = record.get("workflow") or record.get("path")
    if not relative:
        blockers.append(f"{label}: recorded path is missing.")
        return
    path = ROOT / relative
    if not path.is_file():
        blockers.append(f"{label}: {relative} is missing.")
        return
    expected = record.get("workflow_sha256") or record.get("sha256")
    if not expected or not SHA256.fullmatch(expected) or digest(path) != expected:
        blockers.append(f"{label}: hash does not match {relative}.")
    commit = record.get("workflow_commit")
    if commit is not None and not COMMIT.fullmatch(commit):
        blockers.append(f"{label}: transformation commit is invalid.")


def check_biological_object(model, folder: Path, reference_evidence: dict, blockers: list[str]) -> None:
    profile = read_json(folder / model.validation.profile)
    baseline_relative = profile.get("deseq2_evidence")
    if not baseline_relative:
        blockers.append(f"{model.id}: DESeq2 evidence is not named in the validation profile.")
        return
    baseline_path = declared_path(model, folder, baseline_relative, blockers)
    reference_path = declared_path(model, folder, "reference.json", blockers)
    rights_path = declared_path(model, folder, "rights.json", blockers)
    attribution_path = declared_path(model, folder, "attribution.json", blockers)
    provenance_path = declared_path(model, folder, "provenance/transform.json", blockers)
    required = (baseline_path, reference_path, rights_path, attribution_path, provenance_path)
    if any(path is None for path in required):
        return

    baseline = read_json(baseline_path)
    if baseline.get("status") != "pass":
        blockers.append(f"{model.id}: DESeq2 comparison did not pass.")
    for name, minimum in baseline.get("thresholds", {}).items():
        value = baseline.get("metrics", {}).get(name)
        if value is None or value < minimum:
            blockers.append(f"{model.id}: DESeq2 {name} is below {minimum}.")
    input_hashes = baseline.get("inputs", {})
    expected_files = {
        "full_counts_sha256": folder / "expected/source-counts.tsv",
        "pocket_counts_sha256": folder / "pocket/counts.tsv",
    }
    for field, path in expected_files.items():
        if input_hashes.get(field) != digest(path):
            blockers.append(f"{model.id}: DESeq2 {field} does not match the certified input.")
    workflow_hashes = baseline.get("workflow", {})
    workflow_files = {
        "python_sha256": ROOT / "workflows/atlas_deseq2.py",
        "r_sha256": ROOT / "workflows/deseq2_baseline.R",
    }
    for field, path in workflow_files.items():
        if workflow_hashes.get(field) != digest(path):
            blockers.append(f"{model.id}: DESeq2 workflow hash has drifted for {path.name}.")

    attribution = read_json(attribution_path)
    if attribution.get("dataset_id") != model.id:
        blockers.append(f"{model.id}: attribution dataset identifier differs.")
    if attribution.get("source_accession") != model.source.accession:
        blockers.append(f"{model.id}: attribution accession differs.")
    owner_matches = (
        attribution.get("curator") == "Vivaan Patni"
        and attribution.get("repository_account") == "vxxqv"
    )
    if not owner_matches:
        blockers.append(f"{model.id}: curator or repository account differs from release ownership.")

    rights = read_json(rights_path)
    if rights.get("status") != "GREEN" or rights.get("license_id") != "CC-BY-4.0":
        blockers.append(f"{model.id}: distributable rights evidence is incomplete.")
    if rights.get("checked") != model.rights.checked.isoformat():
        blockers.append(f"{model.id}: rights review dates differ.")

    reference = read_json(reference_path)
    object_key = attribution.get("object_key")
    check = reference.get("pocket_check", {})
    if check.get("object_key") != object_key or check.get("status") != "pass":
        blockers.append(f"{model.id}: pocket reference check is missing or belongs to another object.")
    if check.get("matched_features") != check.get("pocket_features") or check.get("unmatched_features"):
        blockers.append(f"{model.id}: pocket identifiers do not fully match the declared annotation.")
    reference_path_global = ROOT / "evidence/reference-verification.json"
    if reference.get("collection_evidence_sha256") != digest(reference_path_global):
        blockers.append(f"{model.id}: collection reference evidence hash has drifted.")
    if object_key not in reference_evidence.get("verified_objects", []):
        blockers.append(f"{model.id}: object is absent from the verified reference set.")

    provenance = read_json(provenance_path)
    if provenance.get("dataset_id") != model.id or provenance.get("object_key") != object_key:
        blockers.append(f"{model.id}: provenance identity differs.")
    check_recorded_file(provenance.get("selection", {}), f"{model.id} selection", blockers)
    check_recorded_file(provenance.get("assembly", {}), f"{model.id} assembly", blockers)
    for name in ("diagnostic", "deseq2", "reference", "intake"):
        check_recorded_file(provenance.get(name, {}), f"{model.id} {name}", blockers)


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify the v1 collection and write its release decision.")
    parser.add_argument(
        "--preflight",
        action="store_true",
        help="Exit successfully when only DOI, independent-trial or post-upload gates remain.",
    )
    args = parser.parse_args()

    metadata = read_json(ROOT / "release/metadata-input.json")
    reference_evidence = read_json(ROOT / "evidence/reference-verification.json")
    audit = read_json(ROOT / "release/collection-audit.json")
    blockers = []
    validated = []

    if reference_evidence.get("status") != "pass":
        blockers.append("Collection reference verification has not passed.")
    if audit.get("status") != "pass":
        blockers.append("Collection overlap and prose audit has not passed.")

    try:
        records = registry(ROOT)
    except ValueError as exc:
        records = {}
        blockers.append(str(exc))
    for model, folder in records.values():
        try:
            result = validate(model, folder)
        except ValueError as exc:
            blockers.append(str(exc))
            continue
        if model.kind == "real" and model.status == "validated" and result["status"] == "pass":
            validated.append(model)
            check_biological_object(model, folder, reference_evidence, blockers)

    pockets = [
        model
        for model in validated
        if model.rights.status == "GREEN" and any(record.tier == "pocket" for record in model.files)
    ]
    if not 12 <= len(validated) <= 20:
        blockers.append(f"Biological collection has {len(validated)} objects; target is 12 to 20.")
    if len(pockets) < 8:
        blockers.append(f"Collection has {len(pockets)} redistributable pockets; at least 8 are required.")
    if len({model.archetype for model in validated}) < 4:
        blockers.append("At least four certified design archetypes are required.")
    for model in validated:
        if model.derivation.commit is None or not COMMIT.fullmatch(model.derivation.commit):
            blockers.append(f"{model.id}: transformation commit is missing or invalid.")
        if model.validation.baseline_version.startswith("diagnostic"):
            blockers.append(f"{model.id}: baseline remains diagnostic-only.")

    required_metadata = (
        "authors",
        "repository_url",
        "code_license",
        "metadata_license",
        "scientific_environment_lock",
    )
    for field in required_metadata:
        if not metadata.get(field):
            blockers.append(f"Missing release evidence: {field}.")
    lock = ROOT / metadata.get("scientific_environment_lock", "")
    if not lock.is_file():
        blockers.append("The scientific environment lock is missing.")

    deferred = [field for field in PUBLICATION_FIELDS if not metadata.get(field)]
    report = {
        "scope": "v1_biological_collection",
        "software_version": metadata["software_version"],
        "decision": "GO" if not blockers and not deferred else "NO-GO",
        "preflight": "PASS" if not blockers else "FAIL",
        "certified_biological_objects": len(validated),
        "source_studies": len({model.source.accession for model in validated}),
        "redistributable_pockets": len(pockets),
        "design_archetypes": len({model.archetype for model in validated}),
        "blockers": blockers,
        "deferred_publication_gates": deferred,
    }
    output = ROOT / "release/certification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2))
    if blockers or (deferred and not args.preflight):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
