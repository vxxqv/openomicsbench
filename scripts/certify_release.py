"""Verify the v1 collection and write a fail-closed release decision."""
import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

from omicsbench.hashing import digest
from omicsbench.registry import registry
from omicsbench.validate import validate


ROOT = Path(__file__).resolve().parents[1]
SHA256 = re.compile(r"^[0-9a-f]{64}$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
ORCID = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$")
PUBLICATION_FIELDS = ("version_doi", "external_user_trial", "post_upload_verification")
INVENTORY_PATH = "release/file-manifest.json"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def valid_orcid(value: str) -> bool:
    if not ORCID.fullmatch(value):
        return False
    compact = value.replace("-", "")
    total = 0
    for digit in compact[:-1]:
        total = (total + int(digit)) * 2
    check = (12 - total % 11) % 11
    expected = "X" if check == 10 else str(check)
    return compact[-1] == expected


def check_author_metadata(metadata: dict, blockers: list[str]) -> None:
    authors = metadata.get("authors", [])
    expected_orcid = "0009-0005-1859-5107"
    if len(authors) != 1:
        blockers.append("Release metadata must contain exactly one author.")
        return
    author = authors[0]
    if author.get("name") != "Vivaan Patni" or author.get("zenodo_name") != "Vivaan Patni":
        blockers.append("Release and Zenodo author names must be Vivaan Patni.")
    if author.get("github") != "vxxqv":
        blockers.append("The GitHub account must be vxxqv.")
    if author.get("orcid") != expected_orcid or not valid_orcid(author.get("orcid", "")):
        blockers.append("The release ORCID is missing or invalid.")

    zenodo = read_json(ROOT / ".zenodo.json")
    creators = zenodo.get("creators", [])
    if len(creators) != 1 or creators[0].get("name") != "Vivaan Patni" or creators[0].get("orcid") != expected_orcid:
        blockers.append("Zenodo creator metadata differs from the approved author record.")

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    required_lines = (
        'version: "1.0.0"',
        'family-names: "Patni"',
        'given-names: "Vivaan"',
        f'orcid: "https://orcid.org/{expected_orcid}"',
    )
    if any(line not in citation for line in required_lines):
        blockers.append("CITATION.cff differs from the approved version or author record.")


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


def tracked_files() -> set[str]:
    git = os.environ.get("OPENOMICSBENCH_GIT", "git")
    output = subprocess.run(
        [git, "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "-z"],
        cwd=ROOT, check=True, capture_output=True,
    ).stdout.decode("utf-8")
    return {value for value in output.split("\0") if value}


def tracked_blob(relative: str) -> bytes:
    git = os.environ.get("OPENOMICSBENCH_GIT", "git")
    return subprocess.run(
        [git, "-c", f"safe.directory={ROOT.as_posix()}", "show", f":{relative}"],
        cwd=ROOT, check=True, capture_output=True,
    ).stdout


def check_release_inventory(metadata: dict, blockers: list[str]) -> None:
    path = ROOT / INVENTORY_PATH
    if not path.is_file():
        blockers.append("The release file inventory is missing.")
        return
    inventory = read_json(path)
    if inventory.get("version") != metadata.get("software_version"):
        blockers.append("The release inventory version differs from the software version.")
    if inventory.get("created") != metadata.get("candidate_date"):
        blockers.append("The release inventory date differs from the candidate date.")
    records = inventory.get("files", [])
    declared = [record.get("path") for record in records]
    if inventory.get("file_count") != len(records) or len(declared) != len(set(declared)):
        blockers.append("The release inventory count is wrong or contains duplicate paths.")
        return
    expected = tracked_files() - {INVENTORY_PATH}
    if set(declared) != expected:
        missing = sorted(expected - set(declared))
        extra = sorted(set(declared) - expected)
        blockers.append(f"The release inventory path set differs: missing={missing}, extra={extra}.")
        return
    for record in records:
        relative = record["path"]
        blob = tracked_blob(relative)
        if len(blob) != record.get("bytes") or hashlib.sha256(blob).hexdigest() != record.get("sha256"):
            blockers.append(f"Release inventory hash or byte count differs: {relative}.")


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
        "candidate_date",
        "repository_url",
        "code_license",
        "metadata_license",
        "scientific_environment_lock",
    )
    for field in required_metadata:
        if not metadata.get(field):
            blockers.append(f"Missing release evidence: {field}.")
    check_author_metadata(metadata, blockers)
    lock = ROOT / metadata.get("scientific_environment_lock", "")
    if not lock.is_file():
        blockers.append("The scientific environment lock is missing.")
    check_release_inventory(metadata, blockers)

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
