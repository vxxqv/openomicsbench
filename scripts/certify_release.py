"""Verify the complete collection and write a fail-closed release decision."""
import argparse
import hashlib
import json
import os
import re
import subprocess
from functools import lru_cache
from pathlib import Path

from omicsbench.hashing import digest
from omicsbench.registry import registry
from omicsbench.validate import validate


ROOT = Path(__file__).resolve().parents[1]
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
ORCID = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$")
ZENODO_DOI = re.compile(r"^10\.5281/zenodo\.\d+$")
PUBLICATION_FIELDS = ("version_doi", "external_user_trial", "post_upload_verification", "pypi_verification")
INVENTORY_PATH = "release/file-manifest.json"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


@lru_cache(maxsize=None)
def historical_digest(commit: str, relative: str) -> str | None:
    if not COMMIT.fullmatch(commit):
        return None
    git = os.environ.get("OPENOMICSBENCH_GIT", "git")
    result = subprocess.run(
        [git, "-c", f"safe.directory={ROOT.as_posix()}", "show", f"{commit}:{relative}"],
        cwd=ROOT, capture_output=True,
    )
    return hashlib.sha256(result.stdout).hexdigest() if result.returncode == 0 else None


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

    release_commit = metadata.get("release_commit")
    git = os.environ.get("OPENOMICSBENCH_GIT", "git")
    citation = (
        subprocess.run(
            [git, "-c", f"safe.directory={ROOT.as_posix()}", "show", f"{release_commit}:CITATION.cff"],
            cwd=ROOT, check=True, capture_output=True, text=True,
        ).stdout
        if release_commit else (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    )
    required_lines = (
        f'version: "{metadata["software_version"]}"',
        'family-names: "Patni"',
        'given-names: "Vivaan"',
        f'orcid: "https://orcid.org/{expected_orcid}"',
    )
    if any(line not in citation for line in required_lines):
        blockers.append("CITATION.cff differs from the approved version or author record.")


def check_publication_metadata(metadata: dict, blockers: list[str]) -> None:
    trial = metadata.get("external_user_trial")
    if trial:
        valid_statuses = {"independent trial", "owner-authorized clean-room substitute"}
        if trial.get("status") not in valid_statuses:
            blockers.append("The quickstart trial status is not recognized.")
        if trial.get("outcome") != "pass" or not COMMIT.fullmatch(trial.get("commit", "")):
            blockers.append("The quickstart trial outcome or commit is invalid.")
        if not trial.get("date") or not trial.get("operating_system") or not trial.get("python"):
            blockers.append("The quickstart trial environment is incomplete.")
        if len(trial.get("commands", [])) < 7:
            blockers.append("The quickstart trial does not record enough release commands.")
        if trial.get("status") == "owner-authorized clean-room substitute" and not trial.get("limitation"):
            blockers.append("The clean-room substitute must state its independence limitation.")

    version_doi = metadata.get("version_doi")
    if version_doi and not ZENODO_DOI.fullmatch(version_doi):
        blockers.append("The Zenodo version DOI is invalid.")

    release_commit = metadata.get("release_commit")
    if release_commit:
        git = os.environ.get("OPENOMICSBENCH_GIT", "git")
        if not COMMIT.fullmatch(release_commit):
            blockers.append("The release commit is invalid.")
        else:
            tag = f"v{metadata['software_version']}"
            tagged_commit = subprocess.run(
                [git, "-c", f"safe.directory={ROOT.as_posix()}", "rev-list", "-n", "1", tag],
                cwd=ROOT, check=True, capture_output=True, text=True,
            ).stdout.strip()
            if tagged_commit != release_commit:
                blockers.append(f"The {tag} tag does not identify the recorded release commit.")

    upload = metadata.get("post_upload_verification")
    if upload:
        if upload.get("status") != "pass" or upload.get("version_doi") != version_doi:
            blockers.append("The post-upload verification does not match the version DOI.")
        if not upload.get("record_url") or not upload.get("verified_at"):
            blockers.append("The post-upload verification record is incomplete.")
        if upload.get("git_commit") != release_commit or upload.get("git_tag") != f"v{metadata['software_version']}":
            blockers.append("The post-upload verification does not identify the recorded release tag.")
        if upload.get("archive_content_match") is not True:
            blockers.append("The GitHub and Zenodo archive contents have not been matched.")
        if upload.get("missing_files") or upload.get("extra_files") or upload.get("changed_files"):
            blockers.append("The post-upload comparison contains file differences.")
        if upload.get("archived_files") != upload.get("inventory_files_verified", 0) + 1:
            blockers.append("The post-upload file counts differ from the self-excluding inventory.")
        for field in ("zenodo_sha256", "github_sha256"):
            if not SHA256.fullmatch(upload.get(field, "")):
                blockers.append(f"The recorded {field} is invalid.")
        for field in ("zenodo_md5", "github_md5"):
            if not MD5.fullmatch(upload.get(field, "")):
                blockers.append(f"The recorded {field} is invalid.")

    pypi = metadata.get("pypi_verification")
    if pypi:
        version = metadata["software_version"]
        if pypi.get("status") != "pass" or pypi.get("version") != version:
            blockers.append("The PyPI verification does not match the release version.")
        if pypi.get("url") != f"https://pypi.org/project/openomicsbench/{version}/":
            blockers.append("The PyPI release URL is missing or incorrect.")
        files = pypi.get("files", [])
        expected = {f"openomicsbench-{version}-py3-none-any.whl", f"openomicsbench-{version}.tar.gz"}
        if {item.get("filename") for item in files} != expected or len(files) != 2:
            blockers.append("PyPI verification must name the wheel and source distribution.")
        for item in files:
            if not isinstance(item.get("bytes"), int) or item["bytes"] <= 0 or not SHA256.fullmatch(item.get("sha256", "")):
                blockers.append("A PyPI file has an invalid size or SHA-256 digest.")
            if item.get("file_set_match") is not True or item.get("normalized_content_match") is not True:
                blockers.append("A published PyPI file has not been matched with the tested build.")


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
    commit = record.get("workflow_commit")
    recorded_digest = historical_digest(commit, relative) if commit else None
    if not expected or not SHA256.fullmatch(expected) or (digest(path) != expected and recorded_digest != expected):
        blockers.append(f"{label}: hash does not match {relative}.")
    if commit is not None and not COMMIT.fullmatch(commit):
        blockers.append(f"{label}: transformation commit is invalid.")


def release_files(metadata: dict) -> set[str]:
    git = os.environ.get("OPENOMICSBENCH_GIT", "git")
    commit = metadata.get("release_commit")
    command = [git, "-c", f"safe.directory={ROOT.as_posix()}"]
    command += ["ls-tree", "-r", "--name-only", "-z", commit] if commit else ["ls-files", "-z"]
    output = subprocess.run(
        command,
        cwd=ROOT, check=True, capture_output=True,
    ).stdout.decode("utf-8")
    return {value for value in output.split("\0") if value}


def release_blob(relative: str, metadata: dict) -> bytes:
    git = os.environ.get("OPENOMICSBENCH_GIT", "git")
    reference = metadata.get("release_commit") or ""
    object_name = f"{reference}:{relative}" if reference else f":{relative}"
    return subprocess.run(
        [git, "-c", f"safe.directory={ROOT.as_posix()}", "show", object_name],
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
    expected = release_files(metadata) - {INVENTORY_PATH}
    if set(declared) != expected:
        missing = sorted(expected - set(declared))
        extra = sorted(set(declared) - expected)
        blockers.append(f"The release inventory path set differs: missing={missing}, extra={extra}.")
        return
    for record in records:
        relative = record["path"]
        blob = release_blob(relative, metadata)
        if len(blob) != record.get("bytes") or hashlib.sha256(blob).hexdigest() != record.get("sha256"):
            blockers.append(f"Release inventory hash or byte count differs: {relative}.")


def check_biological_object(model, folder: Path, blockers: list[str]) -> None:
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
    provenance = read_json(provenance_path)
    original_commit = provenance.get("assembly", {}).get("workflow_commit")
    for field, path in workflow_files.items():
        expected_hash = workflow_hashes.get(field)
        relative = path.relative_to(ROOT).as_posix()
        if expected_hash != digest(path) and expected_hash != historical_digest(original_commit or "", relative):
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
    reference_record = provenance.get("reference", {})
    reference_relative = reference_record.get("path", "evidence/reference-verification.json")
    reference_path_global = ROOT / reference_relative
    if not reference_path_global.is_file():
        blockers.append(f"{model.id}: collection reference evidence is missing: {reference_relative}.")
        return
    object_reference_evidence = read_json(reference_path_global)
    if object_reference_evidence.get("status") != "pass":
        blockers.append(f"{model.id}: collection reference verification has not passed.")
    if reference.get("collection_evidence_sha256") != digest(reference_path_global):
        blockers.append(f"{model.id}: collection reference evidence hash has drifted.")
    if object_key not in object_reference_evidence.get("verified_objects", []):
        blockers.append(f"{model.id}: object is absent from the verified reference set.")

    if provenance.get("dataset_id") != model.id or provenance.get("object_key") != object_key:
        blockers.append(f"{model.id}: provenance identity differs.")
    check_recorded_file(provenance.get("selection", {}), f"{model.id} selection", blockers)
    check_recorded_file(provenance.get("assembly", {}), f"{model.id} assembly", blockers)
    for name in ("diagnostic", "deseq2", "reference", "intake"):
        check_recorded_file(provenance.get(name, {}), f"{model.id} {name}", blockers)


def check_real_single_cell(model, folder: Path, blockers: list[str], source_audit: dict,
                           reference_audit: dict) -> None:
    """Check the pinned source and annotation lineage of a derived cell case."""
    provenance_path = folder / "provenance/source-and-transform.json"
    license_path = folder / "provenance/source-license.txt"
    if not provenance_path.is_file() or not license_path.is_file():
        blockers.append(f"{model.id}: source provenance or license notice is missing.")
        return
    provenance = read_json(provenance_path)
    stem = provenance.get("source_stem")
    sources = {case["stem"]: case for case in source_audit.get("cases", [])}
    references = {case["stem"]: case for case in reference_audit.get("cases", [])}
    if stem not in sources or stem not in references:
        blockers.append(f"{model.id}: source is absent from pinned intake audits.")
        return
    source = sources[stem]
    reference = references[stem]
    if provenance.get("source_commit") != source_audit.get("upstream_commit"):
        blockers.append(f"{model.id}: source revision differs from intake audit.")
    if provenance.get("source_count_sha256") != source["source_files"]["counts"]["sha256"] or model.source.sha256 != provenance.get("source_count_sha256"):
        blockers.append(f"{model.id}: source count hash differs from intake audit.")
    if provenance.get("source_metadata_sha256") != source["source_files"]["metadata"]["sha256"]:
        blockers.append(f"{model.id}: source metadata hash differs from intake audit.")
    if provenance.get("reference_gtf_sha256") != reference_audit.get("annotations", {}).get("ensembl88", {}).get("sha256"):
        blockers.append(f"{model.id}: reference annotation hash differs from intake audit.")
    if sorted(provenance.get("excluded_unmatched_feature_ids", [])) != sorted(reference.get("unmatched_examples", [])) or reference.get("annotation_unmatched", 999) > 25:
        blockers.append(f"{model.id}: excluded feature IDs differ from reference audit.")
    if provenance.get("selected_features") != reference.get("feature_count", 0) - reference.get("annotation_unmatched", 0):
        blockers.append(f"{model.id}: selected feature count differs from reference audit.")
    if provenance.get("workflow_commit") != model.derivation.commit or not COMMIT.fullmatch(model.derivation.commit or ""):
        blockers.append(f"{model.id}: transformation revision is missing or inconsistent.")
    elif digest(ROOT / model.derivation.workflow) != historical_digest(model.derivation.commit, model.derivation.workflow):
        blockers.append(f"{model.id}: transformation workflow differs from pinned revision.")
    if digest(license_path) != digest(ROOT / "workflows/cellbench-source-license.txt"):
        blockers.append(f"{model.id}: upstream license notice differs from pinned copy.")
    if model.rights.status != "GREEN" or model.rights.license != "MIT" or model.reference.compatibility != "verified":
        blockers.append(f"{model.id}: rights or reference compatibility are not certified.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify the collection and write its release decision.")
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
    sequence_validated = []
    single_cell_validated = []
    real_single_cell_validated = []
    interval_validated = []

    if reference_evidence.get("status") != "pass":
        blockers.append("Collection reference verification has not passed.")
    if audit.get("status") != "pass":
        blockers.append("Collection overlap and prose audit has not passed.")

    try:
        records = registry(ROOT)
    except ValueError as exc:
        records = {}
        blockers.append(str(exc))
    major_version = int(metadata["software_version"].split(".")[0])
    source_audit = read_json(ROOT / "curation/v15-cellbench-source-audit.json") if major_version >= 15 else {}
    reference_audit = read_json(ROOT / "curation/v15-cellbench-reference-audit.json") if major_version >= 15 else {}
    for model, folder in records.values():
        try:
            result = validate(model, folder)
        except ValueError as exc:
            blockers.append(str(exc))
            continue
        if model.kind == "real" and model.assay == "bulk_rna_seq" and model.status == "validated" and result["status"] == "pass":
            validated.append(model)
            check_biological_object(model, folder, blockers)
        if model.sequence is not None and model.status == "validated" and result["status"] == "pass":
            sequence_validated.append(model)
        if model.single_cell is not None and model.status == "validated" and result["status"] == "pass":
            single_cell_validated.append(model)
            if model.kind == "real":
                real_single_cell_validated.append(model)
                if major_version >= 15:
                    check_real_single_cell(model, folder, blockers, source_audit, reference_audit)
        if model.intervals is not None and model.status == "validated" and result["status"] == "pass":
            interval_validated.append(model)

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
    if len(sequence_validated) < 4:
        blockers.append("The release requires at least four validated sequence objects.")
    expected_sequence_families = {"sequence_dna", "sequence_rna", "sequence_protein", "short_read_dna", "whole_genome_dna_seq", "shotgun_metagenomics"}
    if not expected_sequence_families.issubset({model.assay for model in sequence_validated}):
        blockers.append("The DNA, RNA, protein, paired-read and DNA-seq sequence families are incomplete.")
    if not single_cell_validated:
        blockers.append("The release requires a validated single-cell object.")
    if major_version >= 15 and len(real_single_cell_validated) < 6:
        blockers.append(f"Version 15 requires six sourced single-cell cases; found {len(real_single_cell_validated)}.")
    if not interval_validated:
        blockers.append("The release requires a validated genomic interval object.")
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
    check_publication_metadata(metadata, blockers)
    lock = ROOT / metadata.get("scientific_environment_lock", "")
    if not lock.is_file():
        blockers.append("The scientific environment lock is missing.")
    check_release_inventory(metadata, blockers)

    deferred = [field for field in PUBLICATION_FIELDS if not metadata.get(field)]
    report = {
        "scope": f"v{metadata['software_version'].split('.')[0]}_collection",
        "software_version": metadata["software_version"],
        "decision": "GO" if not blockers and not deferred else "NO-GO",
        "preflight": "PASS" if not blockers else "FAIL",
        "certified_biological_objects": len(validated),
        "certified_sequence_objects": len(sequence_validated),
        "certified_single_cell_objects": len(single_cell_validated),
        "certified_real_single_cell_objects": len(real_single_cell_validated),
        "certified_interval_objects": len(interval_validated),
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
