# Version 1 release procedure

## Current position

OpenOmicsBench 1.0.0 is published from commit 5ce08c6bf57e2886b667b77e8c33b238390c2e6e. The GitHub release is normal rather than prerelease, and the tag-triggered integrity workflow passed. Zenodo archived the release at DOI 10.5281/zenodo.22679414 with Vivaan Patni and ORCID 0009-0005-1859-5107. The GitHub and Zenodo archives contain the same files, and the certified inventory comparison passes.

## Freeze the candidate

1. Confirm that `release/metadata-input.json` names Vivaan Patni for Zenodo and `vxxqv` for GitHub. Add the author ORCID only when it has been supplied by the owner.
2. Run the unit tests, catalog builder, overlap audit, repository checker and release preflight.
3. Regenerate the release file inventory after every content change.
4. Push the candidate commits and wait for the GitHub integrity workflow to pass.
5. Stop changing scientific data, thresholds, source configurations and workflow code unless a failed check requires a new candidate.

## Quickstart trial

Use Python 3.11 or newer in a fresh checkout and run:

```sh
python -m pip install -r requirements-tested.txt
python -m pip install --no-deps -e .
omicsbench list --assay bulk_rna_seq
omicsbench info rnaseq-002
omicsbench validate rnaseq-002
omicsbench get rnaseq-002 --size pocket
omicsbench verify-cache
```

Record the date, operating system, Python version, exact commit, outcome and any point where the written instructions were insufficient. State whether the tester was independent. Do not record a tester's name or contact details unless they agree. Fix any hidden maintainer step and repeat the trial if a fix changes the user path. If the owner authorizes a release-preparer trial instead, record that decision and do not describe the run as independent.

## Publish GitHub version 1

1. Confirm that the repository is public, the Zenodo connection is active and the release candidate workflow is green.
2. Review the generated archive inventory, catalog, citation file, Zenodo metadata and release notes.
3. Create the annotated tag `v1.0.0` at the approved commit and push that tag.
4. Publish a GitHub release titled `OpenOmicsBench 1.0.0`. Do not mark it as a prerelease. State that further improvements are planned without weakening the version 1 status.
5. Wait for the connected Zenodo deposit to appear and confirm that the creator is Vivaan Patni.

Tag creation and release publication are owner-approved actions. They must use the exact reviewed commit.

## Complete Zenodo verification

1. Copy the exact Zenodo version DOI and record URL into the release metadata. Keep the concept DOI for the all-versions citation.
2. Update the GitHub release notes so the recommended citation points to the exact v1 DOI.
3. Download the files exposed by Zenodo. Compare their SHA-256 hashes and byte counts with the certified GitHub archive and release inventory.
4. Confirm that the DOI resolves, the record says version 1.0.0, the author is Vivaan Patni, the licence is correct and the keywords are present.
5. Record the verification result and run `python scripts/certify_release.py` without `--preflight`. A GO decision is valid only after all three publication fields are present and their evidence has been checked.

## Later corrections

Do not move the `v1.0.0` tag after publication. Correct code or metadata on a new patch version. Keep the v1 Zenodo record immutable apart from permitted metadata corrections, and document any correction in the changelog.
