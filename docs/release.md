# Release procedure

The development archive is a handoff for continued work. Do not tag it v1.0 or upload it as the completed biological collection.

## Establish the publication inputs

Record the project owner, credited authors and any ORCIDs they provide. Confirm the code license and the license for project-written metadata. Choose the GitHub repository and Zenodo account. Populate `release/metadata-input.json` with those reviewed values. A DOI must come from the actual deposit; do not invent one to make the citation file look complete.

The publication metadata generator should create `CITATION.cff` and `.zenodo.json` from the same reviewed inputs. Zenodo gives `.zenodo.json` precedence when both exist. This snapshot retains the unresolved inputs instead of supplying invalid citation files or fictional contributors.

## Complete biological certification

Resolve rights, sample-to-accession mapping and compatible references before distributing each source. Build 12 to 20 meaningful objects, with at least eight GREEN pockets and at least four distinct design archetypes. Do not count multiple tiers of the same object as different datasets or count synthetic fixtures toward this requirement.

Freeze baseline versions, parameters, thresholds and candidate-size policies. Execute the scientific workflow in an isolated environment, record a complete dependency lock and container digest where used, and rebuild several objects from a clean environment. Recompute every declared metric from the certified inputs. Preserve logs for failures as well as successes.

## Review the user path

Ask an independent person to complete the quickstart from a fresh environment. Record their consented feedback, environment, commands and outcome. Fix any hidden maintainer step. Review every dataset card and figure at normal reading size. Check links and narrow-width command rendering in the actual publication surfaces.

## Assemble and publish

Run `python scripts/certify_release.py`. The command must exit successfully before creating a v1.0 release tag. It currently exits with a no-go report. Store the final registry, release manifest, checksums, changelog and citation metadata together after all gates are satisfied.

Reserve the real version DOI, insert it into the reviewed citation inputs, rebuild the metadata and freeze hashes of the final upload files. Publish only after the owner has approved that concrete bundle. Retrieve the deposited files and compare their hashes with the certification inventory. Check that the DOI resolves to the exact version and that the displayed citation names the correct authors.

## Preserve evidence of use

Record actual requests, issues, workshops, external contributions and maintenance decisions. Aggregate usage signals without collecting personal student profiles. Views, downloads and stars can describe interest; they do not establish scientific impact. Write the retrospective from observed events and specific changes needed for v2.
