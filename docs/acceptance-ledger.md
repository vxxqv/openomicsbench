# Acceptance ledger

This ledger maps the version 1 acceptance tests in OmicsBenchPlan_v1.0.docx to current evidence. Verified means the named evidence exists and passes. Partial marks a bounded result with a remaining condition. Open marks work that depends on an event that has not happened.

## 1. Release thesis and quality bar

- **VERIFIED** Project charter defines audience, inclusion criteria, object tiers, and exclusions. README.md and docs/methods.md define the user, bulk RNA-seq scope, object tiers, admission tests and exclusions.
- **VERIFIED** No dataset merges without provenance and rights review. The collection plan links every biological object to dated intake, rights, reference and transformation records; the overlap audit passes.
- **VERIFIED** v1 scope remains transcriptomics plus core infrastructure. The registry contains bulk RNA-seq objects and two infrastructure fixtures; no other assay is represented.
- **VERIFIED** Documentation style rules live in CONTRIBUTING.md and are reviewed like code. CONTRIBUTING.md states the writing rules; CI checks local links, unfinished text, punctuation and attribution traces.

## 2. Repository architecture and data-object contract

- **VERIFIED** Fresh clone can enumerate all datasets from manifests only. registry.py discovers manifests without a hand-maintained dataset list; the clean temporary-repository test passes.
- **VERIFIED** Validator catches duplicate IDs, missing files, invalid hashes, and undeclared roles. Regression tests cover duplicate identifiers and roles, missing and changed files, bad hashes and undeclared paths.
- **VERIFIED** Dataset-specific logic is isolated in adapters/workflows rather than scattered through core code. Source retrieval and biological selection live under workflows; src/omicsbench contains source-neutral models, registry, transfer, count and validation code.

## 3. Manifest schema, identifiers, and provenance

- **VERIFIED** Schema round-trip tests preserve values. test_manifest_round_trip passes against the current strict Pydantic contract.
- **VERIFIED** Provenance report can be generated without reading Python source. The provenance command reads the manifest and transformation record and is exercised by the command-line test.
- **VERIFIED** Help text for fields is written as human guidance, not terse schema jargon. Every public manifest field has a plain-language description in the generated JSON Schema.
- **VERIFIED** Schema changes require tests and migration notes. CONTRIBUTING.md requires both; schema drift and manifest round trips are checked automatically.

## 4. Source discovery, licensing, and intake

- **VERIFIED** Every candidate has an intake record before processing starts. The curation intake directory records every evaluated archive candidate, including excluded and deferred sources.
- **VERIFIED** Rights decision includes evidence and date checked. Every biological object carries rights.json and manifest rights fields with the licence, evidence URL, decision and review date.
- **VERIFIED** Link-only status is visible in catalog and CLI. rnaseq-001 is marked original_link in its manifest and command output; retrieval refuses to package it.
- **VERIFIED** No source is called open simply because a browser download exists. GEO records remain AMBER without a specific redistribution grant; downloadable material is not packaged on access alone.

## 5. Bulk RNA-seq collection design

- **VERIFIED** At least four distinct RNA design archetypes are represented. The 12 certified objects declare 12 design descriptions spanning balanced, paired, blocked factorial, imbalanced and stratified comparisons.
- **VERIFIED** Counts are never overwritten by normalized values. Source integer counts and pocket integer counts remain separate; normalization is computed only during validation.
- **VERIFIED** Strandedness and sample design metadata are explicit. Every sample record declares condition, replicate, batch or blocking value and strandedness; unavailable strandedness is stated as unknown.
- **VERIFIED** Reference/annotation compatibility is checked automatically. verify_references.py checks official annotation files and pocket identifiers; certification rejects missing or unresolved evidence.

## 6. Deterministic pocketing algorithms

- **VERIFIED** All stochastic reducers accept and record a seed. Reducers accept a seed and every manifest derivation records it, including zero for deterministic feature selection.
- **VERIFIED** Paired reads and sample identity are tested. Regression tests cover paired-mate identity, duplicate runs, sample order and design-to-matrix agreement.
- **VERIFIED** Pocket size is justified quantitatively. Each object retains its candidate curve and selects the smallest feature count that meets all four predeclared preservation thresholds.
- **VERIFIED** No spreadsheet-only/manual cleanup is required for regeneration. Intake, selection, DESeq2 comparison, reference checking, assembly, catalog and certification are scripted.

## 7. Expected outputs and quantitative validation

- **PARTIAL** Every dataset has hard and quantitative validation. All 12 distributed biological objects and the synthetic fixture have integrity and metric checks. The deliberately link-only rnaseq-001 record carries no redistributed matrix to score.
- **VERIFIED** Metric definitions include filtering universe and direction. Every validation profile states the feature universe, calculation and inclusive minimum for all four metrics.
- **VERIFIED** Expected outputs rebuild from a pinned baseline. All 12 DESeq2 evidence records were rebuilt under R 4.5.3 and DESeq2 1.50.2 with pinned inputs, hashes and one workflow revision.
- **VERIFIED** Figures are never the only location where expected values exist. Exact results and thresholds live in expected JSON, the catalog JSON and catalog TSV; the SVG is generated from those records.

## 8. Python CLI and API architecture

- **VERIFIED** CLI has unit and end-to-end tests on nano objects. All seven public commands are covered by the 45-test suite, including the generated fixture path.
- **VERIFIED** All downloads verify hashes. Transfers verify declared size and SHA-256 before atomic replacement; interrupted and corrupt transfers have regression tests.
- **VERIFIED** A user can discover, retrieve, inspect, and validate without reading source code. The documented command path passed from a fresh Windows checkout under Python 3.12.14 at commit feb2e1ef0c9df24ecb4b45c0ee17c16baa2680ee.
- **VERIFIED** Help and error messages are concise and human. Every command has focused help and an example; expected user failures return short messages without tracebacks.

## 9. Reproducible workflows and environments

- **PARTIAL** Several datasets regenerate end-to-end from a clean environment. A clean Ubuntu workflow validates all 14 registered records and regenerates the fixture, catalog and figure. It does not redownload and rebuild several biological sources, so this remains partial.
- **VERIFIED** Resolved parameters and versions are preserved. Manifests, transformation records, DESeq2 evidence and runtime/environment-lock.json preserve inputs, parameters, code hashes and software versions.
- **VERIFIED** No notebook-only hidden step is needed for release artifacts. The repository contains scripted workflows for every generated object, catalog, figure, audit and certification output.
- **PARTIAL** Failures retain logs and a short human summary. Rejected pocket and reference candidates retain structured evidence and explanations. A general failure-log contract is not enforced for every workflow.

## 10. CI, tests, and integrity gates

- **VERIFIED** Nano CI stays within a defined time budget. GitHub Actions run 34308797805 completed the five-minute job in 35 seconds, including all 45 tests and the full preflight.
- **VERIFIED** Release tag is blocked by schema/hash/registry failures. The integrity workflow runs on every push and tag and includes schema drift, manifest hashes, registry validation, overlap audit and preflight.
- **VERIFIED** Corruption and sample-order failures have regression tests. Dedicated corruption, missing-file, duplicate, sample-order and design-mismatch tests pass.
- **VERIFIED** Documentation quality is reviewed beyond spelling. The repository checker covers links and unfinished text; the collection audit detects repeated long prose; the fidelity figure received a rendered layout review.

## 11. Human documentation and code-comment policy

- **VERIFIED** Every public command has concise help and one example. list, info, get, validate, provenance, verify-cache and doctor each expose an example in command help.
- **VERIFIED** Every dataset README is usable without opening YAML. Each of the 12 biological cards states its design, file contents, commands, measured result, limits, citation and rights location.
- **VERIFIED** Release notes are manually edited for flow. CHANGELOG.md describes the finished v1 collection by user impact, scientific scope and validation evidence.
- **VERIFIED** No em dash appears in docs/help/release notes. The repository checker scans every Markdown file and the current documentation passes.
- **VERIFIED** Review removes long narration comments. A source-wide comment and docstring review found only concise module descriptions and constraint explanations.

## 12. Figures, catalog, and rendering QA

- **VERIFIED** All core figures regenerate from code and source tables. scripts/build_catalog.py regenerates figures/v1-fidelity.svg directly from manifests and expected-result records; CI checks drift.
- **VERIFIED** Captions are concise and units are present where applicable. The fidelity figure labels each unitless preservation metric, threshold and object; its README caption states the comparison.
- **PARTIAL** Rendered documentation is manually inspected. The full-resolution fidelity preview has been checked for clipping, collisions, labels and value placement. GitHub and Zenodo rendering remain part of publication verification.
- **VERIFIED** No visual hides rights or warning status. The figure is explicitly scoped to the 12 certified biological objects; the adjacent catalog and README state rights and release status.

## 13. Zenodo release engineering, CAS milestones, and final go/no-go

- **PARTIAL** One external user completes quickstart from a fresh environment. The release preparer completed the full quickstart from a fresh checkout at commit feb2e1ef0c9df24ecb4b45c0ee17c16baa2680ee. No independent tester is claimed; the owner explicitly authorized this substitute before publication.
- **VERIFIED** Final uploaded hashes match certified files. The GitHub and Zenodo archives contain the same 282 files; all 281 inventory-controlled files match byte counts and SHA-256 values, with no missing, extra or changed files.
- **VERIFIED** DOI/citation instructions point to exact version. OpenOmicsBench 1.0.0 is archived at DOI 10.5281/zenodo.22679414; the concept DOI remains 10.5281/zenodo.22551734.
- **PARTIAL** v1 retrospective identifies concrete improvements for v2. Development lessons and candidate v2 work are recorded, but post-release user experience and archive behavior cannot be assessed before publication.
