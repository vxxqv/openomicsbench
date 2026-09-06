# Acceptance ledger

Each acceptance item below comes from OmicsBenchPlan_v1.0.docx. Verified describes the narrow evidence stated; it does not imply that v1.0 is ready. Partial and open items remain release work.

## 1. Release thesis and quality bar

- **OPEN** Project charter defines audience, inclusion criteria, object tiers, and exclusions. Not yet demonstrated for the requested v1 biological release.
- **OPEN** No dataset merges without provenance and rights review. Not yet demonstrated for the requested v1 biological release.
- **VERIFIED** v1 scope remains transcriptomics plus core infrastructure. The intake excludes single-cell data; code supports bulk count objects and infrastructure fixtures.
- **VERIFIED** Documentation style rules live in CONTRIBUTING.md and are reviewed like code. CONTRIBUTING.md and documentation lint are present; no external editorial review is claimed.

## 2. Repository architecture and data-object contract

- **VERIFIED** Fresh clone can enumerate all datasets from manifests only. registry.py enumerates local manifests; CLI list passed in a temporary repository.
- **VERIFIED** Validator catches duplicate IDs, missing files, invalid hashes, and undeclared roles. Covered by schema, inventory, duplicate-ID and corruption tests.
- **OPEN** Dataset-specific logic is isolated in adapters/workflows rather than scattered through core code. Not yet demonstrated for the requested v1 biological release.

## 3. Manifest schema, identifiers, and provenance

- **VERIFIED** Schema round-trip tests preserve values. test_manifest_round_trip passed.
- **OPEN** Provenance report can be generated without reading Python source. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Help text for fields is written as human guidance, not terse schema jargon. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Schema changes require tests and migration notes. Not yet demonstrated for the requested v1 biological release.

## 4. Source discovery, licensing, and intake

- **VERIFIED** Every candidate has an intake record before processing starts. Twelve intake records were written before running the Pasilla diagnostic reduction.
- **OPEN** Rights decision includes evidence and date checked. Not yet demonstrated for the requested v1 biological release.
- **VERIFIED** Link-only status is visible in catalog and CLI. rnaseq-001 has original_link status and AMBER rights; get rejects it.
- **VERIFIED** No source is called open simply because a browser download exists. All GEO intake records remain AMBER.

## 5. Bulk RNA-seq collection design

- **OPEN** At least four distinct RNA design archetypes are represented. Not yet demonstrated for the requested v1 biological release.
- **VERIFIED** Counts are never overwritten by normalized values. Raw-count files are preserved; diagnostic transforms are separate arrays.
- **OPEN** Strandedness and sample design metadata are explicit. Not yet demonstrated for the requested v1 biological release.
- **PARTIAL** Reference/annotation compatibility is checked automatically. The contract blocks unresolved references for validated objects; actual biological reference-file compatibility checks are not implemented.

## 6. Deterministic pocketing algorithms

- **VERIFIED** All stochastic reducers accept and record a seed. Feature and read sampling accept seeds; fixture and Pasilla configurations record them.
- **OPEN** Paired reads and sample identity are tested. Not yet demonstrated for the requested v1 biological release.
- **PARTIAL** Pocket size is justified quantitatively. Pasilla has six measured candidates, but the result is diagnostic-only and not a certified pocket.
- **VERIFIED** No spreadsheet-only/manual cleanup is required for regeneration. Both executed local workflows are scripted.

## 7. Expected outputs and quantitative validation

- **PARTIAL** Every dataset has hard and quantitative validation. The fixture passes both. The real-source record is link-only and is not scientifically certified.
- **PARTIAL** Metric definitions include filtering universe and direction. Current diagnostics define both. Full inferential and read-QC profiles remain open.
- **PARTIAL** Expected outputs rebuild from a pinned baseline. The fixture workflow and full-table Pasilla DESeq2 baseline run in the locked local environment. A clean rebuild on another machine remains open.
- **PARTIAL** Figures are never the only location where expected values exist. Diagnostic tables and JSON back the current plots; the full biological figure set remains open.

## 8. Python CLI and API architecture

- **VERIFIED** CLI has unit and end-to-end tests on nano objects. All seven CLI commands passed on the generated fixture.
- **VERIFIED** All downloads verify hashes. transfer verifies size and SHA-256 before atomic replacement; interrupted and corrupt transfers tested.
- **PARTIAL** A user can discover, retrieve, inspect, and validate without reading source code. The fixture path works. The independent-user trial has not occurred.
- **VERIFIED** Help and error messages are concise and human. Command help includes examples and expected failures omit tracebacks by default.

## 9. Reproducible workflows and environments

- **OPEN** Several datasets regenerate end-to-end from a clean environment. Not yet demonstrated for the requested v1 biological release.
- **PARTIAL** Resolved parameters and versions are preserved. Executed Python runs record versions and parameters; R, workflow-engine and container certification remain open.
- **OPEN** No notebook-only hidden step is needed for release artifacts. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Failures retain logs and a short human summary. Not yet demonstrated for the requested v1 biological release.

## 10. CI, tests, and integrity gates

- **PARTIAL** Nano CI stays within a defined time budget. The local 28-test run completed in under ten seconds; a five-minute CI job is authored but has not run remotely.
- **OPEN** Release tag is blocked by schema/hash/registry failures. Not yet demonstrated for the requested v1 biological release.
- **VERIFIED** Corruption and sample-order failures have regression tests. Both failures have passing regression tests.
- **OPEN** Documentation quality is reviewed beyond spelling. Not yet demonstrated for the requested v1 biological release.

## 11. Human documentation and code-comment policy

- **OPEN** Every public command has concise help and one example. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Every dataset README is usable without opening YAML. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Release notes are manually edited for flow. Not yet demonstrated for the requested v1 biological release.
- **VERIFIED** No em dash appears in docs/help/release notes. The repository checker scans documentation; help text has been reviewed.
- **OPEN** Review removes long narration comments. Not yet demonstrated for the requested v1 biological release.

## 12. Figures, catalog, and rendering QA

- **OPEN** All core figures regenerate from code and source tables. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Captions are concise and units are present where applicable. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Rendered documentation is manually inspected. Not yet demonstrated for the requested v1 biological release.
- **OPEN** No visual hides rights or warning status. Not yet demonstrated for the requested v1 biological release.

## 13. Zenodo release engineering, CAS milestones, and final go/no-go

- **OPEN** One external user completes quickstart from a fresh environment. Not yet demonstrated for the requested v1 biological release.
- **OPEN** Final uploaded hashes match certified files. Not yet demonstrated for the requested v1 biological release.
- **OPEN** DOI/citation instructions point to exact version. Not yet demonstrated for the requested v1 biological release.
- **OPEN** v1 retrospective identifies concrete improvements for v2. Not yet demonstrated for the requested v1 biological release.

