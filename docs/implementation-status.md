# Version 1 status

The v1 repository preflight passes. It finds 12 certified biological objects, 12 redistributable pockets and no unresolved scientific or structural blocker. The publication decision remains NO-GO until the independent quickstart trial, exact version DOI and post-upload file check are recorded.

## 1. Scope and quality bar

OpenOmicsBench is limited to bulk RNA-seq count-matrix objects and the software needed to inspect, retrieve and validate them. The README defines the audience and exclusions. Biological data enter the collection only after source, rights, design, reference and quantitative review.

## 2. Repository and object contract

Core models, registry handling, hashing, caching, count operations, validation and command handling are separate modules. Source-specific work stays in versioned workflow scripts and configuration files. Every dataset folder has one manifest that declares its complete file inventory, roles, tiers, byte counts and SHA-256 hashes.

## 3. Manifests and provenance

Pydantic models reject unknown fields and generate the checked JSON Schema. Manifests record source citation, dated rights evidence, organism, taxonomy, ordered samples, reference compatibility, derivation parameters, software versions and validation minimums. The CLI can report provenance without reading Python source. Biological objects also include a detailed JSON transformation record with workflow hashes and commits.

## 4. Intake and rights

Every evaluated source has a dated intake record. Expression Atlas objects use a documented CC BY 4.0 decision with provider attribution. GEO candidates remain AMBER and are not redistributed. E-MTAB-5477 and E-MTAB-7126 retain their rejection evidence so failed candidates are not silently recycled.

## 5. Biological collection

The collection has 12 objects from seven studies, three organisms and five annotation profiles. It covers balanced knockouts, a paired tumour design, blocked factorial comparisons, RNA interference, an imbalanced disease comparison and genotype-specific strata. Sample designs record condition, replicate, blocking value and strandedness. Unknown strandedness is stated rather than inferred.

## 6. Deterministic pocketing

Candidate grids, selection seeds, contrast fields and sample subsets are declared in configuration. All sample columns and integer counts are preserved. The smallest candidate meeting all four predeclared minimums is selected. The paired FASTQ helper remains covered by tests but is outside the count-matrix collection scope.

## 7. Expected results

Each biological object has a full selected-sample count matrix, a pocket count matrix, a validation profile and a compact DESeq2 evidence record. Validation first checks file integrity and exact count equality, then recomputes the diagnostic metrics. DESeq2 evidence independently confirms the full-source versus pocket comparison under R 4.5.3 and DESeq2 1.50.2.

## 8. Command line and Python core

The seven public commands are `list`, `info`, `get`, `validate`, `provenance`, `verify-cache` and `doctor`. Each has concise help and an example. Expected failures return a short message and a nonzero exit code. The cache verifies size and SHA-256 before an atomic replacement and records a completion receipt.

## 9. Workflows and environments

Collection intake, diagnostic selection, DESeq2 comparison, reference verification, object assembly, catalog generation, overlap auditing and release certification are scripted. No notebook or spreadsheet step is required. The tested Conda and R package set is preserved in `runtime/environment-lock.json`. A container is not claimed because no container was used.

## 10. Continuous integration

Forty-five unit and end-to-end tests pass locally. CI rebuilds the synthetic fixture and catalog, checks generated-file drift, runs the overlap audit, validates every registered object, scans documentation and metadata, and runs the release preflight. Pushes and tags use the same workflow.

## 11. Documentation

The README provides the working quickstart and complete v1 object list. The methods and source review explain the reduction, thresholds, source decisions, reference checks and limits. The release procedure identifies the evidence needed before and after publication. Automated review catches broken local links, unfinished placeholders, prohibited punctuation and attribution traces.

## 12. Catalog and figure

The catalog is generated from manifests and object evidence in JSON and tabular forms. The black-and-white fidelity figure is generated from the same records. Its four panels show every exact metric against the relevant predeclared minimum. The rendered preview was checked for clipping, collisions, legibility and faithful values.

## 13. Release decision

The machine-readable certification report has no local blockers. It remains NO-GO because three steps depend on the final publication process: an independent person must complete the quickstart, Zenodo must provide the exact v1 DOI, and the published archive must be downloaded and matched to the certified inventory. These fields are left empty until the events occur.
