# Implementation status against the v1 plan

The published prerelease establishes the tested software core and one real-source diagnostic workflow. Version 1.0.0 is in development. Its collection target is 12 to 20 curated RNA objects with at least eight redistributable pockets; the current certified biological count is zero. Twenty source candidates and one synthetic fixture do not satisfy that target.

The [acceptance ledger](acceptance-ledger.md) preserves each acceptance item from the supplied plan. A partial entry means there is useful implementation but the stated release gate has not been demonstrated.

## 1 Release thesis and quality bar

The README and contribution policy define audience, scope, rights review and the distinction between fixtures and biological pockets. Bulk RNA-seq is the first biological assay. Single-cell material is explicitly excluded from the intake shortlist. No leaderboard, clinical claim or expansion to unrelated assays has been introduced.

## 2 Repository architecture and data contract

The source package separates models, registry, file hashing, downloading, cache management, RNA diagnostics and command handling. Dataset IDs live in manifests. The generated registry is a build output. File declarations include path, role, tier, size and SHA-256. The validator rejects missing and undeclared files and duplicate logical roles. Documentation and provenance may have multiple files under distinct paths.

## 3 Manifest and provenance

One Pydantic model generates the JSON Schema. Nested models reject unknown fields. The contract includes dated rights decisions, organism, taxon, sample order, reference status, source citation and derivation. JSON is used instead of YAML for the first editable representation. This is a deliberate format decision, not an assertion that general YAML works.

Local derivations record a workflow file hash and tested package versions. A public workflow commit remains unset until the repository history and ownership are established. The release checker treats that as an open gate. Schema round trips and required field families have regression coverage. A migration has not been needed yet.

## 4 Discovery rights and intake

Twelve GEO series records have dated intake records and source response hashes. Four are prioritized for deeper work; others need design resolution or are deferred. One single-cell record is excluded. Eight Expression Atlas candidates have official resource and design checks; E-MTAB-8572 and E-MTAB-6866 also have full file hashes and repeated reduction diagnostics. The Expression Atlas licence supports GREEN decisions with attribution. The Pasilla source package and annotation were inspected at a pinned commit. No biological data are republished in this archive.

## 5 Bulk RNA collection

The Pasilla diagnostic uses a real seven-library count table. A link-only manifest documents it. The shortlist includes paired donor, factorial, reference-material and non-human designs, but those archetypes are not yet represented by certified distributed objects. Low-signal, transcript-level and read-QC challenge objects still need selection and full construction. The count parser rejects fractional values rather than treating estimated quantifications as raw counts.

## 6 Deterministic reduction

The general feature selector preserves all sample columns, combines four declared feature lists, uses a seed for the background ordering, and preserves original row order. The E-MTAB-8572 extension anchors 100 high-effect genes before filling the same general order. The miniature fixture rebuilds exactly in the tested environment. The paired-read helper checks synchronization and returns actual selected pair counts. It has not yet been used to certify a real read pocket.

Six Pasilla candidate sizes were evaluated. The 8,000-gene candidate is the smallest tested size that met the preset diagnostic envelope. That conclusion applies to this grid and this baseline. It is not proof of the globally smallest useful object or of inferential DESeq2 fidelity.

Seven E-MTAB-8572 candidates were evaluated with median-ratio normalization. The 2,000-gene candidate was the smallest tested size meeting the same four thresholds, and the complete report reproduced exactly. Its final DESeq2 baseline and exact reference-file identities remain open.

The corresponding E-MTAB-6866 grid first passed at 4,000 of 32,833 genes. This adds a plant knockout design, although it is still a diagnostic candidate until its reference files, attribution package and DESeq2 baseline are certified.

## 7 Expected outputs and quantitative validation

The fixture has a machine-readable baseline, exact sample and shape invariants, a candidate curve and four metric definitions. Validation recomputes metrics after checking integrity. The Pasilla diagnostic retains failures and explains the direction of its effect contrast. The pinned DESeq2 baseline has now run on the full Pasilla count table. Independent envelope review and full-universe top-feature recall remain necessary for biological certification.

## 8 Python CLI and API

The seven planned commands work on the local registry: list, info, get, validate, provenance, verify-cache and doctor. argparse provides the interface without adding an unneeded CLI dependency. Command tests check structured output and nonzero error exits. Cache receipts identify the source manifest hash. Interrupted transfers leave no valid completion marker.

The cache is intended for a local trusted workspace. Concurrent writes are not yet stress-tested, and a remote catalog update protocol is not implemented. Distribution of datasets inside the Python wheel is intentionally avoided; the CLI uses the selected repository root.

## 9 Workflows and environments

The fixture and Pasilla diagnostic are executable Python workflows. Snakemake 9.26.1 ran the fixture from generation through validation. The DESeq2 script checked integer input, sample order and model rank, then fitted the full Pasilla table with DESeq2 1.50.2 under R 4.5.3.

The tested Conda and R package versions are recorded in runtime/environment-lock.json. A container image digest and clean rebuild on another machine remain open. No Docker tag or unverified environment digest is used as evidence.

## 10 CI and integrity

Thirty-nine local tests pass. The checks cover invalid schema families, corruption, paired-read mismatches, sample order, interrupted downloads, seeded and contrast-aware selection, median-ratio normalization, Expression Atlas intake, fixture regeneration and every CLI command. The integrity suite at commit 71ce874 passes on GitHub Actions; the expanded suite will run when these commits reach GitHub. Scheduled remote health checks, rotating biological regeneration and a protected release-tag policy remain future improvements.

## 11 Human documentation

The README, candidate card, contribution guide, methods, source review and release procedure describe concrete behavior and limitations. Documentation is checked for broken local file links, unfinished placeholders and prohibited em-dash punctuation. This catches mechanical errors; it does not replace the required external usability review.

## 12 Figures catalog and rendering

The diagnostic review includes a size-versus-fidelity figure generated from the saved metrics, with explicit threshold lines and pass/fail text. Plot data are retained. The local catalog is generated from manifests. A public documentation site, browser inspection of GitHub and Zenodo rendering, and the complete per-dataset figure set remain open.

## 13 Release and CAS evidence

The collection checker produces a machine-readable no-go report for biological certification. It distinguishes fixture tests from biological evidence and requires the collection, licensing, external trial and deposited-file checks before version 1.0.0 can claim a certified biological collection.

The source intake, corrected confounding example, failed pocket sizes, code changes and test results are useful project evidence. Workshop participation, external feedback, impact metrics and reflections over time cannot be fabricated in a local build. They remain future evidence to collect through actual use.
