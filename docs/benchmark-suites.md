# Benchmark suite reference

The suite commands run the same validation and comparison code used by the individual `validate` and `compare` commands, then collect the outcomes into one report. They are intended for regression testing and continuous integration.

## Validate a collection

```sh
omicsbench suite validate
```

By default, every registered object is checked. Link-only records are reported as skipped because they have no redistributed files to validate. Use `--assay` to select one manifest assay or repeat `--id` to select exact objects.

```sh
omicsbench suite validate --assay sequence_dna
omicsbench suite validate --id sequence-001 --id sequence-004
```

Validation checks each manifest, declared file, byte count, SHA-256 digest and scientific expected result. One corrupt object does not stop the suite from checking the remaining selection.

## Compare differential-expression results

```sh
omicsbench suite compare RESULTS_DIRECTORY
```

The directory must contain one result per selected biological object. The base name is the stable benchmark ID. Accepted names are `<id>.csv`, `<id>.tsv`, `<id>.csv.gz` and `<id>.tsv.gz`. Each table must contain unique `gene_id` values and finite `log2_fold_change` values.

When no IDs are supplied, all comparable bulk RNA-seq objects are required. A missing result is a failed case rather than a silent skip. If more than one accepted file exists for the same ID, the case fails as ambiguous.

Each comparison reports shared and missing genes, effect-rank correlation, top-gene overlap, sign agreement and the thresholds declared by that benchmark. An unexpected gene universe or a missed threshold fails the case.

## Evaluate bulk RNA, single-cell, chromatin, DNA and protein results together

```sh
omicsbench suite evaluate RESULTS_DIRECTORY
```

The unified evaluator chooses a comparator from the benchmark contract. Biological bulk RNA-seq objects accept the same CSV or TSV effects as `suite compare`. The single-cell object accepts `<id>.cells.tsv`, `<id>.cells.csv` or a gzip-compressed equivalent. The chromatin object accepts `<id>.peaks.bed`, `<id>.narrowPeak` or a gzip-compressed equivalent. The synthetic whole-genome DNA object accepts `<id>.vcf` or `<id>.vcf.gz` and reports exact-allele precision, recall, F1 and TP, FP and FN counts. The synthetic protein-identification object accepts `<id>.peptides.tsv` or `<id>.peptides.tsv.gz` with `protein_id` and `peptide` columns.

A selected output is required. The suite fails a missing file, more than one accepted file for the same benchmark, malformed input, reference disagreement or a missed scientific threshold. Without `--id`, all objects with a result comparator are required.

The VCF path accepts VCF 4.x, expands comma-separated ALT alleles and assesses records whose FILTER value is `PASS` or `.`. The current truth object contains A/C/G/T substitutions. It compares contig, one-based position, REF and ALT; it does not compare genotypes or reconcile equivalent haplotype representations. The report states this boundary directly.

Peptide tables may include `q_value`. When present, rows above the benchmark's declared maximum are excluded before scoring. The comparator reports peptide precision, recall and F1, protein recall and the fraction of truth-peptide residues recovered. It assesses exact unmodified protein-peptide pairs from the packaged digest. Spectra, modifications, charge states, protein inference and FDR calibration are outside its scope.

Single-cell tables use `cell_id` and `label` by default. Cluster names do not need to match the truth labels. The comparator reports adjusted Rand index, normalized mutual information, purity and cell coverage. It fails on unexpected cells or missed thresholds. These measurements describe recovery of a synthetic partition, not biological cell-type annotation.

Chromatin peak files use zero-based, half-open coordinates and are checked against the fixture chromosome sizes. The comparator reports interval and unioned base-pair precision, recall and F1 under the benchmark's reciprocal-overlap rule. These measurements describe coordinate recovery against synthetic truth, not biological peak validity.

## Index pipeline outputs

Benchmark-named files are convenient for small tests, but production workflows often place outputs in nested directories and use tool-specific headers. A result manifest connects those existing files to benchmark IDs without copying or renaming them.

```sh
omicsbench suite index results \
  --result rnaseq-002=rna/deseq2.csv \
  --result singlecell-001=single-cell/clusters.csv \
  --result interval-001=chromatin/peaks.narrowPeak \
  --result sequence-005=dna/final.vcf.gz \
  --result sequence-006=protein/identifications.tsv \
  --column rnaseq-002:gene_id=Geneid \
  --column rnaseq-002:log2_fold_change=log2FoldChange \
  --column singlecell-001:cell_id=Barcode \
  --column singlecell-001:label=Cluster \
  --column sequence-006:protein_id=Protein \
  --column sequence-006:peptide=Sequence \
  --column sequence-006:q_value=QValue
```

The command writes `results.json` inside the result directory. Each entry contains one benchmark ID, a normalized relative path, the observed SHA-256 digest and any nondefault column mappings. The contract is published as [result-manifest.schema.json](../schemas/result-manifest.schema.json).

When `results.json` is present, `suite evaluate` and `suite evaluate-matrix` use it instead of the filename convention. They verify the complete manifest before scoring. Absolute paths, parent traversal, duplicate benchmark IDs, unknown fields, missing files, changed hashes, unsupported suffixes and invalid mappings are rejected as input errors.

RNA tables can map `gene_id` and `log2_fold_change`. Single-cell tables can map `cell_id` and `label`. Peptide tables can map `protein_id`, `peptide` and `q_value`. VCF and BED-family fields cannot be remapped because their meaning comes from their format contracts. Omitted mappings retain the standard OpenOmicsBench headers.

The effective mapping is part of the run fingerprint because it changes how bytes are interpreted. Directory names and JSON formatting are excluded. A conventional layout and a manifest layout therefore receive the same fingerprint when their file bytes, benchmark selection, method receipt and effective columns are identical.

An optional `method.json` in the result directory uses this contract:

```json
{
  "schema_version": "1.0",
  "name": "Example workflow",
  "version": "4.2.0",
  "command": "example run --mode strict",
  "container": "registry.example.org/workflow@sha256:...",
  "source_revision": "0123456789abcdef",
  "runtime_seconds": 91.4,
  "peak_memory_mb": 820.0,
  "threads": 4,
  "parameters": {
    "mode": "strict"
  }
}
```

Only `schema_version`, `name` and `version` are required. Resource values must be finite and non-negative, and thread count must be positive.

## Compare several methods

```sh
omicsbench suite matrix \
  --method current=results/current \
  --method candidate=results/candidate
```

Supply at least two unique method names. Each directory follows the same filename rules as `suite compare`, and every method is checked against the same benchmark selection. The matrix contains one case per method and benchmark, per-method mean metrics and a deterministic leaderboard. The leaderboard favours passed cases first; its metric values are descriptive and do not replace the declared pass thresholds.

For mixed bulk RNA, single-cell, chromatin, DNA and protein outputs, use:

```sh
omicsbench suite evaluate-matrix \
  --method current=results/current \
  --method candidate=results/candidate \
  --id rnaseq-002 \
  --id singlecell-001 \
  --id interval-001 \
  --id sequence-005 \
  --id sequence-006
```

Every unified evaluation records SHA-256 hashes for its benchmark manifests and submitted files. A `run_fingerprint` binds those hashes to the package version, effective input contract and validated method receipt while excluding local paths. Identical content therefore has the same fingerprint in different directories or on different machines.

## Apply an acceptance policy

```sh
omicsbench suite gate evaluation.json policy.json
```

Gate policies are strict JSON documents. They can require exact benchmark IDs, overall and per-assay pass rates, a method receipt, runtime and peak-memory ceilings, and metric floors limited to an assay or comparator. Required evidence that is absent fails the gate.

```json
{
  "schema_version": "1.0",
  "name": "release gate",
  "required_benchmarks": ["rnaseq-002", "interval-001", "sequence-005", "sequence-006"],
  "minimum_pass_rate": 1.0,
  "minimum_assay_pass_rate": {
    "bulk_rna_seq": 1.0,
    "chromatin_accessibility": 1.0,
    "whole_genome_dna_seq": 1.0,
    "sequence_protein": 1.0
  },
  "require_method_receipt": true,
  "maximum_runtime_seconds": 3600,
  "maximum_peak_memory_mb": 4096,
  "metric_floors": [
    {"metric": "precision", "minimum": 0.99, "comparator": "exact_snv"},
    {"metric": "interval_f1", "minimum": 0.90, "comparator": "genomic_intervals"},
    {"metric": "protein_recall", "minimum": 1.0, "comparator": "peptide_identification"}
  ]
}
```

## Compare a run history

```sh
omicsbench suite trend \
  --run accepted=evaluation-accepted.json \
  --run candidate=evaluation-candidate.json \
  --json trend.json \
  --html trend.html
```

Run entries are ordered from earliest to latest. Each input must be a fingerprinted `suite evaluate` report. The result preserves method and resource receipts, pass rates, source fingerprints, benchmark status history, current metrics and first-to-latest metric deltas.

The multi-assay leaderboard ranks methods by passed cases, failed cases and method name. Assay summaries, runtime and peak memory remain visible, but metrics from different assays are not pooled.

## Check for regressions

```sh
omicsbench suite regress accepted.json candidate.json --absolute-tolerance 0.01
```

Regression checks accept validation, comparison, matrix, evaluation and multi-assay matrix JSON reports. A baseline case that changes from pass or skipped to fail is a regression. Decreases in effect-rank correlation, top-gene overlap, sign agreement, coverage, variant metrics, clustering metrics, protein recovery, interval overlap or base-pair overlap also fail when they exceed the absolute tolerance. Missing baseline cases fail unless `--allow-missing` is supplied. New candidate cases are recorded without failing the check.

## Report formats

The complete report is always printed as JSON. It contains the operation, selection, summary counts and one result per benchmark.

Use any combination of these output options:

```sh
omicsbench suite validate \
  --json reports/validation.json \
  --markdown reports/validation.md \
  --junit reports/validation.xml \
  --csv reports/validation.csv \
  --html reports/validation.html
```

JSON preserves all calculated details. Markdown is a short review table. JUnit XML records one test case per benchmark and can be consumed by CI systems. CSV provides a flat table for further analysis. HTML is a responsive, self-contained, script-free dashboard with summary cards, rankings, assay breakdowns and a method-by-benchmark matrix when applicable. Report files are written atomically and are not replaced unless `--force` is present.

Evaluation commands can also write one portable evidence crate:

```sh
omicsbench suite evaluate results --id rnaseq-002 --evidence evidence.zip
omicsbench evidence verify evidence.zip
```

The deterministic ZIP contains normalized JSON, CSV and HTML reports, the exact submitted result files, an optional result manifest and method receipt, `checksums.sha256` and an `ro-crate-metadata.json` document using the [RO-Crate 1.3 context](https://www.researchobject.org/ro-crate/specification). Submitted files use stable benchmark-based names inside the crate. Verification checks the complete member inventory, every SHA-256 digest, the evaluation report contract and the required RO-Crate root entities without extracting files.

Exit status 0 means the suite completed and every required case passed. Status 2 means the suite completed with at least one failed case. Status 1 means the request itself was invalid, such as an unknown benchmark ID, an absent results directory or a protected output path.
