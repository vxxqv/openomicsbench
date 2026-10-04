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

## Evaluate RNA and DNA results together

```sh
omicsbench suite evaluate RESULTS_DIRECTORY
```

The unified evaluator chooses a comparator from the benchmark contract. Biological bulk RNA-seq objects accept the same CSV or TSV effects as `suite compare`. The synthetic whole-genome DNA object accepts `<id>.vcf` or `<id>.vcf.gz` and reports exact-allele precision, recall, F1 and TP, FP and FN counts.

A selected output is required. The suite fails a missing file, more than one accepted file for the same benchmark, malformed input, reference disagreement or a missed scientific threshold. Without `--id`, all objects with a result comparator are required.

The VCF path accepts VCF 4.x, expands comma-separated ALT alleles and assesses records whose FILTER value is `PASS` or `.`. The current truth object contains A/C/G/T substitutions. It compares contig, one-based position, REF and ALT; it does not compare genotypes or reconcile equivalent haplotype representations. The report states this boundary directly.

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

For mixed RNA and DNA outputs, use:

```sh
omicsbench suite evaluate-matrix \
  --method current=results/current \
  --method candidate=results/candidate \
  --id rnaseq-002 \
  --id sequence-005
```

The multi-assay leaderboard ranks methods by passed cases, failed cases and method name. Assay summaries, runtime and peak memory remain visible, but metrics from different assays are not pooled.

## Check for regressions

```sh
omicsbench suite regress accepted.json candidate.json --absolute-tolerance 0.01
```

Regression checks accept validation, comparison, matrix, evaluation and multi-assay matrix JSON reports. A baseline case that changes from pass or skipped to fail is a regression. Decreases in effect-rank correlation, top-gene overlap, sign agreement, gene coverage, variant precision, recall or F1 also fail when they exceed the absolute tolerance. Missing baseline cases fail unless `--allow-missing` is supplied. New candidate cases are recorded without failing the check.

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

The deterministic ZIP contains normalized JSON, CSV and HTML reports, the exact submitted result files, optional method receipts, `checksums.sha256` and an `ro-crate-metadata.json` document using the [RO-Crate 1.3 context](https://www.researchobject.org/ro-crate/specification). Verification checks the complete member inventory, every SHA-256 digest, the evaluation report contract and the required RO-Crate root entities without extracting files.

Exit status 0 means the suite completed and every required case passed. Status 2 means the suite completed with at least one failed case. Status 1 means the request itself was invalid, such as an unknown benchmark ID, an absent results directory or a protected output path.
