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

## Report formats

The complete report is always printed as JSON. It contains the operation, selection, summary counts and one result per benchmark.

Use any combination of these output options:

```sh
omicsbench suite validate \
  --json reports/validation.json \
  --markdown reports/validation.md \
  --junit reports/validation.xml
```

JSON preserves all calculated details. Markdown is a short review table. JUnit XML records one test case per benchmark and can be consumed by CI systems. Report files are written atomically and are not replaced unless `--force` is present.

Exit status 0 means the suite completed and every required case passed. Status 2 means the suite completed with at least one failed case. Status 1 means the request itself was invalid, such as an unknown benchmark ID, an absent results directory or a protected output path.
