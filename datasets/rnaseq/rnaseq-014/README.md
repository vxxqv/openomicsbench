# Paired mesenchymal stem-cell strain response

`rnaseq-014` is a compact test object derived from [E-MTAB-7925](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-7925).

## Design

- Organism: Homo sapiens (NCBI taxonomy 9606)
- Contrast: none versus cyclic tensile strain
- Samples: 12
- Pocket: 8,000 of 61,860 source genes
- Reference: Ensembl 107 on GRCh38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `expected/reference-effects.tsv.gz` contains finite full-source DESeq2 effects for result comparison.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-014
omicsbench validate rnaseq-014
omicsbench get rnaseq-014 --size pocket
```

rnaseq-014 passed all four predeclared preservation minimums with DESeq2 1.54.0. All 8,000 pocket identifiers matched Ensembl 107.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: RNA-seq of human mesenchymal stem cells from six donors with and without cyclic tensile strain. Expression Atlas E-MTAB-7925.

rnaseq-014 includes material from E-MTAB-7925 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
