# ELP3 depletion in BT549 breast cancer cells

`rnaseq-007` is a compact test object derived from [E-MTAB-9206](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9206).

## Design

- Organism: Homo sapiens (NCBI taxonomy 9606)
- Contrast: scrambled shRNA versus ELP3 shRNA
- Samples: 6
- Pocket: 8,000 of 58,735 source genes
- Reference: Ensembl 95 on GRCh38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-007
omicsbench validate rnaseq-007
omicsbench get rnaseq-007 --size pocket
```

rnaseq-007 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 8,000 pocket identifiers matched Ensembl 95.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: BT549 cells depleted of ELP3 compared with a scrambled shRNA control. Expression Atlas E-MTAB-9206.

rnaseq-007 includes material from E-MTAB-9206 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
