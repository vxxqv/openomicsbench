# Myofibroblast differentiation versus control

`rnaseq-017` is a compact test object derived from [E-MTAB-3101](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-3101).

## Design

- Organism: Homo sapiens (NCBI taxonomy 9606)
- Contrast: control versus TGF-beta-1 induced myofibroblast differentiation
- Samples: 6
- Pocket: 1,000 of 58,735 source genes
- Reference: Ensembl 95 on GRCh38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-017
omicsbench validate rnaseq-017
omicsbench get rnaseq-017 --size pocket
```

rnaseq-017 passed all four predeclared preservation minimums with DESeq2 1.54.0. All 1,000 pocket identifiers matched Ensembl 95.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Mellone et al. Induction of fibroblast senescence generates a non-fibrogenic myofibroblast phenotype that differentially impacts on cancer prognosis. Expression Atlas E-MTAB-3101.

rnaseq-017 includes material from E-MTAB-3101 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
