# Duchenne muscular dystrophy myoblast comparison

`rnaseq-008` is a compact test object derived from [E-MTAB-10322](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10322).

## Design

- Organism: Mus musculus (NCBI taxonomy 10090)
- Contrast: normal versus Duchenne muscular dystrophy
- Samples: 9
- Pocket: 8,000 of 56,748 source genes
- Reference: Ensembl 107 on GRCm39

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-008
omicsbench validate rnaseq-008
omicsbench get rnaseq-008 --size pocket
```

rnaseq-008 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 8,000 pocket identifiers matched Ensembl 107.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Gosselin et al. Loss of full-length dystrophin expression causes major cell-autonomous abnormalities in proliferating myoblasts. eLife. DOI 10.7554/eLife.75521.

rnaseq-008 includes material from E-MTAB-10322 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
