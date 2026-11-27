# Mild drought response in Arabidopsis ICE163

`rnaseq-016` is a compact test object derived from [E-MTAB-10374](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10374).

## Design

- Organism: Arabidopsis thaliana (NCBI taxonomy 3702)
- Contrast: none versus drought environment
- Samples: 6
- Pocket: 2,000 of 32,833 source genes
- Reference: Ensembl Genomes 53 on TAIR10

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-016
omicsbench validate rnaseq-016
omicsbench get rnaseq-016 --size pocket
```

rnaseq-016 passed all four predeclared preservation minimums with DESeq2 1.54.0. All 2,000 pocket identifiers matched Ensembl Genomes 53.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Chen et al. Distinct cellular strategies determine sensitivity to mild drought of Arabidopsis natural accessions. Expression Atlas E-MTAB-10374.

rnaseq-016 includes material from E-MTAB-10374 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
