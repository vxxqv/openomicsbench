# Phytophthora response in gsnor1 Arabidopsis

`rnaseq-011` is a compact test object derived from [E-MTAB-8845](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8845).

## Design

- Organism: Arabidopsis thaliana (NCBI taxonomy 3702)
- Contrast: none versus Phytophthora parasitica
- Samples: 6
- Pocket: 500 of 32,833 source genes
- Reference: Ensembl Genomes 44 on TAIR10

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-011
omicsbench validate rnaseq-011
omicsbench get rnaseq-011 --size pocket
```

rnaseq-011 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 500 pocket identifiers matched Ensembl Genomes 44.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: RNA sequencing of Arabidopsis thaliana gsnor1 and Col-0 after Phytophthora parasitica infection. Expression Atlas E-MTAB-8845.

rnaseq-011 includes material from E-MTAB-8845 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
