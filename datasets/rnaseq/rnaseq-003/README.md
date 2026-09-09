# AtRsgA knockout in Arabidopsis seedlings

`rnaseq-003` is a compact test object derived from [E-MTAB-6866](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-6866).

## Design

- Organism: Arabidopsis thaliana (NCBI taxonomy 3702)
- Contrast: wild type genotype versus AtRsgA-i -/-
- Samples: 6
- Pocket: 4,000 of 32,833 source genes
- Reference: Ensembl Genomes 44 on TAIR10

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-003
omicsbench validate rnaseq-003
omicsbench get rnaseq-003 --size pocket
```

rnaseq-003 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 4,000 pocket identifiers matched Ensembl Genomes 44.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Janowski et al. AtRsgA from Arabidopsis thaliana is important for maturation of the small subunit of the chloroplast ribosome. Plant Journal. DOI 10.1111/tpj.14040.

rnaseq-003 includes material from E-MTAB-6866 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
