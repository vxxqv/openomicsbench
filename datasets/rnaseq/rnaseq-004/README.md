# Klf1 knockout in mouse erythroid tissue

`rnaseq-004` is a compact test object derived from [E-GEOD-33979](https://www.ebi.ac.uk/gxa/experiments/E-GEOD-33979).

## Design

- Organism: Mus musculus (NCBI taxonomy 10090)
- Contrast: wild type genotype versus Klf1-/-
- Samples: 6
- Pocket: 4,000 of 55,573 source genes
- Reference: Ensembl 97 on GRCm38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-004
omicsbench validate rnaseq-004
omicsbench get rnaseq-004 --size pocket
```

rnaseq-004 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 4,000 pocket identifiers matched Ensembl 97.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Tallack et al. Novel roles for Klf1 in regulating the erythroid transcriptome revealed by mRNA-seq. Genome Research. DOI 10.1101/gr.135707.111.

rnaseq-004 includes material from E-GEOD-33979 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
