# Paired prostate tumour and adjacent tissue

`rnaseq-005` is a compact test object derived from [E-MTAB-567](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-567).

## Design

- Organism: Homo sapiens (NCBI taxonomy 9606)
- Contrast: tumor tissue versus adjacent non-tumor tissue
- Samples: 28
- Pocket: 8,000 of 61,860 source genes
- Reference: Ensembl 107 on GRCh38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-005
omicsbench validate rnaseq-005
omicsbench get rnaseq-005 --size pocket
```

rnaseq-005 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 8,000 pocket identifiers matched Ensembl 107.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Ren et al. RNA-seq analysis of prostate cancer in the Chinese population. Cell Research. DOI 10.1038/cr.2012.30.

rnaseq-005 includes material from E-MTAB-567 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
