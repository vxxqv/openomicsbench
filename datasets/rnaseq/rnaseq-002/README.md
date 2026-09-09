# SLC2A5 knockout in A549 xenografts

`rnaseq-002` is a compact test object derived from [E-MTAB-8572](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8572).

## Design

- Organism: Homo sapiens (NCBI taxonomy 9606)
- Contrast: wild type genotype versus SLC2A5 knockout mediated by CRISPR/Cas9
- Samples: 10
- Pocket: 2,000 of 58,735 source genes
- Reference: Ensembl 95 on GRCh38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-002
omicsbench validate rnaseq-002
omicsbench get rnaseq-002 --size pocket
```

rnaseq-002 passed all four predeclared preservation minimums with DESeq2 1.50.2. All 2,000 pocket identifiers matched Ensembl 95.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: RNA-seq of human lung cancer cell line A549-NC and A549-SLC2A5 knockout xenograft in BALB/c nude mice. Expression Atlas E-MTAB-8572.

rnaseq-002 includes material from E-MTAB-8572 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
