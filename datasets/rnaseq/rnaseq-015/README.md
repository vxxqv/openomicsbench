# FOXK1/2 knockdown in neuronal stem cells

`rnaseq-015` is a compact test object derived from [E-MTAB-9639](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9639).

## Design

- Organism: Homo sapiens (NCBI taxonomy 9606)
- Contrast: siRNA knockdown of FOXK1/2 versus scrambled siRNA
- Samples: 6
- Pocket: 1,000 of 58,735 source genes
- Reference: Ensembl 95 on GRCh38

## Files

- `pocket/counts.tsv` contains the selected integer counts.
- `expected/source-counts.tsv` contains the full source matrix for the declared samples.
- `expected/reference-effects.tsv.gz` contains finite full-source DESeq2 effects for result comparison.
- `samples.json` fixes sample order, condition, replicate and blocking values.
- `expected/validation.json` and `expected/deseq2.json` define the checks and results.
- `reference.json`, `rights.json`, `attribution.json` and `provenance/transform.json` record origin and use.

## Check it

```sh
omicsbench info rnaseq-015
omicsbench validate rnaseq-015
omicsbench get rnaseq-015 --size pocket
```

rnaseq-015 passed all four predeclared preservation minimums with DESeq2 1.54.0. All 1,000 pocket identifiers matched Ensembl 95.

## Limits

- Validation begins with the archive count matrix and does not repeat alignment or feature counting.
- Use the cited source study for the complete experiment or a different contrast.
- Strandedness is unknown because count-matrix validation does not require it.

Source citation: Ji et al. The forkhead transcription factor FOXK2 premarks lineage-specific genes in human embryonic stem cells for activation during differentiation. Nucleic Acids Research. DOI 10.1093/nar/gkaa1281.

rnaseq-015 includes material from E-MTAB-9639 under CC BY 4.0 with provider credit. See `rights.json` and `attribution.json` for the dated record.
