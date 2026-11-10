# RNA-seq input packs and result tables

The 12 biological bulk RNA-seq objects have two useful scales. The default `pocket` tier is a small feature subset selected for declared preservation measures. The `full` tier contains the complete count matrix for the selected samples. Both use the same sample IDs and study design. The full tier is larger; it is not a substitute for FASTQ-level alignment or quantification tests.

## Make an input pack

```sh
omicsbench rnaseq pack quick-pack --id rnaseq-002 --id rnaseq-005
omicsbench rnaseq pack full-pack --id rnaseq-002 --tier full
omicsbench rnaseq verify-pack quick-pack
```

Omit `--id` to include all 12 biological objects. The command creates a new directory and refuses to replace an existing one. Each study directory contains:

| File | Contents |
| --- | --- |
| `counts.tsv` | `gene_id` followed by sample-ID columns; non-negative integer counts |
| `samples.tsv` | `sample`, `condition`, `batch`, `replicate`, `strandedness` in count-column order |
| `contrasts.tsv` | One declared `id`, `variable`, `reference`, `target`, `blocking` row |
| `pack.json` | Source, rights, limitations, design, dimensions and file checksums |

The top-level `pack-index.json` identifies every study and checks each `pack.json`. Expected differential-expression answers are not copied. They remain available in the installed collection. The files are suitable as small matrix-level workflow inputs; external pipelines may need a matching annotation file or other method-specific inputs. OpenOmicsBench does not claim that these packs have passed a run of any named external workflow.

The reference contrast is written as **target minus reference**. For `rnaseq-005`, this is tumour tissue minus adjacent non-tumour tissue, with the paired sample identifier in `batch`. For objects with a declared `~ batch + condition` design, `blocking` is `batch`; otherwise it is empty. This derives from each object's archived DESeq2 evidence, and pack creation refuses an unknown formula or inconsistent condition labels.

The matrix, observation sheet and contrast layout follows the fields described in the [nf-core differentialabundance input guide](https://nf-co.re/differentialabundance/latest/docs/usage/). That guide also explains method-specific requirements. In particular, some Salmon-derived raw gene counts need transcript-length offsets or scaled counts; the provenance of each OpenOmicsBench matrix must be checked before choosing an external method. These packs are not asserted to be drop-in replacements for a pipeline's full input contract.

## Check an external design

```sh
omicsbench rnaseq check-design counts.tsv samples.tsv contrasts.tsv
omicsbench rnaseq check-design counts.tsv samples.tsv contrasts.tsv --id rnaseq-005
```

The independent check reads CSV, TSV and their gzip-compressed forms. It requires unique gene and sample IDs, non-negative integer counts, nonempty libraries, matching matrix and sample-sheet ID sets, and a single contrast target and reference with at least two samples each. With blocking variables, it checks the fixed-effect design matrix for rank and residual degrees of freedom. A mismatch returns exit status 2 with a JSON report; a malformed file returns 1. A reordered sample sheet is reported but passes if the IDs still match. `--id` additionally checks the benchmark's sample annotations and reference contrast. Only two-level fixed-effect contrasts are currently checked; the command does not evaluate mixed models, interactions, gene annotation compatibility or biological validity.

## Use native differential-expression results

```sh
omicsbench rnaseq score-effects rnaseq-002 results/deseq2.csv
omicsbench rnaseq score-effects rnaseq-002 results/limma.tsv --reverse-sign
omicsbench rnaseq normalize-effects rnaseq-002 results/edger.tsv results/rnaseq-002.tsv
```

`score-effects` compares a native table directly with its named benchmark, regardless of the source filename. `normalize-effects` creates a two-column `gene_id` and `log2_fold_change` TSV for the other suite commands. Both accept CSV, TSV and `.gz` files. Recognized effect headers include DESeq2's `log2FoldChange` and edgeR or limma's `logFC`; recognized ID headers include `gene_id`, `Geneid`, `GeneID`, `gene`, `Gene`, `feature_id`, `Unnamed: 0` and an empty first column from an R CSV export. These conventions are documented by [DESeq2](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), [edgeR](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf) and [limma](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).

Use `--gene-column HEADER` or `--effect-column HEADER` when more than one recognized column exists. Empty or `NA` effects are omitted and counted; duplicate gene IDs, nonnumeric values and infinite effects are rejected. No sign is inferred from the filename or method. `--reverse-sign` is an explicit choice for a source result whose target and reference are opposite the benchmark's declared contrast.

The score uses the object's full-source DESeq2 reference, reports shared and missing gene IDs, rank correlation, top-gene overlap and sign agreement, then applies its existing thresholds. It tests agreement with that particular reference. Method differences and biological uncertainty are not resolved by a passing score. A published [SEQanswers discussion of RNA-seq benchmarking](https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/9192-testing-benchmarking-rnaseq-de-tools) gives context for this limitation.
