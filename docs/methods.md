# Methods

## Collection scope

Version 1 contains 12 bulk RNA-seq benchmark objects derived from seven Expression Atlas experiments. The source organisms are human, mouse and Arabidopsis thaliana. Each object represents one declared two-level contrast, with an optional blocking factor or predeclared sample subset. Several contrasts may come from one study when they exercise different valid designs. They keep separate identifiers, sample metadata, validation profiles and provenance.

The collection starts from archive-supplied gene count matrices. It does not repeat read alignment or quantify transcripts from FASTQ. The objects are intended for software tests, tutorials and method checks at the count-matrix stage.

## Source intake

The Expression Atlas adapter retrieves the raw-count and experiment-design resources named by the official experiment page. Downloads use temporary files and are moved into place only after completion. The workflow records the resource URL, byte count and SHA-256 hash.

Count columns must agree with analysed run identifiers in the experiment design. A column marked unanalysed is excluded. Repeated run columns are accepted only when their complete integer count vectors are identical. The parser rejects missing or duplicate gene identifiers, non-integer values, negative counts, inconsistent row widths and empty sample libraries.

Each source configuration names the experimental factor, the two retained factor values, an optional blocking column and any sample subset. These choices are made before candidate pockets are scored. The workflow rejects a two-condition design with fewer than two samples per condition, rank deficiency or no residual degrees of freedom.

## Feature selection

All selected sample columns are retained. No library is downsampled and no count is normalized in the distributed matrices.

The feature order begins with the 100 genes having the largest absolute source contrast effects. Gene identifier settles ties. The remaining positions alternate among four deterministic lists:

1. high variance of `log2(count + 1)`;
2. mean expression closest to the median feature mean;
3. low variance; and
4. a SHA-256 ordering of the seed and gene identifier.

A gene already selected from one list is skipped in the others. Each candidate size is a prefix of this order, then rows are restored to source order. The configuration records the seed, candidate grid, contrast, subset and top-k value. The saved diagnostic report includes passing and failing sizes, so the chosen pocket is the smallest tested candidate that met every threshold rather than an undocumented manual choice.

## Diagnostic envelope

The first validation layer uses median-ratio size factors calculated from genes with positive counts in every sample. Size factors are centred to a geometric mean of one. The diagnostic fits ordinary least squares to `log2(normalized count + 1)`, using condition and the declared blocking factor when present. It is used to screen the candidate grid efficiently. It does not report p-values or claim to reproduce negative-binomial inference.

The four predeclared minimums are:

| Metric | Minimum |
|---|---:|
| Spearman effect correlation | 0.90 |
| Top-50 effect Jaccard overlap | 0.60 |
| Sample-distance correlation | 0.90 |
| Effect-sign agreement | 0.90 |

Spearman preservation is the average-tie rank correlation between full-source and pocket condition effects for genes retained in the pocket. Non-finite or constant inputs are errors.

Top-50 overlap ranks absolute condition effects within the retained feature universe. Gene identifier settles ties. The score is the intersection divided by the union of the two top-50 sets.

Sample-distance preservation is the Pearson correlation between the upper triangles of the full and pocket Euclidean distance matrices after the declared normalization and log transform.

Sign agreement is the fraction of retained genes with matching effect signs when the absolute full-source effect is at least 0.1. An empty eligible set is an error.

## DESeq2 confirmation

Every selected pocket is then compared with its full source matrix using DESeq2 1.50.2 under R 4.5.3. The model is `~ condition` for unblocked designs and `~ batch + condition` where a blocking factor is declared. The baseline script checks sample order, integer counts, model rank and residual degrees of freedom before fitting.

The DESeq2 comparison calculates the same four preservation measures from fitted log2 fold changes and DESeq2 normalized counts. All 12 objects pass all four minimums. Each `expected/deseq2.json` records the formula, contrast direction, feature counts, exact metric values, runtime versions, input hashes, output hashes and hashes of the Python and R workflows. Gene-level DESeq2 tables remain build evidence rather than release files; the compact report contains enough information to verify the certified input and result identity.

## Reference verification

Expression Atlas processing metadata names the genome and gene-model release used for each experiment. The v1 objects resolve to five reference profiles: Ensembl 95 GRCh38, Ensembl 97 GRCm38, Ensembl 107 GRCh38, Ensembl 107 GRCm39 and Ensembl Genomes 44 TAIR10.

For each profile, the workflow pins the official genome and GTF URLs and their checksum records. The GTF is downloaded, checked as a gzip stream, hashed with SHA-256 and parsed for gene identifiers. Every pocket identifier must occur in the matching GTF. The release evidence records matched and unmatched counts for each object. All 12 included objects have complete matches. E-MTAB-7126 is excluded because its current count matrix does not agree with the Ensembl 95 annotation stated by the experiment page.

## Object validation

Before scientific metrics are calculated, the validator compares the actual files with the manifest inventory. It checks every byte count and SHA-256 hash, then checks sample order, count-matrix shape and exact equality of each retained integer count with the source matrix. Normalized values never replace the source counts.

The scientific calculation runs only after these hard checks pass. The result is accepted only if every declared metric meets its minimum. The release preflight separately checks the DESeq2 report, reference record, rights decision, attribution, workflow commits and collection overlap audit.

## Reproducibility and limits

The tested scientific environment uses Python 3.12.14, NumPy 2.3.5, R 4.5.3 and DESeq2 1.50.2. `runtime/environment-lock.json` records the complete Conda and R package set. Source configurations and transformation scripts are versioned, and each biological object names the commits and hashes used to select and assemble it.

The release does not claim that a chosen pocket is the globally smallest possible representation. It is the smallest passing member of its declared grid. Preservation of a count-level contrast does not show that every downstream method, pathway result or biological conclusion is unchanged. Users who need the complete study, alternative contrasts or read-level quality control should return to the cited source archive.
