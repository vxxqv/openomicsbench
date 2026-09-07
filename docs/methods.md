# Reduction and validation methods

## Scope of the current baseline

The working Python baseline is a diagnostic for preservation of count-matrix structure. It converts each sample to counts per million, applies log2(CPM + 1), and fits an ordinary least-squares model with an intercept, a two-level condition and categorical batch indicators. It reports the condition coefficient. It does not estimate negative-binomial dispersion, report p-values or replace DESeq2.

Condition levels are sorted lexically. The reported contrast is the second level minus the first. For the Pasilla labels, that means untreated minus treated. Every output records this direction. Library preparation type is the blocking covariate for this diagnostic; it is not asserted to be the experimental batch.

The model rejects rank-deficient designs and designs with no residual degrees of freedom. At least two samples per condition are required. The supplied plan's example places controls only in batch A and treated samples only in batch B; it therefore cannot identify a treatment effect adjusted for batch. The implementation tests this failure explicitly.

## Feature selection

All sample columns are retained. The selector alternates among four ordered feature lists: high variance of log2(count + 1), expression near the median feature mean, low variance, and a seeded SHA-256 ordering of feature IDs. Already selected genes are skipped. Source row order is restored in the final matrix. This uses no condition labels to rank features.

The procedure provides nested candidate sizes, deterministic tie breaking and an exact seed. It does not claim that low-variance genes are biological housekeeping genes. Zeros and low-expression features are allowed. A future biological profile may require a different documented mixture; changing the algorithm creates a new baseline version.

The Pasilla grid is 500, 1,000, 2,000, 4,000, 8,000 and 12,000 genes. The minimum requirements were saved before this diagnostic run: rank correlation 0.90, top-50 Jaccard 0.60, distance correlation 0.90 and sign agreement 0.90. They are development thresholds. They have not been independently calibrated for a biological release.

E-MTAB-8572 uses a contrast-aware extension. The 100 largest absolute full-source effects are included first, then the four general feature lists fill a single deterministic order. Every candidate size is a prefix of that order before source row order is restored. This prevents a small pocket from losing the strongest source effects by chance while retaining genes from the four general selection categories.

Its diagnostic uses median-ratio size factors calculated from genes with positive counts in every sample. Size factors are centered to a geometric mean of one, and the model fits log2(normalized count + 1). This is close to the normalization stage used by DESeq2, but the diagnostic remains ordinary least squares and produces no inferential statistics. The first passing size in the declared grid was 2,000 genes. The saved report was identical across two complete source retrieval and calculation runs.

## Quantitative definitions

Spearman preservation uses average ranks for tied coefficients on the shared selected feature universe. Constant or non-finite input is an error, not a passing score.

Top-50 Jaccard ranks the absolute condition coefficients within that same shared universe. Feature ID resolves ties. The denominator is the union size. This metric does not measure recovery of top genes omitted from the source universe. A release baseline must add full-universe top-feature recall before claiming preservation of the source's strongest findings.

Distance preservation is Pearson correlation between the strict upper triangles of Euclidean sample-distance matrices after the declared normalization and log2 transform. The source distance uses all source genes; the pocket distance uses its selected genes. This check can fail even when effect rankings remain correlated.

Sign agreement considers shared genes with an absolute source coefficient of at least 0.1. A zero pocket coefficient disagrees with a nonzero source coefficient. The implementation rejects an empty eligible set. It does not assign arbitrary success to an undefined metric.

PCA uses singular-value decomposition after centering each transformed feature across samples. For repeatable plotting, each component is oriented so its largest absolute sample coordinate is positive. Distances, rather than exact PCA axis signs, are the validation target. Nearly equal singular values may rotate components between environments.

## Hard invariants

The validator checks declared paths, byte sizes and SHA-256 hashes before scientific calculations. The baseline sample order must agree with both the manifest and profile. The pocket shape must match the declared shape. Every retained gene must occur in the source and every retained count must equal the corresponding source value.

Count input rejects empty or duplicate IDs, inconsistent row widths, fractional counts, negative values and empty sample libraries. It uses signed 64-bit integer storage. This is a count-matrix contract, not a quantification-estimate contract; fractional Salmon counts need a separate adapter.

## Read sampling

The paired-read helper validates FASTQ structure, equal mate counts and matching IDs. It hashes the seed, sample ID, record ordinal and pair ID to select both mates together at the requested fraction. This samples an expected fraction, not an exact number of pairs. It writes plain FASTQ and returns actual input and selected pair counts. It is tested on miniature records; a real read-level pocket and its QC profile remain to be built.

## Reproducibility record

The exact fixture rebuild was checked with Python 3.12.14 and NumPy 2.3.5. Snakemake 9.26.1 ran the fixture workflow from generation through validation. The full Pasilla table was fitted with DESeq2 1.50.2 under R 4.5.3. The Conda and R package versions are recorded in runtime/environment-lock.json. A clean rebuild on another machine and a container digest remain to be completed.
