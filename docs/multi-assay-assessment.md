# Multi-assay assessment (version 15)

These commands read processed results and produce JSON reports. They do not download data or infer biological truth. Supply a versioned reference with documented origin, assembly or accession, and rights before interpreting a report as a benchmark.

## Single-cell labels and embeddings

```sh
omicsbench cell assess singlecell-002 clusters.tsv --annotation annotations.tsv --embedding embedding.tsv --neighbors 10
```

`clusters.tsv` and `annotations.tsv` need `cell_id` and `label` columns. Cluster names are arbitrary: adjusted Rand index and normalized mutual information measure partition agreement. Annotation names are compared exactly to the declared reference, and missing or extra cells are reported. `embedding.tsv` needs `cell_id` and at least two finite numeric coordinates. Its same-state and same-batch neighbor fractions are shown beside random-label baselines. Batch mixing and state preservation must be considered together; neither figure is a standalone quality score. Cases `singlecell-002` through `singlecell-007` contain sourced CellBench counts and labels across three platforms and three technical plates, all from one study; they are not six independent cohorts. Their source, transformations, reference-identifier exclusions and reuse notice are recorded in each case and in the [intake record](../curation/v15-source-intake.md). The separate `singlecell-001` is a deterministic synthetic demonstration.

An optional AnnData adapter exports explicit TSVs:

```sh
python -m pip install anndata
omicsbench cell export-h5ad analysis.h5ad exported --label-column cell_type --batch-column batch --embedding-key X_umap
```

Review the exported labels and embedding before scoring. The adapter selects named `obs` and `obsm` fields; it does not infer cell types. [AnnData documents `read_h5ad`](https://anndata.readthedocs.io/en/stable/generated/anndata.io.read_h5ad.html).

## Small variants

```sh
omicsbench variant benchmark truth.vcf.gz calls.vcf.gz benchmark-regions.bed.gz --assembly GRCh38 --truth-assembly GRCh38
```

Only PASS or unfiltered literal A/C/G/T alleles are compared. Each allele's full REF span must lie within a supplied BED region. The report separates SNVs, indels and multi-nucleotide variants, and counts excluded outside-region alleles. Assembly names are *user declarations*; the command cannot prove that the files use those assemblies. It compares exact `CHROM/POS/REF/ALT` representations, not haplotypes, genotypes or left-normalized equivalents.

For a scientific germline benchmark, use a matched GIAB VCF/BED/reference set and a representation-aware tool such as [hap.py](https://github.com/Illumina/hap.py):

```sh
hap.py truth.vcf.gz calls.vcf.gz -f benchmark-regions.bed.gz -r reference.fa -o happy-result
```

NIST [publishes matched benchmark calls and regions](https://www.nist.gov/programs-projects/genome-bottle) and [explains why complex calls can sit outside the BED regions](https://www.nist.gov/programs-projects/faqs-genome-bottle). This command's exact-match numbers should not be substituted for hap.py's haplotype-aware numbers.

## Chromatin and peaks

```sh
omicsbench interval assess peaks.narrowPeak --replicate replicate-2.narrowPeak --reference reference-peaks.narrowPeak --qc-json qc.json --idr-json idr.json
```

`qc.json` may contain precomputed `frip` (0–1) and `tss_enrichment` (nonnegative), and must name their `source`. `idr.json` must contain an external `source` and nonnegative `idr_peaks`. OpenOmicsBench imports these values and does **not** calculate FRiP, TSS enrichment or IDR from reads. The report separates peak width, reciprocal replicate overlap, supplied-reference overlap, and width bins. Summit distance appears when both files have valid narrowPeak summits. Use matching assemblies, blacklists, versions and replicate definitions. [ENCODE's ATAC-seq standards](https://www.encodeproject.org/data-standards/atac-seq/atac-encode4/) explain those measurements and IDR's role.

## Processed peptides

```sh
omicsbench peptide assess reference.tsv calls.tsv --mapping accessions.tsv --q 0.001 --q 0.01 --q 0.05
```

Each table requires `protein_id` and `peptide`, with optional `q_value` and positive `quantity`. The optional mapping TSV has `source_id` and `reference_id` columns; one-to-many mappings are reported as ambiguous and excluded. Reports distinguish exact protein-peptide pair recovery, distinct peptide-sequence recovery, protein recovery, q-threshold behavior, and (when both tables contain quantities) the median absolute log2 ratio error on matched pairs. This does not reprocess spectra, calibrate false discovery rate, identify modifications, or solve protein inference. [PRIDE describes processed identification data](https://www.ebi.ac.uk/training/online/courses/pride-quick-tour/what-consitutes-a-pride-dataset/); a PRIDE accession alone is not a scored truth set.

## Case intake before release claims

For each new biological case, record source accession and version, permitted redistribution, file hashes, sample/replicate design, reference and truth type, transformations, and known confounders. Small synthetic fixtures are useful for checking commands but are not counted as independently sourced biological benchmarks. Version 15 ships these processed-result commands with the stated boundaries; it does not certify their use as clinical or biological ground truth.
