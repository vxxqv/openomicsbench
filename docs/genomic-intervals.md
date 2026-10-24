# Genomic interval benchmark

`interval-001` is a deterministic software fixture for peak-calling outputs. It is intended for parser tests, continuous integration, method wiring and teaching. It does not represent a biological sample.

The object contains two synthetic contigs, a 50-base binned bedGraph signal track and 12 narrowPeak truth intervals. The peak set includes short and broad regions, a close doublet and intervals near contig boundaries. Rebuilding the object with `workflows/build_interval_fixture.py` reproduces the manifest and every declared file exactly.

## Accepted files

The interval reader accepts BED3, BED6 and narrowPeak, with optional gzip compression. Every data row in one file must use the same number of tab-separated fields. Coordinates use the [UCSC BED convention](https://genome.ucsc.edu/FAQ/FAQformat.html): the start is zero-based and included, while the end is excluded. Peak intervals must have positive width. When a chromosome-size file is supplied, every contig must be declared and every end coordinate must be within its contig.

narrowPeak is handled as BED6 plus signal value, negative log10 p-value, negative log10 q-value and a summit offset. A summit of `-1` means unavailable. Track and browser headers are rejected because they are display directives rather than peak records.

## QC

```sh
omicsbench interval qc calls.narrowPeak --chrom-sizes chrom.sizes
```

The report gives the detected format, interval and contig counts, summed interval span, unioned coverage, number of overlapping input intervals, minimum, median and maximum width, and whether the input was already sorted. Union coverage merges overlapping and directly adjacent intervals before counting bases.

## Comparison

```sh
omicsbench interval compare interval-001 calls.narrowPeak
```

A truth interval and a submitted interval match when their overlap covers at least 50% of each interval. Interval precision is the fraction of submitted intervals with a match. Interval recall is the fraction of truth intervals with a match. Their harmonic mean is interval F1.

Base-pair precision and recall are calculated from the intersection of the unioned submitted and truth intervals. This prevents overlapping calls from counting the same bases more than once. The separate base-pair score exposes calls that satisfy the interval rule but are consistently too wide or too narrow. When both matched records contain summit positions, the report also gives the median absolute summit distance.

The included acceptance profile requires at least 0.90 for interval precision, recall and F1, and at least 0.85 for base-pair precision, recall and F1. These gates test coordinate recovery against a known synthetic truth set. They do not establish suitable thresholds for a biological experiment.

The broader assay workflow still needs read quality checks, alignment review, duplicate and mitochondrial handling, blacklist filtering, TSS enrichment where applicable, replicate concordance and an assay-appropriate peak caller. The [ENCODE ATAC-seq standards](https://www.encodeproject.org/data-standards/atac-seq/atac-encode4/) give biological quality considerations outside this fixture's scope.
