# Pasilla count-matrix candidate

This record points to the seven-library Pasilla count table maintained in the Bioconductor source package. It contains four untreated and three knockdown libraries with both single-read and paired-end preparations. The broader GEO series includes other experiments and must not be treated as this seven-library object.

The record is link-only. No source counts or reduced biological data are included in this repository. The pinned retrieval recipe is workflows/pasilla-source.json at the repository root. Run `python workflows/pasilla_diagnostic.py --staging staging/pasilla` to retrieve the source into a local staging directory and evaluate candidate feature subsets.

The diagnostic fits condition and library type on log2(CPM + 1). It checks effect rankings, top-50 overlap, sample distances and effect signs. A separate DESeq2 1.50.2 run on the full table is recorded in evidence/pasilla-deseq2-baseline.json. Reference compatibility, data-specific redistribution terms and held-out threshold review remain open.

Source: Brooks et al., Genome Research 2011, PMID 20921232, GEO GSE18508. Credit the Pasilla package authors Wolfgang Huber and Alejandro Reyes for the processed count representation. There is no OpenOmicsBench version DOI yet.
