# Transcriptomics source review

Twelve GEO series records were retrieved on 6 September 2026. A thirteenth candidate, E-MTAB-8572, was reviewed through Expression Atlas and BioStudies on 7 September. Intake files retain the accession, source checksums, sample count, assay type, source publication identifiers and file links. Long study abstracts are not republished. The numerical scorecard remains incomplete wherever sample-level review is still needed.

## First candidates

E-MTAB-8572 is the first candidate with a GREEN redistribution decision. Expression Atlas supplies a 58,735-gene integer count matrix and a ten-run design: five wild-type A549 xenografts and five SLC2A5-knockout A549 xenografts. The checked files are 2.69 MB in total, their sample identifiers agree exactly, and the Atlas methods name Ensembl release 95, HISAT2 and featureCounts. A deterministic reduction run found 2,000 genes to be the smallest tested candidate meeting all four development thresholds. Its report reproduced byte for byte on a second run. The exact assembly and annotation file identities still need to be pinned, so this candidate is not yet release eligible. See [Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8572), [BioStudies](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8572) and the [Expression Atlas licence](https://www.ebi.ac.uk/gxa/licence.html).

Pasilla is the first implemented retrieval recipe because it offers a small processed integer count table and a non-human example. The pinned package table has 14,599 genes and seven libraries. The original package sample annotation includes library preparation type, which is needed for the current diagnostic. The package description cites six GEO accessions while the table has seven libraries. Resolve that mapping before certifying provenance at library level. See [Pasilla](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html) and [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508).

The airway study offers a useful paired cell-line design. Its parent GEO series has 16 samples; the commonly used dexamethasone subset has eight. GEO's supplementary expression matrix is FPKM, so it must not be relabelled as raw counts. Obtain a versioned count representation with its own processing and reference metadata. See [airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html) and [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).

[GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450) offers a twelve-library mammary-gland design with cell population and developmental stage. It is suitable for a carefully selected replicated contrast; the entire factorial design needs a more explicit adapter than the current two-condition baseline.

[GSE49712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49712) provides ten libraries from reference RNA materials with ERCC controls. This could support technical evaluation. Replicate libraries of reference material should not be presented as independent human donors.

## Candidates needing more work

GSE47774 is a large multicentre SEQC record. Select site, platform and material explicitly, and check overlap with other SEQC accessions before counting independent objects. GSE53960 offers rat tissue and age variation. GSE55347 offers toxicogenomics. Both need sample-level contrast and replicate review before processing.

GSE60314 includes genotype, sex, environment and repeated-library structure in Drosophila. Its record mentions different FlyBase reference releases. The reference and biological unit need to be resolved before a pocket is meaningful.

GSE37704 is a mixed-assay superseries. Resolve the RNA-seq subseries before selecting a transcript-level object. GSE50760 needs subject and tissue pairing review. GSE48138 is a T-cell lineage study with RNA-fraction and custom lincRNA annotation issues; it is not the high-replicate yeast benchmark sometimes sought for this collection.

GSE64016 is explicitly single-cell and is excluded from bulk RNA-seq v1. Its intake record is retained so the same unsuitable source is not reconsidered without that context.

## Rights decision

All screened GEO sources remain AMBER. The [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) does not provide blanket unrestricted copying or distribution permission. The Pasilla package declares LGPL, which is useful evidence, but the exact data obligations and historical modENCODE terms need review before extracted count tables are repackaged. Expression Atlas applies CC BY 4.0 to copyrightable material available on its website; E-MTAB-8572 is therefore GREEN with attribution obligations recorded in its intake file. No source counts or biological reductions are included in the development archive.

## Infrastructure evidence

[Pydantic models](https://docs.pydantic.dev/latest/concepts/models/) provide the typed contract and generated JSON Schema. JSON is the canonical editable manifest in this snapshot; general YAML parsing is not implemented. This avoids silently supporting only part of YAML while claiming full support.

[nf-core test datasets](https://github.com/nf-core/test-datasets) and [nf-core RNA-seq](https://nf-co.re/rnaseq) informed the separation between miniature software fixtures and scientific workflow data. This project has not executed the nf-core RNA-seq pipeline.

[Zenodo's software metadata guidance](https://help.zenodo.org/docs/github/describe-software/) confirms that `.zenodo.json` takes precedence over `CITATION.cff` when both are present. Publication metadata therefore needs one reviewed source of truth. Authors, license and DOI remain explicit release inputs.
