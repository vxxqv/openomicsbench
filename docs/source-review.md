# Source review

## Selection result

The v1 biological collection uses seven Expression Atlas studies and defines 12 benchmark objects. The difference between studies and objects is deliberate: two factorial studies support several predeclared contrasts or strata, each with its own identifier, sample design, pocket, validation and provenance.

| Source | Included object | Main design |
|---|---|---|
| [E-MTAB-8572](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8572) | `rnaseq-002` | A549 SLC2A5 knockout xenografts |
| [E-MTAB-6866](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-6866) | `rnaseq-003` | Arabidopsis AtRsgA knockout |
| [E-GEOD-33979](https://www.ebi.ac.uk/gxa/experiments/E-GEOD-33979) | `rnaseq-004` | Mouse Klf1 knockout |
| [E-MTAB-567](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-567) | `rnaseq-005` | Paired prostate tumour and adjacent tissue |
| [E-MTAB-8845](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8845) | `rnaseq-006`, `rnaseq-009` to `rnaseq-011` | Arabidopsis genotype and infection factorial study |
| [E-MTAB-9206](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9206) | `rnaseq-007` | ELP3 RNA interference in BT549 cells |
| [E-MTAB-10322](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10322) | `rnaseq-008`, `rnaseq-012`, `rnaseq-013` | Duchenne muscular dystrophy myoblast comparisons |

Each included source exposes an integer gene count matrix and an experiment design. Sample identifiers were checked across both resources. Unanalysed columns were removed by rule, repeated run columns required exact equality, and all contrast or subset choices were recorded before pocket scoring.

## Scientific review

The collection spans three organisms and five annotation profiles. It includes balanced two-group studies, a paired subject design, a blocked factorial analysis, an imbalanced disease comparison and stratified contrasts. Each selected comparison has enough replication and residual degrees of freedom for its declared model.

Candidate pockets were tested in a fixed feature-count grid. The smallest passing sizes range from 500 to 8,000 genes. Every selected pocket passed the diagnostic envelope in two complete retrieval and calculation runs, then passed a full-source versus pocket comparison with DESeq2 1.50.2. Exact results are stored in the object and summarized in `catalog/collection.tsv`.

The paired E-MTAB-567 object retains 14 tumour and 14 adjacent-tissue samples and uses the individual as a blocking factor. E-MTAB-8845 produces two blocked 12-sample contrasts and two six-sample infection strata. E-MTAB-10322 produces one nine-sample disease comparison and two six-sample genotype comparisons. The overlap audit records the shared samples and confirms that no sample identifiers cross source accessions.

## Reference review

The reference check uses the genome and annotation releases named by Expression Atlas processing metadata. Official Ensembl checksum records were pinned for each genome and GTF. Downloaded GTF files passed gzip testing and SHA-256 hashing, and every selected gene identifier was found in its matching annotation.

E-MTAB-7126 was rejected at this stage. Its current count matrix contains identifiers absent from the Ensembl 95 annotation stated on the experiment page, including after checking the corresponding patch and haplotype GTF. The mismatch record is preserved in `evidence/e-mtab-7126-reference-mismatch.json`.

E-MTAB-5477 was rejected earlier because none of the tested pockets up to 20,000 genes met every predeclared threshold. Keeping the failed diagnostic prevents the study from being reconsidered as if it had passed.

## Rights and attribution

The [Expression Atlas licence](https://www.ebi.ac.uk/gxa/licence.html) applies CC BY 4.0 to copyrightable material on the service and requires appropriate credit. Every distributed object carries a GREEN rights record with the evidence URL and review date, plus an attribution record naming the study, accession, provider and changes made.

The earlier GEO shortlist remains in the intake history. Those records are not part of the v1 distributed collection because the GEO disclaimer does not grant blanket redistribution rights. Pasilla remains a link-only object for the same reason. A downloadable source is not treated as an unrestricted licence.

## Scope limits

These objects begin at archive count matrices. They do not independently verify read alignment, feature counting, strandedness or every biological conclusion in the source publication. A user needing read-level quality control, another contrast or the complete study should use the cited archive record.
