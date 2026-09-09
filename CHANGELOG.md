# Changes

## 1.0.0

Version 1 introduces 12 compact bulk RNA-seq benchmark objects from seven Expression Atlas studies. The collection covers human, mouse and Arabidopsis data across balanced knockouts, RNA interference, paired tumour samples, blocked factorial comparisons, an imbalanced disease design and genotype-specific strata.

Each object includes an integer pocket matrix, the matching full selected-sample matrix, ordered sample metadata, reference and rights records, a transformation record and exact expected results. Pocket sizes were selected from fixed candidate grids as the smallest size meeting four predeclared preservation thresholds. A second DESeq2 comparison under R 4.5.3 and DESeq2 1.50.2 confirms every selected pocket against its full source matrix.

Reference checks use official Ensembl annotation files and require every selected gene identifier to resolve. E-MTAB-7126 is excluded because its current matrix does not match the annotation stated by the source, and E-MTAB-5477 is excluded because no tested pocket met all thresholds. Their rejection evidence remains in the repository.

The command line can list, inspect, retrieve, validate and report provenance without requiring users to read source code. Validation checks declared files and hashes, exact counts, sample order, design structure, reference status and the four quantitative metrics. Forty-five tests cover normal use and failures such as corrupt transfers, duplicate identifiers, mismatched sample order, confounding and paired-read errors.

A generated JSON and TSV catalog provides the complete machine-readable collection summary. The black-and-white fidelity figure reads the same evidence and shows all objects against their thresholds. Continuous integration rebuilds the fixture, catalog and figure, checks generated drift, validates every registered object, audits sample and prose overlap, and runs the v1 preflight.

Source licences, provider credit and project ownership are separate. Redistributed Expression Atlas material remains under CC BY 4.0 with provider attribution. OpenOmicsBench code is Apache-2.0. Vivaan Patni is the release author on Zenodo, and GitHub development is attributed to `vxxqv`.

Further releases can add read-level objects, more assay types and broader external platform testing without changing the v1 evidence or tag.
