# v14 RNA-seq source screening

## Focused four-study pool (27 November 2026)

The v14 pool now prioritizes four non-infection studies with distinct experimental designs. These are intake records, not released benchmark objects. The comparison, sample subset, candidate-size grid, seed and diagnostic thresholds are frozen in the linked source configurations before count-file processing. The official experiment designs were checked for analysed runs and group sizes. Raw-count alignment, study-specific redistribution rights, reference annotation, DESeq2 confirmation and independent use remain open.

| Study | Locked comparison | Selected runs | Design SHA-256 | Intake and configuration |
| --- | --- | ---: | --- | --- |
| E-MTAB-7925 | Cyclic strain versus none at 24 hours, blocked by donor | 12 (6 + 6) | `bc8a148b78f7a798dd4c26030f052d17d70144ef1c385eefadb902fdac058f9f` | `curation/intake/e-mtab-7925.json`, `workflows/e-mtab-7925-source.json` |
| E-MTAB-9639 | FOXK1/2 knockdown versus scrambled siRNA in neuronal stem cells | 6 (3 + 3) | `0aba876967f3e579dac73936897bcd20ecee2b542599e998575138e11a9d0a00` | `curation/intake/e-mtab-9639.json`, `workflows/e-mtab-9639-source.json` |
| E-MTAB-10374 | Drought versus control within the ICE163 accession | 6 (3 + 3) | `4fbb4a5eadf5aadb611463afd1900c9b7b0244f4247cf01500343577ca292f8e` | `curation/intake/e-mtab-10374.json`, `workflows/e-mtab-10374-source.json` |
| E-MTAB-3101 | Myofibroblast differentiation versus control | 6 (3 + 3) | `f973dab0c38aa3a2a392a872616f43ed93437c26e338b3d2ef716e8c798abb92` | `curation/intake/e-mtab-3101.json`, `workflows/e-mtab-3101-source.json` |

The earlier screen below remains as historical intake evidence. It is not the four-study release pool.

The frozen diagnostic and DESeq2 runs selected the smallest passing pockets below. All four passed the predeclared four-metric DESeq2 comparison, and every selected gene ID matched the official Ensembl release recorded by Expression Atlas. These results establish count-matrix and reference compatibility, not biological ground truth. The source licence review found that Expression Atlas applies CC BY 4.0 to its website material; each object needs provider and study attribution in its card.

| Study | Smallest passing pocket | DESeq2 result | Reference check |
| --- | ---: | --- | --- |
| E-MTAB-7925 | 8,000 genes | Pass | Ensembl 107, 8,000/8,000 IDs matched |
| E-MTAB-9639 | 1,000 genes | Pass | Ensembl 95, 1,000/1,000 IDs matched |
| E-MTAB-10374 | 2,000 genes | Pass | Ensembl Genomes 53, 2,000/2,000 IDs matched |
| E-MTAB-3101 | 1,000 genes | Pass | Ensembl 95, 1,000/1,000 IDs matched |

The exact candidate curves, source hashes, DESeq2 outputs and official reference checks are in `evidence/`. Count files remain in staging until the dataset cards, package checks and clean-environment trial are complete.

Checked on 15 November 2026. These are **intake candidates**, not released benchmark objects. The Expression Atlas experiment pages identify four distinct biological designs; the official raw-count and experiment-design download endpoints responded successfully to a metadata request. No candidate has passed the full rights, annotation, contrast, DESeq2, pocket-fidelity or release checks yet.

| Candidate | Why it broadens the collection | Design rows | Design SHA-256 | Source |
| --- | --- | ---: | --- | --- |
| E-MTAB-10451 | Human intestinal organoid Salmonella response at 2.5 and 8 hours; tests a time-stratified infection contrast. | 48 | `0fa4e2769b0ff77a95d5139ec5131b7509112c03e0e15711a3af5e62af7834a2` | [Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10451) |
| E-MTAB-7925 | Human mesenchymal stem cells from six donors with and without cyclic strain; tests donor-aware treatment analysis. | 24 | `bc8a148b78f7a798dd4c26030f052d17d70144ef1c385eefadb902fdac058f9f` | [Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-7925) |
| E-MTAB-9639 | Human embryonic and neuronal stem-cell knockdowns; tests a cell-type-specific perturbation design. | 21 | `0aba876967f3e579dac73936897bcd20ecee2b542599e998575138e11a9d0a00` | [Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9639) |
| E-MTAB-3649 | Human epithelial-cell Listeria infection with six analysed runs; tests a compact infection comparison. | 6 | `82640ec1018984d0ccdcfc95e0491fbae66fd3520c3a3a9feb2d522c47ecde71` | [Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-3649) |

The official endpoints follow `https://www.ebi.ac.uk/gxa/experiments-content/ACCESSION/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts` and `.../ExperimentDesignFile.RnaSeq/experiment-design`. The reported design hashes are over the downloaded experiment-design bytes, so a later change is detectable. The design-row count is not the final number of analysis samples; each candidate still needs the `Analysed` flag and contrast-specific sample subset checked.

Before adding one to `datasets/`, lock a two-level contrast and blocking design; inspect all samples and library types; verify annotation compatibility, study-specific redistribution terms, source citation and reference method; run the archived DESeq2 baseline; apply the predeclared pocket thresholds without retuning after outcomes; then run repository and installed-wheel validation. A failed candidate stays recorded as failed rather than being counted as a benchmark.

## First diagnostic: E-MTAB-3649

The existing `workflows/atlas_diagnostic.py` pipeline downloaded the official counts and design, selected all six analysed runs (three uninfected and three infected), and evaluated the predeclared pocket sizes in `workflows/e-mtab-3649-source.json`. The full source has 58,735 features. The 500-feature candidate **failed** the diagnostic thresholds; the 1,000-feature candidate is the smallest one that passed the diagnostic OLS preservation checks (effect-rank correlation 0.925, top-50 Jaccard 0.613, sample-distance correlation 0.952 and sign concordance 0.941). The staged diagnostic report and source files are outside the released collection. This is an intake result only: DESeq2 confirmation, reference-annotation audit and rights review remain outstanding.
