# v14 RNA-seq source screening

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
