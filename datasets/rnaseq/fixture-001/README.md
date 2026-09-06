# Synthetic RNA infrastructure fixture

Eight generated libraries exercise parsing, sample order, feature reduction, batch design and cache integrity. No patient or experimental data are included.

Run `omicsbench validate fixture-001` from the repository root. Both data tiers contain the smallest candidate that passes the diagnostic profile; they deliberately share content. This fixture is smaller than a normal pocket and is excluded from the biological collection count.

The baseline fits log2(CPM + 1) to condition and batch. It does not calculate differential-expression p-values and is not a DESeq2 substitute. Exact sample order and unchanged integer counts are mandatory. Thresholds and the complete candidate curve are stored under expected.

The generator and seed are recorded in provenance/transform.json. Cite the software version when reporting a test; do not cite this fixture as a biological study.
