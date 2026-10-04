# Changes

## 3.0.0

Version 3 introduces a unified evaluator for bulk RNA-seq effects and small-variant calls. One run can apply the correct comparator to each selected benchmark, require every expected result, preserve assay-specific metrics and return one CI-ready status. Multi-assay matrices compare several methods without combining RNA and DNA measurements into an artificial score.

The synthetic DNA-seq object now has a user-facing VCF comparator. It accepts VCF 4.x and gzip-compressed VCF, expands multi-allelic records, accounts for filtered alleles, verifies submitted reference alleles and reports exact-allele precision, recall, F1, true positives, false positives and false negatives. Its stated scope is single-nucleotide substitutions; it does not claim haplotype-aware or genotype-aware comparison.

An optional `method.json` receipt records a tool name and version together with its command, container, source revision, parameters, runtime, memory and thread count. These details stay attached to evaluation and leaderboard reports.

HTML reports have been rebuilt as responsive, self-contained dashboards with summary cards, status charts, method rankings, assay breakdowns, a method-by-benchmark coverage matrix and compact case metrics. CSV, Markdown, JUnit and regression reports now carry the DNA metrics as well as the existing RNA measures.

Evaluation evidence can be written as a deterministic ZIP. Each crate contains normalized JSON, CSV and HTML reports, the exact submitted result files, optional method receipts, SHA-256 checksums and RO-Crate 1.3 metadata. Verification checks safe paths, the complete checksum inventory, every member digest, the report contract and the RO-Crate root entities without extracting the archive.

## 2.2.0

Benchmark matrices compare two or more analysis methods or parameter sets against the same RNA-seq selection. Results include every method and benchmark case, per-method pass rates and mean preservation metrics, and a deterministic leaderboard without hiding failed or missing inputs.

Regression gates compare an accepted JSON suite report with a candidate report. They fail on new benchmark failures, missing baseline cases or decreases in effect-rank correlation, top-gene overlap, sign agreement and gene coverage beyond a declared tolerance.

Suite commands can now write flat CSV and self-contained HTML alongside JSON, Markdown and JUnit XML. HTML reports contain no external scripts or network dependencies.

Portable bundles package selected benchmarks, their complete declared file inventories and licence material into deterministic ZIP archives. Verification checks safe paths, exact inventory membership, byte counts, SHA-256 digests and embedded dataset manifests without extracting the archive.

## 2.1.0

Benchmark suites can now validate several objects in one run or compare a directory of differential-expression outputs across the biological collection. Exact IDs and assay filters make the selection explicit, missing comparison files fail closed, and a problem in one object does not prevent the remaining cases from running.

Suite results are printed as JSON and can also be written as full JSON, a compact Markdown table or JUnit XML for continuous integration. Reports are written atomically, existing files are protected by default, and exit status 2 distinguishes completed benchmark failures from invalid command input.

## 2.0.0

OpenOmicsBench now handles DNA, RNA and protein sequences in FASTA and FASTQ files, including gzip-compressed input and output. The new command suite provides strict validation, combined QC, length and composition statistics, positional quality, adapter scans, trimming, filtering, deterministic sampling, deduplication, paired-read checks, interleaving, extraction, reverse complements, transcription, translation, ORF discovery, IUPAC motif searches, k-mer counts and deterministic sequence sketches.

Five synthetic sequence benchmarks cover DNA FASTA, RNA FASTA, protein FASTA, paired FASTQ and a miniature DNA-seq truth set. The DNA-seq object includes a reference, 20 read pairs and three known SNVs; validation checks the truth alleles against the reference and confirms alternate-allele support in the reads. Each object has a version 2 manifest, exact hashes, an expected summary profile and a deterministic rebuild. The 12 biological bulk RNA-seq objects and their version 1 evidence remain unchanged.

Protein FASTA inputs can now be summarized by amino-acid composition, average molecular weight, hydropathy, charge, estimated isoelectric point, aromatic fraction and extinction coefficient. Deterministic digestion supports trypsin, Lys-C, Arg-C and chymotrypsin with missed-cleavage and peptide-length controls.

Fourteen assay profiles cover short-read, long-read, reference and protein inputs. Each profile states the expected inputs, the checks available locally and the point where a reference-aware or database-backed workflow becomes necessary. The plan command returns executable local validation and preprocessing steps.

All file-writing sequence commands protect existing outputs unless replacement is explicitly requested. Reports use stable JSON, sampling is seed-controlled and similarity sketches record their algorithm and parameters.

## 1.1.0

OpenOmicsBench can now compare a CSV or TSV table of gene-level log2 fold changes with the full-source DESeq2 reference for any biological object. The report records effect-rank correlation, top-gene overlap, sign agreement, gene coverage, missing and unexpected identifiers, the declared thresholds and a pass or fail decision suitable for automated tests.

Finite DESeq2 reference effects are bundled as deterministic compressed tables. Their source outputs are checked against the existing evidence hashes before packaging. The synthetic fixture remains available for infrastructure tests but is not presented as a DESeq2 comparison reference.

## 1.0.4

Cache verification now rejects completion receipts that name a dataset tier which is not declared by the current manifest. This keeps verification fail-closed when local cache metadata is altered or stale.

## 1.0.3

Version 1.0.3 is the first PyPI release. The publishing workflow now fetches the complete tag history needed to certify the frozen version 1 release before upload.

## 1.0.2

Version 1.0.2 retains the complete bundled collection from 1.0.1 and makes the synthetic fixture rebuild test portable across supported operating systems and Python environments. Its release tag was not uploaded to PyPI because the publishing checkout lacked the historical v1.0.0 tag required by certification.

## 1.0.1

OpenOmicsBench gained standard Python source and wheel distributions. The wheel contains the complete version 1 collection, so discovery, inspection, retrieval and validation work from any directory without a repository checkout. The command still accepts `--root` and `OMICSBENCH_ROOT` for local or extended collections.

## 1.0.0

Version 1 introduces 12 compact bulk RNA-seq benchmark objects from seven Expression Atlas studies. The collection covers human, mouse and Arabidopsis data across balanced knockouts, RNA interference, paired tumour samples, blocked factorial comparisons, an imbalanced disease design and genotype-specific strata.

Each object includes an integer pocket matrix, the matching full selected-sample matrix, ordered sample metadata, reference and rights records, a transformation record and exact expected results. Pocket sizes were selected from fixed candidate grids as the smallest size meeting four predeclared preservation thresholds. A second DESeq2 comparison under R 4.5.3 and DESeq2 1.50.2 confirms every selected pocket against its full source matrix.

Reference checks use official Ensembl annotation files and require every selected gene identifier to resolve. E-MTAB-7126 is excluded because its current matrix does not match the annotation stated by the source, and E-MTAB-5477 is excluded because no tested pocket met all thresholds. Their rejection evidence remains in the repository.

The command line can list, inspect, retrieve, validate and report provenance without requiring users to read source code. Validation checks declared files and hashes, exact counts, sample order, design structure, reference status and the four quantitative metrics. Forty-five tests cover normal use and failures such as corrupt transfers, duplicate identifiers, mismatched sample order, confounding and paired-read errors.

A generated JSON and TSV catalog provides the complete machine-readable collection summary. The black-and-white fidelity figure reads the same evidence and shows all objects against their thresholds. Continuous integration rebuilds the fixture, catalog and figure, checks generated drift, validates every registered object, audits sample and prose overlap, and runs the v1 preflight.

Source licences, provider credit and project ownership are separate. Redistributed Expression Atlas material remains under CC BY 4.0 with provider attribution. OpenOmicsBench code is Apache-2.0. Vivaan Patni is the release author on Zenodo, and GitHub development is attributed to `vxxqv`.

Further releases can add read-level objects, more assay types and broader external platform testing without changing the v1 evidence or tag.

The archived version 1.0.0 release is available at https://doi.org/10.5281/zenodo.22679414.
