# Sequence command reference

All sequence commands begin with `omicsbench seq`. Inputs may be plain text or gzip-compressed FASTA and FASTQ files. File type is detected from content rather than the extension.

## Inspection

`stats INPUT` reports record and letter counts, length range, mean, median, N50, L50, composition, GC, ambiguity and exact duplicates. FASTQ reports also include the Phred+33 range, mean quality and Q20/Q30 fractions.

`validate INPUT` performs the same complete read while enforcing format structure, unique identifiers, a single declared alphabet, printable qualities and equal sequence and quality lengths.

`qc INPUT` adds per-position FASTQ metrics, normalized dinucleotide complexity, overrepresented exact sequences and optional literal adapter counts. Repeat `--adapter` to scan more than one adapter.

## Record operations

`convert INPUT OUTPUT --to fasta` converts FASTQ to FASTA or rewrites FASTA with consistent wrapping. FASTA to FASTQ is rejected because quality values cannot be recovered.

`filter INPUT OUTPUT` accepts `--min-length`, `--max-length`, `--min-gc`, `--max-gc`, `--max-ambiguity` and `--min-mean-quality`. A record is counted once under the first failed rule in the JSON report.

`trim INPUT OUTPUT` accepts fixed `--left` and `--right` removal, a literal `--adapter`, end `--quality` trimming and `--min-length`. Quality trimming requires FASTQ.

`sample INPUT OUTPUT` requires either `--count` or `--fraction`. Selection uses a SHA-256 value derived from the seed, record order and identifier. The same input and seed give the same subset without depending on Python's random state.

`deduplicate INPUT OUTPUT` keeps the first record for each exact sequence by default. Use `--by identifier` to deduplicate names instead.

`extract INPUT OUTPUT` accepts repeated `--id`, an `--ids-file`, and optional `--start` and `--end` coordinates. Coordinates are zero-based and half-open. Omitting identifiers selects every record.

## Paired FASTQ

`pair-check R1 R2` requires the same number of records and matching identifiers in order. Common `/1`, `/2`, `.1` and `.2` mate suffixes are normalized for comparison.

`interleave R1 R2 OUTPUT` runs the pair check before writing alternating mates. `deinterleave INPUT R1 R2` requires an even record count and validates every adjacent pair before writing either final output.

## Sequence transformations

`transform INPUT OUTPUT --operation reverse-complement` supports DNA and RNA. FASTQ qualities are reversed with the sequence.

`--operation transcribe` converts DNA `T` to RNA `U`. `--operation back-transcribe` converts RNA `U` to DNA `T`. `--operation translate` accepts frames `1`, `2`, `3`, `-1`, `-2` and `-3`; ambiguous codons become `X`.

`orfs INPUT` searches all six DNA or RNA frames. Complete ORFs begin with `ATG` and end at `TAA`, `TAG` or `TGA`. Use `--include-partial` to retain starts that reach the end without a stop. Coordinates refer to the original input sequence and use zero-based, half-open intervals.

## Protein analysis

`protein-stats INPUT` reports amino-acid composition, average molecular weight, Kyte-Doolittle mean hydropathy, net charge at a chosen pH, estimated isoelectric point, aromatic fraction and the reduced-protein extinction coefficient at 280 nm. Molecular weight is left unset when a sequence contains an ambiguous residue rather than assigning an invented mass.

`digest INPUT` supports trypsin, Lys-C, Arg-C and chymotrypsin. It records zero-based, half-open peptide coordinates and accepts missed-cleavage and peptide-length limits. Trypsin, Lys-C and Arg-C respect the proline exception used by the command.

## Motifs, k-mers and similarity

`motif INPUT MOTIF` finds overlapping IUPAC nucleotide motifs. `--both-strands` also scans the reverse complement. Each result contains the record, interval and strand.

`kmers INPUT --k K` counts unambiguous nucleotide k-mers. `--canonical` merges each DNA k-mer with its reverse complement. Ambiguous windows are skipped. `--max-distinct` is a deliberate memory guard.

`sketch INPUT` hashes canonical k-mers with BLAKE2b and retains the smallest values. `compare LEFT RIGHT` builds compatible sketches and reports intersection size, Jaccard similarity, containment and Mash-style distance. This is a lightweight similarity screen, not an alignment or taxonomic classification.

## DNA-seq truth object

`sequence-005` contains a 1,000-base synthetic reference, 20 deterministic paired reads and three declared single-nucleotide variants. Collection validation checks pair order, exact sequence summaries, truth-table structure, reference alleles and read support for the alternate-allele context. It is intended for regression tests around DNA-seq file handling and truth-set plumbing. It does not measure performance on repeats, indels, structural variants or realistic instrument errors.

## Output rules

Reports are JSON. Sequence outputs preserve identifiers and descriptions. FASTQ operations preserve qualities unless the result is no longer a read, such as translation. Existing outputs are protected unless `--force` is present. A failed operation removes its temporary file and does not leave a partial new output.
