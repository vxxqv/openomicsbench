"""Argument parsing and dispatch for ``omicsbench seq``."""
from pathlib import Path

from .sequence_tools import (
    compare_sequence_files,
    convert_file,
    create_sketch,
    deduplicate_file,
    filter_file,
    kmer_report,
    motif_report,
    sample_file,
    transform_file,
    trim_file,
    validate_file,
)
from .sequences import summarize
from .sequence_analysis import (
    deinterleave_pairs,
    extract_records,
    find_orfs,
    interleave_pairs,
    pair_report,
    qc_report,
)
from .proteins import digest_report, protein_report


MOLECULES = ["auto", "dna", "rna", "protein"]


def add_sequence_parser(subparsers):
    parser = subparsers.add_parser("seq", help="Inspect and process FASTA or FASTQ sequences.", epilog="Example: omicsbench seq stats reads.fastq.gz")
    actions = parser.add_subparsers(dest="seq_command", required=True)

    stats = actions.add_parser("stats", help="Report composition, length, duplication and quality statistics.")
    stats.add_argument("input", type=Path)
    stats.add_argument("--molecule", choices=MOLECULES, default="auto")
    stats.add_argument("--format", choices=["auto", "fasta", "fastq"], default="auto")

    validate = actions.add_parser("validate", help="Strictly validate structure, symbols, identifiers and qualities.")
    validate.add_argument("input", type=Path)
    validate.add_argument("--molecule", choices=MOLECULES, default="auto")
    validate.add_argument("--format", choices=["auto", "fasta", "fastq"], default="auto")

    convert = actions.add_parser("convert", help="Rewrite FASTA or convert FASTQ to FASTA.")
    convert.add_argument("input", type=Path)
    convert.add_argument("output", type=Path)
    convert.add_argument("--to", choices=["fasta", "fastq"], required=True)
    convert.add_argument("--wrap", type=int, default=80)
    convert.add_argument("--force", action="store_true")

    transform = actions.add_parser("transform", help="Reverse-complement, transcribe or translate sequences.")
    transform.add_argument("input", type=Path)
    transform.add_argument("output", type=Path)
    transform.add_argument("--operation", choices=["reverse-complement", "transcribe", "back-transcribe", "translate"], required=True)
    transform.add_argument("--molecule", choices=MOLECULES, default="auto")
    transform.add_argument("--frame", type=int, choices=[-3, -2, -1, 1, 2, 3], default=1)
    transform.add_argument("--trim-stop", action="store_true")
    transform.add_argument("--force", action="store_true")

    filters = actions.add_parser("filter", help="Filter records by length, composition and mean quality.")
    filters.add_argument("input", type=Path)
    filters.add_argument("output", type=Path)
    filters.add_argument("--molecule", choices=MOLECULES, default="auto")
    filters.add_argument("--min-length", type=int, default=0)
    filters.add_argument("--max-length", type=int)
    filters.add_argument("--min-gc", type=float)
    filters.add_argument("--max-gc", type=float)
    filters.add_argument("--max-ambiguity", type=float, default=1.0)
    filters.add_argument("--min-mean-quality", type=float)
    filters.add_argument("--force", action="store_true")

    trim = actions.add_parser("trim", help="Trim fixed ends, a literal adapter and low-quality ends.")
    trim.add_argument("input", type=Path)
    trim.add_argument("output", type=Path)
    trim.add_argument("--left", type=int, default=0)
    trim.add_argument("--right", type=int, default=0)
    trim.add_argument("--quality", type=int)
    trim.add_argument("--adapter")
    trim.add_argument("--min-length", type=int, default=1)
    trim.add_argument("--force", action="store_true")

    sample = actions.add_parser("sample", help="Select records reproducibly by count or fraction.")
    sample.add_argument("input", type=Path)
    sample.add_argument("output", type=Path)
    selection = sample.add_mutually_exclusive_group(required=True)
    selection.add_argument("--count", type=int)
    selection.add_argument("--fraction", type=float)
    sample.add_argument("--seed", type=int, default=0)
    sample.add_argument("--force", action="store_true")

    deduplicate = actions.add_parser("deduplicate", help="Keep the first record for each identifier or sequence.")
    deduplicate.add_argument("input", type=Path)
    deduplicate.add_argument("output", type=Path)
    deduplicate.add_argument("--by", choices=["sequence", "identifier"], default="sequence")
    deduplicate.add_argument("--force", action="store_true")

    kmers = actions.add_parser("kmers", help="Count nucleotide k-mers.")
    kmers.add_argument("input", type=Path)
    kmers.add_argument("--k", type=int, required=True)
    kmers.add_argument("--canonical", action="store_true")
    kmers.add_argument("--top", type=int, default=20)
    kmers.add_argument("--max-distinct", type=int, default=1_000_000)

    motif = actions.add_parser("motif", help="Find overlapping IUPAC nucleotide motifs.")
    motif.add_argument("input", type=Path)
    motif.add_argument("motif")
    motif.add_argument("--both-strands", action="store_true")

    sketch = actions.add_parser("sketch", help="Build a deterministic bottom-k nucleotide sketch.")
    sketch.add_argument("input", type=Path)
    sketch.add_argument("--output", type=Path)
    sketch.add_argument("--k", type=int, default=21)
    sketch.add_argument("--size", type=int, default=1000)
    sketch.add_argument("--no-canonical", action="store_true")
    sketch.add_argument("--force", action="store_true")

    compare = actions.add_parser("compare", help="Estimate sequence similarity from deterministic k-mer sketches.")
    compare.add_argument("left", type=Path)
    compare.add_argument("right", type=Path)
    compare.add_argument("--k", type=int, default=21)
    compare.add_argument("--size", type=int, default=1000)
    compare.add_argument("--no-canonical", action="store_true")

    qc = actions.add_parser("qc", help="Run a combined sequence, quality, complexity and adapter report.")
    qc.add_argument("input", type=Path)
    qc.add_argument("--molecule", choices=MOLECULES, default="auto")
    qc.add_argument("--max-positions", type=int, default=300)
    qc.add_argument("--overrepresented", type=int, default=20)
    qc.add_argument("--adapter", action="append", default=[])

    pairs = actions.add_parser("pair-check", help="Check paired FASTQ counts and read identifiers.")
    pairs.add_argument("left", type=Path)
    pairs.add_argument("right", type=Path)
    pairs.add_argument("--detail-limit", type=int, default=20)

    interleave = actions.add_parser("interleave", help="Interleave two validated paired FASTQ files.")
    interleave.add_argument("left", type=Path)
    interleave.add_argument("right", type=Path)
    interleave.add_argument("output", type=Path)
    interleave.add_argument("--force", action="store_true")

    deinterleave = actions.add_parser("deinterleave", help="Split an interleaved paired FASTQ file.")
    deinterleave.add_argument("input", type=Path)
    deinterleave.add_argument("left", type=Path)
    deinterleave.add_argument("right", type=Path)
    deinterleave.add_argument("--force", action="store_true")

    orfs = actions.add_parser("orfs", help="Find complete or partial open reading frames in six frames.")
    orfs.add_argument("input", type=Path)
    orfs.add_argument("--min-aa", type=int, default=30)
    orfs.add_argument("--include-partial", action="store_true")
    orfs.add_argument("--max-results", type=int, default=1000)

    extract = actions.add_parser("extract", help="Select records and optionally slice their sequences.")
    extract.add_argument("input", type=Path)
    extract.add_argument("output", type=Path)
    extract.add_argument("--id", action="append", default=[])
    extract.add_argument("--ids-file", type=Path)
    extract.add_argument("--start", type=int)
    extract.add_argument("--end", type=int)
    extract.add_argument("--force", action="store_true")

    protein = actions.add_parser("protein-stats", help="Report amino-acid composition and physicochemical estimates.")
    protein.add_argument("input", type=Path)
    protein.add_argument("--ph", type=float, default=7.0)

    digest = actions.add_parser("digest", help="Perform deterministic in-silico protein digestion.")
    digest.add_argument("input", type=Path)
    digest.add_argument("--enzyme", choices=["trypsin", "lys-c", "arg-c", "chymotrypsin"], default="trypsin")
    digest.add_argument("--missed-cleavages", type=int, default=0)
    digest.add_argument("--min-length", type=int, default=1)
    digest.add_argument("--max-length", type=int)
    return parser


def run_sequence_command(args):
    if args.seq_command == "stats":
        return summarize(args.input, args.molecule, args.format)
    if args.seq_command == "validate":
        return validate_file(args.input, args.molecule, args.format)
    if args.seq_command == "convert":
        return convert_file(args.input, args.output, args.to, args.wrap, args.force)
    if args.seq_command == "transform":
        return transform_file(args.input, args.output, args.operation, args.molecule, args.frame, args.trim_stop, args.force)
    if args.seq_command == "filter":
        return filter_file(args.input, args.output, args.molecule, args.min_length, args.max_length, args.min_gc, args.max_gc, args.max_ambiguity, args.min_mean_quality, args.force)
    if args.seq_command == "trim":
        return trim_file(args.input, args.output, args.left, args.right, args.quality, args.adapter, args.min_length, args.force)
    if args.seq_command == "sample":
        return sample_file(args.input, args.output, args.count, args.fraction, args.seed, args.force)
    if args.seq_command == "deduplicate":
        return deduplicate_file(args.input, args.output, args.by, args.force)
    if args.seq_command == "kmers":
        return kmer_report(args.input, args.k, args.canonical, args.top, args.max_distinct)
    if args.seq_command == "motif":
        return motif_report(args.input, args.motif, args.both_strands)
    if args.seq_command == "sketch":
        return create_sketch(args.input, args.output, args.k, args.size, not args.no_canonical, args.force)
    if args.seq_command == "compare":
        return compare_sequence_files(args.left, args.right, args.k, args.size, not args.no_canonical)
    if args.seq_command == "qc":
        return qc_report(args.input, args.molecule, args.max_positions, args.overrepresented, args.adapter)
    if args.seq_command == "pair-check":
        return pair_report(args.left, args.right, args.detail_limit)
    if args.seq_command == "interleave":
        return interleave_pairs(args.left, args.right, args.output, args.force)
    if args.seq_command == "deinterleave":
        return deinterleave_pairs(args.input, args.left, args.right, args.force)
    if args.seq_command == "orfs":
        return find_orfs(args.input, args.min_aa, args.include_partial, args.max_results)
    if args.seq_command == "extract":
        identifiers = set(args.id)
        if args.ids_file:
            identifiers.update(line.strip() for line in args.ids_file.read_text(encoding="utf-8").splitlines() if line.strip())
        return extract_records(args.input, args.output, identifiers, args.start, args.end, args.force)
    if args.seq_command == "protein-stats":
        return protein_report(args.input, args.ph)
    if args.seq_command == "digest":
        return digest_report(args.input, args.enzyme, args.missed_cleavages, args.min_length, args.max_length)
    raise ValueError(f"unknown sequence command: {args.seq_command}")
