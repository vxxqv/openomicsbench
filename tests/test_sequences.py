import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.sequences import (
    SequenceRecord,
    compare_sketches,
    count_kmers,
    detect_format,
    find_motifs,
    infer_molecule,
    read_records,
    reverse_complement,
    sketch,
    summarize,
    translate,
    trim_record,
    write_records,
)
from omicsbench.sequence_tools import (
    convert_file,
    deduplicate_file,
    filter_file,
    sample_file,
    transform_file,
    trim_file,
)


class SequenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding="utf-8", newline="\n")
        return path

    def test_multiline_fasta(self):
        path = self.write("input.fa", ">a first\nACGT\nNN\n>b\nUUAC\n")
        records = list(read_records(path))
        self.assertEqual(records[0], SequenceRecord("a", "first", "ACGTNN"))
        self.assertEqual(records[1].sequence, "UUAC")

    def test_multiline_fastq(self):
        path = self.write("reads.fastq", "@r1 note\nACGT\nAC\n+\nIIII\nII\n@r2\nNN\n+\n!!\n")
        records = list(read_records(path))
        self.assertEqual(records[0], SequenceRecord("r1", "note", "ACGTAC", "IIIIII"))
        self.assertEqual(records[1].quality, "!!")

    def test_gzip_input_and_output(self):
        source = self.root / "reads.fq.gz"
        with gzip.open(source, "wt", encoding="utf-8", newline="\n") as stream:
            stream.write("@r\nACGT\n+\nIIII\n")
        self.assertEqual(detect_format(source), "fastq")
        target = self.root / "copy.fq.gz"
        self.assertEqual(write_records(read_records(source), target), 1)
        self.assertEqual(list(read_records(target))[0].sequence, "ACGT")

    def test_duplicate_identifiers_rejected(self):
        path = self.write("bad.fa", ">same\nAC\n>same\nGT\n")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            list(read_records(path))

    def test_truncated_fastq_rejected(self):
        path = self.write("bad.fq", "@r\nACGT\n+\nIII\n")
        with self.assertRaisesRegex(ValueError, "lengths differ"):
            list(read_records(path))

    def test_molecule_inference(self):
        self.assertEqual(infer_molecule("ACGTRYN"), "dna")
        self.assertEqual(infer_molecule("ACGURYN"), "rna")
        self.assertEqual(infer_molecule("MKWVTF"), "protein")
        with self.assertRaisesRegex(ValueError, "both T and U"):
            infer_molecule("AUTG")

    def test_iupac_reverse_complement(self):
        self.assertEqual(reverse_complement("ACGTRYMKWSBDHVN"), "NBDHVSWMKRYACGT")
        self.assertEqual(reverse_complement("ACGURYMKWSBDHVN", "rna"), "NBDHVSWMKRYACGU")

    def test_translation_frames(self):
        self.assertEqual(translate("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"), "MAIVMGR*KGAR*")
        self.assertEqual(translate("AATGGCC", 2), "MA")
        self.assertEqual(translate("TGGCCAT", -1), "MA")
        self.assertEqual(translate("ATGNAA"), "MX")

    def test_summary_fasta(self):
        path = self.write("dna.fa", ">a\nACGT\n>b\nGGNN\n>c\nACGT\n")
        report = summarize(path)
        self.assertEqual(report["format"], "fasta")
        self.assertEqual(report["molecule"], "dna")
        self.assertEqual(report["records"], 3)
        self.assertEqual(report["length"]["n50"], 4)
        self.assertEqual(report["duplicate_sequences"], 1)
        self.assertAlmostEqual(report["gc_fraction"], 0.6)

    def test_summary_fastq_quality(self):
        path = self.write("reads.fq", "@a\nAC\n+\nI5\n@b\nGT\n+\n+!\n")
        report = summarize(path, "dna")
        self.assertEqual(report["quality"]["minimum"], 0)
        self.assertEqual(report["quality"]["maximum"], 40)
        self.assertEqual(report["quality"]["q30_fraction"], 0.25)

    def test_write_fasta_wrap(self):
        target = self.root / "out.fa"
        write_records([SequenceRecord("x", "description", "ACGTAC")], target, "fasta", wrap=4)
        self.assertEqual(target.read_text(), ">x description\nACGT\nAC\n")

    def test_fastq_output_requires_quality(self):
        with self.assertRaisesRegex(ValueError, "requires qualities"):
            write_records([SequenceRecord("x", "", "AC")], self.root / "out.fq", "fastq")

    def test_trim_coordinates_adapter_and_quality(self):
        record = SequenceRecord("r", "", "AACCGGADAPT", "!!IIIIIIIIII")
        trimmed = trim_record(record, left=1, quality=20, adapter="ADAPT")
        self.assertEqual(trimmed.sequence, "CCGG")
        self.assertEqual(trimmed.quality, "IIII")

    def test_overlapping_motif_and_reverse_strand(self):
        result = find_motifs([SequenceRecord("x", "", "ATATAT")], "ATA", both_strands=True)
        self.assertEqual([(hit["start"], hit["strand"]) for hit in result["hits"]], [(0, "+"), (2, "+"), (1, "-"), (3, "-")])

    def test_iupac_motif(self):
        result = find_motifs([SequenceRecord("x", "", "AAGAGG")], "ARG")
        self.assertEqual(result["count"], 2)

    def test_canonical_kmers(self):
        counts = count_kmers([SequenceRecord("x", "", "ACGTT")], 3, canonical=True)
        self.assertEqual(counts, {"ACG": 2, "AAC": 1})

    def test_ambiguous_kmers_are_skipped(self):
        counts = count_kmers([SequenceRecord("x", "", "ACNTA")], 2)
        self.assertEqual(counts, {"AC": 1, "TA": 1})

    def test_sketch_is_deterministic(self):
        records = [SequenceRecord("x", "", "ACGTTACGTT")]
        self.assertEqual(sketch(records, k=3, size=4), sketch(records, k=3, size=4))

    def test_sketch_comparison(self):
        left = sketch([SequenceRecord("x", "", "ACGTTACGTT")], k=3, size=20)
        same = compare_sketches(left, left)
        self.assertEqual(same["jaccard"], 1.0)
        self.assertEqual(same["mash_distance"], 0.0)
        right = sketch([SequenceRecord("y", "", "TTTTTAAAAA")], k=3, size=20)
        self.assertLess(compare_sketches(left, right)["jaccard"], 1.0)

    def test_incompatible_sketches_rejected(self):
        left = sketch([SequenceRecord("x", "", "ACGT")], k=2)
        right = sketch([SequenceRecord("x", "", "ACGT")], k=3)
        with self.assertRaisesRegex(ValueError, "differ in k"):
            compare_sketches(left, right)

    def test_convert_fastq_to_fasta(self):
        source = self.write("reads.fq", "@a\nACGT\n+\nIIII\n")
        output = self.root / "reads.fa"
        report = convert_file(source, output, "fasta", 80, False)
        self.assertEqual(report["records"], 1)
        self.assertEqual(output.read_text(), ">a\nACGT\n")

    def test_convert_does_not_invent_qualities(self):
        source = self.write("input.fa", ">a\nACGT\n")
        with self.assertRaisesRegex(ValueError, "invented"):
            convert_file(source, self.root / "output.fq", "fastq", 80, False)

    def test_atomic_output_requires_force(self):
        source = self.write("input.fa", ">a\nACGT\n")
        output = self.write("output.fa", "existing\n")
        with self.assertRaisesRegex(ValueError, "--force"):
            convert_file(source, output, "fasta", 80, False)
        convert_file(source, output, "fasta", 80, True)
        self.assertEqual(output.read_text(), ">a\nACGT\n")

    def test_transform_reverse_complement_preserves_fastq(self):
        source = self.write("reads.fq", "@a\nACGTN\n+\n!5I?+\n")
        output = self.root / "reverse.fq"
        transform_file(source, output, "reverse-complement", "dna", 1, False, False)
        record = list(read_records(output))[0]
        self.assertEqual(record.sequence, "NACGT")
        self.assertEqual(record.quality, "+?I5!")

    def test_transform_translation_drops_quality(self):
        source = self.write("reads.fq", "@a\nATGGCC\n+\nIIIIII\n")
        output = self.root / "protein.fa"
        transform_file(source, output, "translate", "dna", 1, False, False)
        self.assertEqual(list(read_records(output))[0].sequence, "MA")

    def test_filter_reports_first_reason(self):
        source = self.write("reads.fq", "@short\nAC\n+\nII\n@lowq\nACGT\n+\n!!!!\n@keep\nGCGC\n+\nIIII\n")
        output = self.root / "filtered.fq"
        report = filter_file(source, output, "dna", 4, None, 0.5, None, 0, 30, False)
        self.assertEqual(report["kept_records"], 1)
        self.assertEqual(report["removed_by_first_reason"]["short"], 1)
        self.assertEqual(report["removed_by_first_reason"]["low_quality"], 1)

    def test_trim_discards_short_records(self):
        source = self.write("reads.fq", "@a\nAACCGG\n+\nIIIIII\n@b\nAA\n+\nII\n")
        output = self.root / "trimmed.fq"
        report = trim_file(source, output, 1, 1, None, None, 3, False)
        self.assertEqual(report["kept_records"], 1)
        self.assertEqual(report["discarded_short"], 1)

    def test_sampling_is_reproducible_and_ordered(self):
        source = self.write("input.fa", "".join(f">r{i}\nACGT{'A' if i % 2 == 0 else 'C'}\n" for i in range(10)))
        first, second = self.root / "one.fa", self.root / "two.fa"
        sample_file(source, first, 4, None, 17, False)
        sample_file(source, second, 4, None, 17, False)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        identifiers = [record.identifier for record in read_records(first)]
        self.assertEqual(identifiers, sorted(identifiers, key=lambda value: int(value[1:])))

    def test_deduplicate_keeps_first_record(self):
        source = self.write("input.fa", ">a\nACGT\n>b\nACGT\n>c\nGGGG\n")
        output = self.root / "unique.fa"
        report = deduplicate_file(source, output, "sequence", False)
        self.assertEqual(report["removed_records"], 1)
        self.assertEqual([record.identifier for record in read_records(output)], ["a", "c"])

    def test_sequence_cli_commands(self):
        source = self.write("input.fa", ">a\nACGTACGT\n>b\nACGTTCGT\n")
        commands = [
            ["stats", str(source)],
            ["validate", str(source), "--molecule", "dna"],
            ["kmers", str(source), "--k", "3", "--canonical"],
            ["motif", str(source), "ACG", "--both-strands"],
            ["sketch", str(source), "--k", "3", "--size", "10"],
            ["compare", str(source), str(source), "--k", "3", "--size", "10"],
        ]
        for command in commands:
            with self.subTest(command=command):
                result = subprocess.run([sys.executable, "-m", "omicsbench", "seq", *command], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                json.loads(result.stdout)

    def test_sequence_cli_file_pipeline(self):
        source = self.write("reads.fq", "@a\nAACCGG\n+\nIIIIII\n@b\nAACCGT\n+\nIIIIII\n")
        trimmed = self.root / "trimmed.fq"
        sampled = self.root / "sampled.fq"
        commands = [
            ["trim", str(source), str(trimmed), "--left", "1", "--min-length", "4"],
            ["sample", str(trimmed), str(sampled), "--count", "1", "--seed", "9"],
        ]
        for command in commands:
            result = subprocess.run([sys.executable, "-m", "omicsbench", "seq", *command], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            json.loads(result.stdout)
        self.assertEqual(len(list(read_records(sampled))), 1)


if __name__ == "__main__":
    unittest.main()
