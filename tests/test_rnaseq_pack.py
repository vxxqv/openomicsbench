import csv
import json
import tempfile
import unittest
from pathlib import Path

from omicsbench import Collection
from omicsbench.compare import read_effects
from omicsbench.rnaseq_pack import normalize_effects


class RNASeqPackTests(unittest.TestCase):
    def setUp(self):
        self.collection = Collection()

    def test_pack_carries_aligned_ids_and_the_declared_contrast(self):
        with tempfile.TemporaryDirectory() as temporary:
            pack = Path(temporary) / "pack"
            result = self.collection.pack_rnaseq(pack, ["rnaseq-005", "rnaseq-002"])
            self.assertEqual(result["summary"]["passed"], 2)
            self.assertEqual(self.collection.verify_rnaseq_pack(pack)["summary"]["status"], "pass")
            for identifier in ("rnaseq-002", "rnaseq-005"):
                folder = pack / identifier
                with (folder / "counts.tsv").open(encoding="utf-8", newline="") as stream:
                    count_header = next(csv.reader(stream, delimiter="\t"))
                with (folder / "samples.tsv").open(encoding="utf-8", newline="") as stream:
                    sample_rows = list(csv.DictReader(stream, delimiter="\t"))
                with (folder / "contrasts.tsv").open(encoding="utf-8", newline="") as stream:
                    contrasts = list(csv.DictReader(stream, delimiter="\t"))
                self.assertEqual(count_header[1:], [row["sample"] for row in sample_rows])
                self.assertEqual(len(contrasts), 1)
                self.assertEqual({contrasts[0]["reference"], contrasts[0]["target"]},
                                 {row["condition"] for row in sample_rows})
                self.assertEqual(contrasts[0]["blocking"], "batch" if identifier == "rnaseq-005" else "")
                self.assertNotIn("expected", {path.name for path in folder.iterdir()})
            with self.assertRaises(FileExistsError):
                self.collection.pack_rnaseq(pack, ["rnaseq-002"])

    def test_pack_detects_file_and_design_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            pack = Path(temporary) / "pack"
            self.collection.pack_rnaseq(pack, ["rnaseq-002"])
            samples = pack / "rnaseq-002" / "samples.tsv"
            original = samples.read_bytes()
            samples.write_bytes(original.replace(b"wild type genotype", b"wrong condition", 1))
            with self.assertRaisesRegex(ValueError, "checksum"):
                self.collection.verify_rnaseq_pack(pack)
            samples.write_bytes(original)
            contrast = pack / "rnaseq-002" / "contrasts.tsv"
            original = contrast.read_bytes()
            contrast.write_bytes(original.replace(b"wild type genotype", b"wrong condition", 1))
            with self.assertRaisesRegex(ValueError, "checksum"):
                self.collection.verify_rnaseq_pack(pack)
            contrast.write_bytes(original)
            self.assertEqual(self.collection.verify_rnaseq_pack(pack)["summary"]["status"], "pass")

    def test_pack_rejects_unsupported_selection_without_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "bad"
            with self.assertRaisesRegex(ValueError, "only validated biological"):
                self.collection.pack_rnaseq(output, ["sequence-005"])
            self.assertFalse(output.exists())
            with self.assertRaisesRegex(ValueError, "unique"):
                self.collection.pack_rnaseq(output, ["rnaseq-002", "rnaseq-002"])
            self.assertFalse(output.exists())

    def test_native_result_conversion_and_sign_control(self):
        benchmark = self.collection.get("rnaseq-002")
        reference = read_effects(benchmark.file(role="effects", tier="expected"))
        genes = list(reference)[:100]
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "deseq2.csv"
            with source.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("", "log2FoldChange", "padj"))
                writer.writerows((gene, -reference[gene], 0.01) for gene in genes)
                writer.writerow(("ENSG_NOT_FINITE", "NA", "NA"))
            output = Path(temporary) / "rnaseq-002.tsv"
            result = normalize_effects(benchmark, source, output, reverse_sign=True)
            self.assertEqual(result["finite_effects"], 100)
            self.assertEqual(result["missing_effects_dropped"], 1)
            self.assertEqual(result["shared_reference_genes"], 100)
            self.assertEqual(result["columns"], {"gene_id": "", "log2_fold_change": "log2FoldChange"})
            self.assertEqual(read_effects(output), {gene: reference[gene] for gene in genes})
            self.assertTrue(result["sign_reversed"])
            with self.assertRaises(FileExistsError):
                normalize_effects(benchmark, source, output)

    def test_native_result_conversion_rejects_ambiguity_and_bad_rows(self):
        benchmark = self.collection.get("rnaseq-002")
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "edgeR.tsv"
            output = Path(temporary) / "out.tsv"
            source.write_text("gene_id\tGeneid\tlogFC\nA\tA\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ambiguous"):
                normalize_effects(benchmark, source, output)
            self.assertFalse(output.exists())
            result = normalize_effects(benchmark, source, output, gene_column="gene_id")
            self.assertEqual(result["finite_effects"], 1)
            self.assertEqual(read_effects(output), {"A": 1.0})
            output.unlink()
            source.write_text("gene_id\tlogFC\nA\t1\nA\t2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                normalize_effects(benchmark, source, output)
            self.assertFalse(output.exists())
            source.write_text("gene_id\tlogFC\nA\tinf\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "finite"):
                normalize_effects(benchmark, source, output)
            self.assertFalse(output.exists())

    def test_direct_scoring_keeps_contrast_direction_explicit(self):
        benchmark = self.collection.get("rnaseq-002")
        reference = read_effects(benchmark.file(role="effects", tier="expected"))
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "native-limma.tsv"
            with source.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("Geneid", "logFC", "P.Value"))
                writer.writerows((gene, format(-value, ".17g"), "0.01") for gene, value in reference.items())
            wrong = benchmark.score_native_effects(source)
            self.assertEqual(wrong["status"], "fail")
            self.assertEqual(wrong["conversion"]["sign_reversed"], False)
            correct = benchmark.score_native_effects(source, reverse_sign=True)
            self.assertEqual(correct["status"], "pass")
            self.assertEqual(correct["comparison"]["genes"]["shared"], len(reference))
            self.assertEqual(correct["comparison"]["input"], str(source))


if __name__ == "__main__":
    unittest.main()
