import csv
import tempfile
import unittest
from pathlib import Path

from omicsbench import Collection


class RNASeqDesignTests(unittest.TestCase):
    def setUp(self):
        self.collection = Collection()

    def test_packs_have_identifiable_designs(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary) / "pack"
            self.collection.pack_rnaseq(folder)
            for row in self.collection.list(assay="bulk_rna_seq", status="validated", kind="real"):
                identifier = row["id"]
                case = folder / identifier
                with self.subTest(id=identifier):
                    report = self.collection.check_rnaseq_design(
                        case / "counts.tsv", case / "samples.tsv", case / "contrasts.tsv", identifier)
                    self.assertEqual(report["summary"]["status"], "pass")
                    self.assertEqual(report["contrasts"][0]["rank"], report["contrasts"][0]["columns"])
                    self.assertTrue(report["sample_order_matches"])

    def test_sample_mismatch_and_confounded_design_fail_visibly(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "counts.tsv"
            samples = root / "samples.tsv"
            contrasts = root / "contrasts.tsv"
            counts.write_text("gene_id\ta\tb\tc\td\nG1\t1\t2\t3\t4\nG2\t2\t3\t4\t5\n", encoding="utf-8")
            samples.write_text("sample\tcondition\tbatch\na\tcontrol\tone\nb\tcontrol\tone\nc\ttreated\ttwo\nd\ttreated\ttwo\n", encoding="utf-8")
            contrasts.write_text("id\tvariable\treference\ttarget\tblocking\ntreat\tcondition\tcontrol\ttreated\tbatch\n", encoding="utf-8")
            report = self.collection.check_rnaseq_design(counts, samples, contrasts)
            self.assertEqual(report["summary"]["status"], "fail")
            self.assertIn("confounded", report["contrasts"][0]["issues"][0])
            samples.write_text(samples.read_text(encoding="utf-8").replace("\nd\ttreated", "\ne\ttreated"), encoding="utf-8")
            report = self.collection.check_rnaseq_design(counts, samples, contrasts)
            self.assertEqual(report["missing_sample_ids"], ["d"])
            self.assertEqual(report["extra_sample_ids"], ["e"])
            self.assertEqual(report["summary"]["status"], "fail")

    def test_matrix_rejects_duplicate_ids_and_fractional_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "counts.tsv"
            samples = root / "samples.tsv"
            contrasts = root / "contrasts.tsv"
            samples.write_text("sample\tcondition\na\tcontrol\nb\tcontrol\nc\ttreated\nd\ttreated\n", encoding="utf-8")
            contrasts.write_text("id\tvariable\treference\ttarget\nx\tcondition\tcontrol\ttreated\n", encoding="utf-8")
            counts.write_text("gene_id\ta\ta\tc\td\nG1\t1\t2\t3\t4\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicated"):
                self.collection.check_rnaseq_design(counts, samples, contrasts)
            counts.write_text("gene_id\ta\tb\tc\td\nG1\t1\t2.5\t3\t4\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-negative integers"):
                self.collection.check_rnaseq_design(counts, samples, contrasts)


if __name__ == "__main__":
    unittest.main()
