import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench import Collection


def write_counts(path, genes, samples, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])
        writer.writerows(([gene, *row] for gene, row in zip(genes, rows)))


class RNASeqChecksTests(unittest.TestCase):
    def setUp(self):
        self.benchmark = Collection().get("rnaseq-002")

    def test_export_verifies_both_tiers(self):
        with tempfile.TemporaryDirectory() as temporary:
            for tier in ("pocket", "full"):
                with self.subTest(tier=tier):
                    directory = Path(temporary) / tier
                    self.benchmark.export_rnaseq(directory, tier)
                    report = self.benchmark.verify_rnaseq_export(directory)
                    self.assertEqual(report["status"], "pass")
                    self.assertEqual(report["count_comparison"]["mismatched_cells"], 0)

    def test_reordered_matrix_passes_and_changed_count_fails(self):
        data = self.benchmark.rnaseq()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reordered.tsv"
            write_counts(path, reversed(data.genes), reversed(data.samples),
                         (reversed(row.tolist()) for row in reversed(data.counts)))
            self.assertEqual(self.benchmark.compare_counts(path)["status"], "pass")
            with path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.reader(stream, delimiter="\t"))
            rows[1][1] = str(int(rows[1][1]) + 1)
            write_counts(path, [row[0] for row in rows[1:]], rows[0][1:], [row[1:] for row in rows[1:]])
            report = self.benchmark.compare_counts(path, detail_limit=1)
            self.assertEqual(report["status"], "fail")
            self.assertEqual(report["mismatched_cells"], 1)
            self.assertEqual(len(report["examples"]), 1)
            self.assertEqual(report["missing_gene_count"], 0)

    def test_missing_ids_and_export_tampering_fail(self):
        data = self.benchmark.rnaseq()
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "export"
            self.benchmark.export_rnaseq(directory)
            metadata = directory / "samples.tsv"
            metadata.write_text(metadata.read_text(encoding="utf-8").replace("wild type genotype", "other", 1), encoding="utf-8")
            report = self.benchmark.verify_rnaseq_export(directory)
            self.assertEqual(report["status"], "fail")
            self.assertTrue(any("samples.tsv" in issue for issue in report["issues"]))
            details = directory / "export.json"
            original = json.loads(details.read_text(encoding="utf-8"))
            original["source"]["accession"] = "wrong"
            details.write_text(json.dumps(original), encoding="utf-8")
            self.assertTrue(any("source" in issue for issue in self.benchmark.verify_rnaseq_export(directory)["issues"]))
            path = Path(temporary) / "partial.tsv"
            write_counts(path, data.genes[1:], data.samples[1:], data.counts[1:, 1:].tolist())
            compared = self.benchmark.compare_counts(path, detail_limit=0)
            self.assertEqual(compared["status"], "fail")
            self.assertEqual(compared["missing_gene_count"], 1)
            self.assertEqual(compared["missing_sample_count"], 1)
            self.assertEqual(compared["missing_genes"], [])
            with self.assertRaisesRegex(ValueError, "detail_limit"):
                self.benchmark.compare_counts(path, detail_limit=-1)

    def test_cli_returns_failure_status_for_changed_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "export"
            self.benchmark.export_rnaseq(directory)
            command = [sys.executable, "-m", "omicsbench", "rnaseq", "check-counts", "rnaseq-002", str(directory / "counts.tsv")]
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
            path = directory / "counts.tsv"
            with path.open("ab") as stream:
                stream.write(b"broken\n")
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 1)
            self.assertEqual(subprocess.run([sys.executable, "-m", "omicsbench", "rnaseq", "verify-export", "rnaseq-002", str(directory)], capture_output=True).returncode, 2)


if __name__ == "__main__":
    unittest.main()
