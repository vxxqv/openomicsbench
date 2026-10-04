import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.evaluate import evaluate_matrix, evaluate_suite, read_method_receipt


class EvaluateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]
        cls.effects = cls.collection / "datasets/rnaseq/rnaseq-002/expected/reference-effects.tsv.gz"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.a = self.root / "method-a"
        self.b = self.root / "method-b"
        self.a.mkdir()
        self.b.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def add_truth(self, directory: Path):
        shutil.copyfile(self.effects, directory / "rnaseq-002.tsv.gz")
        (directory / "sequence-005.vcf").write_text(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "synthetic_reference\t150\t.\tC\tA\t60\tPASS\t.\n"
            "synthetic_reference\t500\t.\tA\tC\t60\tPASS\t.\n"
            "synthetic_reference\t850\t.\tC\tA\t60\tPASS\t.\n",
            encoding="utf-8",
            newline="\n",
        )

    def test_mixed_assay_evaluation(self):
        self.add_truth(self.a)
        report = evaluate_suite(self.collection, self.a, ["rnaseq-002", "sequence-005"])
        self.assertEqual(report["summary"]["status"], "pass")
        self.assertEqual(report["summary"]["assays"]["bulk_rna_seq"]["passed"], 1)
        self.assertEqual(report["summary"]["assays"]["whole_genome_dna_seq"]["passed"], 1)
        self.assertEqual([item["comparator"] for item in report["results"]], ["differential_expression", "exact_snv"])

    def test_missing_results_fail_closed(self):
        report = evaluate_suite(self.collection, self.a, ["rnaseq-002", "sequence-005"])
        self.assertEqual(report["summary"]["failed"], 2)
        self.assertIn("missing result file", report["results"][0]["reason"])

    def test_method_receipt_is_validated(self):
        receipt = {
            "schema_version": "1.0",
            "name": "Example caller",
            "version": "4.2",
            "runtime_seconds": 12.5,
            "peak_memory_mb": 256.0,
            "threads": 4,
            "parameters": {"mode": "strict"},
        }
        (self.a / "method.json").write_text(json.dumps(receipt), encoding="utf-8")
        loaded = read_method_receipt(self.a)
        self.assertEqual(loaded["runtime_seconds"], 12.5)
        receipt["runtime_seconds"] = -1
        (self.a / "method.json").write_text(json.dumps(receipt), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "invalid method receipt"):
            read_method_receipt(self.a)

    def test_multi_assay_matrix_ranks_complete_method_first(self):
        self.add_truth(self.a)
        shutil.copyfile(self.effects, self.b / "rnaseq-002.tsv.gz")
        report = evaluate_matrix(
            self.collection,
            [("complete", self.a), ("rna-only", self.b)],
            ["rnaseq-002", "sequence-005"],
        )
        self.assertEqual(report["summary"]["benchmarks"], 2)
        self.assertEqual([row["name"] for row in report["leaderboard"]], ["complete", "rna-only"])
        self.assertIn("not pooled", report["ranking"])

    def test_cli_writes_mixed_report_and_uses_exit_two_for_failure(self):
        self.add_truth(self.a)
        output = self.root / "evaluation.json"
        passed = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "evaluate", str(self.a),
                "--id", "rnaseq-002", "--id", "sequence-005", "--json", str(output),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(passed.returncode, 0, passed.stderr)
        self.assertEqual(json.loads(output.read_text())["summary"]["passed"], 2)
        failed = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "evaluate", str(self.b),
                "--id", "sequence-005",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(failed.returncode, 2, failed.stderr)


if __name__ == "__main__":
    unittest.main()
