import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.matrix import compare_matrix, parse_method


class MatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]
        cls.reference = cls.collection / "datasets" / "rnaseq" / "rnaseq-002" / "expected" / "reference-effects.tsv.gz"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.a = self.root / "method-a"
        self.b = self.root / "method-b"
        self.a.mkdir()
        self.b.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def test_parse_method(self):
        self.assertEqual(parse_method("deseq2=results/deseq2"), ("deseq2", Path("results/deseq2")))
        with self.assertRaisesRegex(ValueError, "NAME=RESULTS_DIRECTORY"):
            parse_method("deseq2")

    def test_matrix_requires_two_unique_methods(self):
        with self.assertRaisesRegex(ValueError, "at least two"):
            compare_matrix(self.collection, [("one", self.a)], ["rnaseq-002"])
        with self.assertRaisesRegex(ValueError, "unique"):
            compare_matrix(self.collection, [("one", self.a), ("one", self.b)], ["rnaseq-002"])

    def test_matrix_ranks_complete_method_before_missing_method(self):
        (self.a / "rnaseq-002.tsv.gz").write_bytes(self.reference.read_bytes())
        report = compare_matrix(self.collection, [("complete", self.a), ("missing", self.b)], ["rnaseq-002"])
        self.assertEqual(report["summary"]["methods"], 2)
        self.assertEqual(report["summary"]["benchmarks"], 1)
        self.assertEqual(report["summary"]["failed"], 1)
        self.assertEqual([row["name"] for row in report["leaderboard"]], ["complete", "missing"])
        self.assertEqual(report["leaderboard"][0]["mean_metrics"]["spearman_logfc"], 1.0)

    def test_matrix_cli_returns_scientific_failure_status(self):
        (self.a / "rnaseq-002.tsv.gz").write_bytes(self.reference.read_bytes())
        completed = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "matrix",
                "--method", f"complete={self.a}", "--method", f"missing={self.b}", "--id", "rnaseq-002",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 2, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)["summary"]["failed"], 1)


if __name__ == "__main__":
    unittest.main()
