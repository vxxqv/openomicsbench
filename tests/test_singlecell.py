import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.registry import lookup
from omicsbench.singlecell import clustering_metrics, compare_labels, read_matrix_market, summarize_matrix
from omicsbench.validate import validate


ROOT = Path(__file__).resolve().parents[1]


class SingleCellTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model, cls.folder = lookup(ROOT, "singlecell-001")

    def test_fixture_validates_with_exact_qc_summary(self):
        report = validate(self.model, self.folder)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["summary"]["cells"], 90)
        self.assertEqual(report["summary"]["features"], 60)
        self.assertEqual(report["summary"]["nonzero"], 3996)
        self.assertEqual(report["summary"]["zero_count_cells"], 0)

    def test_matrix_summary_checks_companion_dimensions(self):
        with tempfile.TemporaryDirectory() as value:
            barcodes = Path(value) / "barcodes.tsv"
            barcodes.write_text("one\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "dimensions do not match"):
                summarize_matrix(
                    self.folder / "nano/matrix.mtx",
                    barcodes,
                    self.folder / "nano/features.tsv",
                    feature_columns=3,
                )

    def test_matrix_parser_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as value:
            matrix = Path(value) / "matrix.mtx"
            matrix.write_text(
                "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 3\n1 1 4\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate matrix coordinate"):
                read_matrix_market(matrix)

    def test_clustering_metrics_ignore_label_names(self):
        metrics = clustering_metrics(["a", "a", "b", "b"], ["2", "2", "1", "1"])
        self.assertEqual(metrics["adjusted_rand_index"], 1.0)
        self.assertEqual(metrics["normalized_mutual_information"], 1.0)
        self.assertEqual(metrics["purity"], 1.0)
        trivial = clustering_metrics(["one", "one"], ["renamed", "renamed"])
        self.assertEqual(trivial["adjusted_rand_index"], 1.0)

    def test_exact_labels_pass_and_partial_labels_fail(self):
        truth = self.folder / "expected/cell-labels.tsv"
        passed = compare_labels(self.model, self.folder, truth)
        self.assertEqual(passed["status"], "pass")
        self.assertEqual(passed["metrics"]["coverage"], 1.0)
        with tempfile.TemporaryDirectory() as value:
            partial = Path(value) / "partial.tsv"
            partial.write_text("cell_id\tlabel\nCELL-001\tstate-1\n", encoding="utf-8")
            failed = compare_labels(self.model, self.folder, partial)
            self.assertEqual(failed["status"], "fail")
            self.assertLess(failed["metrics"]["coverage"], 0.02)

    def test_fixture_rebuild_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as value:
            target = Path(value)
            run = subprocess.run(
                [sys.executable, str(ROOT / "workflows/build_singlecell_fixture.py"), str(target)],
                capture_output=True,
                text=True,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            rebuilt = target / "datasets/singlecell/singlecell-001"
            for relative in ("manifest.json", "nano/matrix.mtx", "nano/barcodes.tsv", "nano/features.tsv", "expected/cell-labels.tsv", "expected/validation.json"):
                self.assertEqual((rebuilt / relative).read_bytes(), (self.folder / relative).read_bytes(), relative)

    def test_cli_qc_and_compare(self):
        qc = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "cell", "qc",
                str(self.folder / "nano/matrix.mtx"),
                str(self.folder / "nano/barcodes.tsv"),
                str(self.folder / "nano/features.tsv"),
                "--feature-columns", "3",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(qc.returncode, 0, qc.stderr)
        self.assertEqual(json.loads(qc.stdout)["cells"], 90)
        compared = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "--root", str(ROOT), "cell", "compare",
                "singlecell-001", str(self.folder / "expected/cell-labels.tsv"),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(compared.returncode, 0, compared.stderr)
        self.assertEqual(json.loads(compared.stdout)["status"], "pass")


if __name__ == "__main__":
    unittest.main()
