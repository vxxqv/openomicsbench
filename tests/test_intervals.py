import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.intervals import Interval, compare_intervals, interval_metrics, read_chrom_sizes, read_intervals, summarize_intervals
from omicsbench.registry import lookup


ROOT = Path(__file__).resolve().parents[1]
CLI = [sys.executable, "-c", "from omicsbench.cli import main; main()"]


class IntervalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_summary_and_strict_reference_checks(self):
        truth = ROOT / "datasets/intervals/interval-001/expected/truth-peaks.narrowPeak"
        sizes = ROOT / "datasets/intervals/interval-001/nano/chrom.sizes"
        report = summarize_intervals(truth, sizes)
        self.assertEqual(report["intervals"], 12)
        self.assertEqual(report["format"], "narrowPeak")
        invalid = self.directory / "invalid.bed"
        invalid.write_text("chrSynthetic1\t0\t12001\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "outside the declared reference"):
            read_intervals(invalid, read_chrom_sizes(sizes))

    def test_metrics_use_reciprocal_and_union_overlap(self):
        truth = [Interval("chr1", 0, 100), Interval("chr1", 200, 300)]
        calls = [Interval("chr1", 10, 90), Interval("chr1", 400, 500)]
        metrics = interval_metrics(truth, calls, 0.5)
        self.assertEqual(metrics["interval_precision"], 0.5)
        self.assertEqual(metrics["interval_recall"], 0.5)
        self.assertEqual(metrics["intersection_bp"], 80)
        self.assertAlmostEqual(metrics["base_recall"], 0.4)

    def test_truth_passes_and_incomplete_calls_fail(self):
        model, folder = lookup(ROOT, "interval-001")
        truth = folder / "expected/truth-peaks.narrowPeak"
        self.assertEqual(compare_intervals(model, folder, truth)["status"], "pass")
        incomplete = self.directory / "calls.bed"
        incomplete.write_text("chrSynthetic1\t500\t700\n", encoding="utf-8")
        failed = compare_intervals(model, folder, incomplete)
        self.assertEqual(failed["status"], "fail")
        self.assertTrue(any("interval_recall" in reason for reason in failed["reasons"]))

    def test_cli_qc_and_compare(self):
        truth = ROOT / "datasets/intervals/interval-001/expected/truth-peaks.narrowPeak"
        sizes = ROOT / "datasets/intervals/interval-001/nano/chrom.sizes"
        qc = subprocess.run(CLI + ["interval", "qc", str(truth), "--chrom-sizes", str(sizes)], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(qc.returncode, 0, qc.stderr)
        self.assertEqual(json.loads(qc.stdout)["intervals"], 12)
        scored = subprocess.run(CLI + ["interval", "compare", "interval-001", str(truth)], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(scored.returncode, 0, scored.stderr)
        self.assertEqual(json.loads(scored.stdout)["status"], "pass")

    def test_builder_is_reproducible(self):
        result = subprocess.run([sys.executable, str(ROOT / "workflows/build_interval_fixture.py"), str(self.directory)], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        rebuilt = self.directory / "datasets/intervals/interval-001"
        original = ROOT / "datasets/intervals/interval-001"
        original_files = {path.relative_to(original).as_posix(): path.read_bytes() for path in original.rglob("*") if path.is_file()}
        rebuilt_files = {path.relative_to(rebuilt).as_posix(): path.read_bytes() for path in rebuilt.rglob("*") if path.is_file()}
        self.assertEqual(original_files, rebuilt_files)


if __name__ == "__main__":
    unittest.main()
