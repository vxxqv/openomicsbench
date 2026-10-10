import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.evaluate import evaluate_suite
from omicsbench.trends import parse_run, trend_reports


ROOT = Path(__file__).resolve().parents[1]


class TrendReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        truth = ROOT / "datasets/sequences/sequence-006/expected/peptides.tsv"
        self.reports = []
        for name in ("baseline", "candidate"):
            results = self.directory / name
            results.mkdir()
            shutil.copyfile(truth, results / "sequence-006.peptides.tsv")
            (results / "method.json").write_text(json.dumps({
                "schema_version": "1.0", "name": "search", "version": name,
                "runtime_seconds": 10.0 if name == "baseline" else 8.0,
                "peak_memory_mb": 128.0 if name == "baseline" else 96.0,
            }), encoding="utf-8")
            report = evaluate_suite(ROOT, results, ["sequence-006"])
            path = self.directory / f"{name}.json"
            path.write_text(json.dumps(report), encoding="utf-8")
            self.reports.append(path)

    def tearDown(self):
        self.temp.cleanup()

    def test_history_preserves_runs_metrics_resources_and_fingerprints(self):
        report = trend_reports([("baseline", self.reports[0]), ("candidate", self.reports[1])])
        self.assertEqual(report["summary"]["status"], "pass")
        self.assertEqual(report["summary"]["runs"], 2)
        self.assertEqual([run["name"] for run in report["runs"]], ["baseline", "candidate"])
        self.assertEqual(report["runs"][1]["method"]["runtime_seconds"], 8.0)
        self.assertEqual(len(report["results"][0]["history"]), 2)
        self.assertEqual(report["results"][0]["metric_deltas"]["f1"], 0.0)
        self.assertRegex(report["run_fingerprint"], r"^sha256:[0-9a-f]{64}$")

    def test_latest_missing_benchmark_fails_trend(self):
        candidate = json.loads(self.reports[1].read_text())
        candidate["results"] = []
        candidate["summary"] = {"status": "pass", "total": 0, "passed": 0, "failed": 0, "skipped": 0, "assays": {}}
        self.reports[1].write_text(json.dumps(candidate), encoding="utf-8")
        report = trend_reports([("baseline", self.reports[0]), ("candidate", self.reports[1])])
        self.assertEqual(report["summary"]["status"], "fail")
        self.assertEqual(report["results"][0]["history"][-1]["status"], "missing")

    def test_parser_and_duplicate_names_fail_closed(self):
        self.assertEqual(parse_run(f"v1={self.reports[0]}")[0], "v1")
        with self.assertRaisesRegex(ValueError, "NAME=PATH"):
            parse_run("bad")
        with self.assertRaisesRegex(ValueError, "unique"):
            trend_reports([("same", self.reports[0]), ("same", self.reports[1])])

    def test_cli_writes_json_csv_and_offline_html(self):
        json_path = self.directory / "trend.json"
        csv_path = self.directory / "trend.csv"
        html_path = self.directory / "trend.html"
        completed = subprocess.run([
            sys.executable, "-m", "omicsbench", "suite", "trend",
            "--run", f"baseline={self.reports[0]}", "--run", f"candidate={self.reports[1]}",
            "--json", str(json_path), "--csv", str(csv_path), "--html", str(html_path),
        ], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(json_path.read_text())["summary"]["runs"], 2)
        with csv_path.open(encoding="utf-8") as stream:
            self.assertIn("delta_f1", next(csv.reader(stream)))
        html = html_path.read_text(encoding="utf-8")
        self.assertIn("Run history", html)
        self.assertNotIn("<script", html)


if __name__ == "__main__":
    unittest.main()
