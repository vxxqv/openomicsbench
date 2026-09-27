import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.regression import compare_reports, read_report


def comparison(value=1.0, status="pass"):
    return {
        "operation": "compare",
        "summary": {"status": status, "total": 1, "passed": int(status == "pass"), "failed": int(status == "fail"), "skipped": 0},
        "results": [{
            "id": "rnaseq-002",
            "status": status,
            "details": {
                "metrics": {"spearman_logfc": value, "top_k_jaccard": value, "sign_concordance": value},
                "genes": {"coverage": value},
            },
        }],
    }


class RegressionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.baseline = self.root / "baseline.json"
        self.candidate = self.root / "candidate.json"

    def tearDown(self):
        self.temp.cleanup()

    def write(self, path, payload):
        path.write_text(json.dumps(payload), encoding="utf-8")

    def test_detects_metric_drop_and_status_failure(self):
        self.write(self.baseline, comparison())
        self.write(self.candidate, comparison(0.9, "fail"))
        report = compare_reports(self.baseline, self.candidate)
        self.assertEqual(report["summary"]["status"], "fail")
        self.assertIn("status changed", report["results"][0]["reason"])
        self.assertIn("spearman_logfc decreased", report["results"][0]["reason"])

    def test_tolerance_allows_small_metric_drift(self):
        self.write(self.baseline, comparison())
        self.write(self.candidate, comparison(0.999))
        report = compare_reports(self.baseline, self.candidate, absolute_tolerance=0.01)
        self.assertEqual(report["summary"]["status"], "pass")

    def test_missing_case_fails_closed(self):
        self.write(self.baseline, comparison())
        empty = {"operation": "compare", "summary": {}, "results": []}
        self.write(self.candidate, empty)
        self.assertEqual(compare_reports(self.baseline, self.candidate)["summary"]["failed"], 1)
        self.assertEqual(compare_reports(self.baseline, self.candidate, allow_missing=True)["summary"]["failed"], 0)

    def test_invalid_report_is_rejected(self):
        self.candidate.write_text('{"operation":"unknown","results":[]}', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            read_report(self.candidate)

    def test_cli_exit_status_and_junit(self):
        self.write(self.baseline, comparison())
        self.write(self.candidate, comparison(0.8))
        junit = self.root / "regression.xml"
        completed = subprocess.run(
            [sys.executable, "-m", "omicsbench", "suite", "regress", str(self.baseline), str(self.candidate), "--junit", str(junit)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 2, completed.stderr)
        self.assertTrue(junit.is_file())


if __name__ == "__main__":
    unittest.main()
