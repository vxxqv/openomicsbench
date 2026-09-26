import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from omicsbench.suite import compare_suite, junit_report, markdown_report, validate_suite, write_reports


class SuiteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_validation_suite_reports_pass_and_skip(self):
        report = validate_suite(self.collection, dataset_ids=["sequence-001", "rnaseq-001"])
        self.assertEqual(report["summary"], {"status": "pass", "total": 2, "passed": 1, "failed": 0, "skipped": 1})
        self.assertEqual([item["status"] for item in report["results"]], ["skipped", "pass"])

    def test_validation_suite_rejects_unknown_ids(self):
        with self.assertRaisesRegex(ValueError, "unknown dataset IDs"):
            validate_suite(self.collection, dataset_ids=["missing"])

    def test_compare_suite_passes_reference_effects(self):
        source = self.collection / "datasets" / "rnaseq" / "rnaseq-002" / "expected" / "reference-effects.tsv.gz"
        target = self.root / "rnaseq-002.tsv.gz"
        target.write_bytes(source.read_bytes())
        report = compare_suite(self.collection, self.root, ["rnaseq-002"])
        self.assertEqual(report["summary"]["status"], "pass")
        self.assertEqual(report["results"][0]["details"]["metrics"]["spearman_logfc"], 1.0)

    def test_compare_suite_marks_missing_result_as_failure(self):
        report = compare_suite(self.collection, self.root, ["rnaseq-002"])
        self.assertEqual(report["summary"]["failed"], 1)
        self.assertIn("missing rnaseq-002", report["results"][0]["reason"])

    def test_reports_are_machine_and_human_readable(self):
        report = validate_suite(self.collection, dataset_ids=["sequence-001"])
        markdown = markdown_report(report)
        self.assertIn("`sequence-001`", markdown)
        xml = junit_report(report)
        root = ET.fromstring(xml)
        self.assertEqual(root.attrib["tests"], "1")
        json_path = self.root / "reports" / "suite.json"
        markdown_path = self.root / "reports" / "suite.md"
        junit_path = self.root / "reports" / "suite.xml"
        paths = write_reports(report, json_path, markdown_path, junit_path)
        self.assertEqual(json.loads(json_path.read_text())["summary"]["status"], "pass")
        self.assertTrue(markdown_path.is_file())
        self.assertTrue(junit_path.is_file())
        self.assertEqual(paths["junit"], str(junit_path))
        with self.assertRaisesRegex(ValueError, "--force"):
            write_reports(report, json_path=json_path)

    def test_cli_suite_exit_codes_and_reports(self):
        report_path = self.root / "suite.xml"
        passed = subprocess.run(
            [sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "validate", "--id", "sequence-001", "--junit", str(report_path)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(passed.returncode, 0, passed.stderr)
        self.assertTrue(report_path.is_file())
        missing = subprocess.run(
            [sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "compare", str(self.root), "--id", "rnaseq-002"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(missing.returncode, 2, missing.stderr)
        self.assertEqual(json.loads(missing.stdout)["summary"]["failed"], 1)


if __name__ == "__main__":
    unittest.main()
