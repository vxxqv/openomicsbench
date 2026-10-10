import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.evaluate import evaluate_suite
from omicsbench.policy import evaluate_policy, read_policy


ROOT = Path(__file__).resolve().parents[1]


class GatePolicyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        self.results = self.directory / "results"
        self.results.mkdir()
        shutil.copyfile(
            ROOT / "datasets/sequences/sequence-006/expected/peptides.tsv",
            self.results / "sequence-006.peptides.tsv",
        )
        (self.results / "method.json").write_text(json.dumps({
            "schema_version": "1.0",
            "name": "test search",
            "version": "1.0",
            "runtime_seconds": 8.5,
            "peak_memory_mb": 96.0,
            "threads": 2,
        }), encoding="utf-8")
        report = evaluate_suite(ROOT, self.results, ["sequence-006"])
        self.report = self.directory / "evaluation.json"
        self.report.write_text(json.dumps(report), encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def write_policy(self, **changes) -> Path:
        policy = {
            "schema_version": "1.0",
            "name": "release gate",
            "required_benchmarks": ["sequence-006"],
            "minimum_pass_rate": 1.0,
            "minimum_assay_pass_rate": {"sequence_protein": 1.0},
            "require_method_receipt": True,
            "maximum_runtime_seconds": 10.0,
            "maximum_peak_memory_mb": 128.0,
            "metric_floors": [
                {"metric": "precision", "minimum": 1.0, "comparator": "peptide_identification"},
                {"metric": "protein_recall", "minimum": 1.0, "assay": "sequence_protein"},
            ],
        }
        policy.update(changes)
        path = self.directory / "policy.json"
        path.write_text(json.dumps(policy), encoding="utf-8")
        return path

    def test_complete_evaluation_passes_policy(self):
        gate = evaluate_policy(self.report, self.write_policy())
        self.assertEqual(gate["summary"]["status"], "pass")
        self.assertEqual(gate["source_run_fingerprint"], json.loads(self.report.read_text())["run_fingerprint"])
        self.assertEqual(gate["summary"]["failed"], 0)

    def test_resource_and_required_benchmark_fail_closed(self):
        gate = evaluate_policy(self.report, self.write_policy(
            required_benchmarks=["sequence-006", "rnaseq-002"],
            maximum_peak_memory_mb=64.0,
        ))
        self.assertEqual(gate["summary"]["status"], "fail")
        failures = {item["id"] for item in gate["results"] if item["status"] == "fail"}
        self.assertEqual(failures, {"benchmark:rnaseq-002", "peak-memory-mb"})

    def test_missing_resource_measurement_fails_a_declared_limit(self):
        report = json.loads(self.report.read_text())
        report["method"]["runtime_seconds"] = None
        self.report.write_text(json.dumps(report), encoding="utf-8")
        gate = evaluate_policy(self.report, self.write_policy())
        failure = next(item for item in gate["results"] if item["id"] == "runtime-seconds")
        self.assertEqual(failure["status"], "fail")
        self.assertIn("missing", failure["reason"])

    def test_invalid_or_duplicate_policy_scopes_are_rejected(self):
        path = self.write_policy(metric_floors=[
            {"metric": "precision", "minimum": 0.9},
            {"metric": "precision", "minimum": 0.8},
        ])
        with self.assertRaisesRegex(ValueError, "scopes must be unique"):
            read_policy(path)

    def test_cli_returns_two_and_writes_report_for_failed_gate(self):
        policy = self.write_policy(maximum_runtime_seconds=1.0)
        output = self.directory / "gate.json"
        completed = subprocess.run(
            [sys.executable, "-m", "omicsbench", "suite", "gate", str(self.report), str(policy), "--json", str(output)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 2, completed.stderr)
        self.assertEqual(json.loads(output.read_text())["summary"]["status"], "fail")


if __name__ == "__main__":
    unittest.main()
