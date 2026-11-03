import json
import shutil
import tempfile
import unittest
from pathlib import Path

from omicsbench.demo import create_demo
from omicsbench.evidence import verify_evidence_crate
from omicsbench.registry import default_root
from omicsbench.results import create_result_manifest
from omicsbench.workspace import grade_workspace, prepare_workspace, verify_workspace


IDS = ("rnaseq-002", "singlecell-001", "interval-001", "sequence-005", "sequence-006", "sequence-007")


class WorkspaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = default_root()

    def test_prepare_and_grade_six_assays(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            workspace = base / "work"
            prepared = prepare_workspace(self.root, workspace, list(IDS))
            self.assertEqual(prepared["summary"]["passed"], 6)
            manifest = json.loads((workspace / "workspace.json").read_text(encoding="utf-8"))
            self.assertEqual([item["id"] for item in manifest["cases"]], list(IDS))
            self.assertFalse(any("/expected/" in path.as_posix() for path in (workspace / "inputs").rglob("*")))
            self.assertEqual(verify_workspace(self.root, workspace)["summary"]["status"], "pass")
            create_demo(self.root, base / "demo")
            for source in (base / "demo/reference-results").iterdir():
                if source.name != "results.json":
                    shutil.copyfile(source, workspace / "results" / source.name)
            graded = grade_workspace(self.root, workspace)
            self.assertEqual(graded["summary"]["passed"], 6)
            self.assertTrue((workspace / "report.html").is_file())
            self.assertEqual(verify_evidence_crate(workspace / "evidence.zip")["status"], "pass")
            with self.assertRaisesRegex(ValueError, "already exists"):
                grade_workspace(self.root, workspace)
            self.assertEqual(grade_workspace(self.root, workspace, force=True)["summary"]["passed"], 6)

    def test_changed_input_and_extra_file_fail_verification(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["sequence-007"])
            path = workspace / "inputs/sequence-007/nano/reference.fasta"
            original = path.read_bytes()
            path.write_bytes(original + b"A")
            with self.assertRaisesRegex(ValueError, "checksum"):
                verify_workspace(self.root, workspace)
            path.write_bytes(original)
            (workspace / "inputs/sequence-007/extra.txt").write_text("unexpected", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "file set differs"):
                verify_workspace(self.root, workspace)

    def test_missing_result_produces_a_failure_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["sequence-007"])
            graded = grade_workspace(self.root, workspace)
            self.assertEqual(graded["summary"]["failed"], 1)
            self.assertEqual(verify_evidence_crate(workspace / "evidence.zip")["status"], "pass")

    def test_native_result_path_uses_the_result_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["sequence-007"])
            result = workspace / "results/taxonomy/calls.tsv"
            result.parent.mkdir()
            shutil.copyfile(self.root / "datasets/sequences/sequence-007/expected/read-taxonomy.tsv", result)
            create_result_manifest(self.root, workspace / "results", ["sequence-007=taxonomy/calls.tsv"], [])
            self.assertEqual(grade_workspace(self.root, workspace)["summary"]["passed"], 1)

    def test_changed_case_provenance_fails_verification(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["sequence-007"])
            path = workspace / "workspace.json"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["cases"][0]["source_citation"] = "changed"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "provenance"):
                verify_workspace(self.root, workspace)

    def test_invalid_selection_does_not_create_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "work"
            with self.assertRaisesRegex(ValueError, "no validated result comparator"):
                prepare_workspace(self.root, output, ["sequence-001"])
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
