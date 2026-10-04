import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

from omicsbench.evidence import create_evidence_crate, verify_evidence_crate
from omicsbench.evaluate import evaluate_suite


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]
        cls.effects = cls.collection / "datasets/rnaseq/rnaseq-002/expected/reference-effects.tsv.gz"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.results = self.root / "results"
        self.results.mkdir()
        shutil.copyfile(self.effects, self.results / "rnaseq-002.tsv.gz")

    def tearDown(self):
        self.temp.cleanup()

    def test_crate_is_deterministic_and_verifiable(self):
        report = evaluate_suite(self.collection, self.results, ["rnaseq-002"])
        first = self.root / "first.zip"
        second = self.root / "second.zip"
        created = create_evidence_crate(report, first)
        create_evidence_crate(report, second)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        verified = verify_evidence_crate(first)
        self.assertEqual(verified["status"], "pass")
        self.assertEqual(created["sha256"], verified["sha256"])
        with zipfile.ZipFile(first) as archive:
            self.assertIn("ro-crate-metadata.json", archive.namelist())
            self.assertIn("inputs/rnaseq-002.tsv.gz", archive.namelist())
            self.assertNotIn(str(self.root), archive.read("report.json").decode("utf-8"))

    def test_tampered_crate_fails_verification(self):
        report = evaluate_suite(self.collection, self.results, ["rnaseq-002"])
        original = self.root / "original.zip"
        changed = self.root / "changed.zip"
        create_evidence_crate(report, original)
        with zipfile.ZipFile(original) as source, zipfile.ZipFile(changed, "w") as target:
            for name in source.namelist():
                data = source.read(name)
                if name == "report.csv":
                    data += b"changed"
                target.writestr(name, data)
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            verify_evidence_crate(changed)

    def test_cli_creates_and_verifies_evidence(self):
        archive = self.root / "evidence.zip"
        created = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "evaluate", str(self.results),
                "--id", "rnaseq-002", "--evidence", str(archive),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(created.returncode, 0, created.stderr)
        self.assertTrue(archive.is_file())
        verified = subprocess.run(
            [sys.executable, "-m", "omicsbench", "evidence", "verify", str(archive)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(verified.returncode, 0, verified.stderr)
        self.assertIn('"status": "pass"', verified.stdout)


if __name__ == "__main__":
    unittest.main()
