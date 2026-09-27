import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

from omicsbench.bundle import create_bundle, verify_bundle


class BundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_bundle_is_deterministic_and_verified(self):
        first = self.root / "first.zip"
        second = self.root / "second.zip"
        one = create_bundle(self.collection, first, ["sequence-001"])
        two = create_bundle(self.collection, second, ["sequence-001"])
        self.assertEqual(one["sha256"], two["sha256"])
        self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(one["datasets"], ["sequence-001"])
        self.assertEqual(verify_bundle(first)["status"], "pass")

    def test_bundle_rejects_changed_content(self):
        source = self.root / "source.zip"
        broken = self.root / "broken.zip"
        create_bundle(self.collection, source, ["sequence-001"])
        with zipfile.ZipFile(source) as original, zipfile.ZipFile(broken, "w") as target:
            for info in original.infolist():
                data = original.read(info.filename)
                if info.filename.endswith("sequences.fasta"):
                    data += b"A"
                target.writestr(info, data)
        with self.assertRaisesRegex(ValueError, "size or SHA-256 differs"):
            verify_bundle(broken)

    def test_bundle_rejects_unknown_dataset(self):
        with self.assertRaisesRegex(ValueError, "unknown dataset"):
            create_bundle(self.collection, self.root / "x.zip", ["missing"])

    def test_bundle_cli_create_and_verify(self):
        archive = self.root / "portable.zip"
        created = subprocess.run(
            [sys.executable, "-m", "omicsbench", "--root", str(self.collection), "bundle", "create", str(archive), "--id", "sequence-001"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(created.returncode, 0, created.stderr)
        self.assertEqual(json.loads(created.stdout)["dataset_count"], 1)
        verified = subprocess.run(
            [sys.executable, "-m", "omicsbench", "bundle", "verify", str(archive)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(verified.returncode, 0, verified.stderr)
        self.assertEqual(json.loads(verified.stdout)["status"], "pass")


if __name__ == "__main__":
    unittest.main()
