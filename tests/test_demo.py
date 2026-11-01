import json
import tempfile
import unittest
from pathlib import Path

from omicsbench.demo import create_demo
from omicsbench.evidence import verify_evidence_crate
from omicsbench.registry import default_root


class DemoTests(unittest.TestCase):
    def test_reference_walkthrough_and_existing_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "walkthrough"
            result = create_demo(default_root(), destination)
            self.assertEqual(result["summary"]["passed"], 6)
            self.assertEqual(result["summary"]["failed"], 0)
            report = json.loads((destination / "report.json").read_text(encoding="utf-8"))
            self.assertEqual(len(report["results"]), 6)
            self.assertIn("does not measure an external method", report["example_scope"])
            self.assertIn("does not measure an external method", (destination / "report.html").read_text(encoding="utf-8"))
            self.assertTrue((destination / "report.html").is_file())
            self.assertEqual(verify_evidence_crate(destination / "evidence.zip")["status"], "pass")
            with self.assertRaisesRegex(ValueError, "already exists"):
                create_demo(default_root(), destination)
            self.assertTrue((destination / "report.json").is_file())


if __name__ == "__main__":
    unittest.main()
