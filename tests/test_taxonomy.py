import json
import tempfile
import unittest
from pathlib import Path

from omicsbench.evaluate import evaluate_suite
from omicsbench.evidence import create_evidence_crate, verify_evidence_crate
from omicsbench.registry import default_root, lookup
from omicsbench.results import create_result_manifest
from omicsbench.taxonomy import classify_exact, compare_taxonomy, read_assignments
from omicsbench.validate import validate


class TaxonomyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = default_root()
        cls.model, cls.folder = lookup(cls.root, "sequence-007")
        cls.truth = cls.folder / "expected/read-taxonomy.tsv"

    def test_fixture_and_exact_classifier(self):
        report = validate(self.model, self.folder)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["taxonomy_truth"]["counts"], {"0": 12, "1": 12, "2": 12, "3": 12})
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "calls.tsv"
            classified = classify_exact(self.folder / "nano/reference.fasta", self.folder / "nano/reads_R1.fastq", self.folder / "nano/reads_R2.fastq", output)
            self.assertEqual((classified["assigned"], classified["unclassified"]), (36, 12))
            self.assertEqual(output.read_bytes(), self.truth.read_bytes())
            with self.assertRaisesRegex(ValueError, "already exists"):
                classify_exact(self.folder / "nano/reference.fasta", self.folder / "nano/reads_R1.fastq", self.folder / "nano/reads_R2.fastq", output)

    def test_kraken_and_failure_cases(self):
        truth = read_assignments(self.truth)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "calls.kraken"
            output.write_text("".join(f"{'C' if code else 'U'}\t{name}\t{code}\t50|50\t0:1\n" for name, code in truth.items()), encoding="utf-8")
            self.assertEqual(compare_taxonomy(self.model, self.folder, output)["status"], "pass")
            lines = output.read_text(encoding="utf-8").splitlines()
            lines[0] = lines[0].replace("\t1\t", "\t2\t")
            output.write_text("\n".join(lines) + "\n", encoding="utf-8")
            self.assertLess(compare_taxonomy(self.model, self.folder, output)["metrics"]["accuracy"], 1.0)
            output.write_text("\n".join(lines[1:]) + "\n", encoding="utf-8")
            failed = compare_taxonomy(self.model, self.folder, output)
            self.assertEqual(failed["status"], "fail")
            self.assertLess(failed["metrics"]["coverage"], 1)
            output.write_text(lines[1] + "\n" + lines[1] + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate read ID"):
                read_assignments(output)

    def test_unified_evaluation_and_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            result = directory / "sequence-007.taxonomy.tsv"
            result.write_bytes(self.truth.read_bytes())
            create_result_manifest(self.root, directory, [f"sequence-007={result.name}"], [])
            report = evaluate_suite(self.root, directory, ["sequence-007"])
            self.assertEqual(report["summary"]["passed"], 1)
            self.assertEqual(report["results"][0]["details"]["metrics"]["macro_f1"], 1.0)
            archive = directory / "evidence.zip"
            create_evidence_crate(report, archive)
            self.assertEqual(verify_evidence_crate(archive)["status"], "pass")


if __name__ == "__main__":
    unittest.main()
