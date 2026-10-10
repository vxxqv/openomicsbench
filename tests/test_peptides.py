import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.peptides import compare_peptides, read_peptides
from omicsbench.registry import lookup


ROOT = Path(__file__).resolve().parents[1]


class PeptideComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model, cls.folder = lookup(ROOT, "sequence-006")
        cls.truth = cls.folder / "expected/peptides.tsv"

    def result_table(self, directory: Path, rows: list[dict], include_q_value: bool = False) -> Path:
        path = directory / "sequence-006.peptides.tsv"
        fields = ["protein_id", "peptide"] + (["q_value"] if include_q_value else [])
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def truth_rows(self) -> list[dict]:
        with self.truth.open(encoding="utf-8", newline="") as stream:
            return [{"protein_id": row["protein_id"], "peptide": row["peptide"]} for row in csv.DictReader(stream, delimiter="\t")]

    def test_exact_truth_passes_with_complete_metrics(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = self.result_table(Path(temporary), self.truth_rows())
            report = compare_peptides(self.model, self.folder, result)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["metrics"]["f1"], 1.0)
        self.assertEqual(report["metrics"]["protein_recall"], 1.0)
        self.assertEqual(report["metrics"]["residue_coverage"], 1.0)

    def test_false_identification_and_missing_truth_fail(self):
        rows = self.truth_rows()[1:]
        rows.append({"protein_id": "proteome_alpha", "peptide": "AAAAAAAA"})
        with tempfile.TemporaryDirectory() as temporary:
            report = compare_peptides(self.model, self.folder, self.result_table(Path(temporary), rows))
        self.assertEqual(report["status"], "fail")
        self.assertEqual(report["counts"]["false_positive"], 1)
        self.assertEqual(report["counts"]["false_negative"], 1)

    def test_q_value_filter_is_applied_before_scoring(self):
        rows = [{**row, "q_value": "0.001"} for row in self.truth_rows()]
        rows.append({"protein_id": "proteome_alpha", "peptide": "AAAAAAAA", "q_value": "0.2"})
        with tempfile.TemporaryDirectory() as temporary:
            result = self.result_table(Path(temporary), rows, include_q_value=True)
            report = compare_peptides(self.model, self.folder, result)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["input_summary"]["filtered_identifications"], 1)

    def test_malformed_and_duplicate_identifications_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            duplicate = self.truth_rows()[0]
            path = self.result_table(directory, [duplicate, duplicate])
            with self.assertRaisesRegex(ValueError, "duplicate identification"):
                read_peptides(path, 0.01)
            self.result_table(directory, [{"protein_id": "proteome_alpha", "peptide": "PEP*IDE"}])
            with self.assertRaisesRegex(ValueError, "canonical"):
                read_peptides(path, 0.01)

    def test_cli_uses_exit_two_for_scientific_failure(self):
        rows = self.truth_rows()[1:]
        with tempfile.TemporaryDirectory() as temporary:
            path = self.result_table(Path(temporary), rows)
            completed = subprocess.run(
                [sys.executable, "-m", "omicsbench", "--root", str(ROOT), "peptide", "compare", "sequence-006", str(path)],
                capture_output=True,
                text=True,
            )
        self.assertEqual(completed.returncode, 2, completed.stderr)


if __name__ == "__main__":
    unittest.main()
