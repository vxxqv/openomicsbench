import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.proteins import (
    cleavage_sites,
    digest_sequence,
    isoelectric_point,
    net_charge,
    protein_record_stats,
)


class ProteinTests(unittest.TestCase):
    def test_known_small_peptide_mass(self):
        report = protein_record_stats("x", "ACD")
        self.assertAlmostEqual(report["molecular_weight_da"], 307.32148, places=5)
        self.assertEqual(report["composition"], {"A": 1, "C": 1, "D": 1})

    def test_ambiguous_mass_is_not_invented(self):
        report = protein_record_stats("x", "ABX")
        self.assertIsNone(report["molecular_weight_da"])
        self.assertEqual(report["unknown_mass_residues"], ["B", "X"])

    def test_charge_changes_with_ph(self):
        sequence = "DEKKR"
        self.assertGreater(net_charge(sequence, 2), net_charge(sequence, 12))
        self.assertGreaterEqual(isoelectric_point(sequence), 0)
        self.assertLessEqual(isoelectric_point(sequence), 14)

    def test_hydropathy_aromaticity_and_extinction(self):
        report = protein_record_stats("x", "FWYI")
        self.assertEqual(report["aromatic_fraction"], 0.75)
        self.assertEqual(report["extinction_coefficient_reduced_280nm"], 6990)
        self.assertGreater(report["mean_hydropathy"], 0)

    def test_trypsin_proline_exception(self):
        self.assertEqual(cleavage_sites("AKRPQK", "trypsin"), [0, 2, 6])

    def test_digest_with_missed_cleavage(self):
        peptides = digest_sequence("p", "AKRTAK", "trypsin", missed_cleavages=1)
        values = {(item["sequence"], item["missed_cleavages"]) for item in peptides}
        self.assertIn(("AK", 0), values)
        self.assertIn(("AKR", 1), values)
        self.assertIn(("RTAK", 1), values)

    def test_digest_length_filter(self):
        peptides = digest_sequence("p", "AKRTAK", "trypsin", minimum_length=2, maximum_length=3)
        self.assertEqual([item["sequence"] for item in peptides], ["AK", "TAK"])

    def test_cli_protein_commands(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "proteins.fa"
            path.write_text(">p\nAKRTAK\n", encoding="utf-8")
            for command in (["protein-stats", str(path), "--ph", "7.4"], ["digest", str(path), "--enzyme", "trypsin", "--missed-cleavages", "1"]):
                with self.subTest(command=command):
                    result = subprocess.run([sys.executable, "-m", "omicsbench", "seq", *command], capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    json.loads(result.stdout)


if __name__ == "__main__":
    unittest.main()
