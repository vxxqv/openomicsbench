import json
import subprocess
import sys
import unittest
from pathlib import Path

from omicsbench.assays import PROFILES, build_plan, get_profile, list_profiles


class AssayProfileTests(unittest.TestCase):
    def test_profiles_cover_major_sequence_families(self):
        required = {
            "whole-genome", "exome", "targeted-panel", "bulk-rna-seq", "single-cell-rna-seq",
            "atac-seq", "chip-seq", "amplicon", "shotgun-metagenomics", "long-read-dna",
            "long-read-rna", "reference", "rna-fasta", "protein-fasta",
        }
        self.assertEqual(set(PROFILES), required)
        for profile in PROFILES.values():
            self.assertIn(profile["molecule"], {"dna", "rna", "protein"})
            self.assertTrue(profile["inputs"])
            self.assertTrue(profile["checks"])
            self.assertTrue(profile["local_scope"])
            self.assertTrue(profile["downstream"])

    def test_list_and_lookup(self):
        self.assertEqual(len(list_profiles()), 14)
        self.assertEqual(get_profile("protein-fasta")["molecule"], "protein")
        with self.assertRaisesRegex(ValueError, "unknown assay"):
            get_profile("unknown")

    def test_single_input_plan(self):
        report = build_plan("long-read-dna", [Path("reads.fastq.gz")], Path("results"))
        self.assertEqual([step["id"] for step in report["steps"]], ["validate-sequences", "qc-sequences"])
        self.assertIn("--molecule dna", report["steps"][0]["command"])

    def test_paired_short_read_plan(self):
        report = build_plan("whole-genome", [Path("r1.fq.gz"), Path("r2.fq.gz")], Path("results"), ["AGATCGGAAGAGC"])
        identifiers = [step["id"] for step in report["steps"]]
        self.assertEqual(identifiers, ["validate-r1", "qc-r1", "validate-r2", "qc-r2", "pair-check", "filter-r1", "filter-r2"])
        self.assertIn("--adapter AGATCGGAAGAGC", report["steps"][1]["command"])
        self.assertIn("results", report["steps"][-1]["command"])

    def test_plan_rejects_unsupported_input_count(self):
        with self.assertRaisesRegex(ValueError, "at least one"):
            build_plan("reference", [], Path("out"))
        with self.assertRaisesRegex(ValueError, "one paired"):
            build_plan("reference", [Path("a"), Path("b"), Path("c")], Path("out"))

    def test_assay_cli(self):
        commands = [
            ["assay", "list"],
            ["assay", "info", "bulk-rna-seq"],
            ["assay", "plan", "whole-genome", "r1.fq.gz", "r2.fq.gz", "--adapter", "AGATCGGAAGAGC"],
        ]
        for command in commands:
            with self.subTest(command=command):
                result = subprocess.run([sys.executable, "-m", "omicsbench", *command], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                json.loads(result.stdout)


if __name__ == "__main__":
    unittest.main()
