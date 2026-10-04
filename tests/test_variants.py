import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.registry import lookup
from omicsbench.variants import compare_variants, read_vcf


class VariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]
        cls.model, cls.folder = lookup(cls.collection, "sequence-005")

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write_vcf(self, body: str, name: str = "calls.vcf") -> Path:
        path = self.root / name
        content = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + body
        if name.endswith(".gz"):
            with gzip.open(path, "wt", encoding="utf-8", newline="") as stream:
                stream.write(content)
        else:
            path.write_text(content, encoding="utf-8", newline="\n")
        return path

    def truth_body(self) -> str:
        return "synthetic_reference\t150\t.\tC\tA\t60\tPASS\t.\nsynthetic_reference\t500\t.\tA\tC\t60\tPASS\t.\nsynthetic_reference\t850\t.\tC\tA\t60\tPASS\t.\n"

    def test_exact_truth_passes(self):
        report = compare_variants(self.model, self.folder, self.write_vcf(self.truth_body()))
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["counts"]["true_positive"], 3)
        self.assertEqual(report["metrics"], {"precision": 1.0, "recall": 1.0, "f1": 1.0})

    def test_false_positive_and_negative_fail(self):
        body = "synthetic_reference\t150\t.\tC\tA\t60\tPASS\t.\nsynthetic_reference\t151\t.\tT\tG\t60\tPASS\t.\n"
        report = compare_variants(self.model, self.folder, self.write_vcf(body))
        self.assertEqual(report["status"], "fail")
        self.assertEqual(report["counts"]["false_positive"], 1)
        self.assertEqual(report["counts"]["false_negative"], 2)

    def test_multiallelic_and_filtered_records_are_accounted_for(self):
        body = "synthetic_reference\t150\t.\tC\tA,G\t60\tPASS\t.\nsynthetic_reference\t500\t.\tA\tC\t60\tLowQual\t.\n"
        alleles, summary = read_vcf(self.write_vcf(body, "calls.vcf.gz"))
        self.assertEqual(len(alleles), 2)
        self.assertEqual(summary["filtered_alleles"], 1)

    def test_reference_mismatch_is_invalid_input(self):
        body = "synthetic_reference\t150\t.\tG\tA\t60\tPASS\t.\n"
        with self.assertRaisesRegex(ValueError, "disagrees with the reference"):
            compare_variants(self.model, self.folder, self.write_vcf(body))

    def test_rejects_non_snv_and_bad_header(self):
        with self.assertRaisesRegex(ValueError, "single-nucleotide"):
            read_vcf(self.write_vcf("synthetic_reference\t150\t.\tC\tCA\t60\tPASS\t.\n"))
        bad = self.root / "bad.vcf"
        bad.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "fileformat"):
            read_vcf(bad)

    def test_cli_exit_status(self):
        path = self.write_vcf(self.truth_body())
        completed = subprocess.run(
            [sys.executable, "-m", "omicsbench", "--root", str(self.collection), "variant", "compare", "sequence-005", str(path)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)["status"], "pass")


if __name__ == "__main__":
    unittest.main()
