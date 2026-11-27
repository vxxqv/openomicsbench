import tempfile
import unittest
from pathlib import Path

from omicsbench import Collection


class RNASeqAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.collection = Collection()

    def test_count_qc_reports_sample_metrics_and_empty_library(self):
        with tempfile.TemporaryDirectory() as temporary:
            counts = Path(temporary) / "counts.tsv"
            counts.write_text("gene_id\ta\tb\nMT-X\t2\t0\nG2\t3\t0\nG3\t0\t0\n", encoding="utf-8")
            report = self.collection.qc_rnaseq(counts)
            self.assertEqual(report["summary"]["status"], "fail")
            self.assertEqual(report["sample_qc"][0]["library_size"], 5)
            self.assertEqual(report["sample_qc"][0]["detected_features"], 2)
            self.assertAlmostEqual(report["sample_qc"][0]["mitochondrial_fraction"], 0.4)
            self.assertIsNone(report["sample_qc"][1]["largest_feature_fraction"])
            self.assertIsNone(report["sample_qc"][1]["mitochondrial_fraction"])

    def test_de_assessment_counts_adjusted_calls_and_reference_direction(self):
        benchmark = self.collection.get("rnaseq-005")
        reference = benchmark.file(role="effects", tier="expected")
        import gzip
        with gzip.open(reference, "rt", encoding="utf-8") as stream:
            next(stream)
            first = next(stream).strip().split("\t")
            second = next(stream).strip().split("\t")
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.csv"
            result.write_text(
                "gene_id,log2FoldChange,padj\n"
                f"{first[0]},{first[1]},0.01\n"
                f"{second[0]},{-float(second[1])},0.5\n"
                "UNLISTED,NA,NA\n", encoding="utf-8")
            report = benchmark.assess_de(result)
            self.assertEqual(report["rows"], 3)
            self.assertEqual(report["tested_rows"], 2)
            self.assertEqual(report["significant_rows"], 1)
            self.assertEqual(report["missing_adjusted_p_values"], 1)
            self.assertEqual(report["shared_finite_effects"], 2)
            self.assertEqual(report["effect_direction_agreement"], 0.5)
            self.assertIsNone(report["effect_rank_spearman"])
            self.assertIsNone(report["independent_truth"])
            truth = Path(temporary) / "truth.csv"
            truth.write_text(
                "gene_id,is_de\n"
                f"{first[0]},1\n{second[0]},0\n", encoding="utf-8")
            scored = benchmark.assess_de(result, truth=truth)
            self.assertEqual(scored["independent_truth"]["true_positive"], 1)
            self.assertEqual(scored["independent_truth"]["true_negative"], 1)
            self.assertEqual(scored["independent_truth"]["f1"], 1.0)
            truth.write_text("gene_id,is_de\nUNLISTED,1\n" + f"{second[0]},0\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "absent or have no adjusted p-value"):
                benchmark.assess_de(result, truth=truth)
            with self.assertRaisesRegex(ValueError, "between zero and one"):
                benchmark.assess_de(result, alpha=1)
            result.write_text(result.read_text(encoding="utf-8").replace(",0.5", ",1.2"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "adjusted p-value must be"):
                benchmark.assess_de(result)

    def test_de_rank_agreement_on_three_genes(self):
        import gzip
        benchmark = self.collection.get("rnaseq-005")
        with gzip.open(benchmark.file(role="effects", tier="expected"), "rt", encoding="utf-8") as stream:
            next(stream)
            rows = [next(stream).strip().split("\t") for _ in range(3)]
        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "results.tsv"
            result.write_text(
                "gene_id\tlog2FoldChange\tpadj\n" +
                "".join(f"{gene}\t{effect}\t0.5\n" for gene, effect, *_ in rows),
                encoding="utf-8",
            )
            self.assertAlmostEqual(benchmark.assess_de(result)["effect_rank_spearman"], 1.0)


if __name__ == "__main__":
    unittest.main()
