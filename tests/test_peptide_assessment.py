import tempfile
import unittest
from pathlib import Path

from omicsbench.peptide_assessment import assess_peptides


class PeptideAssessmentTests(unittest.TestCase):
    def test_identification_thresholds_mapping_and_quantity(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            truth, calls, mapping = (base / name for name in ('truth.tsv', 'calls.tsv', 'map.tsv'))
            truth.write_text('protein_id\tpeptide\tq_value\tquantity\n'
                             'P1\tACDEFG\t0.001\t100\nP2\tLMNPQR\t0.001\t200\n', encoding='utf-8')
            calls.write_text('protein_id\tpeptide\tq_value\tquantity\n'
                             'ALIAS\tACDEFG\t0.001\t100\nP2\tLMNPQR\t0.02\t100\n'
                             'BAD\tGGGGGK\t0.03\t100\n', encoding='utf-8')
            mapping.write_text('source_id\treference_id\nALIAS\tP1\n', encoding='utf-8')
            result = assess_peptides(truth, calls, mapping_table=mapping)
            self.assertEqual(result['thresholds']['0.01']['identification']['matched_pairs'], 1)
            self.assertEqual(result['thresholds']['0.05']['identification']['matched_pairs'], 2)
            self.assertEqual(result['thresholds']['0.05']['identification']['called_pairs'], 3)
            self.assertEqual(result['quantification']['matched_pairs'], 2)
            self.assertEqual(result['quantification']['median_absolute_log2_ratio_error'], 0.5)
            mapping.write_text('source_id\treference_id\nALIAS\tP1\nALIAS\tP2\n', encoding='utf-8')
            result = assess_peptides(truth, calls, mapping_table=mapping)
            self.assertEqual(result['thresholds']['0.01']['ambiguous_accessions_excluded'], 1)
            self.assertEqual(result['thresholds']['0.01']['identification']['matched_pairs'], 0)


if __name__ == '__main__':
    unittest.main()
