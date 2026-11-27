import tempfile
import unittest
from pathlib import Path

from omicsbench.variant_benchmark import benchmark_variants


HEADER = '##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


class VariantBenchmarkTests(unittest.TestCase):
    def test_regions_strata_and_assembly(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            truth, calls, regions = (base / name for name in ('truth.vcf', 'calls.vcf', 'regions.bed'))
            truth.write_text(HEADER + 'chr1\t2\t.\tA\tG\t.\tPASS\t.\n'
                             'chr1\t5\t.\tAT\tA\t.\tPASS\t.\n'
                             'chr1\t20\t.\tC\tT\t.\tPASS\t.\n', encoding='utf-8')
            calls.write_text(HEADER + 'chr1\t2\t.\tA\tG\t.\tPASS\t.\n'
                             'chr1\t5\t.\tAT\tA\t.\tPASS\t.\n'
                             'chr1\t20\t.\tC\tA\t.\tPASS\t.\n'
                             'chr1\t10\t.\tC\tA\t.\tPASS\t.\n', encoding='utf-8')
            regions.write_text('chr1\t0\t6\n', encoding='utf-8')
            result = benchmark_variants(truth, calls, regions, assembly='GRCh38', truth_assembly='GRCh38')
            self.assertEqual(result['overall']['f1'], 1)
            self.assertEqual(result['by_type']['SNV']['true_positive'], 1)
            self.assertEqual(result['by_type']['indel']['true_positive'], 1)
            self.assertEqual(result['outside_regions'], {'truth': 1, 'calls': 2})
            with self.assertRaisesRegex(ValueError, 'assemblies must match'):
                benchmark_variants(truth, calls, regions, assembly='GRCh37', truth_assembly='GRCh38')
            regions.write_text('chr1\t0\t5\n', encoding='utf-8')
            result = benchmark_variants(truth, calls, regions, assembly='GRCh38', truth_assembly='GRCh38')
            self.assertEqual(result['by_type']['indel']['truth'], 0)
            self.assertEqual(result['outside_regions']['truth'], 2)


if __name__ == '__main__':
    unittest.main()
