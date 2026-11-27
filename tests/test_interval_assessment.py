import json
import tempfile
import unittest
from pathlib import Path

from omicsbench.interval_assessment import assess_intervals


ROOT = Path(__file__).resolve().parents[1]
PEAKS = ROOT / 'datasets/intervals/interval-001/expected/truth-peaks.narrowPeak'


class IntervalAssessmentTests(unittest.TestCase):
    def test_reference_replicate_and_imported_qc_are_separate(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            qc = base / 'qc.json'
            idr = base / 'idr.json'
            qc.write_text(json.dumps({'source': 'external pipeline', 'frip': 0.4,
                                      'tss_enrichment': 8.0}), encoding='utf-8')
            idr.write_text(json.dumps({'source': 'external IDR run', 'idr_peaks': 4}), encoding='utf-8')
            report = assess_intervals(PEAKS, replicate=PEAKS, reference=PEAKS,
                                      qc_json=qc, idr_json=idr)
            self.assertEqual(report['reference_recovery']['interval_f1'], 1)
            self.assertEqual(report['replicate_overlap']['interval_f1'], 1)
            self.assertEqual(report['supplied_qc']['frip'], 0.4)
            self.assertEqual(report['supplied_idr']['idr_peaks'], 4)
            self.assertTrue(report['width_sensitivity'])
            qc.write_text(json.dumps({'source': 'external pipeline', 'frip': 1.5}), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'frip'):
                assess_intervals(PEAKS, qc_json=qc)


if __name__ == '__main__':
    unittest.main()
