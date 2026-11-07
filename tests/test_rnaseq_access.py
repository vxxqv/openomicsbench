import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from omicsbench import Collection
from omicsbench.hashing import digest
from omicsbench.rnaseq import read_counts


class RNASeqAccessTests(unittest.TestCase):
    def setUp(self):
        self.collection = Collection()

    def test_all_biological_objects_have_aligned_tiers(self):
        records = self.collection.list(assay="bulk_rna_seq", status="validated", kind="real")
        self.assertEqual(len(records), 12)
        for record in records:
            benchmark = self.collection.get(record["id"])
            for tier in ("pocket", "full"):
                with self.subTest(id=benchmark.id, tier=tier):
                    data = benchmark.rnaseq(tier)
                    self.assertEqual(data.counts.shape, (len(data.genes), len(data.samples)))
                    self.assertEqual(data.samples, tuple(row["sample_id"] for row in data.metadata))
                    self.assertTrue(np.issubdtype(data.counts.dtype, np.integer))
                    self.assertFalse(data.counts.flags.writeable)
                    self.assertTrue(data.source["accession"])

    def test_export_is_aligned_and_hash_recorded(self):
        benchmark = self.collection.get("rnaseq-002")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "export"
            result = benchmark.export_rnaseq(output, "pocket")
            genes, samples, counts = read_counts(output / "counts.tsv")
            with (output / "samples.tsv").open(newline="", encoding="utf-8") as stream:
                metadata = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(samples, [row["sample_id"] for row in metadata])
            self.assertEqual(genes, list(benchmark.rnaseq().genes))
            np.testing.assert_array_equal(counts, benchmark.rnaseq().counts)
            self.assertEqual(result, json.loads((output / "export.json").read_text(encoding="utf-8")))
            for name, record in result["files"].items():
                self.assertEqual(record["sha256"], digest(output / name))
                self.assertEqual(record["bytes"], (output / name).stat().st_size)
            with self.assertRaises(FileExistsError):
                benchmark.export_rnaseq(output)

    def test_rejects_unsupported_and_modified_inputs(self):
        with self.assertRaisesRegex(ValueError, "only validated biological"):
            self.collection.get("rnaseq-001").rnaseq()
        with self.assertRaisesRegex(ValueError, "only validated biological"):
            self.collection.get("singlecell-001").rnaseq()
        with self.assertRaisesRegex(ValueError, "tier must"):
            self.collection.get("rnaseq-002").rnaseq("nano")
        with tempfile.TemporaryDirectory() as temporary:
            source = self.collection.get("rnaseq-002")
            folder = Path(temporary) / "datasets/rnaseq/rnaseq-002"
            shutil.copytree(source._folder, folder)
            local = Collection(temporary).get("rnaseq-002")
            path = folder / "pocket/counts.tsv"
            path.write_bytes(path.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "missing or changed"):
                local.rnaseq()


if __name__ == "__main__":
    unittest.main()
