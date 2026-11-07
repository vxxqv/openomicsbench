import shutil
import tempfile
import unittest
from pathlib import Path

from omicsbench import Collection


class ClientTests(unittest.TestCase):
    def test_collection_search_and_verified_file(self):
        collection = Collection()
        records = collection.list(assay="bulk_rna_seq", status="validated", kind="real")
        self.assertEqual(len(records), 12)
        benchmark = collection.get("rnaseq-002")
        self.assertEqual(benchmark.assay, "bulk_rna_seq")
        self.assertEqual(benchmark.manifest["id"], "rnaseq-002")
        self.assertTrue(benchmark.file(role="samples", tier="metadata").is_file())
        with self.assertRaisesRegex(ValueError, "matches 0 files"):
            benchmark.file(role="not-a-role")
        self.assertEqual(benchmark.validate()["status"], "pass")

    def test_python_workflow_from_workspace_to_evidence(self):
        collection = Collection()
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            collection.prepare(workspace, ["sequence-007"])
            self.assertEqual(collection.verify_workspace(workspace)["summary"]["passed"], 1)
            result = workspace / "results/sequence-007.taxonomy.tsv"
            shutil.copyfile(collection.get("sequence-007").file(role="effects", tier="expected"), result)
            self.assertEqual(collection.discover(workspace)["summary"]["passed"], 1)
            self.assertEqual(collection.discover(workspace, write=True)["result_manifest"]["entries"], 1)
            self.assertEqual(collection.evaluate(workspace / "results", ["sequence-007"])["summary"]["passed"], 1)
            self.assertEqual(collection.grade(workspace)["summary"]["passed"], 1)

    def test_file_integrity_check_rejects_changed_content(self):
        collection = Collection()
        with tempfile.TemporaryDirectory() as temporary:
            source = collection.get("sequence-007")
            root = Path(temporary)
            folder = root / "datasets/sequences/sequence-007"
            shutil.copytree(source._folder, folder)
            local = Collection(root).get("sequence-007")
            self.assertTrue(local.file(role="effects", tier="expected").is_file())
            truth = local.file(role="effects", tier="expected")
            truth.write_bytes(truth.read_bytes() + b"extra")
            with self.assertRaisesRegex(ValueError, "missing or changed"):
                local.file(role="effects", tier="expected")


if __name__ == "__main__":
    unittest.main()
