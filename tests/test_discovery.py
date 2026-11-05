import csv
import gzip
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from omicsbench.discovery import discover_workspace_results
from omicsbench.evidence import verify_evidence_crate
from omicsbench.registry import default_root
from omicsbench.workspace import grade_workspace, prepare_workspace


class DiscoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = default_root()

    def test_named_native_outputs_across_six_assays(self):
        ids = ["rnaseq-002", "singlecell-001", "interval-001", "sequence-005", "sequence-006", "sequence-007"]
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ids)
            results = workspace / "results"
            sources = {
                "rnaseq-002": ("datasets/rnaseq/rnaseq-002/expected/reference-effects.tsv.gz", "rnaseq-002/effects.tsv.gz"),
                "singlecell-001": ("datasets/singlecell/singlecell-001/expected/cell-labels.tsv", "singlecell-001_clustering.tsv"),
                "interval-001": ("datasets/intervals/interval-001/expected/truth-peaks.narrowPeak", "interval-001_peaks.narrowPeak"),
                "sequence-006": ("datasets/sequences/sequence-006/expected/peptides.tsv", "sequence-006_peptides.tsv"),
                "sequence-007": ("datasets/sequences/sequence-007/expected/read-taxonomy.tsv", "sequence-007_taxonomy.tsv"),
            }
            for source, target in sources.values():
                destination = results / target
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(self.root / source, destination)
            (results / "sequence-005.vcf").write_text(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "synthetic_reference\t150\t.\tC\tA\t60\tPASS\t.\n"
                "synthetic_reference\t500\t.\tA\tC\t60\tPASS\t.\n"
                "synthetic_reference\t850\t.\tC\tA\t60\tPASS\t.\n",
                encoding="utf-8", newline="\n",
            )
            preview = discover_workspace_results(self.root, workspace)
            self.assertEqual(preview["summary"]["passed"], 6)
            self.assertFalse((results / "results.json").exists())
            written = discover_workspace_results(self.root, workspace, write=True)
            self.assertEqual(written["result_manifest"]["entries"], 6)
            self.assertEqual(grade_workspace(self.root, workspace)["summary"]["passed"], 6)
            self.assertEqual(verify_evidence_crate(workspace / "evidence.zip")["status"], "pass")

    def test_header_aliases_and_explicit_path(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["rnaseq-002"])
            source = self.root / "datasets/rnaseq/rnaseq-002/expected/reference-effects.tsv.gz"
            result = workspace / "results/native/effects.csv"
            result.parent.mkdir()
            with gzip.open(source, "rt", encoding="utf-8", newline="") as stream, result.open("w", encoding="utf-8", newline="") as target:
                reader = csv.DictReader(stream, delimiter="\t")
                writer = csv.DictWriter(target, fieldnames=["Geneid", "log2FoldChange"])
                writer.writeheader()
                for row in reader:
                    writer.writerow({"Geneid": row["gene_id"], "log2FoldChange": row["log2_fold_change"]})
            report = discover_workspace_results(self.root, workspace, ["rnaseq-002=native/effects.csv"], write=True)
            self.assertEqual(report["summary"]["status"], "pass")
            entry = json.loads((workspace / "results/results.json").read_text(encoding="utf-8"))["results"][0]
            self.assertEqual(entry["columns"]["gene_id"], "Geneid")
            self.assertEqual(entry["columns"]["log2_fold_change"], "log2FoldChange")
            self.assertEqual(grade_workspace(self.root, workspace)["summary"]["status"], "pass")

    def test_ambiguity_and_missing_columns_fail_without_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["singlecell-001"])
            results = workspace / "results"
            (results / "singlecell-001_one.tsv").write_text("Barcode\tCluster\nx\ta\n", encoding="utf-8")
            (results / "singlecell-001_two.tsv").write_text("Barcode\tCluster\nx\ta\n", encoding="utf-8")
            report = discover_workspace_results(self.root, workspace)
            self.assertIn("multiple matching files", report["results"][0]["reason"])
            with self.assertRaisesRegex(ValueError, "no manifest was written"):
                discover_workspace_results(self.root, workspace, write=True)
            self.assertFalse((results / "results.json").exists())
            (results / "singlecell-001_two.tsv").unlink()
            (results / "singlecell-001_one.tsv").write_text("unknown\tCluster\nx\ta\n", encoding="utf-8")
            report = discover_workspace_results(self.root, workspace)
            self.assertIn("no cell_id column", report["results"][0]["reason"])

    def test_explicit_column_resolves_ambiguity_and_hash_detects_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "work"
            prepare_workspace(self.root, workspace, ["singlecell-001"])
            result = workspace / "results/singlecell-001_labels.tsv"
            source = self.root / "datasets/singlecell/singlecell-001/expected/cell-labels.tsv"
            with source.open("r", encoding="utf-8", newline="") as stream, result.open("w", encoding="utf-8", newline="") as target:
                reader = csv.DictReader(stream, delimiter="\t")
                writer = csv.DictWriter(target, fieldnames=["Barcode", "cell_id", "Cluster"], delimiter="\t")
                writer.writeheader()
                for row in reader:
                    writer.writerow({"Barcode": row["cell_id"], "cell_id": row["cell_id"], "Cluster": row["label"]})
            self.assertIn("ambiguous cell_id", discover_workspace_results(self.root, workspace)["results"][0]["reason"])
            discover_workspace_results(self.root, workspace, column_values=["singlecell-001:cell_id=Barcode"], write=True)
            with result.open("a", encoding="utf-8") as stream:
                stream.write("extra\tcell\tlabel\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                grade_workspace(self.root, workspace)


if __name__ == "__main__":
    unittest.main()
