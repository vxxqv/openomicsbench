import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from omicsbench.cell_assessment import assess_cells
from omicsbench.registry import lookup


ROOT = Path(__file__).resolve().parents[1]


class CellAssessmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model, cls.folder = lookup(ROOT, "singlecell-001")
        cls.truth = cls.folder / "expected/cell-labels.tsv"

    def test_clustering_annotation_and_embedding_are_separate(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            clusters = directory / "clusters.tsv"
            embedding = directory / "embedding.tsv"
            with self.truth.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            with clusters.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t")
                writer.writerow(("cell_id", "label"))
                for row in reversed(rows):
                    writer.writerow((row["cell_id"], row["label"].replace("state-", "cluster-")))
            with embedding.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t")
                writer.writerow(("cell_id", "x", "y"))
                for index, row in enumerate(reversed(rows)):
                    writer.writerow((row["cell_id"], int(row["label"][-1]) * 10 + index / 1000, index / 100))
            result = assess_cells(self.model, self.folder, clusters, annotation=self.truth,
                                  embedding=embedding, neighbors=5)
            self.assertEqual(result["status"], "pass")
            self.assertEqual(result["clustering"]["metrics"]["adjusted_rand_index"], 1)
            self.assertEqual(result["annotation"]["exact_label_accuracy"], 1)
            self.assertIn("mean_same_batch_neighbor_fraction", result["embedding"])
            cli = subprocess.run([sys.executable, "-m", "omicsbench", "--root", str(ROOT),
                                  "cell", "assess", "singlecell-001", str(clusters),
                                  "--annotation", str(self.truth), "--embedding", str(embedding),
                                  "--neighbors", "5"], capture_output=True, text=True)
            self.assertEqual(cli.returncode, 0, cli.stderr)
            self.assertEqual(json.loads(cli.stdout)["status"], "pass")

    def test_mismatched_embedding_ids_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            embedding = Path(temporary) / "embedding.tsv"
            embedding.write_text("cell_id\tx\ty\nCELL-001\t0\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exactly match"):
                assess_cells(self.model, self.folder, self.truth, embedding=embedding)

    def test_incorrect_annotation_is_reported_without_renaming(self):
        with tempfile.TemporaryDirectory() as temporary:
            annotation = Path(temporary) / "annotation.tsv"
            with self.truth.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            with annotation.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t")
                writer.writerow(("cell_id", "label"))
                for row in rows:
                    writer.writerow((row["cell_id"], row["label"].replace("state-", "cluster-")))
            result = assess_cells(self.model, self.folder, self.truth, annotation=annotation)
            self.assertEqual(result["clustering"]["status"], "pass")
            self.assertEqual(result["annotation"]["exact_label_accuracy"], 0)
            self.assertEqual(result["status"], "fail")


if __name__ == "__main__":
    unittest.main()
