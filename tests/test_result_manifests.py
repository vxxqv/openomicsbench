import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

from omicsbench.evaluate import evaluate_matrix, evaluate_suite
from omicsbench.evidence import create_evidence_crate, verify_evidence_crate
from omicsbench.results import ResultManifest, create_result_manifest, read_result_manifest


class ResultManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.collection = Path(__file__).resolve().parents[1]
        cls.effects = cls.collection / "datasets/rnaseq/rnaseq-002/expected/reference-effects.tsv.gz"
        cls.peptides = cls.collection / "datasets/sequences/sequence-006/expected/peptides.tsv"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.results = self.root / "pipeline-output"
        (self.results / "rna").mkdir(parents=True)
        (self.results / "dna").mkdir()
        (self.results / "protein").mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def write_native_outputs(self):
        with gzip.open(self.effects, "rt", encoding="utf-8", newline="") as source, (self.results / "rna/effects.csv").open("w", encoding="utf-8", newline="") as target:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(target, fieldnames=["Geneid", "log2FoldChange"])
            writer.writeheader()
            for row in reader:
                writer.writerow({"Geneid": row["gene_id"], "log2FoldChange": row["log2_fold_change"]})
        (self.results / "dna/calls.vcf").write_text(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "synthetic_reference\t150\t.\tC\tA\t60\tPASS\t.\n"
            "synthetic_reference\t500\t.\tA\tC\t60\tPASS\t.\n"
            "synthetic_reference\t850\t.\tC\tA\t60\tPASS\t.\n",
            encoding="utf-8",
            newline="\n",
        )
        with self.peptides.open("r", encoding="utf-8", newline="") as source, (self.results / "protein/identifications.tsv").open("w", encoding="utf-8", newline="") as target:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(target, fieldnames=["Protein", "Sequence", "QValue"], delimiter="\t")
            writer.writeheader()
            for row in reader:
                writer.writerow({"Protein": row["protein_id"], "Sequence": row["peptide"], "QValue": "0.001"})

    def create_native_manifest(self):
        return create_result_manifest(
            self.collection,
            self.results,
            [
                "rnaseq-002=rna/effects.csv",
                "sequence-005=dna/calls.vcf",
                "sequence-006=protein/identifications.tsv",
            ],
            [
                "rnaseq-002:gene_id=Geneid",
                "rnaseq-002:log2_fold_change=log2FoldChange",
                "sequence-006:protein_id=Protein",
                "sequence-006:peptide=Sequence",
                "sequence-006:q_value=QValue",
            ],
        )

    def test_nested_native_outputs_are_hash_locked_and_scored(self):
        self.write_native_outputs()
        indexed = self.create_native_manifest()
        self.assertEqual(indexed["summary"]["passed"], 3)
        manifest = json.loads((self.results / "results.json").read_text(encoding="utf-8"))
        self.assertEqual([item["benchmark_id"] for item in manifest["results"]], ["rnaseq-002", "sequence-005", "sequence-006"])
        report = evaluate_suite(self.collection, self.results, ["rnaseq-002", "sequence-005", "sequence-006"])
        self.assertEqual(report["summary"]["status"], "pass")
        self.assertEqual(report["fingerprint_contract"], "openomicsbench-evaluation-v2")
        self.assertEqual(report["result_manifest"]["entries"], 3)
        self.assertEqual(report["results"][0]["details"]["input_columns"]["gene_id"], "Geneid")
        self.assertEqual(report["results"][2]["details"]["input_summary"]["columns"]["q_value"], "QValue")

    def test_changed_result_fails_before_scoring(self):
        self.write_native_outputs()
        self.create_native_manifest()
        with (self.results / "dna/calls.vcf").open("a", encoding="utf-8") as stream:
            stream.write("synthetic_reference\t900\t.\tA\tG\t60\tPASS\t.\n")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            evaluate_suite(self.collection, self.results, ["sequence-005"])

    def test_manifest_rejects_traversal_duplicates_and_unknown_columns(self):
        invalid = {
            "schema_version": "1.0",
            "results": [
                {"benchmark_id": "rnaseq-002", "path": "../outside.tsv", "sha256": "0" * 64},
            ],
        }
        with self.assertRaises(ValueError):
            ResultManifest.model_validate(invalid)
        invalid["results"] = [
            {"benchmark_id": "rnaseq-002", "path": "one.tsv", "sha256": "0" * 64},
            {"benchmark_id": "rnaseq-002", "path": "two.tsv", "sha256": "1" * 64},
        ]
        with self.assertRaisesRegex(ValueError, "unique"):
            ResultManifest.model_validate(invalid)
        self.write_native_outputs()
        with self.assertRaisesRegex(ValueError, "unsupported column"):
            create_result_manifest(
                self.collection,
                self.results,
                ["rnaseq-002=rna/effects.csv"],
                ["rnaseq-002:not_a_field=value"],
            )

    def test_manifest_and_conventional_layout_share_content_fingerprint(self):
        conventional = self.root / "conventional"
        manifested = self.root / "manifested"
        conventional.mkdir()
        manifested.mkdir()
        shutil.copyfile(self.effects, conventional / "rnaseq-002.tsv.gz")
        shutil.copyfile(self.effects, manifested / "native.tsv.gz")
        create_result_manifest(self.collection, manifested, ["rnaseq-002=native.tsv.gz"], [])
        first = evaluate_suite(self.collection, conventional, ["rnaseq-002"])
        second = evaluate_suite(self.collection, manifested, ["rnaseq-002"])
        self.assertEqual(first["run_fingerprint"], second["run_fingerprint"])

    def test_evidence_contains_manifest_and_stable_benchmark_names(self):
        self.write_native_outputs()
        self.create_native_manifest()
        report = evaluate_suite(self.collection, self.results, ["rnaseq-002", "sequence-005", "sequence-006"])
        archive = self.root / "evidence.zip"
        create_evidence_crate(report, archive)
        verified = verify_evidence_crate(archive)
        self.assertEqual(verified["status"], "pass")
        with zipfile.ZipFile(archive) as value:
            self.assertIn("results.json", value.namelist())
            self.assertIn("inputs/rnaseq-002.csv", value.namelist())
            self.assertIn("inputs/sequence-005.vcf", value.namelist())
            self.assertIn("inputs/sequence-006.tsv", value.namelist())

    def test_matrix_evidence_preserves_each_method_manifest(self):
        methods = []
        for name in ("first", "second"):
            directory = self.root / name
            directory.mkdir()
            shutil.copyfile(self.effects, directory / "native.tsv.gz")
            create_result_manifest(self.collection, directory, ["rnaseq-002=native.tsv.gz"], [])
            methods.append((name, directory))
        report = evaluate_matrix(self.collection, methods, ["rnaseq-002"])
        archive = self.root / "matrix-evidence.zip"
        create_evidence_crate(report, archive)
        with zipfile.ZipFile(archive) as value:
            self.assertIn("methods/first/results.json", value.namelist())
            self.assertIn("methods/second/results.json", value.namelist())
            self.assertIn("inputs/first/rnaseq-002.tsv.gz", value.namelist())
            self.assertIn("inputs/second/rnaseq-002.tsv.gz", value.namelist())

    def test_cli_indexes_and_evaluates_native_output(self):
        self.write_native_outputs()
        indexed = subprocess.run(
            [
                sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "index", str(self.results),
                "--result", "rnaseq-002=rna/effects.csv",
                "--column", "rnaseq-002:gene_id=Geneid",
                "--column", "rnaseq-002:log2_fold_change=log2FoldChange",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(indexed.returncode, 0, indexed.stderr)
        self.assertEqual(json.loads(indexed.stdout)["operation"], "index_results")
        evaluated = subprocess.run(
            [sys.executable, "-m", "omicsbench", "--root", str(self.collection), "suite", "evaluate", str(self.results), "--id", "rnaseq-002"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(evaluated.returncode, 0, evaluated.stderr)
        self.assertEqual(json.loads(evaluated.stdout)["summary"]["passed"], 1)

    def test_manifest_loader_rejects_extra_fields(self):
        self.write_native_outputs()
        self.create_native_manifest()
        value = json.loads((self.results / "results.json").read_text(encoding="utf-8"))
        value["unexpected"] = True
        (self.results / "results.json").write_text(json.dumps(value), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "invalid result manifest"):
            read_result_manifest(self.results)


if __name__ == "__main__":
    unittest.main()
