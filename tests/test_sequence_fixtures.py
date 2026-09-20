import copy
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from omicsbench.hashing import digest
from omicsbench.models import Dataset
from omicsbench.registry import lookup, registry
from omicsbench.validate import validate
from workflows.build_sequence_fixtures import build


ROOT = Path(__file__).resolve().parents[1]


class SequenceFixtureTests(unittest.TestCase):
    def test_registered_sequence_fixtures(self):
        records = registry(ROOT)
        identifiers = [identifier for identifier in records if identifier.startswith("sequence-")]
        self.assertEqual(identifiers, ["sequence-001", "sequence-002", "sequence-003", "sequence-004"])
        expected = {
            "sequence-001": ("sequence_dna", "dna", "fasta", False),
            "sequence-002": ("sequence_rna", "rna", "fasta", False),
            "sequence-003": ("sequence_protein", "protein", "fasta", False),
            "sequence-004": ("short_read_dna", "dna", "fastq", True),
        }
        for identifier in identifiers:
            model, _ = records[identifier]
            self.assertEqual((model.assay, model.sequence.molecule, model.sequence.format, model.sequence.paired), expected[identifier])

    def test_all_sequence_fixtures_validate(self):
        for identifier in ["sequence-001", "sequence-002", "sequence-003", "sequence-004"]:
            model, folder = lookup(ROOT, identifier)
            with self.subTest(identifier=identifier):
                self.assertEqual(validate(model, folder)["status"], "pass")

    def test_sequence_summary_drift_is_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary)
            shutil.copytree(ROOT / "datasets", target / "datasets")
            model, folder = lookup(target, "sequence-001")
            path = folder / "nano/sequences.fasta"
            path.write_text(path.read_text().replace("GGGGCCCCAAAATTTT", "GGGGCCCCAAAATTTA", 1))
            record = next(item for item in model.files if item.path == "nano/sequences.fasta")
            record.bytes = path.stat().st_size
            record.sha256 = digest(path)
            with self.assertRaisesRegex(ValueError, "summary differs"):
                validate(model, folder)

    def test_sequence_contract_gates(self):
        model, _ = lookup(ROOT, "sequence-001")
        data = model.model_dump(mode="json")
        for mutation in ("schema", "missing_sequence", "paired_fasta", "quality"):
            value = copy.deepcopy(data)
            if mutation == "schema":
                value["schema_version"] = "1.0"
            elif mutation == "missing_sequence":
                value["sequence"] = None
            elif mutation == "paired_fasta":
                value["sequence"]["paired"] = True
            else:
                value["sequence"]["quality_encoding"] = "phred33"
            with self.subTest(mutation=mutation), self.assertRaises(ValidationError):
                Dataset.model_validate(value)

    def test_fixture_rebuild_is_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            build(destination)
            for original in sorted((ROOT / "datasets/sequences").rglob("*")):
                if not original.is_file():
                    continue
                relative = original.relative_to(ROOT)
                self.assertEqual(digest(original), digest(destination / relative), str(relative))

    def test_cli_lists_and_validates_sequence_objects(self):
        listed = subprocess.run([sys.executable, "-m", "omicsbench", "--root", str(ROOT), "list", "--assay", "sequence_dna"], capture_output=True, text=True)
        self.assertEqual(listed.returncode, 0, listed.stderr)
        self.assertEqual([item["id"] for item in json.loads(listed.stdout)], ["sequence-001"])
        checked = subprocess.run([sys.executable, "-m", "omicsbench", "--root", str(ROOT), "validate", "sequence-004"], capture_output=True, text=True)
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertEqual(json.loads(checked.stdout)["pairing"]["pairs"], 4)


if __name__ == "__main__":
    unittest.main()
