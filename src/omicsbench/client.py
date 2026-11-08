"""Small Python interface to the installed benchmark collection."""

from __future__ import annotations

from pathlib import Path

from .discovery import discover_workspace_results
from .evaluate import evaluate_suite
from .hashing import contained, digest
from .registry import default_root, lookup, registry
from .validate import validate
from .workspace import grade_workspace, prepare_workspace, verify_workspace


class Benchmark:
    def __init__(self, model, folder: Path):
        self._model = model
        self._folder = folder

    @property
    def id(self) -> str:
        return self._model.id

    @property
    def assay(self) -> str:
        return self._model.assay

    @property
    def title(self) -> str:
        return self._model.title

    @property
    def manifest(self) -> dict:
        return self._model.model_dump(mode="json")

    @property
    def files(self) -> list[dict]:
        return [item.model_dump(mode="json") for item in self._model.files]

    def file(self, *, role: str, tier: str | None = None) -> Path:
        """Return one declared file after checking its size and SHA-256 digest."""
        matches = [item for item in self._model.files if item.role == role and (tier is None or item.tier == tier)]
        if len(matches) != 1:
            raise ValueError(f"{self.id}: role {role!r} at tier {tier!r} matches {len(matches)} files; choose an exact role and tier")
        item = matches[0]
        path = contained(self._folder, item.path)
        if not path.is_file() or path.stat().st_size != item.bytes or digest(path) != item.sha256:
            raise ValueError(f"{self.id}: declared file is missing or changed: {item.path}")
        return path

    def validate(self) -> dict:
        return validate(self._model, self._folder)

    def rnaseq(self, tier: str = "pocket"):
        """Load aligned biological RNA counts and sample metadata."""
        from .rnaseq_access import load_rnaseq
        return load_rnaseq(self, tier)

    def export_rnaseq(self, destination: str | Path, tier: str = "pocket") -> dict:
        """Export aligned counts and metadata to a new analysis directory."""
        from .rnaseq_access import export_rnaseq
        return export_rnaseq(self, destination, tier)

    def compare_counts(self, path: str | Path, tier: str = "pocket", detail_limit: int = 20) -> dict:
        """Check an integer count matrix against a bundled RNA-seq tier."""
        from .rnaseq_access import compare_counts
        return compare_counts(self, path, tier, detail_limit)

    def verify_rnaseq_export(self, directory: str | Path) -> dict:
        """Check an exported directory, including its counts and metadata."""
        from .rnaseq_access import verify_export
        return verify_export(self, directory)


class Collection:
    """Access local benchmarks and score pipeline results without a CLI subprocess."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root if root is not None else default_root()).resolve()
        if not (self.root / "datasets").is_dir():
            raise ValueError(f"{self.root}: benchmark collection is missing")

    def list(self, *, assay: str | None = None, status: str | None = None,
             kind: str | None = None) -> list[dict]:
        records = registry(self.root)
        return [
            {"id": model.id, "title": model.title, "assay": model.assay, "kind": model.kind,
             "status": model.status, "limitations": model.limitations}
            for model, _ in records.values()
            if (assay is None or model.assay == assay) and (status is None or model.status == status)
            and (kind is None or model.kind == kind)
        ]

    def get(self, identifier: str) -> Benchmark:
        model, folder = lookup(self.root, identifier)
        return Benchmark(model, folder)

    def prepare(self, destination: str | Path, identifiers: list[str]) -> dict:
        return prepare_workspace(self.root, Path(destination), identifiers)

    def verify_workspace(self, directory: str | Path) -> dict:
        return verify_workspace(self.root, Path(directory))

    def discover(self, directory: str | Path, *, results: list[str] | None = None,
                 columns: list[str] | None = None, write: bool = False, force: bool = False) -> dict:
        return discover_workspace_results(self.root, Path(directory), results, columns, write, force)

    def grade(self, directory: str | Path, *, force: bool = False) -> dict:
        return grade_workspace(self.root, Path(directory), force)

    def evaluate(self, directory: str | Path, identifiers: list[str] | None = None,
                 *, detail_limit: int = 20) -> dict:
        return evaluate_suite(self.root, Path(directory), identifiers, detail_limit)
