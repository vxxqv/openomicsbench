"""The dataset contract; JSON Schema is generated from these models."""
from datetime import date
from pathlib import PurePosixPath
from typing import Literal
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True, allow_inf_nan=False)

class Source(StrictModel):
    repository: str = Field(min_length=1, description="Archive or synthetic generator that supplied the object.")
    accession: str = Field(min_length=1)
    url: str = Field(pattern=r"^https://")
    citation: str = Field(min_length=1)
    retrieved: date
    sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")

class Rights(StrictModel):
    status: Literal["GREEN", "AMBER", "RED", "RECHECK"]
    license: str = Field(min_length=1)
    evidence: str = Field(min_length=1)
    checked: date

class Reference(StrictModel):
    genome: str = Field(min_length=1)
    annotation: str = Field(min_length=1)
    namespace: str = Field(min_length=1)
    compatibility: Literal["verified", "unresolved", "not_applicable"]

class Sample(StrictModel):
    sample_id: str = Field(min_length=1)
    condition: str = Field(min_length=1)
    replicate: int = Field(gt=0, strict=True)
    batch: str = Field(min_length=1)
    strandedness: Literal["forward", "reverse", "unstranded", "unknown", "not_applicable"]

class Derivation(StrictModel):
    workflow: str = Field(min_length=1)
    commit: str | None = Field(pattern=r"^[0-9a-f]{40}$")
    seed: int = Field(ge=0, strict=True)
    algorithm: str = Field(min_length=1)
    parameters: dict[str, str | int | float | bool | list[int]]
    versions: dict[str, str] = Field(min_length=1)

def safe_relative(value: str) -> str:
    p = PurePosixPath(value)
    if not value or p.is_absolute() or ".." in p.parts or "\\" in value or ":" in value or str(p) != value:
        raise ValueError("Use a normalized relative path inside the dataset directory")
    return value

class File(StrictModel):
    path: str
    tier: Literal["nano", "pocket", "expected", "metadata"]
    role: Literal["raw_counts", "samples", "fastq_r1", "fastq_r2", "quantification", "tx2gene", "baseline", "metrics", "figure", "provenance", "documentation", "license", "reference"]
    sample_id: str | None = None
    media_type: str = Field(min_length=1)
    bytes: int = Field(ge=0, strict=True)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    url: str | None = Field(default=None, pattern=r"^https://")
    _path = field_validator("path")(safe_relative)

class Metric(StrictModel):
    name: Literal["spearman_logfc", "top_k_jaccard", "distance_correlation", "sign_concordance"]
    minimum: float = Field(ge=-1, le=1)
    definition: str = Field(min_length=20)

class Validation(StrictModel):
    profile: str
    baseline_version: str = Field(min_length=1)
    metrics: list[Metric] = Field(min_length=1)
    _path = field_validator("profile")(safe_relative)

class Dataset(StrictModel):
    schema_version: Literal["1.0"]
    id: str = Field(pattern=r"^(rnaseq|fixture)-[0-9]{3}$")
    title: str = Field(min_length=8)
    release: str = Field(pattern=r"^[0-9]+\.[0-9]+\.[0-9]+(?:[.-][A-Za-z0-9.]+)?$")
    assay: Literal["bulk_rna_seq"]
    kind: Literal["real", "synthetic_fixture"]
    status: Literal["candidate", "validated", "original_link"]
    archetype: str = Field(min_length=1)
    organism: str = Field(min_length=1)
    taxon_id: int | None = Field(gt=0, strict=True)
    source: Source
    rights: Rights
    reference: Reference
    samples: list[Sample] = Field(min_length=2)
    derivation: Derivation
    files: list[File]
    validation: Validation | None
    limitations: list[str] = Field(min_length=1)

    @model_validator(mode="after")
    def consistent(self):
        if self.kind == "real" and self.taxon_id is None:
            raise ValueError("Biological objects require an NCBI taxonomy identifier")
        ids = [s.sample_id for s in self.samples]
        if len(ids) != len(set(ids)):
            raise ValueError("Sample IDs must be unique and ordered")
        paths = [f.path for f in self.files]
        if len(paths) != len(set(paths)):
            raise ValueError("A file path is declared more than once")
        keys = [(f.tier, f.role, f.sample_id) for f in self.files if f.role not in {"figure", "provenance", "documentation", "license"}]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate logical file role within a tier and sample")
        if any(f.sample_id is not None and f.sample_id not in ids for f in self.files):
            raise ValueError("File sample_id is absent from the sample design")
        if self.status == "original_link" and any(f.tier in {"nano", "pocket"} for f in self.files):
            raise ValueError("Link-only objects cannot distribute data tiers")
        if self.rights.status != "GREEN" and any(f.tier in {"nano", "pocket"} for f in self.files):
            raise ValueError("Distributed data need a GREEN rights decision")
        if self.status == "validated":
            if not self.validation or not self.files or self.reference.compatibility == "unresolved":
                raise ValueError("Validated objects need files, validation and resolved references")
        if self.validation and self.validation.profile not in paths:
            raise ValueError("Validation profile must be a declared file")
        return self
