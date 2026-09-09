"""The dataset contract; JSON Schema is generated from these models."""
from datetime import date
from pathlib import PurePosixPath
from typing import Literal
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True, allow_inf_nan=False)

class Source(StrictModel):
    repository: str = Field(min_length=1, description="Archive or synthetic generator that supplied the object.")
    accession: str = Field(min_length=1, description="Stable source accession or fixture identifier.")
    url: str = Field(pattern=r"^https://", description="Public HTTPS page for the source record or generator.")
    citation: str = Field(min_length=1, description="Citation that a user should retain when reusing this object.")
    retrieved: date = Field(description="Date when the source material was retrieved or generated.")
    sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$", description="SHA-256 of the original source file when that file can be redistributed.")

class Rights(StrictModel):
    status: Literal["GREEN", "AMBER", "RED", "RECHECK"] = Field(description="Redistribution decision under the documented intake policy.")
    license: str = Field(min_length=1, description="Licence name or the reason redistribution remains restricted.")
    evidence: str = Field(min_length=1, description="Source and reasoning that support the rights decision.")
    checked: date = Field(description="Date when the rights evidence was reviewed.")

class Reference(StrictModel):
    genome: str = Field(min_length=1, description="Genome assembly used for the source quantification.")
    annotation: str = Field(min_length=1, description="Gene annotation release used for the source quantification.")
    namespace: str = Field(min_length=1, description="Identifier system used in the count matrix.")
    compatibility: Literal["verified", "unresolved", "not_applicable"] = Field(description="Result of checking matrix identifiers against the declared annotation.")

class Sample(StrictModel):
    sample_id: str = Field(min_length=1, description="Stable sample or run label in the order used by the count matrix.")
    condition: str = Field(min_length=1, description="Biological group used by the declared comparison.")
    replicate: int = Field(gt=0, strict=True, description="One-based replicate number within the relevant design group.")
    batch: str = Field(min_length=1, description="Blocking or batch value used by the declared model.")
    strandedness: Literal["forward", "reverse", "unstranded", "unknown", "not_applicable"] = Field(description="Library strandedness reported by the source, or unknown when the archive does not state it.")

class Derivation(StrictModel):
    workflow: str = Field(min_length=1, description="Repository path of the workflow that produced the object.")
    commit: str | None = Field(pattern=r"^[0-9a-f]{40}$", description="Exact source-control commit used for the derivation when available.")
    seed: int = Field(ge=0, strict=True, description="Recorded random seed, including zero for deterministic methods that do not draw random values.")
    algorithm: str = Field(min_length=1, description="Short name of the selection or generation method.")
    parameters: dict[str, str | int | float | bool | list[int]] = Field(description="Resolved inputs that affect the generated files.")
    versions: dict[str, str] = Field(min_length=1, description="Software and data versions needed to interpret the derivation.")

def safe_relative(value: str) -> str:
    p = PurePosixPath(value)
    if not value or p.is_absolute() or ".." in p.parts or "\\" in value or ":" in value or str(p) != value:
        raise ValueError("Use a normalized relative path inside the dataset directory")
    return value

class File(StrictModel):
    path: str = Field(description="Normalized path inside the dataset directory.")
    tier: Literal["nano", "pocket", "expected", "metadata"] = Field(description="Size or function class used for retrieval and validation.")
    role: Literal["raw_counts", "samples", "fastq_r1", "fastq_r2", "quantification", "tx2gene", "baseline", "metrics", "figure", "provenance", "documentation", "license", "reference"] = Field(description="File purpose used by the validator and command line tools.")
    sample_id: str | None = Field(default=None, description="Sample linked to this file when the file is sample-specific.")
    media_type: str = Field(min_length=1, description="Internet media type for the file content.")
    bytes: int = Field(ge=0, strict=True, description="Expected file size in bytes.")
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$", description="Expected lowercase SHA-256 digest.")
    url: str | None = Field(default=None, pattern=r"^https://", description="Optional public download URL for remotely retrieved content.")
    _path = field_validator("path")(safe_relative)

class Metric(StrictModel):
    name: Literal["spearman_logfc", "top_k_jaccard", "distance_correlation", "sign_concordance"] = Field(description="Supported preservation measure recomputed by validation.")
    minimum: float = Field(ge=-1, le=1, description="Inclusive pass threshold for the metric.")
    definition: str = Field(min_length=20, description="Human-readable feature universe, calculation and preferred direction.")

class Validation(StrictModel):
    profile: str = Field(description="Path to the expected validation profile inside the dataset directory.")
    baseline_version: str = Field(min_length=1, description="Version label for the expected values and calculation rules.")
    metrics: list[Metric] = Field(min_length=1, description="Quantitative minimums that every validated object must meet.")
    _path = field_validator("profile")(safe_relative)

class Dataset(StrictModel):
    schema_version: Literal["1.0"] = Field(description="Version of the OpenOmicsBench dataset contract.")
    id: str = Field(pattern=r"^(rnaseq|fixture)-[0-9]{3}$", description="Stable collection identifier used by the command line interface.")
    title: str = Field(min_length=8, description="Specific title that distinguishes this comparison from other objects.")
    release: str = Field(pattern=r"^[0-9]+\.[0-9]+\.[0-9]+(?:[.-][A-Za-z0-9.]+)?$", description="First project release that carries this manifest state.")
    assay: Literal["bulk_rna_seq"] = Field(description="Assay family represented by the object.")
    kind: Literal["real", "synthetic_fixture"] = Field(description="Whether the object comes from a biological source or a software fixture.")
    status: Literal["candidate", "validated", "original_link"] = Field(description="Current review and distribution state.")
    archetype: str = Field(min_length=1, description="Short description of the experimental design exercised by the object.")
    organism: str = Field(min_length=1, description="Scientific organism name, or synthetic for a generated fixture.")
    taxon_id: int | None = Field(gt=0, strict=True, description="NCBI taxonomy identifier for a biological object.")
    source: Source = Field(description="Origin, citation and retrieval record.")
    rights: Rights = Field(description="Dated redistribution decision and supporting evidence.")
    reference: Reference = Field(description="Genome, annotation and identifier compatibility record.")
    samples: list[Sample] = Field(min_length=2, description="Ordered experimental design matching the matrix columns or read files.")
    derivation: Derivation = Field(description="Workflow, code revision and resolved parameters used to create the object.")
    files: list[File] = Field(description="Complete declared inventory for the dataset directory.")
    validation: Validation | None = Field(description="Expected-result profile for a validated object, or null when validation does not apply.")
    limitations: list[str] = Field(min_length=1, description="Known constraints that affect interpretation or reuse.")

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
