"""High-level, atomic operations for FASTA and FASTQ files."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Iterable, Iterator

from .sequences import (
    SequenceRecord,
    ambiguity_fraction,
    compare_sketches,
    count_kmers,
    detect_format,
    find_motifs,
    gc_fraction,
    phred_scores,
    read_records,
    resolve_molecule,
    reverse_complement,
    sketch,
    summarize,
    translate,
    trim_record,
    write_records,
)


def _temporary_for(path: Path) -> Path:
    path = Path(path)
    suffix = ".tmp.gz" if path.suffix.lower() == ".gz" else ".tmp"
    return path.with_name(path.name + suffix)


def atomic_write(records: Iterable[SequenceRecord], path: Path, file_format: str = "auto", wrap: int = 80, force: bool = False) -> int:
    path = Path(path)
    if path.exists() and not force:
        raise ValueError(f"{path}: output already exists; use --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_for(path)
    try:
        count = write_records(records, temporary, file_format, wrap)
        os.replace(temporary, path)
        return count
    finally:
        temporary.unlink(missing_ok=True)


def prepare_output(path: Path, force: bool) -> None:
    if Path(path).exists():
        if not force:
            raise ValueError(f"{path}: output already exists; pass --force to replace it")


def validate_file(path: Path, molecule: str = "auto", file_format: str = "auto") -> dict:
    report = summarize(Path(path), molecule, file_format)
    return {"status": "pass", **report}


def convert_file(source: Path, destination: Path, output_format: str, wrap: int, force: bool) -> dict:
    input_format = detect_format(source)
    if output_format == "fastq" and input_format != "fastq":
        raise ValueError("FASTA cannot be converted to FASTQ without invented quality scores")
    prepare_output(destination, force)
    count = atomic_write(read_records(source), destination, output_format, wrap, force)
    return {"input": str(source), "output": str(destination), "input_format": input_format, "output_format": output_format, "records": count}


def transform_file(source: Path, destination: Path, operation: str, molecule: str, frame: int, trim_stop: bool, force: bool) -> dict:
    selected_molecule, records = resolve_molecule(read_records(source), molecule)
    if operation == "reverse-complement":
        if selected_molecule == "protein":
            raise ValueError("reverse complement requires DNA or RNA")
        transformed = (
            SequenceRecord(record.identifier, record.description, reverse_complement(record.sequence, selected_molecule), record.quality[::-1] if record.quality is not None else None)
            for record in records
        )
        output_format = detect_format(source)
    elif operation == "transcribe":
        if selected_molecule != "dna":
            raise ValueError("transcription requires DNA input")
        transformed = (SequenceRecord(record.identifier, record.description, record.sequence.replace("T", "U")) for record in records)
        output_format = "fasta"
    elif operation == "back-transcribe":
        if selected_molecule != "rna":
            raise ValueError("back-transcription requires RNA input")
        transformed = (SequenceRecord(record.identifier, record.description, record.sequence.replace("U", "T")) for record in records)
        output_format = "fasta"
    elif operation == "translate":
        if selected_molecule not in {"dna", "rna"}:
            raise ValueError("translation requires DNA or RNA")
        transformed = (
            SequenceRecord(record.identifier, f"{record.description} frame={frame}".strip(), translate(record.sequence, frame, trim_stop))
            for record in records
        )
        output_format = "fasta"
    else:
        raise ValueError(f"unknown transformation: {operation}")
    prepare_output(destination, force)
    count = atomic_write(transformed, destination, output_format, force=force)
    return {"input": str(source), "output": str(destination), "operation": operation, "molecule": selected_molecule, "records": count}


def filter_file(
    source: Path,
    destination: Path,
    molecule: str,
    minimum_length: int,
    maximum_length: int | None,
    minimum_gc: float | None,
    maximum_gc: float | None,
    maximum_ambiguity: float,
    minimum_mean_quality: float | None,
    force: bool,
) -> dict:
    if minimum_length < 0 or maximum_length is not None and maximum_length < minimum_length:
        raise ValueError("length bounds are invalid")
    for value in (minimum_gc, maximum_gc, maximum_ambiguity):
        if value is not None and not 0 <= value <= 1:
            raise ValueError("GC and ambiguity fractions must be between 0 and 1")
    selected_molecule, records = resolve_molecule(read_records(source), molecule)
    if selected_molecule == "protein" and (minimum_gc is not None or maximum_gc is not None):
        raise ValueError("GC filters require DNA or RNA")
    kept: list[SequenceRecord] = []
    reasons = {"short": 0, "long": 0, "low_gc": 0, "high_gc": 0, "ambiguous": 0, "low_quality": 0}
    for record in records:
        length = len(record.sequence)
        reason = None
        value_gc = gc_fraction(record.sequence)
        value_ambiguity = ambiguity_fraction(record.sequence, selected_molecule)
        mean_quality = sum(phred_scores(record.quality)) / length if record.quality is not None else None
        if length < minimum_length:
            reason = "short"
        elif maximum_length is not None and length > maximum_length:
            reason = "long"
        elif minimum_gc is not None and value_gc < minimum_gc:
            reason = "low_gc"
        elif maximum_gc is not None and value_gc > maximum_gc:
            reason = "high_gc"
        elif value_ambiguity > maximum_ambiguity:
            reason = "ambiguous"
        elif minimum_mean_quality is not None:
            if mean_quality is None:
                raise ValueError("mean-quality filtering requires FASTQ input")
            if mean_quality < minimum_mean_quality:
                reason = "low_quality"
        if reason:
            reasons[reason] += 1
        else:
            kept.append(record)
    prepare_output(destination, force)
    count = atomic_write(kept, destination, detect_format(source), force=force)
    return {"input": str(source), "output": str(destination), "molecule": selected_molecule, "input_records": len(records), "kept_records": count, "removed_records": len(records) - count, "removed_by_first_reason": reasons}


def trim_file(source: Path, destination: Path, left: int, right: int, quality: int | None, adapter: str | None, minimum_length: int, force: bool) -> dict:
    if minimum_length < 0 or quality is not None and not 0 <= quality <= 93:
        raise ValueError("minimum length or quality threshold is invalid")
    input_records = list(read_records(source))
    trimmed = [trim_record(record, left, right, quality, adapter) for record in input_records]
    kept = [record for record in trimmed if len(record.sequence) >= minimum_length]
    prepare_output(destination, force)
    count = atomic_write(kept, destination, detect_format(source), force=force)
    return {
        "input": str(source), "output": str(destination), "input_records": len(input_records), "kept_records": count,
        "discarded_short": len(trimmed) - count,
        "bases_removed": sum(len(before.sequence) - len(after.sequence) for before, after in zip(input_records, trimmed)),
    }


def _record_hash(record: SequenceRecord, ordinal: int, seed: int) -> int:
    value = f"{seed}:{ordinal}:{record.identifier}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(value).digest()[:8], "big")


def sample_file(source: Path, destination: Path, count: int | None, fraction: float | None, seed: int, force: bool) -> dict:
    if (count is None) == (fraction is None):
        raise ValueError("choose exactly one of --count or --fraction")
    if count is not None and count < 1:
        raise ValueError("count must be at least 1")
    if fraction is not None and not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    records = list(read_records(source))
    if count is not None:
        if count > len(records):
            raise ValueError("sample count exceeds the number of records")
        ranked = sorted((_record_hash(record, index, seed), index, record) for index, record in enumerate(records))[:count]
        selected = [item[2] for item in sorted(ranked, key=lambda item: item[1])]
    else:
        threshold = fraction * 2**64
        selected = [record for index, record in enumerate(records) if _record_hash(record, index, seed) < threshold]
    prepare_output(destination, force)
    written = atomic_write(selected, destination, detect_format(source), force=force)
    return {"input": str(source), "output": str(destination), "seed": seed, "input_records": len(records), "selected_records": written, "count": count, "fraction": fraction}


def deduplicate_file(source: Path, destination: Path, by: str, force: bool) -> dict:
    seen: set[str] = set()
    kept: list[SequenceRecord] = []
    total = 0
    for record in read_records(source):
        total += 1
        key = record.sequence if by == "sequence" else record.identifier
        if key in seen:
            continue
        seen.add(key)
        kept.append(record)
    prepare_output(destination, force)
    written = atomic_write(kept, destination, detect_format(source), force=force)
    return {"input": str(source), "output": str(destination), "key": by, "input_records": total, "kept_records": written, "removed_records": total - written}


def kmer_report(source: Path, k: int, canonical: bool, top: int, maximum_distinct: int) -> dict:
    if top < 0:
        raise ValueError("top must be zero or greater")
    counts = count_kmers(read_records(source), k, canonical, maximum_distinct)
    total = sum(counts.values())
    return {"input": str(source), "k": k, "canonical": canonical, "total_kmers": total, "distinct_kmers": len(counts), "top": [{"kmer": value, "count": count} for value, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:top]]}


def motif_report(source: Path, motif: str, both_strands: bool) -> dict:
    return {"input": str(source), **find_motifs(read_records(source), motif, both_strands)}


def create_sketch(source: Path, destination: Path | None, k: int, size: int, canonical: bool, force: bool) -> dict:
    result = {"input": str(source), **sketch(read_records(source), k, size, canonical)}
    if destination is not None:
        prepare_output(destination, force)
        temporary = _temporary_for(destination)
        try:
            temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
    return result


def compare_sequence_files(left: Path, right: Path, k: int, size: int, canonical: bool) -> dict:
    a = sketch(read_records(left), k, size, canonical)
    b = sketch(read_records(right), k, size, canonical)
    return {"left_path": str(left), "right_path": str(right), "k": k, "sketch_size": size, **compare_sketches(a, b)}
