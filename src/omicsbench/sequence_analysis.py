"""Sequence-level QC, paired-read checks and small deterministic analyses."""
from __future__ import annotations

import math
from collections import Counter
from itertools import zip_longest
from pathlib import Path

from .sequence_tools import atomic_write, prepare_output
from .sequences import (
    CODONS,
    SequenceRecord,
    detect_format,
    phred_scores,
    read_records,
    resolve_molecule,
    reverse_complement,
    summarize,
)


def read_pair_key(identifier: str) -> str:
    for suffix in ("/1", "/2", ".1", ".2"):
        if identifier.endswith(suffix):
            return identifier[:-2]
    return identifier


def pair_report(left: Path, right: Path, detail_limit: int = 20) -> dict:
    if detail_limit < 0:
        raise ValueError("detail limit must be zero or greater")
    mismatches = []
    pairs = 0
    for ordinal, (first, second) in enumerate(zip_longest(read_records(left), read_records(right)), start=1):
        if first is None or second is None:
            mismatches.append({"pair": ordinal, "reason": "record count differs"})
            continue
        pairs += 1
        if first.quality is None or second.quality is None:
            mismatches.append({"pair": ordinal, "reason": "paired input must be FASTQ"})
        elif read_pair_key(first.identifier) != read_pair_key(second.identifier):
            mismatches.append({"pair": ordinal, "left": first.identifier, "right": second.identifier, "reason": "read identifiers differ"})
    return {"status": "pass" if not mismatches else "fail", "left": str(left), "right": str(right), "pairs": pairs, "mismatches": len(mismatches), "mismatch_examples": mismatches[:detail_limit], "examples_truncated": len(mismatches) > detail_limit}


def interleave_pairs(left: Path, right: Path, destination: Path, force: bool) -> dict:
    report = pair_report(left, right)
    if report["status"] != "pass":
        raise ValueError("paired FASTQ files differ; run seq pair-check for details")

    def records():
        for first, second in zip(read_records(left), read_records(right)):
            yield first
            yield second

    prepare_output(destination, force)
    written = atomic_write(records(), destination, "fastq", force=force)
    return {**report, "output": str(destination), "records": written}


def deinterleave_pairs(source: Path, left: Path, right: Path, force: bool) -> dict:
    prepare_output(left, force)
    prepare_output(right, force)
    records = list(read_records(source))
    if len(records) % 2:
        raise ValueError("interleaved FASTQ contains an odd number of records")
    first, second = records[::2], records[1::2]
    for ordinal, (a, b) in enumerate(zip(first, second), start=1):
        if a.quality is None or b.quality is None or read_pair_key(a.identifier) != read_pair_key(b.identifier):
            raise ValueError(f"interleaved FASTQ pair {ordinal} is invalid")
    atomic_write(first, left, "fastq", force=force)
    try:
        atomic_write(second, right, "fastq", force=force)
    except Exception:
        left.unlink(missing_ok=True)
        raise
    return {"input": str(source), "left_output": str(left), "right_output": str(right), "pairs": len(first)}


def positional_quality(records: list[SequenceRecord], maximum_positions: int = 300) -> list[dict]:
    if maximum_positions < 1:
        raise ValueError("maximum positions must be at least 1")
    totals = [0] * maximum_positions
    observations = [0] * maximum_positions
    q20 = [0] * maximum_positions
    q30 = [0] * maximum_positions
    bases = [Counter() for _ in range(maximum_positions)]
    for record in records:
        if record.quality is None:
            raise ValueError("positional quality requires FASTQ input")
        for index, (base, score) in enumerate(zip(record.sequence, phred_scores(record.quality))):
            if index >= maximum_positions:
                break
            totals[index] += score
            observations[index] += 1
            q20[index] += score >= 20
            q30[index] += score >= 30
            bases[index][base] += 1
    return [
        {
            "position": index + 1,
            "observations": count,
            "mean_quality": totals[index] / count,
            "q20_fraction": q20[index] / count,
            "q30_fraction": q30[index] / count,
            "bases": dict(sorted(bases[index].items())),
        }
        for index, count in enumerate(observations) if count
    ]


def shannon_complexity(sequence: str, k: int = 2) -> float:
    if k < 1:
        raise ValueError("complexity k must be at least 1")
    if len(sequence) < k:
        return 0.0
    counts = Counter(sequence[index:index + k] for index in range(len(sequence) - k + 1))
    total = sum(counts.values())
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    maximum = math.log2(min(4**k, total))
    return entropy / maximum if maximum else 0.0


def qc_report(path: Path, molecule: str = "auto", maximum_positions: int = 300, overrepresented: int = 20, adapters: list[str] | None = None) -> dict:
    selected_molecule, records = resolve_molecule(read_records(path), molecule)
    summary = summarize(path, selected_molecule)
    sequences = Counter(record.sequence for record in records)
    complexity = [shannon_complexity(record.sequence.replace("U", "T")) for record in records if selected_molecule in {"dna", "rna"}]
    adapter_counts = {}
    for adapter in adapters or []:
        value = adapter.upper()
        adapter_counts[value] = sum(value in record.sequence for record in records)
    report = {
        "summary": summary,
        "low_complexity_fraction": sum(value < 0.5 for value in complexity) / len(complexity) if complexity else None,
        "mean_complexity": sum(complexity) / len(complexity) if complexity else None,
        "overrepresented_sequences": [
            {"sequence": sequence, "count": count, "fraction": count / len(records)}
            for sequence, count in sorted(sequences.items(), key=lambda item: (-item[1], item[0]))[:overrepresented]
        ],
        "adapter_hits": adapter_counts,
    }
    if detect_format(path) == "fastq":
        report["per_position"] = positional_quality(records, maximum_positions)
    return report


def find_orfs(path: Path, minimum_amino_acids: int = 30, include_partial: bool = False, maximum_results: int = 1000) -> dict:
    if minimum_amino_acids < 1 or maximum_results < 1:
        raise ValueError("minimum amino acids and maximum results must be at least 1")
    molecule, records = resolve_molecule(read_records(path), "auto")
    if molecule not in {"dna", "rna"}:
        raise ValueError("ORF discovery requires DNA or RNA")
    results = []
    starts = {"ATG"}
    stops = {"TAA", "TAG", "TGA"}
    for record in records:
        forward = record.sequence.replace("U", "T")
        length = len(forward)
        for strand, sequence in (("+", forward), ("-", reverse_complement(forward, "dna"))):
            for offset in range(3):
                open_starts: list[int] = []
                for position in range(offset, len(sequence) - 2, 3):
                    codon = sequence[position:position + 3]
                    if codon in starts:
                        open_starts.append(position)
                    if codon not in stops or not open_starts:
                        continue
                    for start in open_starts:
                        end = position + 3
                        amino_acids = (end - start) // 3 - 1
                        if amino_acids >= minimum_amino_acids:
                            left, right = (start, end) if strand == "+" else (length - end, length - start)
                            protein = "".join(CODONS.get(sequence[index:index + 3], "X") for index in range(start, end, 3))
                            results.append({"id": record.identifier, "start": left, "end": right, "strand": strand, "frame": offset + 1 if strand == "+" else -(offset + 1), "amino_acids": amino_acids, "complete": True, "protein": protein})
                    open_starts = []
                if include_partial:
                    end = len(sequence) - (len(sequence) - offset) % 3
                    for start in open_starts:
                        amino_acids = (end - start) // 3
                        if amino_acids >= minimum_amino_acids:
                            left, right = (start, end) if strand == "+" else (length - end, length - start)
                            protein = "".join(CODONS.get(sequence[index:index + 3], "X") for index in range(start, end, 3))
                            results.append({"id": record.identifier, "start": left, "end": right, "strand": strand, "frame": offset + 1 if strand == "+" else -(offset + 1), "amino_acids": amino_acids, "complete": False, "protein": protein})
    results.sort(key=lambda item: (-item["amino_acids"], item["id"], item["start"], item["frame"]))
    return {"input": str(path), "coordinate_system": "zero-based half-open on the input sequence", "minimum_amino_acids": minimum_amino_acids, "orfs": results[:maximum_results], "count": len(results), "truncated": len(results) > maximum_results}


def extract_records(source: Path, destination: Path, identifiers: set[str], start: int | None, end: int | None, force: bool) -> dict:
    if start is not None and start < 0 or end is not None and end < 0 or start is not None and end is not None and end < start:
        raise ValueError("slice coordinates are invalid")
    selected = []
    found = set()
    for record in read_records(source):
        if identifiers and record.identifier not in identifiers:
            continue
        found.add(record.identifier)
        left, right = start or 0, end if end is not None else len(record.sequence)
        selected.append(SequenceRecord(record.identifier, record.description, record.sequence[left:right], record.quality[left:right] if record.quality is not None else None))
    missing = sorted(identifiers - found)
    if missing:
        raise ValueError(f"sequence identifiers not found: {', '.join(missing[:20])}")
    prepare_output(destination, force)
    written = atomic_write(selected, destination, detect_format(source), force=force)
    return {"input": str(source), "output": str(destination), "records": written, "start": start, "end": end}
