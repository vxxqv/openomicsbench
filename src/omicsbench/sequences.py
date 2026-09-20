"""Streaming FASTA and FASTQ operations used by the sequence command suite."""
from __future__ import annotations

import gzip
import hashlib
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Iterator, Literal, TextIO


Format = Literal["fasta", "fastq"]
Molecule = Literal["dna", "rna", "protein"]

DNA = frozenset("ACGTRYSWKMBDHVN")
RNA = frozenset("ACGURYSWKMBDHVN")
PROTEIN = frozenset("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
COMPLEMENTS = str.maketrans(
    "ACGTRYMKWSBDHVNUacgtrymkwsbdhvnu",
    "TGCAYRKMWSVHDBNAtgcayrkmwsvhdbna",
)
RNA_COMPLEMENTS = str.maketrans(
    "ACGURYMKWSBDHVNTacgurymkwsbdhvnt",
    "UGCAYRKMWSVHDBNAugcayrkmwsvhdbna",
)
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "R": "[AG]", "Y": "[CTU]", "S": "[GC]", "W": "[ATU]",
    "K": "[GTU]", "M": "[AC]", "B": "[CGTU]", "D": "[AGTU]",
    "H": "[ACTU]", "V": "[ACG]", "N": "[ACGTU]",
}

CODONS = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    identifier: str
    description: str
    sequence: str
    quality: str | None = None

    @property
    def format(self) -> Format:
        return "fastq" if self.quality is not None else "fasta"


def _open_text(path: Path, mode: str = "rt") -> TextIO:
    if "b" in mode:
        raise ValueError("sequence files are opened as text")
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    return opener(path, mode, encoding="utf-8", newline="")


def detect_format(path: Path) -> Format:
    with _open_text(path) as stream:
        for line in stream:
            marker = line.strip()[:1]
            if marker == ">":
                return "fasta"
            if marker == "@":
                return "fastq"
            if marker:
                break
    raise ValueError(f"{path}: expected FASTA or FASTQ content")


def _header(text: str, path: Path, line: int) -> tuple[str, str]:
    value = text[1:].strip()
    if not value:
        raise ValueError(f"{path}:{line}: sequence identifier is empty")
    parts = value.split(maxsplit=1)
    return parts[0], parts[1] if len(parts) == 2 else ""


def read_records(path: Path, file_format: str = "auto") -> Iterator[SequenceRecord]:
    path = Path(path)
    selected: Format = detect_format(path) if file_format == "auto" else file_format  # type: ignore[assignment]
    if selected not in {"fasta", "fastq"}:
        raise ValueError("format must be auto, fasta or fastq")
    with _open_text(path) as stream:
        if selected == "fasta":
            identifier = description = None
            sequence: list[str] = []
            seen: set[str] = set()
            for line_number, raw in enumerate(stream, start=1):
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if identifier is not None:
                        if not sequence:
                            raise ValueError(f"{path}:{line_number}: {identifier} has no sequence")
                        yield SequenceRecord(identifier, description or "", "".join(sequence).upper())
                    identifier, description = _header(line, path, line_number)
                    if identifier in seen:
                        raise ValueError(f"{path}:{line_number}: duplicate sequence identifier {identifier}")
                    seen.add(identifier)
                    sequence = []
                elif identifier is None:
                    raise ValueError(f"{path}:{line_number}: sequence appears before a FASTA header")
                elif any(character.isspace() for character in line):
                    raise ValueError(f"{path}:{line_number}: whitespace inside sequence line")
                else:
                    sequence.append(line)
            if identifier is None:
                raise ValueError(f"{path}: no FASTA records")
            if not sequence:
                raise ValueError(f"{path}: {identifier} has no sequence")
            yield SequenceRecord(identifier, description or "", "".join(sequence).upper())
            return

        seen: set[str] = set()
        line_number = 0
        while True:
            for raw in stream:
                line_number += 1
                if raw.strip():
                    break
            else:
                return
            header_line = raw.rstrip("\r\n")
            if not header_line.startswith("@"):
                raise ValueError(f"{path}:{line_number}: FASTQ record must start with @")
            identifier, description = _header(header_line, path, line_number)
            if identifier in seen:
                raise ValueError(f"{path}:{line_number}: duplicate sequence identifier {identifier}")
            seen.add(identifier)
            sequence_lines: list[str] = []
            for raw in stream:
                line_number += 1
                value = raw.rstrip("\r\n")
                if value.startswith("+"):
                    break
                if not value or any(character.isspace() for character in value):
                    raise ValueError(f"{path}:{line_number}: invalid FASTQ sequence line")
                sequence_lines.append(value)
            else:
                raise ValueError(f"{path}: truncated FASTQ record {identifier}")
            sequence = "".join(sequence_lines).upper()
            if not sequence:
                raise ValueError(f"{path}: {identifier} has no sequence")
            qualities: list[str] = []
            quality_length = 0
            for raw in stream:
                line_number += 1
                value = raw.rstrip("\r\n")
                qualities.append(value)
                quality_length += len(value)
                if quality_length >= len(sequence):
                    break
            quality = "".join(qualities)
            if len(quality) != len(sequence):
                raise ValueError(f"{path}: {identifier} sequence and quality lengths differ")
            if any(ord(character) < 33 or ord(character) > 126 for character in quality):
                raise ValueError(f"{path}: {identifier} has quality characters outside printable Phred+33 range")
            yield SequenceRecord(identifier, description, sequence, quality)


def infer_molecule(sequence: str) -> Molecule:
    letters = set(sequence.upper()) - {"-", "."}
    if not letters:
        raise ValueError("cannot infer molecule from an empty or gap-only sequence")
    if "U" in letters and "T" in letters:
        raise ValueError("sequence contains both T and U; specify or correct the molecule type")
    if letters <= RNA and "U" in letters:
        return "rna"
    if letters <= DNA:
        return "dna"
    if letters <= PROTEIN:
        return "protein"
    invalid = "".join(sorted(letters - (DNA | RNA | PROTEIN)))
    raise ValueError(f"sequence contains unsupported symbols: {invalid}")


def validate_letters(sequence: str, molecule: Molecule, allow_gaps: bool = False) -> None:
    allowed = {"dna": DNA, "rna": RNA, "protein": PROTEIN}[molecule]
    if allow_gaps:
        allowed = allowed | {"-", "."}
    invalid = set(sequence.upper()) - allowed
    if invalid:
        raise ValueError(f"{molecule} sequence contains unsupported symbols: {''.join(sorted(invalid))}")


def resolve_molecule(records: Iterable[SequenceRecord], requested: str) -> tuple[Molecule, list[SequenceRecord]]:
    materialized = list(records)
    if not materialized:
        raise ValueError("sequence file contains no records")
    if requested != "auto":
        molecule: Molecule = requested  # type: ignore[assignment]
    else:
        inferred = {infer_molecule(record.sequence) for record in materialized}
        if len(inferred) != 1:
            raise ValueError(f"records do not share one molecule type: {', '.join(sorted(inferred))}")
        molecule = inferred.pop()
    for record in materialized:
        validate_letters(record.sequence, molecule)
    return molecule, materialized


def reverse_complement(sequence: str, molecule: Molecule = "dna") -> str:
    if molecule not in {"dna", "rna"}:
        raise ValueError("reverse complement requires DNA or RNA")
    validate_letters(sequence, molecule)
    table = RNA_COMPLEMENTS if molecule == "rna" else COMPLEMENTS
    return sequence.translate(table)[::-1]


def translate(sequence: str, frame: int = 1, trim_stop: bool = False) -> str:
    if frame not in {1, 2, 3, -1, -2, -3}:
        raise ValueError("translation frame must be one of 1, 2, 3, -1, -2 or -3")
    molecule = "rna" if "U" in sequence.upper() and "T" not in sequence.upper() else "dna"
    validate_letters(sequence, molecule)
    value = sequence.upper().replace("U", "T")
    if frame < 0:
        value = reverse_complement(value, "dna")
    offset = abs(frame) - 1
    protein = "".join(CODONS.get(value[index:index + 3], "X") for index in range(offset, len(value) - 2, 3))
    return protein.rstrip("*") if trim_stop else protein


def phred_scores(quality: str) -> list[int]:
    return [ord(character) - 33 for character in quality]


def gc_fraction(sequence: str) -> float:
    bases = Counter(sequence.upper())
    denominator = sum(bases[symbol] for symbol in "ACGTU")
    return (bases["G"] + bases["C"]) / denominator if denominator else 0.0


def ambiguity_fraction(sequence: str, molecule: Molecule) -> float:
    canonical = set("ACGT") if molecule == "dna" else set("ACGU") if molecule == "rna" else set("ACDEFGHIKLMNPQRSTVWY")
    return sum(character not in canonical for character in sequence.upper()) / len(sequence)


def n50(lengths: list[int]) -> tuple[int, int]:
    threshold = sum(lengths) / 2
    cumulative = 0
    for index, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= threshold:
            return length, index
    return 0, 0


def summarize(path: Path, molecule: str = "auto", file_format: str = "auto") -> dict:
    selected_format = detect_format(path) if file_format == "auto" else file_format
    selected_molecule, records = resolve_molecule(read_records(path, selected_format), molecule)
    lengths = [len(record.sequence) for record in records]
    composition = Counter(character for record in records for character in record.sequence)
    total = sum(lengths)
    n50_value, l50_value = n50(lengths)
    digests = Counter(hashlib.sha256(record.sequence.encode("ascii")).digest() for record in records)
    result = {
        "path": str(path),
        "format": selected_format,
        "molecule": selected_molecule,
        "records": len(records),
        "letters": total,
        "length": {
            "minimum": min(lengths), "maximum": max(lengths), "mean": total / len(lengths),
            "median": median(lengths), "n50": n50_value, "l50": l50_value,
        },
        "composition": dict(sorted(composition.items())),
        "gc_fraction": gc_fraction("".join(record.sequence for record in records)) if selected_molecule in {"dna", "rna"} else None,
        "ambiguity_fraction": sum(ambiguity_fraction(record.sequence, selected_molecule) * len(record.sequence) for record in records) / total,
        "duplicate_sequences": sum(count - 1 for count in digests.values()),
    }
    if selected_format == "fastq":
        scores = [score for record in records for score in phred_scores(record.quality or "")]
        result["quality"] = {
            "encoding": "phred33", "minimum": min(scores), "maximum": max(scores),
            "mean": sum(scores) / len(scores),
            "q20_fraction": sum(score >= 20 for score in scores) / len(scores),
            "q30_fraction": sum(score >= 30 for score in scores) / len(scores),
        }
    return result


def write_records(records: Iterable[SequenceRecord], path: Path, file_format: str = "auto", wrap: int = 80) -> int:
    if wrap < 0:
        raise ValueError("wrap must be zero or greater")
    count = 0
    selected = None if file_format == "auto" else file_format
    with _open_text(Path(path), "wt") as stream:
        for record in records:
            output_format = selected or record.format
            if output_format == "fastq" and record.quality is None:
                raise ValueError(f"{record.identifier}: FASTQ output requires qualities")
            header = record.identifier + (f" {record.description}" if record.description else "")
            if output_format == "fasta":
                stream.write(f">{header}\n")
                width = wrap or len(record.sequence)
                for index in range(0, len(record.sequence), width):
                    stream.write(record.sequence[index:index + width] + "\n")
            elif output_format == "fastq":
                stream.write(f"@{header}\n{record.sequence}\n+\n{record.quality}\n")
            else:
                raise ValueError("output format must be auto, fasta or fastq")
            count += 1
    return count


def trim_record(record: SequenceRecord, left: int = 0, right: int = 0, quality: int | None = None, adapter: str | None = None) -> SequenceRecord:
    if left < 0 or right < 0:
        raise ValueError("trim lengths must be zero or greater")
    start, end = min(left, len(record.sequence)), max(min(len(record.sequence) - right, len(record.sequence)), 0)
    if end < start:
        end = start
    sequence = record.sequence[start:end]
    qualities = record.quality[start:end] if record.quality is not None else None
    if adapter:
        position = sequence.find(adapter.upper())
        if position >= 0:
            sequence = sequence[:position]
            qualities = qualities[:position] if qualities is not None else None
    if quality is not None:
        if qualities is None:
            raise ValueError(f"{record.identifier}: quality trimming requires FASTQ input")
        scores = phred_scores(qualities)
        while scores and scores[-1] < quality:
            scores.pop(); sequence = sequence[:-1]; qualities = qualities[:-1]
        while scores and scores[0] < quality:
            scores.pop(0); sequence = sequence[1:]; qualities = qualities[1:]
    return SequenceRecord(record.identifier, record.description, sequence, qualities)


def motif_pattern(motif: str) -> re.Pattern[str]:
    if not motif:
        raise ValueError("motif cannot be empty")
    try:
        body = "".join(IUPAC[character] for character in motif.upper())
    except KeyError as exc:
        raise ValueError(f"motif contains unsupported IUPAC symbol {exc.args[0]}") from None
    return re.compile(f"(?=({body}))", re.IGNORECASE)


def find_motifs(records: Iterable[SequenceRecord], motif: str, both_strands: bool = False) -> dict:
    pattern = motif_pattern(motif)
    reverse = reverse_complement(motif.upper().replace("U", "T"), "dna") if both_strands else None
    reverse_pattern = motif_pattern(reverse) if reverse and reverse != motif.upper() else None
    hits = []
    for record in records:
        for match in pattern.finditer(record.sequence):
            hits.append({"id": record.identifier, "start": match.start(), "end": match.start() + len(motif), "strand": "+"})
        if reverse_pattern:
            for match in reverse_pattern.finditer(record.sequence):
                hits.append({"id": record.identifier, "start": match.start(), "end": match.start() + len(motif), "strand": "-"})
    return {"motif": motif.upper(), "coordinate_system": "zero-based half-open", "hits": hits, "count": len(hits)}


def iter_kmers(sequence: str, k: int, canonical: bool = False) -> Iterator[str]:
    if k < 1:
        raise ValueError("k must be at least 1")
    for index in range(len(sequence) - k + 1):
        value = sequence[index:index + k]
        if set(value) <= set("ACGT"):
            if canonical:
                reverse = reverse_complement(value, "dna")
                value = min(value, reverse)
            yield value


def count_kmers(records: Iterable[SequenceRecord], k: int, canonical: bool = False, maximum_distinct: int = 1_000_000) -> Counter[str]:
    counts: Counter[str] = Counter()
    for record in records:
        counts.update(iter_kmers(record.sequence.replace("U", "T"), k, canonical))
        if len(counts) > maximum_distinct:
            raise ValueError(f"more than {maximum_distinct} distinct k-mers; increase --max-distinct deliberately")
    return counts


def sketch(records: Iterable[SequenceRecord], k: int = 21, size: int = 1000, canonical: bool = True) -> dict:
    if size < 1:
        raise ValueError("sketch size must be at least 1")
    hashes = {
        int.from_bytes(hashlib.blake2b(value.encode("ascii"), digest_size=8).digest(), "big")
        for record in records for value in iter_kmers(record.sequence.replace("U", "T"), k, canonical)
    }
    selected = sorted(hashes)[:size]
    return {"algorithm": "bottom-k-blake2b-64-v1", "k": k, "canonical": canonical, "requested_size": size, "observed_kmers": len(hashes), "hashes": [f"{value:016x}" for value in selected]}


def compare_sketches(left: dict, right: dict) -> dict:
    for key in ("algorithm", "k", "canonical"):
        if left.get(key) != right.get(key):
            raise ValueError(f"sketches differ in {key}")
    a, b = set(left["hashes"]), set(right["hashes"])
    union = a | b
    intersection = a & b
    jaccard = len(intersection) / len(union) if union else 1.0
    containment = len(intersection) / min(len(a), len(b)) if a and b else 0.0
    k = int(left["k"])
    mash = -math.log(2 * jaccard / (1 + jaccard)) / k if jaccard > 0 else None
    return {"shared": len(intersection), "left": len(a), "right": len(b), "jaccard": jaccard, "containment": containment, "mash_distance": mash}
