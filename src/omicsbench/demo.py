"""Create a small, honest example of the multi-assay evaluation workflow."""

import shutil
from pathlib import Path

from .evaluate import evaluate_suite
from .evidence import create_evidence_crate
from .results import create_result_manifest
from .suite import write_reports


EXAMPLES = (
    ("rnaseq-002", "rnaseq/rnaseq-002/expected/reference-effects.tsv.gz", "rnaseq-002.tsv.gz"),
    ("singlecell-001", "singlecell/singlecell-001/expected/cell-labels.tsv", "singlecell-001.cells.tsv"),
    ("interval-001", "intervals/interval-001/expected/truth-peaks.narrowPeak", "interval-001.narrowPeak"),
    ("sequence-006", "sequences/sequence-006/expected/peptides.tsv", "sequence-006.peptides.tsv"),
)


def create_demo(root: Path, output: Path) -> dict:
    """Copy bundled reference answers and score them as a first-run example.

    A passing self-comparison demonstrates the interfaces, not the performance
    of an external analysis method.
    """
    output = Path(output).absolute()
    if output.exists():
        raise ValueError(f"{output}: output already exists; choose a new directory")
    if not output.parent.is_dir():
        raise ValueError(f"{output.parent}: parent directory does not exist")
    root = Path(root)
    sources = [(identifier, root / "datasets" / relative, name) for identifier, relative, name in EXAMPLES]
    missing = [str(path) for _, path, _ in sources if not path.is_file()]
    variants = root / "datasets/sequences/sequence-005/expected/variants.tsv"
    if not variants.is_file():
        missing.append(str(variants))
    if missing:
        raise ValueError(f"bundled example files are missing: {', '.join(missing)}")

    output.mkdir()
    try:
        results = output / "reference-results"
        results.mkdir()
        mappings = []
        for identifier, source, name in sources:
            shutil.copyfile(source, results / name)
            mappings.append(f"{identifier}={name}")
        rows = [line.split("\t") for line in variants.read_text(encoding="utf-8").splitlines()]
        if not rows or rows[0] != ["position", "ref", "alt"]:
            raise ValueError(f"{variants}: unexpected variant truth columns")
        calls = ["##fileformat=VCFv4.3", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        for row in rows[1:]:
            if len(row) != 3:
                raise ValueError(f"{variants}: malformed variant truth row")
            position, ref, alt = row
            calls.append(f"synthetic_reference\t{position}\t.\t{ref}\t{alt}\t60\tPASS\t.")
        (results / "sequence-005.vcf").write_text("\n".join(calls) + "\n", encoding="utf-8", newline="\n")
        mappings.append("sequence-005=sequence-005.vcf")
        create_result_manifest(root, results, mappings, [])
        report = evaluate_suite(root, results, [item.split("=", 1)[0] for item in mappings])
        if report["summary"]["status"] != "pass":
            raise ValueError("the bundled reference example did not pass evaluation")
        write_reports(report, json_path=output / "report.json", html_path=output / "report.html")
        create_evidence_crate(report, output / "evidence.zip")
        return {
            "operation": "demo",
            "scope": "Bundled reference answers scored against their own truth sets; this does not assess an external method.",
            "summary": report["summary"],
            "results": str(results),
            "report_json": str(output / "report.json"),
            "report_html": str(output / "report.html"),
            "evidence": str(output / "evidence.zip"),
        }
    except Exception:
        shutil.rmtree(output)
        raise
