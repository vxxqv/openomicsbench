"""Small commands with structured output and concise failure messages."""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from . import __version__
from .registry import default_root, registry, lookup
from .cache import get, verify_cache
from .compare import compare
from .sequence_cli import add_sequence_parser, run_sequence_command
from .assays import build_plan, get_profile, list_profiles
from .suite import compare_suite, validate_suite, write_reports
from .matrix import compare_matrix, parse_method
from .regression import compare_reports
from .bundle import create_bundle, verify_bundle
from .validate import validate
from .variants import compare_variants
from .peptides import compare_peptides
from .evaluate import evaluate_matrix, evaluate_suite
from .evidence import create_evidence_crate, verify_evidence_crate
from .policy import evaluate_policy
from .trends import parse_run, trend_reports
from .results import create_result_manifest
from .singlecell import compare_labels, summarize_matrix
from .intervals import compare_intervals, summarize_intervals
from .taxonomy import classify_exact, compare_taxonomy
from .demo import create_demo
from .workspace import prepare_workspace, verify_workspace, grade_workspace
from .discovery import discover_workspace_results
from .client import Collection
from .rnaseq_pack import normalize_effects

ASSAYS = ["bulk_rna_seq", "single_cell_rna_seq", "chromatin_accessibility", "sequence_dna", "sequence_rna", "sequence_protein", "short_read_dna", "whole_genome_dna_seq", "shotgun_metagenomics"]

def main():
    p = argparse.ArgumentParser(
        description="Find, verify and process compact omics benchmarks.",
        epilog=("Example: omicsbench info rnaseq-002\n"
                "If useful, please star the repo: https://github.com/vxxqv/openomicsbench"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--root",type=Path,default=default_root(),help="Collection directory. Defaults to OMICSBENCH_ROOT, the current checkout or the installed collection.")
    p.add_argument("--cache",type=Path,default=Path(os.environ.get("OMICSBENCH_CACHE",".omicsbench-cache")),help="Local cache directory.")
    p.add_argument("--debug",action="store_true",help="Show a traceback when a command fails.")
    p.add_argument("--version",action="version",version=__version__)
    sub = p.add_subparsers(dest="command",required=True)
    examples = {"list":"list --assay bulk_rna_seq","info":"info fixture-001","get":"get fixture-001 --size nano","validate":"validate fixture-001","compare":"compare rnaseq-002 results.csv","provenance":"provenance fixture-001","verify-cache":"verify-cache","doctor":"doctor"}
    for name,example in examples.items():
        command = sub.add_parser(name,help={"list":"List local manifests.","info":"Show files, rights and limitations.","get":"Cache an exact tier after checking hashes.","validate":"Check declared files and recompute metrics.","compare":"Compare differential-expression effects with the full reference.","provenance":"Show source and transformations.","verify-cache":"Verify completed cache entries.","doctor":"Check local runtime and optional tools."}[name],epilog=f"Example: omicsbench {example}")
        if name in {"info","get","validate","compare","provenance"}:
            command.add_argument("id",help="Stable dataset identifier.")
        if name == "get":
            command.add_argument("--size",choices=["nano","pocket","expected"],default="nano")
        if name == "list":
            command.add_argument("--assay", choices=ASSAYS)
        if name == "compare":
            command.add_argument("results",type=Path,help="CSV or TSV file with gene_id and log2_fold_change columns.")
            command.add_argument("--detail-limit",type=int,default=20,help="Maximum missing and unexpected gene examples to return.")
    add_sequence_parser(sub)
    rnaseq = sub.add_parser("rnaseq", help="Export and check biological bulk RNA-seq counts and sample metadata.")
    rnaseq_sub = rnaseq.add_subparsers(dest="rnaseq_command", required=True)
    rnaseq_export = rnaseq_sub.add_parser("export", help="Write counts.tsv, samples.tsv and export.json to a new directory.")
    rnaseq_export.add_argument("id", help="Validated biological bulk RNA-seq benchmark ID.")
    rnaseq_export.add_argument("output", type=Path, help="New export directory; existing directories are preserved.")
    rnaseq_export.add_argument("--tier", choices=["pocket", "full"], default="pocket")
    rnaseq_check = rnaseq_sub.add_parser("check-counts", help="Compare integer counts by gene and sample ID against a bundled tier.")
    rnaseq_check.add_argument("id", help="Validated biological bulk RNA-seq benchmark ID.")
    rnaseq_check.add_argument("counts", type=Path, help="Feature-by-sample TSV count matrix.")
    rnaseq_check.add_argument("--tier", choices=["pocket", "full"], default="pocket")
    rnaseq_check.add_argument("--detail-limit", type=int, default=20)
    rnaseq_verify = rnaseq_sub.add_parser("verify-export", help="Check an exported directory's hashes, counts and sample metadata.")
    rnaseq_verify.add_argument("id", help="Validated biological bulk RNA-seq benchmark ID.")
    rnaseq_verify.add_argument("directory", type=Path, help="Directory created by rnaseq export.")
    rnaseq_pack = rnaseq_sub.add_parser("pack", help="Export checked count, sample and contrast files for one or more studies.")
    rnaseq_pack.add_argument("output", type=Path, help="New pack directory.")
    rnaseq_pack.add_argument("--id", action="append", help="Biological RNA-seq benchmark ID; repeat to select several. Default: all 12.")
    rnaseq_pack.add_argument("--tier", choices=["pocket", "full"], default="pocket")
    rnaseq_pack_verify = rnaseq_sub.add_parser("verify-pack", help="Recheck all pack files, IDs and contrasts.")
    rnaseq_pack_verify.add_argument("directory", type=Path)
    rnaseq_normalize = rnaseq_sub.add_parser("normalize-effects", help="Convert native differential-expression tables to benchmark columns.")
    rnaseq_normalize.add_argument("id", help="Biological bulk RNA-seq benchmark ID.")
    rnaseq_normalize.add_argument("source", type=Path, help="CSV or TSV table, optionally gzip-compressed.")
    rnaseq_normalize.add_argument("output", type=Path, help="New canonical TSV result file.")
    rnaseq_normalize.add_argument("--gene-column", help="Exact source header for gene IDs when auto-detection is ambiguous.")
    rnaseq_normalize.add_argument("--effect-column", help="Exact source header for log2 effects when auto-detection is ambiguous.")
    rnaseq_normalize.add_argument("--reverse-sign", action="store_true", help="Explicitly reverse the source effect direction.")
    rnaseq_design = rnaseq_sub.add_parser("check-design", help="Check matrix, sample and contrast ID alignment and design rank.")
    rnaseq_design.add_argument("counts", type=Path)
    rnaseq_design.add_argument("samples", type=Path)
    rnaseq_design.add_argument("contrasts", type=Path)
    rnaseq_design.add_argument("--id", help="Also compare sample and contrast design with a biological benchmark.")
    rnaseq_score = rnaseq_sub.add_parser("score-effects", help="Score a native DESeq2, edgeR or limma table without renaming it.")
    rnaseq_score.add_argument("id", help="Biological bulk RNA-seq benchmark ID.")
    rnaseq_score.add_argument("source", type=Path, help="CSV or TSV table, optionally gzip-compressed.")
    rnaseq_score.add_argument("--gene-column")
    rnaseq_score.add_argument("--effect-column")
    rnaseq_score.add_argument("--reverse-sign", action="store_true")
    rnaseq_score.add_argument("--detail-limit", type=int, default=20)
    demo = sub.add_parser("demo", help="Run a six-assay example using bundled reference answers.")
    demo.add_argument("output", type=Path, help="New directory for example results, reports and evidence.")
    variant = sub.add_parser("variant", help="Compare VCF calls with a bundled small-variant truth set.")
    variant_sub = variant.add_subparsers(dest="variant_command", required=True)
    variant_compare = variant_sub.add_parser("compare", help="Score exact SNV alleles in a VCF or VCF.GZ file.")
    variant_compare.add_argument("id", help="Benchmark ID with a declared variant truth set.")
    variant_compare.add_argument("results", type=Path, help="VCF or VCF.GZ call set.")
    variant_compare.add_argument("--detail-limit", type=int, default=20, help="Maximum false-positive and false-negative examples to return.")
    peptide = sub.add_parser("peptide", help="Compare protein-peptide identifications with a bundled digest truth set.")
    peptide_sub = peptide.add_subparsers(dest="peptide_command", required=True)
    peptide_compare = peptide_sub.add_parser("compare", help="Score exact protein-peptide pairs from a TSV table.")
    peptide_compare.add_argument("id", help="Benchmark ID with a declared peptide truth set.")
    peptide_compare.add_argument("results", type=Path, help="TSV or TSV.GZ table with protein_id and peptide columns.")
    peptide_compare.add_argument("--detail-limit", type=int, default=20, help="Maximum false-positive and false-negative examples to return.")
    cell = sub.add_parser("cell", help="Inspect a feature-barcode matrix or score single-cell cluster labels.")
    cell_sub = cell.add_subparsers(dest="cell_command", required=True)
    cell_qc = cell_sub.add_parser("qc", help="Summarize a 10x-style Matrix Market feature-barcode matrix.")
    cell_qc.add_argument("matrix", type=Path)
    cell_qc.add_argument("barcodes", type=Path)
    cell_qc.add_argument("features", type=Path)
    cell_qc.add_argument("--mitochondrial-prefix", default="MT-")
    cell_qc.add_argument("--feature-columns", type=int, choices=[2, 3])
    cell_compare = cell_sub.add_parser("compare", help="Score cell cluster labels against a bundled truth set.")
    cell_compare.add_argument("id", help="Single-cell benchmark ID.")
    cell_compare.add_argument("results", type=Path, help="TSV, CSV, or gzip-compressed table with cell_id and label columns.")
    cell_compare.add_argument("--detail-limit", type=int, default=20)
    interval = sub.add_parser("interval", help="Inspect BED-family peaks or score them against bundled truth.")
    interval_sub = interval.add_subparsers(dest="interval_command", required=True)
    interval_qc = interval_sub.add_parser("qc", help="Validate and summarize a BED or narrowPeak file.")
    interval_qc.add_argument("peaks", type=Path)
    interval_qc.add_argument("--chrom-sizes", type=Path)
    interval_compare = interval_sub.add_parser("compare", help="Score peaks using reciprocal interval and base-pair overlap.")
    interval_compare.add_argument("id", help="Benchmark ID with declared genomic interval truth.")
    interval_compare.add_argument("results", type=Path, help="BED, BED.GZ, narrowPeak or narrowPeak.GZ call set.")
    taxon = sub.add_parser("taxon", help="Classify against mock references or score read-level taxonomy assignments.")
    taxon_sub = taxon.add_subparsers(dest="taxon_command", required=True)
    taxon_classify = taxon_sub.add_parser("classify", help="Exact-match paired reads against project-local mock references.")
    taxon_classify.add_argument("reference", type=Path, help="FASTA with mock_taxon_<code> sequence IDs.")
    taxon_classify.add_argument("read1", type=Path, help="R1 FASTQ with /1 read IDs.")
    taxon_classify.add_argument("read2", type=Path, help="R2 FASTQ with /2 read IDs.")
    taxon_classify.add_argument("output", type=Path, help="New canonical assignment TSV.")
    taxon_classify.add_argument("--force", action="store_true")
    taxon_compare = taxon_sub.add_parser("compare", help="Score TSV or Kraken 2 per-read assignments against bundled truth.")
    taxon_compare.add_argument("id", help="Benchmark ID with a read-level taxonomy truth set.")
    taxon_compare.add_argument("results", type=Path, help="Assignment TSV or standard five-column Kraken output, optionally gzip-compressed.")
    taxon_compare.add_argument("--detail-limit", type=int, default=20)
    assay = sub.add_parser("assay", help="Inspect sequencing assay profiles and build local command plans.")
    assay_sub = assay.add_subparsers(dest="assay_command", required=True)
    assay_sub.add_parser("list", help="List supported assay profiles.")
    assay_info = assay_sub.add_parser("info", help="Show the inputs, checks and scope for one profile.")
    assay_info.add_argument("profile")
    assay_plan = assay_sub.add_parser("plan", help="Build an executable local validation and preprocessing plan.")
    assay_plan.add_argument("profile")
    assay_plan.add_argument("input", type=Path, nargs="+")
    assay_plan.add_argument("--output-directory", type=Path, default=Path("omicsbench-output"))
    assay_plan.add_argument("--adapter", action="append", default=[])
    suite = sub.add_parser("suite", help="Run collection-wide checks and write CI reports.")
    suite_sub = suite.add_subparsers(dest="suite_command", required=True)
    suite_validate = suite_sub.add_parser("validate", help="Validate several bundled benchmarks in one run.")
    suite_validate.add_argument("--assay", choices=ASSAYS)
    suite_compare = suite_sub.add_parser("compare", help="Compare a directory of differential-expression results.")
    suite_compare.add_argument("results", type=Path, help="Directory containing files named <dataset-id>.csv or .tsv, optionally gzip-compressed.")
    suite_compare.add_argument("--detail-limit", type=int, default=20, help="Maximum missing and unexpected gene examples per benchmark.")
    suite_evaluate = suite_sub.add_parser("evaluate", help="Evaluate one result directory across RNA, DNA and protein benchmarks.")
    suite_evaluate.add_argument("results", type=Path, help="Directory containing benchmark-named RNA tables, VCF files and peptide tables.")
    suite_evaluate.add_argument("--detail-limit", type=int, default=20, help="Maximum mismatch examples per benchmark.")
    suite_index = suite_sub.add_parser("index", help="Write a hash-locked manifest for pipeline result files.")
    suite_index.add_argument("results", type=Path, help="Directory containing the result files and the new manifest.")
    suite_index.add_argument("--result", action="append", required=True, help="Benchmark and relative file path as ID=PATH. Repeat for every result.")
    suite_index.add_argument("--column", action="append", default=[], help="Source column mapping as ID:FIELD=HEADER. Repeat when native headers differ.")
    suite_matrix = suite_sub.add_parser("matrix", help="Compare several analysis methods across the same benchmarks.")
    suite_matrix.add_argument("--method", action="append", required=True, help="Method and result directory as NAME=PATH. Repeat for every method.")
    suite_matrix.add_argument("--detail-limit", type=int, default=20, help="Maximum missing and unexpected gene examples per benchmark.")
    suite_evaluate_matrix = suite_sub.add_parser("evaluate-matrix", help="Compare several methods across RNA, DNA and protein benchmarks.")
    suite_evaluate_matrix.add_argument("--method", action="append", required=True, help="Method and result directory as NAME=PATH. Repeat for every method.")
    suite_evaluate_matrix.add_argument("--detail-limit", type=int, default=20, help="Maximum mismatch examples per benchmark.")
    suite_regress = suite_sub.add_parser("regress", help="Fail when a candidate suite report regresses from a baseline report.")
    suite_regress.add_argument("baseline", type=Path, help="Previously accepted JSON suite report.")
    suite_regress.add_argument("candidate", type=Path, help="Candidate JSON suite report to check.")
    suite_regress.add_argument("--absolute-tolerance", type=float, default=0.0, help="Largest permitted absolute decrease in a tracked metric.")
    suite_regress.add_argument("--allow-missing", action="store_true", help="Do not fail when a baseline case is absent from the candidate.")
    suite_gate = suite_sub.add_parser("gate", help="Apply a versioned quality and resource policy to an evaluation report.")
    suite_gate.add_argument("report", type=Path, help="JSON report produced by suite evaluate.")
    suite_gate.add_argument("policy", type=Path, help="JSON gate policy with coverage, metric and resource requirements.")
    suite_trend = suite_sub.add_parser("trend", help="Combine fingerprinted evaluation reports into a run history.")
    suite_trend.add_argument("--run", action="append", required=True, help="Run label and evaluation report as NAME=PATH. Repeat in chronological order.")
    suite_prepare = suite_sub.add_parser("prepare", help="Copy selected benchmark inputs into a checked workspace for your own results.")
    suite_prepare.add_argument("output", type=Path, help="New workspace directory.")
    suite_prepare.add_argument("--id", action="append", required=True, help="Validated benchmark to include. Repeat for several cases.")
    suite_verify = suite_sub.add_parser("verify-workspace", help="Check every staged input against the bundled benchmark inventory.")
    suite_verify.add_argument("workspace", type=Path)
    suite_grade = suite_sub.add_parser("grade", help="Score the workspace results and write reports and evidence.")
    suite_grade.add_argument("workspace", type=Path)
    suite_grade.add_argument("--force", action="store_true", help="Replace an earlier grading report after results change.")
    suite_discover = suite_sub.add_parser("discover", help="Find workspace results and check their column contracts before grading.")
    suite_discover.add_argument("workspace", type=Path)
    suite_discover.add_argument("--result", action="append", default=[], help="Use an explicit ID=RELATIVE_PATH for a result with a different name.")
    suite_discover.add_argument("--column", action="append", default=[], help="Resolve a header as ID:FIELD=HEADER.")
    suite_discover.add_argument("--write", action="store_true", help="Write a hash-locked results.json when every case is ready.")
    suite_discover.add_argument("--force", action="store_true", help="Replace an existing results.json with --write.")
    for command in (suite_validate, suite_compare, suite_evaluate, suite_matrix, suite_evaluate_matrix):
        command.add_argument("--id", action="append", default=[], help="Benchmark ID to include. Repeat to select several; omit to run all eligible benchmarks.")
    for command in (suite_validate, suite_compare, suite_evaluate, suite_index, suite_matrix, suite_evaluate_matrix, suite_regress, suite_gate, suite_trend):
        command.add_argument("--json", type=Path, help="Write the complete report as JSON.")
        command.add_argument("--markdown", type=Path, help="Write a concise Markdown report.")
        command.add_argument("--junit", type=Path, help="Write a JUnit XML report for CI systems.")
        command.add_argument("--csv", type=Path, help="Write a flat CSV report for analysis or plotting.")
        command.add_argument("--html", type=Path, help="Write a self-contained HTML report.")
        command.add_argument("--force", action="store_true", help="Replace existing report files.")
    bundle = sub.add_parser("bundle", help="Create or verify a portable benchmark bundle.")
    bundle_sub = bundle.add_subparsers(dest="bundle_command", required=True)
    bundle_create = bundle_sub.add_parser("create", help="Create a deterministic ZIP containing selected benchmarks.")
    bundle_create.add_argument("output", type=Path)
    bundle_create.add_argument("--id", action="append", default=[], help="Benchmark ID to include. Repeat to select several.")
    bundle_create.add_argument("--assay", choices=ASSAYS)
    bundle_create.add_argument("--force", action="store_true", help="Replace an existing bundle.")
    bundle_verify = bundle_sub.add_parser("verify", help="Verify a bundle inventory, hashes and dataset manifests.")
    bundle_verify.add_argument("archive", type=Path)
    evidence = sub.add_parser("evidence", help="Verify a portable evaluation evidence crate.")
    evidence_sub = evidence.add_subparsers(dest="evidence_command", required=True)
    evidence_verify = evidence_sub.add_parser("verify", help="Verify checksums, inventory, report and RO-Crate metadata.")
    evidence_verify.add_argument("archive", type=Path)
    for command in (suite_evaluate, suite_evaluate_matrix):
        command.add_argument("--evidence", type=Path, help="Write a deterministic ZIP with the reports, submitted results and RO-Crate metadata.")
    args = p.parse_args()
    try:
        if args.command == "demo":
            out = create_demo(args.root, args.output)
        elif args.command == "rnaseq":
            collection = Collection(args.root)
            if args.rnaseq_command == "pack":
                out = collection.pack_rnaseq(args.output, args.id, args.tier)
            elif args.rnaseq_command == "verify-pack":
                out = collection.verify_rnaseq_pack(args.directory)
            elif args.rnaseq_command == "normalize-effects":
                out = normalize_effects(collection.get(args.id), args.source, args.output,
                                        gene_column=args.gene_column, effect_column=args.effect_column,
                                        reverse_sign=args.reverse_sign)
            elif args.rnaseq_command == "check-design":
                out = collection.check_rnaseq_design(args.counts, args.samples, args.contrasts, args.id)
            elif args.rnaseq_command == "score-effects":
                out = collection.get(args.id).score_native_effects(
                    args.source, gene_column=args.gene_column, effect_column=args.effect_column,
                    reverse_sign=args.reverse_sign, detail_limit=args.detail_limit)
            elif args.rnaseq_command == "export":
                benchmark = collection.get(args.id)
                out = benchmark.export_rnaseq(args.output, args.tier)
            elif args.rnaseq_command == "check-counts":
                benchmark = collection.get(args.id)
                out = benchmark.compare_counts(args.counts, args.tier, args.detail_limit)
            else:
                benchmark = collection.get(args.id)
                out = benchmark.verify_rnaseq_export(args.directory)
        elif args.command == "seq":
            out = run_sequence_command(args)
        elif args.command == "variant":
            model, folder = lookup(args.root, args.id)
            out = compare_variants(model, folder, args.results, args.detail_limit)
        elif args.command == "peptide":
            model, folder = lookup(args.root, args.id)
            out = compare_peptides(model, folder, args.results, args.detail_limit)
        elif args.command == "cell":
            if args.cell_command == "qc":
                out = summarize_matrix(args.matrix, args.barcodes, args.features, args.mitochondrial_prefix, args.feature_columns)
            else:
                model, folder = lookup(args.root, args.id)
                out = compare_labels(model, folder, args.results, args.detail_limit)
        elif args.command == "interval":
            if args.interval_command == "qc":
                out = summarize_intervals(args.peaks, args.chrom_sizes)
            else:
                model, folder = lookup(args.root, args.id)
                out = compare_intervals(model, folder, args.results)
        elif args.command == "taxon":
            if args.taxon_command == "classify":
                out = classify_exact(args.reference, args.read1, args.read2, args.output, args.force)
            else:
                model, folder = lookup(args.root, args.id)
                out = compare_taxonomy(model, folder, args.results, args.detail_limit)
        elif args.command == "evidence":
            out = verify_evidence_crate(args.archive)
        elif args.command == "bundle":
            if args.bundle_command == "create":
                out = create_bundle(args.root, args.output, args.id, args.assay, args.force)
            else:
                out = verify_bundle(args.archive)
        elif args.command == "assay":
            if args.assay_command == "list":
                out = list_profiles()
            elif args.assay_command == "info":
                out = get_profile(args.profile)
            else:
                out = build_plan(args.profile, args.input, args.output_directory, args.adapter)
        elif args.command == "suite" and args.suite_command == "prepare":
            out = prepare_workspace(args.root, args.output, args.id)
        elif args.command == "suite" and args.suite_command == "verify-workspace":
            out = verify_workspace(args.root, args.workspace)
        elif args.command == "suite" and args.suite_command == "grade":
            out = grade_workspace(args.root, args.workspace, args.force)
        elif args.command == "suite" and args.suite_command == "discover":
            out = discover_workspace_results(args.root, args.workspace, args.result, args.column, args.write, args.force)
        elif args.command == "suite":
            if args.suite_command == "index":
                out = create_result_manifest(args.root, args.results, args.result, args.column, args.force)
            elif args.suite_command == "validate":
                out = validate_suite(args.root, args.assay, args.id)
            elif args.suite_command == "compare":
                out = compare_suite(args.root, args.results, args.id, args.detail_limit)
            elif args.suite_command == "evaluate":
                out = evaluate_suite(args.root, args.results, args.id, args.detail_limit)
            elif args.suite_command == "matrix":
                out = compare_matrix(args.root, [parse_method(value) for value in args.method], args.id, args.detail_limit)
            elif args.suite_command == "evaluate-matrix":
                out = evaluate_matrix(args.root, [parse_method(value) for value in args.method], args.id, args.detail_limit)
            elif args.suite_command == "regress":
                out = compare_reports(args.baseline, args.candidate, args.absolute_tolerance, args.allow_missing)
            elif args.suite_command == "gate":
                out = evaluate_policy(args.report, args.policy)
            else:
                out = trend_reports([parse_run(value) for value in args.run])
            reports = write_reports(
                out,
                json_path=args.json,
                markdown_path=args.markdown,
                junit_path=args.junit,
                force=args.force,
                csv_path=args.csv,
                html_path=args.html,
            )
            if getattr(args, "evidence", None) is not None:
                reports["evidence"] = create_evidence_crate(out, args.evidence, args.force)
            out["reports"] = reports
        elif args.command == "list":
            out = [{"id":m.id,"title":m.title,"kind":m.kind,"status":m.status,"rights":m.rights.status,"tiers":sorted({f.tier for f in m.files})} for m,_ in registry(args.root).values() if not args.assay or m.assay==args.assay]
        elif args.command == "doctor":
            out = {"python":sys.version.split()[0],"version":__version__,"root_exists":(args.root/"datasets").is_dir(),"cache":str(args.cache.resolve()),"optional_tools":{x:shutil.which(x) for x in ["Rscript","snakemake","docker"]},"network":"not probed; discovery works offline"}
        elif args.command == "verify-cache":
            out = verify_cache(args.root,args.cache)
        else:
            model,folder = lookup(args.root,args.id)
            if args.command == "info":
                out = model.model_dump(mode="json")
            elif args.command == "provenance":
                out = {"id":model.id,"source":model.source.model_dump(mode="json"),"rights":model.rights.model_dump(mode="json"),"derivation":model.derivation.model_dump(mode="json")}
            elif args.command == "get":
                out = {"cache_directory":str(get(args.root,args.cache,args.id,args.size).resolve())}
            elif args.command == "compare":
                out = compare(model,folder,args.results,args.detail_limit)
            else:
                out = validate(model,folder)
        print(json.dumps(out,indent=2,allow_nan=False))
        if args.command == "compare" and out["status"] == "fail":
            sys.exit(2)
        if args.command == "variant" and out["status"] == "fail":
            sys.exit(2)
        if args.command == "peptide" and out["status"] == "fail":
            sys.exit(2)
        if args.command == "cell" and args.cell_command == "compare" and out["status"] == "fail":
            sys.exit(2)
        if args.command == "interval" and args.interval_command == "compare" and out["status"] == "fail":
            sys.exit(2)
        if args.command == "taxon" and args.taxon_command == "compare" and out["status"] == "fail":
            sys.exit(2)
        if args.command == "suite" and out["summary"]["status"] == "fail":
            sys.exit(2)
        if args.command == "rnaseq" and (out.get("status") == "fail" or out.get("summary", {}).get("status") == "fail"):
            sys.exit(2)
    except (ValueError,OSError,KeyError) as exc:
        if args.debug:
            raise
        print(f"omicsbench: {exc}",file=sys.stderr)
        sys.exit(1)
