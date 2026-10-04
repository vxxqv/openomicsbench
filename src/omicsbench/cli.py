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
from .evaluate import evaluate_matrix, evaluate_suite

def main():
    p = argparse.ArgumentParser(description="Find, verify and process compact omics benchmarks.",epilog="Example: omicsbench info rnaseq-002")
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
            command.add_argument("--assay",choices=["bulk_rna_seq","sequence_dna","sequence_rna","sequence_protein","short_read_dna","whole_genome_dna_seq"])
        if name == "compare":
            command.add_argument("results",type=Path,help="CSV or TSV file with gene_id and log2_fold_change columns.")
            command.add_argument("--detail-limit",type=int,default=20,help="Maximum missing and unexpected gene examples to return.")
    add_sequence_parser(sub)
    variant = sub.add_parser("variant", help="Compare VCF calls with a bundled small-variant truth set.")
    variant_sub = variant.add_subparsers(dest="variant_command", required=True)
    variant_compare = variant_sub.add_parser("compare", help="Score exact SNV alleles in a VCF or VCF.GZ file.")
    variant_compare.add_argument("id", help="Benchmark ID with a declared variant truth set.")
    variant_compare.add_argument("results", type=Path, help="VCF or VCF.GZ call set.")
    variant_compare.add_argument("--detail-limit", type=int, default=20, help="Maximum false-positive and false-negative examples to return.")
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
    suite_validate.add_argument("--assay", choices=["bulk_rna_seq", "sequence_dna", "sequence_rna", "sequence_protein", "short_read_dna", "whole_genome_dna_seq"])
    suite_compare = suite_sub.add_parser("compare", help="Compare a directory of differential-expression results.")
    suite_compare.add_argument("results", type=Path, help="Directory containing files named <dataset-id>.csv or .tsv, optionally gzip-compressed.")
    suite_compare.add_argument("--detail-limit", type=int, default=20, help="Maximum missing and unexpected gene examples per benchmark.")
    suite_evaluate = suite_sub.add_parser("evaluate", help="Evaluate one result directory across RNA and DNA benchmarks.")
    suite_evaluate.add_argument("results", type=Path, help="Directory containing benchmark-named RNA tables and VCF files.")
    suite_evaluate.add_argument("--detail-limit", type=int, default=20, help="Maximum mismatch examples per benchmark.")
    suite_matrix = suite_sub.add_parser("matrix", help="Compare several analysis methods across the same benchmarks.")
    suite_matrix.add_argument("--method", action="append", required=True, help="Method and result directory as NAME=PATH. Repeat for every method.")
    suite_matrix.add_argument("--detail-limit", type=int, default=20, help="Maximum missing and unexpected gene examples per benchmark.")
    suite_evaluate_matrix = suite_sub.add_parser("evaluate-matrix", help="Compare several methods across RNA and DNA benchmarks.")
    suite_evaluate_matrix.add_argument("--method", action="append", required=True, help="Method and result directory as NAME=PATH. Repeat for every method.")
    suite_evaluate_matrix.add_argument("--detail-limit", type=int, default=20, help="Maximum mismatch examples per benchmark.")
    suite_regress = suite_sub.add_parser("regress", help="Fail when a candidate suite report regresses from a baseline report.")
    suite_regress.add_argument("baseline", type=Path, help="Previously accepted JSON suite report.")
    suite_regress.add_argument("candidate", type=Path, help="Candidate JSON suite report to check.")
    suite_regress.add_argument("--absolute-tolerance", type=float, default=0.0, help="Largest permitted absolute decrease in a tracked metric.")
    suite_regress.add_argument("--allow-missing", action="store_true", help="Do not fail when a baseline case is absent from the candidate.")
    for command in (suite_validate, suite_compare, suite_evaluate, suite_matrix, suite_evaluate_matrix):
        command.add_argument("--id", action="append", default=[], help="Benchmark ID to include. Repeat to select several; omit to run all eligible benchmarks.")
    for command in (suite_validate, suite_compare, suite_evaluate, suite_matrix, suite_evaluate_matrix, suite_regress):
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
    bundle_create.add_argument("--assay", choices=["bulk_rna_seq", "sequence_dna", "sequence_rna", "sequence_protein", "short_read_dna", "whole_genome_dna_seq"])
    bundle_create.add_argument("--force", action="store_true", help="Replace an existing bundle.")
    bundle_verify = bundle_sub.add_parser("verify", help="Verify a bundle inventory, hashes and dataset manifests.")
    bundle_verify.add_argument("archive", type=Path)
    args = p.parse_args()
    try:
        if args.command == "seq":
            out = run_sequence_command(args)
        elif args.command == "variant":
            model, folder = lookup(args.root, args.id)
            out = compare_variants(model, folder, args.results, args.detail_limit)
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
        elif args.command == "suite":
            if args.suite_command == "validate":
                out = validate_suite(args.root, args.assay, args.id)
            elif args.suite_command == "compare":
                out = compare_suite(args.root, args.results, args.id, args.detail_limit)
            elif args.suite_command == "evaluate":
                out = evaluate_suite(args.root, args.results, args.id, args.detail_limit)
            elif args.suite_command == "matrix":
                out = compare_matrix(args.root, [parse_method(value) for value in args.method], args.id, args.detail_limit)
            elif args.suite_command == "evaluate-matrix":
                out = evaluate_matrix(args.root, [parse_method(value) for value in args.method], args.id, args.detail_limit)
            else:
                out = compare_reports(args.baseline, args.candidate, args.absolute_tolerance, args.allow_missing)
            out["reports"] = write_reports(
                out,
                json_path=args.json,
                markdown_path=args.markdown,
                junit_path=args.junit,
                force=args.force,
                csv_path=args.csv,
                html_path=args.html,
            )
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
        if args.command == "suite" and out["summary"]["status"] == "fail":
            sys.exit(2)
    except (ValueError,OSError,KeyError) as exc:
        if args.debug:
            raise
        print(f"omicsbench: {exc}",file=sys.stderr)
        sys.exit(1)
