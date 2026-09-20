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
from .validate import validate

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
    args = p.parse_args()
    try:
        if args.command == "seq":
            out = run_sequence_command(args)
        elif args.command == "assay":
            if args.assay_command == "list":
                out = list_profiles()
            elif args.assay_command == "info":
                out = get_profile(args.profile)
            else:
                out = build_plan(args.profile, args.input, args.output_directory, args.adapter)
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
    except (ValueError,OSError,KeyError) as exc:
        if args.debug:
            raise
        print(f"omicsbench: {exc}",file=sys.stderr)
        sys.exit(1)
