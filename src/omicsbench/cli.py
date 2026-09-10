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
from .validate import validate

def main():
    p = argparse.ArgumentParser(description="Find and verify transcriptomics objects.",epilog="Example: omicsbench info rnaseq-002")
    p.add_argument("--root",type=Path,default=default_root(),help="Collection directory. Defaults to OMICSBENCH_ROOT, the current checkout or the installed collection.")
    p.add_argument("--cache",type=Path,default=Path(os.environ.get("OMICSBENCH_CACHE",".omicsbench-cache")),help="Local cache directory.")
    p.add_argument("--debug",action="store_true",help="Show a traceback when a command fails.")
    p.add_argument("--version",action="version",version=__version__)
    sub = p.add_subparsers(dest="command",required=True)
    examples = {"list":"list --assay bulk_rna_seq","info":"info fixture-001","get":"get fixture-001 --size nano","validate":"validate fixture-001","provenance":"provenance fixture-001","verify-cache":"verify-cache","doctor":"doctor"}
    for name,example in examples.items():
        command = sub.add_parser(name,help={"list":"List local manifests.","info":"Show files, rights and limitations.","get":"Cache an exact tier after checking hashes.","validate":"Check declared files and recompute metrics.","provenance":"Show source and transformations.","verify-cache":"Verify completed cache entries.","doctor":"Check local runtime and optional tools."}[name],epilog=f"Example: omicsbench {example}")
        if name in {"info","get","validate","provenance"}:
            command.add_argument("id",help="Stable dataset identifier.")
        if name == "get":
            command.add_argument("--size",choices=["nano","pocket","expected"],default="nano")
        if name == "list":
            command.add_argument("--assay",choices=["bulk_rna_seq"])
    args = p.parse_args()
    try:
        if args.command == "list":
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
            else:
                out = validate(model,folder)
        print(json.dumps(out,indent=2,allow_nan=False))
    except (ValueError,OSError,KeyError) as exc:
        if args.debug:
            raise
        print(f"omicsbench: {exc}",file=sys.stderr)
        sys.exit(1)
