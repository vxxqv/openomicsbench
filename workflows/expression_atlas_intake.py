"""Stage and verify one Expression Atlas bulk RNA-seq experiment."""
import argparse
import json
from pathlib import Path

from omicsbench.expression_atlas import stage


parser = argparse.ArgumentParser(description="Stage raw counts and design metadata from Expression Atlas.")
parser.add_argument("accession", help="Expression Atlas experiment accession, for example E-MTAB-8572")
parser.add_argument("--staging", type=Path, required=True, help="Directory for source files; keep it outside release data.")
args = parser.parse_args()
print(json.dumps(stage(args.accession, args.staging), indent=2))
