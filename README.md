# OpenOmicsBench

OpenOmicsBench is being built for people who need compact, traceable RNA-seq data for software tests, teaching and method development. Each published object must explain where its files came from, what changed during reduction, and which checks it passes.

The published 0.1.0.dev0 prerelease is the working baseline. Version 1 is now in development. The repository contains the Python core, one generated infrastructure fixture, a pinned Pasilla retrieval and diagnostic workflow, twelve source-intake records, and the collection requirements. Version 1 work will add certified biological pockets as source rights and reference compatibility are resolved.

## Quickstart

Use Python 3.11 or newer from an isolated environment. From this directory:

```sh
python -m pip install -e .
omicsbench list --assay bulk_rna_seq
omicsbench info fixture-001
omicsbench get fixture-001 --size nano
omicsbench validate fixture-001
omicsbench verify-cache
```

The fixture is generated data. It tests software behavior; it is not evidence about gene regulation. The `get` command copies exact declared files and verifies their SHA-256 hashes. It never normalizes the counts.

Run from another directory with `omicsbench --root /path/to/openomicsbench list`. Set `OMICSBENCH_ROOT` and `OMICSBENCH_CACHE` if you prefer environment variables. Dataset files belong to this checkout and are not silently bundled into the installed Python package.

## What works

The core discovers manifests without a network connection, checks the typed schema, validates file inventories, retrieves tiers into a verified cache, reports provenance, and recomputes count-based diagnostics. The tests exercise corruption, invalid metadata, interrupted downloads, confounded designs, sample order and deterministic rebuilds.

The Pasilla record is link-only. To stage its pinned source and reproduce the diagnostic size curve:

```sh
python workflows/pasilla_diagnostic.py --staging staging/pasilla
```

The source table is less than 0.5 MB. This workflow uses the original processed integer counts, preserves all seven libraries, and evaluates six feature subsets. It writes third-party data only to the selected staging directory. Do not add staging data to a release while the rights decision remains AMBER.

The full Pasilla table has also been fitted with DESeq2 1.50.2 under R 4.5.3. The repository keeps the script and a compact evidence record; the source table and gene-level results stay outside the release files while the rights review is open.

## Project documents

- [Implementation and release status](docs/implementation-status.md) follows all thirteen sections of the supplied plan.
- [Scientific methods](docs/methods.md) defines the current diagnostic and its limits.
- [Source review](docs/source-review.md) explains selection decisions and links to evidence.
- [Contributing](CONTRIBUTING.md) describes intake, review and writing standards.
- [Release procedure](docs/release.md) identifies the remaining publication gates.
- [Changes](CHANGELOG.md) records the published prerelease and current v1 development.

## Local verification

```sh
python -m unittest discover -s tests -v
python scripts/check_repository.py
```

Expansion into a certified biological collection requires reference compatibility checks on real sources, independent usability feedback and post-upload hash verification. These requirements remain visible so later releases can be judged against recorded evidence.
