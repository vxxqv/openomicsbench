# Portable benchmark bundles

A bundle is a deterministic ZIP containing selected benchmark objects and the licence files needed to reuse them outside an OpenOmicsBench installation.

## Create a bundle

Select objects by stable ID:

```sh
omicsbench bundle create selected.zip --id rnaseq-002 --id sequence-005
```

Or select one assay family:

```sh
omicsbench bundle create proteins.zip --assay sequence_protein
```

Each archive contains `bundle.json`, the selected dataset manifests, every file declared by those manifests, `CITATION.cff`, `LICENSE`, `LICENSE-METADATA` and `NOTICE`. The bundle inventory records the byte count and SHA-256 digest of every included file. Paths and ZIP timestamps are normalized, so the same package version and selection produce identical archive bytes.

An existing destination is protected unless `--force` is supplied.

## Verify a bundle

```sh
omicsbench bundle verify selected.zip
```

Verification reads the archive without extracting it. It rejects unsafe or duplicate paths, symbolic links, unlisted files, missing files, altered byte counts, altered SHA-256 digests and dataset manifests that do not match the bundle dataset list. A successful result reports the package version, dataset IDs, file count, archive size and archive SHA-256.

The command verifies the archive as created. It does not download sources or update a bundle to a newer OpenOmicsBench release.
