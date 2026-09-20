import json
from .hashing import digest, contained
from .rnaseq import read_counts, diagnostic, preservation
from .sequences import read_records, reverse_complement, summarize
from .sequence_analysis import pair_report


def _sequence_summary(path, molecule):
    report = summarize(path, molecule)
    report.pop("path", None)
    return report


def _variant_truth(model, folder, profile):
    reference_name = profile.get("reference")
    truth_name = profile.get("variant_truth")
    if not reference_name and not truth_name:
        return None
    declared = {item.path for item in model.files}
    if not reference_name or not truth_name or reference_name not in declared or truth_name not in declared:
        raise ValueError(f"{model.id}: DNA-seq truth paths are incomplete or undeclared")
    references = list(read_records(contained(folder, reference_name)))
    if len(references) != 1:
        raise ValueError(f"{model.id}: DNA-seq truth requires one reference sequence")
    reference = references[0].sequence
    lines = contained(folder, truth_name).read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "position\tref\talt":
        raise ValueError(f"{model.id}: variant truth header is invalid")
    variants = []
    seen = set()
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(f"{model.id}: variant truth row is invalid")
        position_text, ref, alt = fields
        try:
            position = int(position_text)
        except ValueError as error:
            raise ValueError(f"{model.id}: variant position is not an integer") from error
        if position in seen or position < 1 or position > len(reference):
            raise ValueError(f"{model.id}: variant position is duplicated or outside the reference")
        if ref not in "ACGT" or alt not in "ACGT" or ref == alt or reference[position - 1] != ref:
            raise ValueError(f"{model.id}: variant truth disagrees with the reference")
        seen.add(position)
        variants.append((position, ref, alt))
    if len(variants) != profile.get("variant_count"):
        raise ValueError(f"{model.id}: variant count differs from the expected profile")
    reads = []
    for item in profile["files"]:
        for record in read_records(contained(folder, item["path"])):
            reads.extend((record.sequence, reverse_complement(record.sequence, "dna")))
    flank = profile.get("context_bases", 10)
    supported = 0
    for position, _, alt in variants:
        center = position - 1
        left = max(0, center - flank)
        right = min(len(reference), center + flank + 1)
        context = reference[left:center] + alt + reference[center + 1:right]
        if not any(context in read for read in reads):
            raise ValueError(f"{model.id}: variant at position {position} is absent from the reads")
        supported += 1
    return {"reference": reference_name, "truth": truth_name, "variants": len(variants), "read_supported": supported}


def validate_sequence(model, folder, profile):
    records = profile.get("files")
    if not isinstance(records, list) or not records:
        raise ValueError(f"{model.id}: sequence validation profile has no files")
    observed = []
    for record in records:
        relative = record.get("path")
        if not relative or relative not in {item.path for item in model.files}:
            raise ValueError(f"{model.id}: sequence validation file is undeclared")
        report = _sequence_summary(contained(folder, relative), model.sequence.molecule)
        if report != record.get("expected"):
            raise ValueError(f"{model.id}/{relative}: sequence summary differs from the expected profile")
        observed.append({"path": relative, "summary": report})
    pairing = None
    if model.sequence.paired:
        if len(records) != 2:
            raise ValueError(f"{model.id}: paired sequence profile requires two files")
        pairing = pair_report(contained(folder, records[0]["path"]), contained(folder, records[1]["path"]))
        if pairing["status"] != "pass" or pairing["pairs"] != profile.get("pairs"):
            raise ValueError(f"{model.id}: paired-read validation differs from the expected profile")
    truth = _variant_truth(model, folder, profile)
    return {"id": model.id, "status": "pass", "kind": model.kind, "sequence": model.sequence.model_dump(mode="json"), "files": observed, "pairing": pairing, "variant_truth": truth}

def validate(model, folder):
    errors = []
    declared = {f.path for f in model.files}
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and p.name != "manifest.json"}
    if actual != declared:
        errors.append(f"File inventory differs: missing={sorted(declared-actual)}, undeclared={sorted(actual-declared)}")
    for f in model.files:
        path = contained(folder,f.path)
        if not path.is_file() or path.stat().st_size != f.bytes or digest(path) != f.sha256:
            errors.append(f"{model.id}/{f.path}: checksum or size mismatch")
    if errors:
        raise ValueError("; ".join(errors))
    if model.status == "original_link":
        return {"id":model.id,"status":"link_only","scientific_validation":"not_run","rights":model.rights.status}
    if model.validation is None:
        raise ValueError(f"{model.id}: no scientific validation profile")
    profile = json.loads(contained(folder, model.validation.profile).read_text(encoding="utf-8"))
    if model.sequence is not None:
        return validate_sequence(model, folder, profile)
    baseline_path = profile["baseline_counts"]
    if baseline_path not in declared:
        raise ValueError("Baseline counts must be declared in the manifest")
    genes,ids,counts = read_counts(contained(folder,baseline_path))
    sample_order = [s.sample_id for s in model.samples]
    if ids != sample_order or ids != profile["sample_ids"]:
        raise ValueError("Baseline sample order differs from manifest or validation profile")
    normalization = profile.get("normalization", "cpm")
    full = diagnostic(counts,model.samples,normalization=normalization)
    results = {}
    for f in model.files:
        if f.role != "raw_counts" or f.tier not in {"nano","pocket"}:
            continue
        g,s,c = read_counts(contained(folder,f.path))
        if s != ids or list(c.shape) != profile["shapes"][f.tier]:
            raise ValueError(f"{f.path}: sample order or shape differs from baseline")
        idx = {gene:i for i,gene in enumerate(genes)}
        if any(gene not in idx for gene in g) or not (c == counts[[idx[gene] for gene in g]]).all():
            raise ValueError(f"{f.path}: source counts changed during feature reduction")
        metrics = preservation(full,diagnostic(c,model.samples,normalization=normalization),genes,g,profile["top_k"])
        for metric in model.validation.metrics:
            if metrics[metric.name] < metric.minimum:
                raise ValueError(f"{f.path}: {metric.name}={metrics[metric.name]:.5f} is below {metric.minimum}")
        results[f.tier] = metrics
    if not results:
        raise ValueError(f"{model.id}: no count tier was validated")
    return {"id":model.id,"status":"pass","kind":model.kind,"metrics":results}
