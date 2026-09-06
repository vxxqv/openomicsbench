import json
from .hashing import digest, contained
from .rnaseq import read_counts, diagnostic, preservation

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
    baseline_path = profile["baseline_counts"]
    if baseline_path not in declared:
        raise ValueError("Baseline counts must be declared in the manifest")
    genes,ids,counts = read_counts(contained(folder,baseline_path))
    sample_order = [s.sample_id for s in model.samples]
    if ids != sample_order or ids != profile["sample_ids"]:
        raise ValueError("Baseline sample order differs from manifest or validation profile")
    full = diagnostic(counts,model.samples)
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
        metrics = preservation(full,diagnostic(c,model.samples),genes,g,profile["top_k"])
        for metric in model.validation.metrics:
            if metrics[metric.name] < metric.minimum:
                raise ValueError(f"{f.path}: {metric.name}={metrics[metric.name]:.5f} is below {metric.minimum}")
        results[f.tier] = metrics
    if not results:
        raise ValueError(f"{model.id}: no count tier was validated")
    return {"id":model.id,"status":"pass","kind":model.kind,"metrics":results}
