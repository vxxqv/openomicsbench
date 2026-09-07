"""Deterministic count and paired-read operations."""
import csv
import gzip
import hashlib
from itertools import zip_longest
from pathlib import Path
import numpy as np

def read_counts(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    if len(rows) < 2 or len(rows[0]) < 3:
        raise ValueError(f"{path}: expected gene IDs and at least two sample columns")
    samples, genes = rows[0][1:], [r[0] for r in rows[1:]]
    if len(samples) != len(set(samples)) or len(genes) != len(set(genes)) or any(not x for x in genes + samples):
        raise ValueError(f"{path}: missing or duplicate sample/feature IDs")
    if any(len(r) != len(rows[0]) for r in rows[1:]):
        raise ValueError(f"{path}: inconsistent row width")
    if any(not v.isascii() or not v.isdigit() for r in rows[1:] for v in r[1:]):
        raise ValueError(f"{path}: raw counts must be non-negative integers")
    counts = np.array([r[1:] for r in rows[1:]], dtype=np.int64)
    if np.any(counts < 0) or np.any(counts.sum(axis=0, dtype=np.float64) <= 0):
        raise ValueError(f"{path}: negative counts or empty sample library")
    return genes, samples, counts

def check_design(samples):
    conditions = sorted({s.condition for s in samples})
    batches = sorted({s.batch for s in samples})
    if len(conditions) != 2:
        raise ValueError("The baseline requires exactly two conditions")
    if any(sum(s.condition == c for s in samples) < 2 for c in conditions):
        raise ValueError("Each condition needs at least two biological samples")
    x = np.array([[1, int(s.condition == conditions[1])] + [int(s.batch == b) for b in batches[1:]] for s in samples], dtype=float)
    if np.linalg.matrix_rank(x) < x.shape[1] or x.shape[0] <= x.shape[1]:
        raise ValueError("Condition and batch design is confounded or has no residual degrees of freedom")
    return x, conditions

def median_ratio_size_factors(counts):
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or counts.shape[1] < 2 or np.any(counts < 0):
        raise ValueError("Size-factor input must be a non-negative feature-by-sample matrix")
    eligible = np.all(counts > 0, axis=1)
    if np.count_nonzero(eligible) < 2:
        raise ValueError("Median-ratio normalization needs at least two features positive in every sample")
    geometric = np.exp(np.mean(np.log(counts[eligible]), axis=1))
    factors = np.median(counts[eligible] / geometric[:, None], axis=0)
    if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
        raise ValueError("Median-ratio normalization produced invalid size factors")
    return factors / np.exp(np.mean(np.log(factors)))


def diagnostic(counts, samples, normalization="cpm"):
    x, conditions = check_design(samples)
    libs = counts.sum(axis=0, dtype=np.float64)
    if normalization == "cpm":
        factors = libs / 1e6
    elif normalization == "median_ratio":
        factors = median_ratio_size_factors(counts)
    else:
        raise ValueError("Normalization must be cpm or median_ratio")
    transformed = np.log2(counts / factors + 1)
    effects = np.linalg.lstsq(x, transformed.T, rcond=None)[0][1]
    centered = transformed.T - transformed.mean(axis=1)
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    coords = u[:, :2] * s[:2]
    for col in range(coords.shape[1]):
        if coords[np.argmax(np.abs(coords[:, col])), col] < 0:
            coords[:, col] *= -1
    total = float(np.sum(s ** 2))
    explained = s[:2] ** 2 / total if total else np.zeros(2)
    dist = np.sqrt(((transformed.T[:, None, :] - transformed.T[None, :, :]) ** 2).sum(axis=2))
    return {"effects":effects,"distances":dist[np.triu_indices(len(samples),1)],"pca":coords,"explained":explained,"libraries":libs,"size_factors":factors,"normalization":normalization,"contrast":f"{conditions[1]} - {conditions[0]}"}

def ranked(values):
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Metric inputs must be finite")
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = (start + stop - 1) / 2
        start = stop
    return ranks

def correlation(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or len(a) < 3 or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Correlation requires matched finite arrays of length at least three")
    if np.std(a) == 0 or np.std(b) == 0:
        raise ValueError("Correlation is undefined for constant inputs")
    return float(np.corrcoef(a,b)[0,1])

def preservation(full, pocket, full_genes, pocket_genes, k=50):
    lookup = {g:i for i,g in enumerate(full_genes)}
    indices = [lookup[g] for g in pocket_genes]
    a,b = full["effects"][indices], pocket["effects"]
    if k > len(a) or k < 1:
        raise ValueError("top-k must fit inside the declared shared feature universe")
    # Gene IDs settle ties without depending on row order.
    top = lambda values: set(sorted(range(len(values)), key=lambda i:(-abs(values[i]),pocket_genes[i]))[:k])
    aa,bb = top(a),top(b)
    eligible = np.abs(a) >= 0.1
    if not np.any(eligible):
        raise ValueError("No effects satisfy the predeclared sign-concordance filter")
    return {"spearman_logfc":correlation(ranked(a),ranked(b)),"top_k_jaccard":len(aa & bb)/len(aa | bb),"distance_correlation":correlation(full["distances"],pocket["distances"]),"sign_concordance":float(np.mean(np.sign(a[eligible]) == np.sign(b[eligible])))}

def selection_order(genes, counts, seed):
    mean = counts.mean(axis=1)
    variance = np.var(np.log2(counts + 1),axis=1)
    priorities = [sorted(range(len(genes)), key=lambda i:(-variance[i],genes[i])), sorted(range(len(genes)),key=lambda i:(abs(mean[i]-float(np.median(mean))),genes[i])), sorted(range(len(genes)),key=lambda i:(variance[i],genes[i])), sorted(range(len(genes)),key=lambda i:hashlib.sha256(f"{seed}:{genes[i]}".encode()).digest())]
    picked = set()
    ordered = []
    cursors = [0]*4
    while len(picked) < len(genes):
        for group,order in enumerate(priorities):
            while order[cursors[group]] in picked:
                cursors[group] += 1
            index = order[cursors[group]]
            picked.add(index)
            ordered.append(index)
            if len(picked) == len(genes):
                break
    return ordered


def contrast_selection_order(genes, counts, effects, seed, anchors=100):
    effects = np.asarray(effects, dtype=float)
    if effects.shape != (len(genes),) or not np.all(np.isfinite(effects)):
        raise ValueError("Contrast effects must match the finite feature vector")
    if not 1 <= anchors <= len(genes):
        raise ValueError("Anchor count must fit inside the source feature count")
    effect_order = sorted(range(len(genes)), key=lambda i:(-abs(effects[i]),genes[i]))
    first = effect_order[:anchors]
    anchored = set(first)
    return first + [index for index in selection_order(genes,counts,seed) if index not in anchored]


def select_genes(genes, counts, size, seed):
    if not 3 <= size <= len(genes):
        raise ValueError("Candidate size must be between three and the source feature count")
    return sorted(selection_order(genes,counts,seed)[:size])

def fastq_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt",encoding="ascii") as stream:
        while True:
            head = stream.readline()
            if not head:
                break
            seq,plus,qual = [stream.readline() for _ in range(3)]
            if not head.startswith("@") or not plus.startswith("+") or not seq or len(seq.rstrip()) != len(qual.rstrip()):
                raise ValueError(f"{path}: malformed FASTQ record")
            yield head,seq,plus,qual

def read_id(record):
    name = record[0].split()[0][1:]
    return name[:-2] if name.endswith(("/1","/2")) else name

def paired_sample(r1,r2,out1,out2,fraction,seed,sample_id):
    if not 0 < fraction <= 1:
        raise ValueError("Read sampling fraction must be in (0,1]")
    selected = total = 0
    with open(out1,"w",encoding="ascii",newline="\n") as a, open(out2,"w",encoding="ascii",newline="\n") as b:
        for i,(left,right) in enumerate(zip_longest(fastq_records(r1),fastq_records(r2))):
            if left is None or right is None or read_id(left) != read_id(right):
                raise ValueError("Paired FASTQ files have different record counts or read IDs")
            total += 1
            value = int.from_bytes(hashlib.sha256(f"{seed}:{sample_id}:{i}:{read_id(left)}".encode()).digest()[:8],"big")
            if value < fraction * 2**64:
                a.writelines(left); b.writelines(right); selected += 1
    return {"input_pairs":total,"selected_pairs":selected,"seed":seed,"fraction":fraction,"algorithm":"sha256-pair-ordinal-v1"}
