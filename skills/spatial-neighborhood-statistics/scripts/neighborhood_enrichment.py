#!/usr/bin/env python3
"""Cell-type neighbourhood enrichment by label permutation, per sample, with a patient-level summary.

For each sample (slide, ROI, field of view) separately, build a spatial neighbour graph from cell or
spot coordinates, count focal-to-neighbour adjacencies for every pair of cell types, and compare each
count with its distribution when the cell-type labels are shuffled over the fixed positions (the
histoCAT neighbourhood permutation test [Schapiro2017]; z-score as in Squidpy [Palla2022]). Then
summarise across patients with one effect per patient, the unit of replication [Hurlbert1984].

Per sample
  C(A,B)       number of directed adjacencies i -> j with i a focal cell of type A and j one of
               its neighbours of type B. Every cell is focal unless --edge-buffer or --focal-col
               says otherwise. With a symmetric graph and every cell focal C(A,B) = C(B,A), so
               pairs are reported unordered (an A-A edge then counts twice); with a directed kNN
               graph, --edge-buffer or --focal-col they are ordered, in every sample alike.
               C(A,B) / n_A is the mean number of B neighbours per A cell; it gives the same
               permutation p-values because n_A does not change under permutation.
  null         labels permuted uniformly among the cells of the sample, or only within each
               stratum (--strata-col and/or --block-size); positions, graph and focal set fixed
               (random labelling of a fixed set of positions).
  expected     exact mean of C(A,B) under that null: the sum over focal edges i -> j of
               P(label_i = A, label_j = B), with P = n_sA (n_sB - [A = B]) / (N_s (N_s - 1)) when
               i and j share stratum s, and (n_sA / N_s)(n_tB / N_t) when they lie in strata s
               and t. Squidpy estimates this mean from the permutations; the two agree up to
               Monte Carlo error.
  sd, z        standard deviation (n - 1) of the permuted counts; z = (observed - expected) / sd
               [Palla2022]. z grows with the number of cells, so do not compare it across samples.
  log2_ratio   log2((observed + 0.5) / (expected + 0.5)); the 0.5 pseudocount is a convention.
  p values     p_enrich = (1 + #{perm >= observed}) / (n_perm + 1), p_deplete likewise with <=
               (one-sided tests for interaction and avoidance [Schapiro2017], with the +1 in the
               numerator so that p is never zero [Phipson2010]); p_two_sided =
               min(1, 2 min(p_enrich, p_deplete)); q = BH [Benjamini1995] or BY [Benjamini2001]
               across the pairs tested in that sample.

Across samples
  The effect of a patient is the mean log2_ratio over that patient's samples (each sample weighted
  equally). Per pair: median effect, direction counts and an exact two-sided sign test over
  patients [Dixon1946], whose smallest attainable p is 2 ** (1 - n) for n patients with a non-zero
  effect; q across pairs. Without --patient-col every sample counts as a separate patient.
  The sign test treats each patient's direction as a fair coin under the null. With small expected
  counts the permutation distribution is discrete and skewed, so it is only approximately fair;
  median_expected (the median over samples of the expected count) shows when that applies.

Pure Python: run time grows with cells x neighbours x permutations.

Usage examples
  python3 neighborhood_enrichment.py cells.tsv --x-col x --y-col y --type-col cell_type \\
      --sample-col slide --patient-col patient --um-per-unit 0.5 --graph knn --k 6 \\
      --max-dist 25 --n-perm 999 --seed 1 --tsv-out pairs.tsv --summary-tsv summary.tsv > result.json
  python3 neighborhood_enrichment.py cells.csv --type-col phenotype --sample-col roi \\
      --graph radius --radius 20 --strata-col compartment --edge-buffer 20 \\
      --pairs CD8T:Tumor,Treg:CD8T --out result.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from collections import Counter, defaultdict
from operator import add

PSEUDOCOUNT = 0.5
DIST_TOL = 1e-9      # relative tolerance: a pair exactly at the radius or cap counts despite rounding
TRUE_VALUES = {"1", "true", "t", "yes", "y"}
FALSE_VALUES = {"0", "false", "f", "no", "n"}
PAIR_COLUMNS = ["sample", "patient", "type_a", "type_b", "n_a", "n_b", "observed", "expected", "sd",
                "z", "log2_ratio", "p_enrich", "p_deplete", "p_two_sided", "q", "call"]
SUMMARY_COLUMNS = ["type_a", "type_b", "n_samples", "n_units", "median_expected", "median_log2_ratio",
                   "n_pos", "n_neg", "n_zero", "sign_test_p", "sign_test_p_floor", "q",
                   "n_samples_enriched", "n_samples_depleted", "unit_effects"]


class InputError(Exception):
    """Invalid input or options; main() prints the message and exits with code 2."""


def sig(x, digits=6):
    """Round to significant digits for output; None stays None."""
    if x is None:
        return None
    return float(f"{x:.{digits}g}")


# --------------------------------------------------------------------------------------- input

def read_table(path, sep=None):
    """Return (header, rows) of a TSV or CSV file; the delimiter follows the extension unless given."""
    if sep is None:
        sep = "," if path.lower().endswith(".csv") else "\t"
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rows = [r for r in csv.reader(fh, delimiter=sep) if any(c.strip() for c in r)]
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc}") from exc
    if not rows:
        raise InputError(f"{path} is empty")
    if len(rows) < 2:
        raise InputError(f"{path} has a header but no data rows")
    return [h.strip() for h in rows[0]], rows[1:]


def parse_samples(header, rows, x_col="x", y_col="y", type_col="cell_type", sample_col=None,
                  patient_col=None, strata_col=None, focal_col=None, um_per_unit=None):
    """Group cells by sample: {sample: {x, y, type, stratum, focal, patient}} in input order."""
    wanted = {"--x-col": x_col, "--y-col": y_col, "--type-col": type_col, "--sample-col": sample_col,
              "--patient-col": patient_col, "--strata-col": strata_col, "--focal-col": focal_col}
    pos = {}
    for opt, col in wanted.items():
        if col is None:
            continue
        if col not in header:
            raise InputError(f"column {col!r} given by {opt} is not in the header: {', '.join(header)}")
        pos[opt] = header.index(col)
    scale = 1.0 if um_per_unit is None else float(um_per_unit)
    samples = {}
    for lineno, row in enumerate(rows, start=2):
        if len(row) < len(header):
            raise InputError(f"line {lineno} has {len(row)} fields but the header has {len(header)}")

        def val(opt):
            return row[pos[opt]].strip()

        try:
            x, y = float(val("--x-col")), float(val("--y-col"))
        except ValueError:
            raise InputError(f"line {lineno}: coordinates must be numbers, got "
                             f"{x_col}={val('--x-col')!r}, {y_col}={val('--y-col')!r}") from None
        if not (math.isfinite(x) and math.isfinite(y)):
            raise InputError(f"line {lineno}: coordinates must be finite numbers")
        label = val("--type-col")
        if not label:
            raise InputError(f"line {lineno}: empty cell type in column {type_col!r}")
        sample = val("--sample-col") if sample_col else "all"
        if not sample:
            raise InputError(f"line {lineno}: empty sample id in column {sample_col!r}")
        rec = samples.setdefault(sample, {"x": [], "y": [], "type": [], "stratum": [], "focal": [],
                                          "patients": set()})
        if patient_col:
            patient = val("--patient-col")
            if not patient:
                raise InputError(f"line {lineno}: empty patient id in column {patient_col!r}")
            rec["patients"].add(patient)
        if strata_col:
            stratum = val("--strata-col")
            if not stratum:
                raise InputError(f"line {lineno}: empty stratum in column {strata_col!r}")
            rec["stratum"].append(stratum)
        if focal_col:
            flag = val("--focal-col").lower()
            if flag in TRUE_VALUES:
                rec["focal"].append(True)
            elif flag in FALSE_VALUES:
                rec["focal"].append(False)
            else:
                raise InputError(f"line {lineno}: column {focal_col!r} must be 1/0 or true/false, got {flag!r}")
        rec["x"].append(x * scale)
        rec["y"].append(y * scale)
        rec["type"].append(label)
    for name, rec in samples.items():
        patients = rec.pop("patients")
        if patient_col and len(patients) != 1:
            raise InputError(f"sample {name!r} maps to {len(patients)} patients "
                             f"({', '.join(sorted(patients))}); each sample must come from one patient")
        rec["patient"] = next(iter(patients)) if patient_col else None
        if not strata_col:
            rec["stratum"] = None
        if not focal_col:
            rec["focal"] = None
    return samples


def parse_pairs(text):
    """'A:B,C:D' -> [('A', 'B'), ('C', 'D')]."""
    pairs = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            raise InputError(f"--pairs entry {item!r} must look like TypeA:TypeB")
        pairs.append((parts[0].strip(), parts[1].strip()))
    if not pairs:
        raise InputError("--pairs is empty")
    return pairs


# --------------------------------------------------------------------------------------- graph

def _ring(cx, cy, r):
    """Grid cells at Chebyshev distance r from (cx, cy)."""
    if r == 0:
        yield (cx, cy)
        return
    for dx in range(-r, r + 1):
        yield (cx + dx, cy - r)
        yield (cx + dx, cy + r)
    for dy in range(-r + 1, r):
        yield (cx - r, cy + dy)
        yield (cx + r, cy + dy)


def _buckets(xs, ys, h):
    buckets = defaultdict(list)
    for i, (x, y) in enumerate(zip(xs, ys)):
        buckets[(math.floor(x / h), math.floor(y / h))].append(i)
    return buckets


def knn_lists(xs, ys, k):
    """The k nearest other points of every point, nearest first; ties go to the lower index."""
    n = len(xs)
    if n < 2 or k < 1:
        return [[] for _ in range(n)]
    k = min(k, n - 1)
    width, height = max(xs) - min(xs), max(ys) - min(ys)
    if width > 0 and height > 0:
        h = math.sqrt(width * height * k / n)
    else:
        extent = max(width, height)
        h = extent * k / n if extent > 0 else 1.0
    buckets = _buckets(xs, ys, h)
    out = []
    for i in range(n):
        xi, yi = xs[i], ys[i]
        cx, cy = math.floor(xi / h), math.floor(yi / h)
        cand, r = [], 0
        while True:
            if (2 * r + 1) ** 2 > 4 * n:      # scanning would cost more than a full pass
                cand = sorted(((xs[j] - xi) ** 2 + (ys[j] - yi) ** 2, j) for j in range(n) if j != i)
                break
            for key in _ring(cx, cy, r):
                for j in buckets.get(key, ()):
                    if j != i:
                        cand.append(((xs[j] - xi) ** 2 + (ys[j] - yi) ** 2, j))
            if len(cand) >= k:
                cand.sort()
                # every point not yet scanned lies at least r * h away (shrunk slightly so that
                # rounding in the bucket assignment cannot hide a closer point)
                if cand[k - 1][0] < (r * h * (1.0 - 1e-9)) ** 2:
                    break
            r += 1
        out.append([j for _, j in cand[:k]])
    return out


def radius_lists(xs, ys, radius):
    """All other points within `radius` (inclusive, relative tolerance DIST_TOL) of every point."""
    if radius <= 0:
        raise InputError("--radius must be positive")
    h = radius * (1.0 + 2 * DIST_TOL)   # a pair within the radius can never sit two buckets apart
    buckets = _buckets(xs, ys, h)
    r2 = (radius * (1.0 + DIST_TOL)) ** 2
    out = []
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        cx, cy = math.floor(xi / h), math.floor(yi / h)
        nb = []
        for gx in (cx - 1, cx, cx + 1):
            for gy in (cy - 1, cy, cy + 1):
                for j in buckets.get((gx, gy), ()):
                    if j != i and (xs[j] - xi) ** 2 + (ys[j] - yi) ** 2 <= r2:
                        nb.append(j)
        nb.sort()
        out.append(nb)
    return out


def build_graph(xs, ys, graph="knn", k=6, knn_mode="union", radius=None, max_dist=None):
    """Neighbour lists, whether the graph is symmetric, and how many edges the distance cap removed."""
    if graph == "knn":
        nbrs = knn_lists(xs, ys, k)
        symmetric = knn_mode == "union"
        if symmetric:
            sets = [set(lst) for lst in nbrs]
            for i, lst in enumerate(nbrs):
                for j in lst:
                    sets[j].add(i)
            nbrs = [sorted(s) for s in sets]
    elif graph == "radius":
        nbrs, symmetric = radius_lists(xs, ys, radius), True
    else:
        raise InputError(f"unknown graph {graph!r}; use knn or radius")
    pruned = 0
    if max_dist is not None:
        cap = (max_dist * (1.0 + DIST_TOL)) ** 2
        capped = []
        for i, lst in enumerate(nbrs):
            keep = [j for j in lst if (xs[j] - xs[i]) ** 2 + (ys[j] - ys[i]) ** 2 <= cap]
            pruned += len(lst) - len(keep)
            capped.append(keep)
        nbrs = capped
    return nbrs, symmetric, pruned


def convex_hull(points):
    """Andrew's monotone chain; returns hull vertices counter-clockwise."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _segment_distance(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    length2 = dx * dx + dy * dy
    if length2 == 0:
        return math.hypot(p[0] - a[0], p[1] - a[1])
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2))
    return math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy))


def focal_by_buffer(xs, ys, buffer, window="bbox"):
    """True for cells at least `buffer` from the edge of the sample's bounding box or convex hull."""
    if window == "bbox":
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        return [min(x - x0, x1 - x, y - y0, y1 - y) >= buffer for x, y in zip(xs, ys)]
    hull = convex_hull(list(zip(xs, ys)))
    if len(hull) < 3:
        return [False] * len(xs)
    edges = list(zip(hull, hull[1:] + hull[:1]))
    return [min(_segment_distance((x, y), a, b) for a, b in edges) >= buffer for x, y in zip(xs, ys)]


# --------------------------------------------------------------------------------------- statistics

def adjust_p(pvalues, method="bh"):
    """Benjamini-Hochberg [Benjamini1995] or Benjamini-Yekutieli [Benjamini2001] adjusted p-values."""
    m = len(pvalues)
    if m == 0:
        return []
    factor = sum(1.0 / i for i in range(1, m + 1)) if method == "by" else 1.0
    order = sorted(range(m), key=lambda i: pvalues[i])
    q, running = [0.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m * factor / rank)
        q[i] = min(1.0, running)
    return q


def sign_test(n_pos, n_neg):
    """Exact two-sided sign test [Dixon1946]; zeros are excluded before calling."""
    n = n_pos + n_neg
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, i) for i in range(min(n_pos, n_neg) + 1)) / 2 ** n
    return min(1.0, 2.0 * tail)


def sign_test_floor(n):
    """Smallest two-sided sign-test p-value attainable with n non-zero units."""
    return 1.0 if n < 1 else min(1.0, 2.0 ** (1 - n))


def permutation_counts(labels, nbrs, focal, strata, codes, n_types, n_perm, rng):
    """Observed counts, exact null means and permutation tallies for the requested pair codes.

    labels: int label per cell; nbrs: neighbour lists; focal: bool per cell; strata: hashable key
    per cell; codes: list of a * n_types + b.
    Returns {code: (observed, expected, sd, n_ge, n_le, permutation_mean)}.
    """
    n = len(labels)
    groups = {}
    for i, s in enumerate(strata):
        groups.setdefault(s, []).append(i)
    group_list = list(groups.values())
    order = [i for g in group_list for i in g]
    new_index = [0] * n
    for new, old in enumerate(order):
        new_index[old] = new
    stratum_of = [0] * n
    for g_id, g in enumerate(group_list):
        for i in g:
            stratum_of[new_index[i]] = g_id
    src, dst = [], []
    for i in range(n):
        if focal[i]:
            for j in nbrs[i]:
                src.append(new_index[i])
                dst.append(new_index[j])
    base = [labels[i] for i in order]

    def count(lab):
        scaled = [v * n_types for v in lab]
        return Counter(map(add, map(scaled.__getitem__, src), map(lab.__getitem__, dst)))

    observed_all = count(base)
    observed = {c: observed_all.get(c, 0) for c in codes}

    # exact expectation under the (stratified) permutation null
    sizes = [len(g) for g in group_list]
    comp = [Counter(labels[i] for i in g) for g in group_list]
    edge_strata = Counter(zip(map(stratum_of.__getitem__, src), map(stratum_of.__getitem__, dst)))
    expected = {}
    for c in codes:
        a, b = divmod(c, n_types)
        total = 0.0
        for (s, t), m in edge_strata.items():
            if s == t:
                ns = sizes[s]
                if ns > 1:
                    total += m * comp[s][a] * (comp[s][b] - (1 if a == b else 0)) / (ns * (ns - 1))
            else:
                total += m * (comp[s][a] / sizes[s]) * (comp[t][b] / sizes[t])
        expected[c] = total

    blocks, start = [], 0
    for g in group_list:
        block = base[start:start + len(g)]
        blocks.append((block, len(set(block)) > 1))
        start += len(g)
    s1 = dict.fromkeys(codes, 0)
    s2 = dict.fromkeys(codes, 0)
    ge = dict.fromkeys(codes, 0)
    le = dict.fromkeys(codes, 0)
    for _ in range(n_perm):
        perm = []
        for block, mixed in blocks:
            if mixed:
                block = block[:]
                rng.shuffle(block)
            perm.extend(block)
        cnt = count(perm)
        for c in codes:
            v = cnt.get(c, 0)
            s1[c] += v
            s2[c] += v * v
            if v >= observed[c]:
                ge[c] += 1
            if v <= observed[c]:
                le[c] += 1
    out = {}
    for c in codes:
        if n_perm > 1:
            var = (s2[c] - s1[c] * s1[c] / n_perm) / (n_perm - 1)
            sd = math.sqrt(var) if var > 1e-12 else 0.0
        else:
            sd = 0.0
        out[c] = (observed[c], expected[c], sd, ge[c], le[c], s1[c] / n_perm)
    return out


def analyze_sample(xs, ys, types, *, sample="all", patient=None, graph="knn", k=6,
                   knn_mode="union", radius=None, max_dist=None, max_dist_factor=None, focal=None,
                   edge_buffer=None, window="bbox", strata=None, block_size=None, min_cells=10,
                   pairs=None, n_perm=999, seed=1, fdr="bh", alpha=0.05, force_ordered=False):
    """Run the neighbourhood permutation test on one sample. Returns {'info', 'rows', 'notes'}.

    force_ordered: report ordered pairs even when this sample's graph is symmetric with every cell
    focal, so that several samples analysed together share one pair definition.
    """
    n = len(xs)
    notes = []
    info = {"sample": sample, "patient": patient, "n_cells": n}
    counts = Counter(types)
    info["cell_types"] = dict(sorted(counts.items()))
    if n < 3:
        notes.append(f"{sample}: fewer than 3 cells; sample skipped")
        info["skipped"] = True
        return {"info": info, "rows": [], "notes": notes}

    nn1 = knn_lists(xs, ys, 1)
    nn_dist = [math.hypot(xs[i] - xs[l[0]], ys[i] - ys[l[0]]) for i, l in enumerate(nn1) if l]
    median_nn = statistics.median(nn_dist) if nn_dist else 0.0
    duplicates = sum(1 for d in nn_dist if d == 0)
    cap = max_dist
    if max_dist_factor is not None:
        if median_nn <= 0:
            raise InputError(f"sample {sample!r}: the median nearest-neighbour distance is 0 ({duplicates} "
                             "cells share coordinates), so --max-dist-factor is undefined; check for "
                             "merged slides or duplicated rows, or give --max-dist")
        cap = max_dist_factor * median_nn
    k_used = min(k, n - 1) if graph == "knn" else None
    if graph == "knn" and k_used < k:
        notes.append(f"{sample}: only {n} cells, k reduced to {k_used}")
    nbrs, symmetric, pruned = build_graph(xs, ys, graph, k, knn_mode, radius, cap)
    lengths = [math.hypot(xs[i] - xs[j], ys[i] - ys[j]) for i, lst in enumerate(nbrs) for j in lst]
    degrees = [len(lst) for lst in nbrs]

    if focal is None and edge_buffer is not None:
        focal = focal_by_buffer(xs, ys, edge_buffer, window)
    if focal is None:
        focal = [True] * n
    n_focal = sum(focal)
    ordered = bool(force_ordered) or (not symmetric) or n_focal < n

    if strata is None:
        strata = [""] * n
    if block_size is not None:
        x0, y0 = min(xs), min(ys)
        strata = [(s, math.floor((x - x0) / block_size), math.floor((y - y0) / block_size))
                  for s, x, y in zip(strata, xs, ys)]
    stratum_types = defaultdict(set)
    stratum_size = Counter(strata)
    for s, t in zip(strata, types):
        stratum_types[s].add(t)
    single_type_cells = sum(stratum_size[s] for s, ts in stratum_types.items() if len(ts) == 1)

    info["graph"] = {
        "type": graph, "k": k_used, "knn_mode": knn_mode if graph == "knn" else None,
        "radius": radius if graph == "radius" else None, "max_dist": sig(cap) if cap is not None else None,
        "symmetric": symmetric, "directed_edges": sum(degrees),
        "mean_degree": sig(sum(degrees) / n), "isolated_cells": sum(1 for d in degrees if d == 0),
        "edges_removed_by_max_dist": pruned, "median_nn_distance": sig(median_nn),
        "median_edge_length": sig(statistics.median(lengths)) if lengths else None,
        "max_edge_length": sig(max(lengths)) if lengths else None,
        "cells_with_duplicate_coordinates": duplicates,
    }
    info["n_focal"] = n_focal
    info["pairs_ordered"] = ordered
    info["n_strata"] = len(stratum_size)
    info["cells_in_single_type_strata"] = single_type_cells

    if duplicates:
        notes.append(f"{sample}: {duplicates} cells share exact coordinates with another cell; if several "
                     "slides were merged into one coordinate frame, give --sample-col")
    if info["graph"]["isolated_cells"] > 0.1 * n:
        notes.append(f"{sample}: {info['graph']['isolated_cells']} of {n} cells have no neighbour; the "
                     "radius or distance cap is small for this cell density")
    if lengths and cap is None and graph == "knn" and max(lengths) > 5 * median_nn > 0:
        notes.append(f"{sample}: longest kNN edge is {max(lengths):.3g} versus a median nearest-neighbour "
                     f"distance of {median_nn:.3g}; edges may bridge gaps in the tissue, consider --max-dist")
    if n_focal == 0:
        notes.append(f"{sample}: no focal cells (edge buffer larger than the tissue); sample skipped")
        info["skipped"] = True
        return {"info": info, "rows": [], "notes": notes}
    if single_type_cells > 0.5 * n:
        notes.append(f"{sample}: {single_type_cells} of {n} cells sit in strata that contain a single cell "
                     "type; the restricted null has little to permute and tests lose power")

    names = sorted(counts)
    index = {t: i for i, t in enumerate(names)}
    eligible = [t for t in names if counts[t] >= min_cells]
    excluded = [t for t in names if counts[t] < min_cells]
    info["types_below_min_cells"] = excluded
    if pairs:
        tested, missing = [], []
        for a, b in pairs:
            if not ordered and b < a:
                a, b = b, a
            if (a, b) in tested or (a, b) in missing:
                continue
            (tested if a in eligible and b in eligible else missing).append((a, b))
        if missing:
            notes.append(f"{sample}: requested pairs not testable here (type absent or below "
                         f"--min-cells): {', '.join(a + ':' + b for a, b in missing)}")
    elif ordered:
        tested = [(a, b) for a in eligible for b in eligible]
    else:
        tested = [(a, b) for i, a in enumerate(eligible) for b in eligible[i:]]
    info["n_tests"] = len(tested)
    info["p_floor_two_sided"] = sig(min(1.0, 2.0 / (n_perm + 1)))
    if not tested:
        notes.append(f"{sample}: no testable pair (needs types with at least {min_cells} cells)")
        return {"info": info, "rows": [], "notes": notes}
    if 2.0 / (n_perm + 1) > alpha / len(tested):
        notes.append(f"{sample}: with {n_perm} permutations the smallest two-sided p is "
                     f"{2.0 / (n_perm + 1):.2g}, above alpha / tests = {alpha / len(tested):.2g}; "
                     "increase --n-perm before relying on adjusted p-values")

    n_types = len(names)
    codes = [index[a] * n_types + index[b] for a, b in tested]
    labels = [index[t] for t in types]
    rng = random.Random(f"{seed}:{sample}")
    stats = permutation_counts(labels, nbrs, focal, strata, codes, n_types, n_perm, rng)
    rows = []
    for (a, b), c in zip(tested, codes):
        obs, exp, sd, n_ge, n_le, _ = stats[c]
        p_enrich = (1 + n_ge) / (n_perm + 1)
        p_deplete = (1 + n_le) / (n_perm + 1)
        rows.append({
            "sample": sample, "patient": patient, "type_a": a, "type_b": b,
            "n_a": counts[a], "n_b": counts[b], "observed": obs, "expected": exp, "sd": sd,
            "z": (obs - exp) / sd if sd > 0 else None,
            "log2_ratio": math.log2((obs + PSEUDOCOUNT) / (exp + PSEUDOCOUNT)),
            "p_enrich": p_enrich, "p_deplete": p_deplete,
            "p_two_sided": min(1.0, 2.0 * min(p_enrich, p_deplete)),
        })
    for row, q in zip(rows, adjust_p([r["p_two_sided"] for r in rows], fdr)):
        row["q"] = q
        if q < alpha:
            row["call"] = "enriched" if row["observed"] > row["expected"] else "depleted"
        else:
            row["call"] = "ns"
    return {"info": info, "rows": rows, "notes": notes}


def summarize_across(rows, alpha=0.05, fdr="bh"):
    """One effect per patient (mean log2_ratio over its samples), then a sign test per pair."""
    by_pair = defaultdict(list)
    for r in rows:
        by_pair[(r["type_a"], r["type_b"])].append(r)
    out = []
    for (a, b), rs in by_pair.items():
        units = defaultdict(list)
        for r in rs:
            units[r["patient"] if r["patient"] is not None else r["sample"]].append(r["log2_ratio"])
        effects = {u: statistics.fmean(v) for u, v in units.items()}
        values = list(effects.values())
        n_pos = sum(1 for v in values if v > 0)
        n_neg = sum(1 for v in values if v < 0)
        out.append({
            "type_a": a, "type_b": b, "n_samples": len(rs), "n_units": len(values),
            "median_expected": statistics.median(r["expected"] for r in rs),
            "median_log2_ratio": statistics.median(values), "n_pos": n_pos, "n_neg": n_neg,
            "n_zero": len(values) - n_pos - n_neg, "sign_test_p": sign_test(n_pos, n_neg),
            "sign_test_p_floor": sign_test_floor(n_pos + n_neg),
            "n_samples_enriched": sum(1 for r in rs if r["call"] == "enriched"),
            "n_samples_depleted": sum(1 for r in rs if r["call"] == "depleted"),
            "unit_effects": effects,
        })
    for item, q in zip(out, adjust_p([o["sign_test_p"] for o in out], fdr)):
        item["q"] = q
    return out


# --------------------------------------------------------------------------------------- output

def _rounded(row):
    return {k: (sig(v) if isinstance(v, float) else v) for k, v in row.items()}


def _tsv_value(v):
    if v is None:
        return "NA"
    if isinstance(v, dict):
        return ";".join(f"{k}:{sig(x):+g}" for k, x in v.items())
    if isinstance(v, float):
        return f"{sig(v):g}"
    return str(v)


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(_tsv_value(r.get(c)) for c in columns) + "\n")


def run(args):
    """Validate options, analyse every sample and return the result dictionary."""
    if args.graph == "radius" and (args.radius is None or args.radius <= 0):
        raise InputError("--graph radius needs a positive --radius")
    if args.graph == "knn" and args.radius is not None:
        raise InputError("--radius applies only to --graph radius (the graph is knn by default); add "
                         "--graph radius, or drop --radius and cap kNN edges with --max-dist")
    if args.graph == "radius" and args.k is not None:
        raise InputError("--k applies only to --graph knn; drop it or use --graph knn")
    if args.graph == "knn" and args.k is None:
        args.k = 6
    if args.graph == "knn" and args.k < 1:
        raise InputError("--k must be at least 1")
    if args.max_dist is not None and args.max_dist_factor is not None:
        raise InputError("give --max-dist or --max-dist-factor, not both")
    if args.edge_buffer is not None and args.focal_col:
        raise InputError("give --edge-buffer or --focal-col, not both")
    for opt in ("max_dist", "max_dist_factor", "block_size", "um_per_unit"):
        value = getattr(args, opt)
        if value is not None and value <= 0:
            raise InputError(f"--{opt.replace('_', '-')} must be positive")
    if args.edge_buffer is not None and args.edge_buffer < 0:
        raise InputError("--edge-buffer must not be negative")
    if args.n_perm < 1:
        raise InputError("--n-perm must be at least 1")
    if args.min_cells < 1:
        raise InputError("--min-cells must be at least 1")
    if not 0 < args.alpha < 1:
        raise InputError("--alpha must lie between 0 and 1")
    sep = {"tab": "\t", "comma": ","}.get(args.sep)
    header, rows = read_table(args.table, sep)
    samples = parse_samples(header, rows, args.x_col, args.y_col, args.type_col, args.sample_col,
                            args.patient_col, args.strata_col, args.focal_col, args.um_per_unit)
    pairs = parse_pairs(args.pairs) if args.pairs else None
    if pairs:
        known = {t for rec in samples.values() for t in rec["type"]}
        unknown = sorted({t for p in pairs for t in p} - known)
        if unknown:
            raise InputError(f"--pairs names cell types absent from the table: {', '.join(unknown)}; "
                             f"known types: {', '.join(sorted(known))}")
    unit = "um" if args.um_per_unit is not None else "coordinate units"
    result = {"distance_unit": unit, "parameters": {k: v for k, v in vars(args).items()},
              "samples": [], "pairs": [], "across_samples": [], "notes": []}
    if args.um_per_unit is None:
        result["notes"].append("distances are in coordinate units; give --um-per-unit to state the "
                               "neighbourhood in micrometres")
    # one pair definition for all samples: ordered whenever the design can make C(A,B) != C(B,A)
    force_ordered = ((args.graph == "knn" and args.knn_mode == "directed")
                     or args.edge_buffer is not None or bool(args.focal_col))
    for name, rec in samples.items():
        res = analyze_sample(rec["x"], rec["y"], rec["type"], sample=name, patient=rec["patient"],
                             graph=args.graph, k=args.k, knn_mode=args.knn_mode, radius=args.radius,
                             max_dist=args.max_dist, max_dist_factor=args.max_dist_factor,
                             focal=rec["focal"], edge_buffer=args.edge_buffer, window=args.window,
                             strata=rec["stratum"], block_size=args.block_size,
                             min_cells=args.min_cells, pairs=pairs, n_perm=args.n_perm,
                             seed=args.seed, fdr=args.fdr, alpha=args.alpha,
                             force_ordered=force_ordered)
        result["samples"].append(res["info"])
        result["pairs"].extend(res["rows"])
        result["notes"].extend(res["notes"])
    result["across_samples"] = summarize_across(result["pairs"], args.alpha, args.fdr)
    units = {r["patient"] if r["patient"] is not None else r["sample"] for r in result["pairs"]}
    if len(samples) > 1 and not args.patient_col:
        result["notes"].append("no --patient-col: every sample is treated as a separate patient; if "
                               "several samples come from one patient, the across-sample summary "
                               "overstates the replication")
    if result["pairs"] and len(units) < 2:
        result["notes"].append("one patient (or sample) only: per-sample tests describe this tissue; "
                               "they say nothing about other patients")
    elif result["pairs"]:
        result["notes"].append(f"across-sample summary over {len(units)} units: the smallest attainable "
                               f"two-sided sign-test p is {sign_test_floor(len(units)):.3g}")
    return result


def build_parser():
    p = argparse.ArgumentParser(
        description="Cell-type neighbourhood enrichment by label permutation within each sample, "
                    "plus a patient-level summary. Coordinates are multiplied by --um-per-unit on "
                    "input, so every distance option is in micrometres when it is given.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Exit code 2 with a message naming the column, line or option when the input is wrong.")
    p.add_argument("table", help="TSV or CSV, one row per cell or spot")
    p.add_argument("--x-col", default="x", help="x coordinate column (default x)")
    p.add_argument("--y-col", default="y", help="y coordinate column (default y)")
    p.add_argument("--type-col", default="cell_type", help="cell-type label column (default cell_type)")
    p.add_argument("--sample-col", help="slide / ROI / section column; graphs never cross samples")
    p.add_argument("--patient-col", help="patient column; the across-sample summary uses one effect per patient")
    p.add_argument("--sep", choices=["tab", "comma"], help="delimiter (default: comma for .csv, else tab)")
    p.add_argument("--um-per-unit", type=float, help="micrometres per coordinate unit (pixel size)")
    p.add_argument("--graph", choices=["knn", "radius"], default="knn", help="neighbour graph (default knn)")
    p.add_argument("--k", type=int, help="neighbours per cell for --graph knn (default 6)")
    p.add_argument("--knn-mode", choices=["union", "directed"], default="union",
                   help="union: symmetric (i~j if either is among the other's k nearest); directed: "
                        "each cell's own k nearest (default union)")
    p.add_argument("--radius", type=float, help="neighbour radius for --graph radius")
    p.add_argument("--max-dist", type=float, help="drop edges longer than this distance")
    p.add_argument("--max-dist-factor", type=float,
                   help="drop edges longer than this multiple of the sample's median nearest-neighbour distance")
    p.add_argument("--edge-buffer", type=float,
                   help="only cells at least this far from the window edge are focal (neighbours still count)")
    p.add_argument("--window", choices=["bbox", "hull"], default="bbox",
                   help="window for --edge-buffer: bounding box or convex hull of the sample's cells")
    p.add_argument("--focal-col", help="column of 1/0 flags marking focal cells (e.g. from a tissue mask)")
    p.add_argument("--strata-col", help="permute labels only within levels of this column (e.g. compartment)")
    p.add_argument("--block-size", type=float, help="also permute only within square blocks of this side")
    p.add_argument("--min-cells", type=int, default=10,
                   help="a type is tested in a sample only with at least this many cells (default 10)")
    p.add_argument("--pairs", help="pre-specified pairs to test, e.g. CD8T:Tumor,Treg:CD8T (default all)")
    p.add_argument("--n-perm", type=int, default=999, help="label permutations per sample (default 999)")
    p.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    p.add_argument("--fdr", choices=["bh", "by"], default="bh", help="multiple-testing adjustment (default bh)")
    p.add_argument("--alpha", type=float, default=0.05, help="level for the enriched/depleted call (default 0.05)")
    p.add_argument("--out", help="write the JSON result here instead of stdout")
    p.add_argument("--tsv-out", help="write the per-sample pair table (TSV)")
    p.add_argument("--summary-tsv", help="write the across-sample summary (TSV)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        result = run(args)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    payload = {
        "distance_unit": result["distance_unit"],
        "parameters": result["parameters"],
        "samples": result["samples"],
        "across_samples": [_rounded(r) for r in result["across_samples"]],
        "pairs": [_rounded(r) for r in result["pairs"]],
        "notes": result["notes"],
    }
    for r in payload["across_samples"]:
        r["unit_effects"] = {u: sig(v) for u, v in r["unit_effects"].items()}
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    try:
        if args.out:
            with open(args.out, "w", encoding="utf-8") as fh:
                fh.write(text + "\n")
        else:
            print(text)
        if args.tsv_out:
            write_tsv(args.tsv_out, PAIR_COLUMNS, result["pairs"])
        if args.summary_tsv:
            write_tsv(args.summary_tsv, SUMMARY_COLUMNS, result["across_samples"])
    except OSError as exc:
        print(f"error: cannot write output: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
