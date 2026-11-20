#!/usr/bin/env python3
r"""Pairwise overlap between T-cell receptor repertoires: Jaccard, Morisita-Horn, shared fractions.

Reads the same clonotype tables, with the same filters (productive rows, one receptor chain,
clonotype key, count or cell column), as repertoire_metrics.py, which must sit in the same
folder: this script loads it for reading, filtering and downsampling. Two or more samples
come from a sample column or from one file per sample.

For samples A and B with clonotype counts a_i and b_i, totals N_A and N_B, and proportions
p_i = a_i / N_A, q_i = b_i / N_B:
  jaccard            shared clonotypes / clonotypes observed in A or B         [Jaccard1912]
  morisita_horn      2 sum p_i q_i / (sum p_i^2 + sum q_i^2)                   [Horn1966]
                     (abundance-weighted; dominated by the largest clonotypes)
  shared_fraction_a  share of N_A in clonotypes also observed in B; shared_fraction_b likewise
Jaccard and the shared fractions change with depth (Jaccard in either direction; a shared
fraction falls as the other sample gets shallower), so pairs of unequal depth are compared
with --downsample: each sample is subsampled without replacement to the same depth, the
overlap recomputed, and this repeated --iterations times from --seed (mean, SD, 2.5 / 97.5
percentiles). "--downsample min" takes the smallest depth among the samples named in the
requested pairs (all samples when --pairs is not given). Source keys resolve in
references/sources.md of the tcr-repertoire-metrics skill.

Examples:
  python3 repertoire_overlap.py --input tumour_blood.tsv --sample-col sample \
      --clone-col junction_aa,v_call --count-col duplicate_count --chain TRB
  python3 repertoire_overlap.py --input pre.tsv post.tsv --clone-col cdr3_aa --count-col count \
      --downsample min --iterations 50 --seed 3 --format tsv
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import os
import random
import sys


def _load_metrics_module():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "repertoire_metrics.py")
    if not os.path.isfile(path):
        sys.exit("error: repertoire_overlap.py needs repertoire_metrics.py in the same folder")
    spec = importlib.util.spec_from_file_location("repertoire_metrics_shared", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rm = _load_metrics_module()
InputError = rm.InputError
OVERLAP_KEYS = ("richness_a", "richness_b", "shared_clonotypes", "jaccard", "morisita_horn",
                "shared_fraction_a", "shared_fraction_b")


def overlap_metrics(a, b):
    """Overlap of two {clonotype: count} mappings (zero counts ignored)."""
    keys_a = {k for k, v in a.items() if v > 0}
    keys_b = {k for k, v in b.items() if v > 0}
    if not keys_a or not keys_b:
        raise InputError("both samples need at least one clonotype with a positive count")
    shared = keys_a & keys_b
    n_a = float(sum(a[k] for k in keys_a))
    n_b = float(sum(b[k] for k in keys_b))
    sum_pq = sum(a[k] * b[k] for k in shared) / (n_a * n_b)
    sum_p2 = sum(a[k] ** 2 for k in keys_a) / n_a ** 2
    sum_q2 = sum(b[k] ** 2 for k in keys_b) / n_b ** 2
    return {
        "richness_a": len(keys_a),
        "richness_b": len(keys_b),
        "shared_clonotypes": len(shared),
        "jaccard": len(shared) / len(keys_a | keys_b),
        "morisita_horn": 2.0 * sum_pq / (sum_p2 + sum_q2),
        "shared_fraction_a": sum(a[k] for k in shared) / n_a,
        "shared_fraction_b": sum(b[k] for k in shared) / n_b,
    }


def parse_pairs(text, names):
    """'A:B,A:C' -> [('A', 'B'), ('A', 'C')]; default all unordered pairs in input order."""
    if not text:
        return list(itertools.combinations(names, 2))
    pairs = []
    for item in text.split(","):
        parts = [p.strip() for p in item.split(":")]
        if len(parts) != 2 or not all(parts) or parts[0] == parts[1]:
            raise InputError(f"--pairs entries look like A:B with two different samples, not '{item}'")
        for p in parts:
            if p not in names:
                raise InputError(f"--pairs names sample '{p}', which is not in the input ({', '.join(names)})")
        if (parts[0], parts[1]) not in pairs:   # a repeated pair is reported once
            pairs.append((parts[0], parts[1]))
    return pairs


def analyse_overlap(samples, qc, notes, pairs=None, downsample=None, iterations=20, seed=1):
    names = list(samples)
    if len(names) < 2:
        raise InputError(f"overlap needs at least two samples; found {len(names)} ({', '.join(names)})")
    if iterations < 1:
        raise InputError("--iterations must be at least 1")
    pair_list = parse_pairs(pairs, names)
    in_pairs = {n for p in pair_list for n in p}
    # the common depth (and the integer-count check) concerns only the samples being compared
    target = rm.resolve_depth(downsample, {n: samples[n] for n in names if n in in_pairs})
    notes = list(notes)
    depths = {n: sum(samples[n].values()) for n in names}
    runs = {}
    if target is not None:
        eligible = [n for n in names if depths[n] >= target]
        for n in names:
            if n in in_pairs and n not in eligible:
                notes.append(f"sample '{n}' (depth {depths[n]}) is below the downsampling depth {target}; its pairs are not downsampled")
        rngs = {n: random.Random(f"{seed}:{n}") for n in eligible}
        keys = {n: list(samples[n]) for n in eligible}
        runs = {p: [] for p in pair_list if p[0] in rngs and p[1] in rngs}
        needed = sorted({n for p in runs for n in p}, key=names.index)
        for _ in range(iterations):
            subs = {}
            for n in needed:
                drawn = rm.subsample([samples[n][k] for k in keys[n]], target, rngs[n])
                subs[n] = {k: c for k, c in zip(keys[n], drawn) if c}
            for p in runs:
                runs[p].append(overlap_metrics(subs[p[0]], subs[p[1]]))
    results = []
    for a, b in pair_list:
        entry = {"sample_a": a, "sample_b": b, "depth_a": depths[a], "depth_b": depths[b],
                 "observed": overlap_metrics(samples[a], samples[b])}
        ratio = max(depths[a], depths[b]) / min(depths[a], depths[b])
        if (a, b) in runs:
            entry["downsampled"] = {"depth": target, "iterations": iterations,
                                    "metrics": {k: rm.summarize([r[k] for r in runs[(a, b)]]) for k in OVERLAP_KEYS}}
        elif ratio > 2:
            notes.append(f"{a} vs {b}: depths differ {ratio:.1f}-fold; Jaccard and shared fractions are "
                         "depth-dependent, so use --downsample before comparing pairs")
        results.append(entry)
    return {"tool": "repertoire_overlap.py",
            "settings": {"downsample_depth": target, "iterations": iterations if target is not None else None,
                         "seed": seed if target is not None else None},
            "samples": [{"sample": n, "depth": depths[n], "clonotypes": len(samples[n]),
                         "chain": qc[n].get("chain")} for n in names],
            "pairs": results, "notes": notes}


def to_tsv(result):
    downsampled = any("downsampled" in e for e in result["pairs"])
    header = ["sample_a", "sample_b", "depth_a", "depth_b"] + list(OVERLAP_KEYS)
    if downsampled:
        header += ["ds_depth", "ds_iterations"] + [f"ds_{k}_{s}" for k in OVERLAP_KEYS for s in ("mean", "sd")]
    lines = ["\t".join(header)]
    for e in result["pairs"]:
        row = [e["sample_a"], e["sample_b"], e["depth_a"], e["depth_b"]] + [e["observed"][k] for k in OVERLAP_KEYS]
        if downsampled:
            ds = e.get("downsampled")
            row += [ds["depth"] if ds else None, ds["iterations"] if ds else None]
            for k in OVERLAP_KEYS:
                row += [ds["metrics"][k]["mean"] if ds else None, ds["metrics"][k]["sd"] if ds else None]
        lines.append("\t".join(rm._cell(v) for v in row))
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Pairwise Jaccard, Morisita-Horn and shared-count fractions between TCR repertoires, "
                    "with optional seeded downsampling to a common depth.",
        epilog="Exit code 2 means the input or the arguments need fixing; the message says which.")
    rm.add_input_arguments(parser)
    parser.add_argument("--pairs", help="restrict to these pairs, e.g. tumour:blood,tumour:normal (default: all pairs)")
    args = parser.parse_args(argv)
    try:
        samples, qc, notes = rm.load_from_args(args)
        result = analyse_overlap(samples, qc, notes, pairs=args.pairs, downsample=args.downsample,
                                 iterations=args.iterations, seed=args.seed)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.format == "tsv":
        for note in result["notes"]:
            print(f"note: {note}", file=sys.stderr)
        text = to_tsv(result)
    else:
        text = json.dumps(rm._round(result), indent=2) + "\n"
    rm.write_output(text, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
