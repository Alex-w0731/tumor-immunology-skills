#!/usr/bin/env python3
"""Compare per-sample cell-type composition between two groups, with the sample as the unit.

For each cell type the per-sample value (default: cells of that type / all cells of the
sample, or of a declared parent population with --within) is compared between two groups:

  unpaired  Wilcoxon rank-sum (Mann-Whitney) test [Wilcoxon1945] [Mann1947].
  paired    (--pair-col, e.g. the same patient before and after treatment) Wilcoxon
            signed-rank test on within-pair differences; zero differences are dropped
            from the test as in [Wilcoxon1945].

Exact p-values: the permutation distribution of the rank statistic is enumerated with
midranks for ties, conditional on the observed tie pattern, and the two-sided p-value is
P(|T - E[T]| >= |t_obs - E[T]|). This is used when the comparison has at most --exact-max
samples (unpaired) or non-zero pairs (paired); above that, a normal approximation with tie
and continuity corrections is used. p_floor is the smallest two-sided p-value the exact test
could return for that design and tie pattern: if it exceeds your alpha, the test cannot
reject, whatever the data.

Effect sizes, always case minus control:
  unpaired  rank_biserial = (U_case - U_control) / (n_case * n_control)  [Kerby2014]
            hl_shift = median of all case - control differences          [Hodges1963]
  paired    rank_biserial = (R+ - R-) / (R+ + R-), zeros dropped           [Kerby2014]
            hl_shift = median of the Walsh averages of all differences,
            zeros included                                                [Hodges1963]
Benjamini-Hochberg adjustment [Benjamini1995] is applied across all tested cell types.

With --transform alr --reference TYPE, each other cell type is analysed as
log2((n_type + 0.5) / (n_reference + 0.5)), so a change is read relative to the reference
cell type rather than to the whole sample; the 0.5 offset follows the logit option of
propeller [Phipson2022]. Rank tests are unchanged by logit or arcsine transforms of a
proportion, so those are not offered here.

Input: TSV or CSV (delimiter taken from the header line). Either one row per cell (sample
and cell-type columns) or one row per sample and cell type with a count column
(--count-col), for example the metadata table written by pseudobulk_aggregate.py. Group,
pair and patient labels must be constant within a sample. A sample with no cells of a type
contributes a value of 0 for that type; it is not treated as missing.

Usage:
  python3 composition_test.py cells.tsv --sample-col sample_id --celltype-col cell_type \\
      --group-col response --case R --control NR --patient-col patient_id > composition.json
  python3 composition_test.py counts.tsv --count-col n_cells --group-col timepoint \\
      --case post --control pre --pair-col patient_id --tsv-out composition.tsv

Exit codes: 0 success, 2 input or design problem (the message names the column or value).
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import sys
from math import comb

DEFAULT_EXACT_MAX = 50
DIGITS = 12  # values are rounded before ranking so floating-point noise cannot split real ties


class InputError(ValueError):
    """An input or design problem that the user must fix (exit code 2)."""


# ----------------------------------------------------------------------------- statistics

def midranks(values):
    """1-based ranks in input order; tied values share the average of their ranks."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def tie_sizes(values):
    """Sizes of the groups of tied values (only groups larger than one)."""
    counts = {}
    for v in values:
        counts[v] = counts.get(v, 0) + 1
    return [c for c in counts.values() if c > 1]


def median(values):
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise ValueError("median of an empty list")
    mid = n // 2
    return ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) / 2.0


def _doubled(ranks):
    # Midranks are multiples of 0.5, so doubling gives exact integer scores.
    return [int(round(2 * r)) for r in ranks]


def _two_sided(dist, observed, mean2, total):
    """Two-sided p and p_floor from a {statistic: count} null distribution."""
    deviation = abs(observed - mean2)
    extreme = sum(c for t, c in dist.items() if abs(t - mean2) >= deviation)
    max_dev = max(abs(t - mean2) for t in dist)
    floor = sum(c for t, c in dist.items() if abs(t - mean2) == max_dev)
    return min(1.0, extreme / total), floor / total


def ranksum_exact(x, y):
    """Exact two-sided permutation p-value and p_floor of the Wilcoxon rank-sum test.

    Subset-sum dynamic programme over doubled midranks; the distribution of |T - E[T]| is
    the same whichever group is enumerated, so the smaller one is used.
    """
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        raise ValueError("both groups need at least one value")
    scores = _doubled(midranks(list(x) + list(y)))
    big_n = n1 + n2
    k = min(n1, n2)
    observed = sum(scores[:n1]) if k == n1 else sum(scores[n1:])
    layers = [dict() for _ in range(k + 1)]
    layers[0][0] = 1
    for idx, s in enumerate(scores):
        for j in range(min(k, idx + 1), 0, -1):
            current = layers[j]
            for t, c in layers[j - 1].items():
                current[t + s] = current.get(t + s, 0) + c
    # E[sum of k doubled midranks] = k * (N + 1), because the midranks always sum to N(N+1)/2.
    return _two_sided(layers[k], observed, k * (big_n + 1), comb(big_n, k))


def ranksum_normal(x, y):
    """Normal approximation to the rank-sum test with tie and continuity corrections."""
    n1, n2 = len(x), len(y)
    pooled = list(x) + list(y)
    big_n = n1 + n2
    ranks = midranks(pooled)
    t = sum(ranks[:n1])
    mean = n1 * (big_n + 1) / 2.0
    tie_term = sum(s ** 3 - s for s in tie_sizes(pooled)) / (big_n * (big_n - 1))
    var = n1 * n2 / 12.0 * ((big_n + 1) - tie_term)
    if var <= 0:
        return 1.0
    z = max(0.0, abs(t - mean) - 0.5) / math.sqrt(var)
    return min(1.0, math.erfc(z / math.sqrt(2.0)))


def signedrank_exact(diffs):
    """Exact two-sided p-value and p_floor of the signed-rank test; diffs must be non-zero."""
    n = len(diffs)
    scores = _doubled(midranks([abs(d) for d in diffs]))
    observed = sum(s for s, d in zip(scores, diffs) if d > 0)
    dist = {0: 1}
    for s in scores:
        nxt = dict(dist)
        for t, c in dist.items():
            nxt[t + s] = nxt.get(t + s, 0) + c
        dist = nxt
    # Doubled midranks sum to n(n+1), so E[W+] in doubled units is n(n+1)/2, an integer.
    return _two_sided(dist, observed, sum(scores) // 2, 2 ** n)


def signedrank_normal(diffs):
    """Normal approximation to the signed-rank test with tie and continuity corrections."""
    n = len(diffs)
    absd = [abs(d) for d in diffs]
    ranks = midranks(absd)
    w = sum(r for r, d in zip(ranks, diffs) if d > 0)
    mean = n * (n + 1) / 4.0
    var = n * (n + 1) * (2 * n + 1) / 24.0 - sum(s ** 3 - s for s in tie_sizes(absd)) / 48.0
    if var <= 0:
        return 1.0
    z = max(0.0, abs(w - mean) - 0.5) / math.sqrt(var)
    return min(1.0, math.erfc(z / math.sqrt(2.0)))


def bh_adjust(pvalues):
    """Benjamini-Hochberg step-up adjusted p-values, returned in input order."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adjusted = [1.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        adjusted[i] = min(1.0, running)
    return adjusted


def walsh_averages(values):
    n = len(values)
    return [(values[i] + values[j]) / 2.0 for i in range(n) for j in range(i, n)]


def _use_exact(method, size, exact_max):
    return method == "exact" or (method == "auto" and size <= exact_max)


def compare_unpaired(case, control, method="auto", exact_max=DEFAULT_EXACT_MAX):
    """Rank-sum comparison of two independent groups of per-sample values."""
    case = [round(v, DIGITS) for v in case]
    control = [round(v, DIGITS) for v in control]
    n1, n2 = len(case), len(control)
    if n1 == 0 or n2 == 0:
        raise InputError("each group needs at least one sample")
    u_case = sum(1.0 if a > b else 0.5 if a == b else 0.0 for a in case for b in control)
    exact = _use_exact(method, n1 + n2, exact_max)
    if exact:
        p, floor = ranksum_exact(case, control)
    else:
        p, floor = ranksum_normal(case, control), None
    return {
        "n_case": n1,
        "n_control": n2,
        "median_case": median(case),
        "median_control": median(control),
        "hl_shift": median([a - b for a in case for b in control]),
        "rank_biserial": 2.0 * u_case / (n1 * n2) - 1.0,
        "statistic": "U_case",
        "statistic_value": u_case,
        "p_value": p,
        "p_floor": floor,
        "method": "exact" if exact else "normal",
    }


def compare_paired(case, control, method="auto", exact_max=DEFAULT_EXACT_MAX):
    """Signed-rank comparison of matched values (case[i] and control[i] share a pair)."""
    if len(case) != len(control) or not case:
        raise InputError("a paired comparison needs at least one complete pair")
    case = [round(v, DIGITS) for v in case]
    control = [round(v, DIGITS) for v in control]
    diffs = [round(a - b, DIGITS) for a, b in zip(case, control)]
    nonzero = [d for d in diffs if d != 0]
    up = sum(1 for d in diffs if d > 0)
    down = sum(1 for d in diffs if d < 0)
    result = {
        "n_pairs": len(diffs),
        "n_up": up,
        "n_down": down,
        "n_zero": len(diffs) - up - down,
        "median_case": median(case),
        "median_control": median(control),
        "median_difference": median(diffs),
        "hl_shift": median(walsh_averages(diffs)),
        "statistic": "W_plus",
    }
    if not nonzero:
        result.update(rank_biserial=0.0, statistic_value=0.0, p_value=1.0, p_floor=1.0, method="exact")
        return result
    ranks = midranks([abs(d) for d in nonzero])
    r_plus = sum(r for r, d in zip(ranks, nonzero) if d > 0)
    r_minus = sum(r for r, d in zip(ranks, nonzero) if d < 0)
    exact = _use_exact(method, len(nonzero), exact_max)
    if exact:
        p, floor = signedrank_exact(nonzero)
    else:
        p, floor = signedrank_normal(nonzero), None
    result.update(
        rank_biserial=(r_plus - r_minus) / (r_plus + r_minus),
        statistic_value=r_plus,
        p_value=p,
        p_floor=floor,
        method="exact" if exact else "normal",
    )
    return result


# ----------------------------------------------------------------------------- input

def read_rows(path):
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            text = fh.read()
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc}") from exc
    first = text.splitlines()[0] if text.strip() else ""
    if not first.strip():
        raise InputError(f"{path} is empty")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t" if "\t" in first else ",")
    rows = []
    for r in reader:
        extra = r.get(None)
        named = [v for k, v in r.items() if k is not None]
        if not extra and not any((v or "").strip() for v in named):
            continue  # blank line
        if extra is not None:
            raise InputError(f"{path}, line {reader.line_num}: more fields than the header has "
                             "(a stray delimiter or an unquoted comma in a value?)")
        if any(v is None for v in named):
            raise InputError(f"{path}, line {reader.line_num}: fewer fields than the header has")
        rows.append(r)
    fields = [f.strip() for f in (reader.fieldnames or [])]
    if not rows:
        raise InputError(f"{path} has a header but no data rows")
    rows = [{(k or "").strip(): v for k, v in r.items()} for r in rows]
    return fields, rows


def parse_count(raw, where):
    try:
        value = float(raw)
    except ValueError:
        raise InputError(f"{where}: {raw!r} is not a number") from None
    if not math.isfinite(value) or value < 0 or value != int(value):
        raise InputError(f"{where}: {raw!r} must be a non-negative whole number of cells")
    return int(value)


def collect(rows, fields, sample_col, celltype_col, group_col, count_col=None,
            pair_col=None, patient_col=None):
    """Return ({sample: {cell_type: n_cells}}, {sample: {"group", "pair", "patient"}})."""
    needed = [c for c in (sample_col, celltype_col, group_col, count_col, pair_col, patient_col) if c]
    missing = [c for c in needed if c not in fields]
    if missing:
        raise InputError(f"column(s) not found: {', '.join(missing)}; available: {', '.join(fields)}")
    label_cols = [("group", group_col)]
    if pair_col:
        label_cols.append(("pair", pair_col))
    if patient_col:
        label_cols.append(("patient", patient_col))
    counts, labels, seen = {}, {}, set()
    for line, row in enumerate(rows, start=2):
        sample = (row.get(sample_col) or "").strip()
        celltype = (row.get(celltype_col) or "").strip()
        if not sample or not celltype:
            raise InputError(f"line {line}: empty value in '{sample_col}' or '{celltype_col}'")
        lab = {}
        for name, col in label_cols:
            value = (row.get(col) or "").strip()
            if not value:
                raise InputError(f"line {line}: empty value in '{col}'")
            lab[name] = value
        previous = labels.setdefault(sample, lab)
        if previous != lab:
            name = next(k for k in lab if lab[k] != previous[k])
            raise InputError(
                f"sample '{sample}' has more than one value in '{dict(label_cols)[name]}' "
                f"({previous[name]!r} and {lab[name]!r}); labels must be constant within a sample")
        if count_col:
            if (sample, celltype) in seen:
                raise InputError(f"line {line}: sample '{sample}' and cell type '{celltype}' appear "
                                 "more than once; aggregate to one row per sample and cell type first")
            seen.add((sample, celltype))
            n = parse_count((row.get(count_col) or "").strip(), f"line {line}, column '{count_col}'")
        else:
            n = 1
        per_sample = counts.setdefault(sample, {})
        per_sample[celltype] = per_sample.get(celltype, 0) + n
    return counts, labels


def build_design(labels, case, control, paired=False, use_patient=False, exclude=()):
    unknown = sorted(set(exclude) - set(labels))
    if unknown:
        raise InputError(f"samples to exclude not found: {', '.join(unknown)}")
    groups = sorted({lab["group"] for lab in labels.values()})
    if case == control:
        raise InputError("--case and --control must be different groups")
    for g in (case, control):
        if g not in groups:
            raise InputError(f"group '{g}' not found; groups present: {', '.join(groups)}")
    active = {s: lab for s, lab in labels.items()
              if lab["group"] in (case, control) and s not in set(exclude)}
    other = sorted(s for s, lab in labels.items() if lab["group"] not in (case, control))
    design = {"other_group_samples": other, "excluded_samples": sorted(set(exclude))}
    if paired:
        by_pair = {}
        for s, lab in sorted(active.items()):
            by_pair.setdefault(lab["pair"], {case: [], control: []})[lab["group"]].append(s)
        pairs, incomplete = [], []
        for pid, members in sorted(by_pair.items()):
            for g in (case, control):
                if len(members[g]) > 1:
                    raise InputError(
                        f"pair '{pid}' has {len(members[g])} samples in group '{g}' "
                        f"({', '.join(members[g])}); keep one per group or merge them before a paired test")
            if members[case] and members[control]:
                pairs.append((pid, members[case][0], members[control][0]))
            else:
                incomplete.append(pid)
        if not pairs:
            raise InputError("no complete pairs: each pair needs one sample in each group")
        design.update(mode="paired", pairs=pairs, incomplete_pairs=incomplete)
        return design
    case_samples = sorted(s for s, lab in active.items() if lab["group"] == case)
    control_samples = sorted(s for s, lab in active.items() if lab["group"] == control)
    if not case_samples or not control_samples:
        raise InputError("after exclusions each group needs at least one sample")
    if use_patient:
        by_patient = {}
        for s in case_samples + control_samples:
            by_patient.setdefault(active[s]["patient"], []).append(s)
        for pid, members in sorted(by_patient.items()):
            if len(members) > 1:
                member_groups = sorted({active[s]["group"] for s in members})
                if len(member_groups) > 1:
                    raise InputError(
                        f"patient '{pid}' has samples in both groups ({', '.join(members)}): this is a "
                        "paired design; rerun with --pair-col set to the patient column")
                raise InputError(
                    f"patient '{pid}' contributes {len(members)} samples to group '{member_groups[0]}' "
                    f"({', '.join(members)}); the patient is the unit of replication: keep one sample "
                    "per patient by a prespecified rule, or merge them, then rerun")
    design.update(mode="unpaired", case_samples=case_samples, control_samples=control_samples)
    return design


# ----------------------------------------------------------------------------- analysis

def analyse(counts, labels, case, control, *, paired=False, use_patient=False, within=None,
            transform="proportion", reference=None, method="auto", exact_max=DEFAULT_EXACT_MAX,
            min_sample_cells=0, exclude=()):
    all_types = sorted({ct for per in counts.values() for ct in per})
    if within:
        absent = [t for t in within if t not in all_types]
        if absent:
            raise InputError(f"--within cell type(s) not in the data: {', '.join(absent)}; "
                             f"present: {', '.join(all_types)}")
        denominator = list(dict.fromkeys(within))
    else:
        denominator = all_types
    if transform == "alr":
        if not reference:
            raise InputError("--transform alr needs --reference (a cell type expected not to change)")
        if reference not in all_types:
            raise InputError(f"reference cell type '{reference}' not in the data; present: {', '.join(all_types)}")
        tested_types = [t for t in denominator if t != reference]
    else:
        tested_types = list(denominator)
    totals = {s: sum(counts[s].get(t, 0) for t in denominator) for s in counts}
    low = sorted(s for s in totals if totals[s] < max(1, min_sample_cells) and s not in set(exclude)
                 and labels[s]["group"] in (case, control))
    design = build_design(labels, case, control, paired, use_patient, tuple(exclude) + tuple(low))
    design["low_total_samples_excluded"] = low

    def value(s, t):
        n = counts[s].get(t, 0)
        if transform == "alr":
            return math.log2((n + 0.5) / (counts[s].get(reference, 0) + 0.5))
        return n / totals[s]

    if design["mode"] == "paired":
        used = [s for _, a, b in design["pairs"] for s in (a, b)]
    else:
        used = design["case_samples"] + design["control_samples"]
    empty_types = [t for t in tested_types if sum(counts[s].get(t, 0) for s in used) == 0]
    results = []
    for t in tested_types:
        if t in empty_types:
            continue
        if design["mode"] == "paired":
            res = compare_paired([value(a, t) for _, a, _ in design["pairs"]],
                                 [value(b, t) for _, _, b in design["pairs"]], method, exact_max)
        else:
            res = compare_unpaired([value(s, t) for s in design["case_samples"]],
                                   [value(s, t) for s in design["control_samples"]], method, exact_max)
        res["cell_type"] = t
        results.append(res)
    for res, q in zip(results, bh_adjust([r["p_value"] for r in results])):
        res["q_bh"] = q
    results.sort(key=lambda r: (r["p_value"], r["cell_type"]))

    notes = []
    floors = [r["p_floor"] for r in results if r["p_floor"] is not None]
    if floors and min(floors) > 0.05:
        notes.append(f"With this design no cell type can reach a two-sided exact p below 0.05 (smallest "
                     f"attainable p = {min(floors):.4g}). Report effect sizes and direction counts; a large "
                     "p-value here is not evidence of no change.")
    if transform == "proportion":
        notes.append("Proportions are relative: a rise in one cell type lowers the proportions of the others. "
                     "Read decreases with care, or rerun with --transform alr and a reference cell type.")
    if design["mode"] == "paired" and design["incomplete_pairs"]:
        notes.append("Left out of the paired test (no sample in one of the groups): "
                     + ", ".join(design["incomplete_pairs"]))
    if design["mode"] == "unpaired" and min(len(design["case_samples"]), len(design["control_samples"])) < 2:
        notes.append("One group has a single sample, so between-sample variation in that group is unknown.")
    if low:
        notes.append("Excluded for low total cell count: " + ", ".join(low))
    if empty_types:
        notes.append("Not tested (no cells in any compared sample): " + ", ".join(empty_types))
    if design["other_group_samples"]:
        notes.append(f"{len(design['other_group_samples'])} sample(s) in other groups were ignored.")

    sample_rows = []
    for s in used:
        row = {"sample": s, "group": labels[s]["group"], "total_cells": totals[s]}
        if "pair" in labels[s]:
            row["pair"] = labels[s]["pair"]
        row["values"] = {t: value(s, t) for t in tested_types if t not in empty_types}
        sample_rows.append(row)
    design_out = {k: v for k, v in design.items() if k not in ("case_samples", "control_samples")}
    if design["mode"] == "paired":
        design_out["pairs"] = [{"pair": p, "case_sample": a, "control_sample": b} for p, a, b in design["pairs"]]
        design_out["n_pairs"] = len(design["pairs"])
    else:
        design_out["n_case"] = len(design["case_samples"])
        design_out["n_control"] = len(design["control_samples"])
    design_out.update(case=case, control=control, transform=transform, reference=reference,
                      denominator="all cell types" if not within else ",".join(denominator),
                      fdr_family=f"Benjamini-Hochberg across {len(results)} tested cell types")
    return {"design": design_out, "results": results, "samples": sample_rows, "notes": notes}


# ----------------------------------------------------------------------------- output

def _clean(obj):
    if isinstance(obj, float):
        return round(obj, 10)
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    return obj


def _fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def results_tsv(report):
    paired = report["design"]["mode"] == "paired"
    cols = ["cell_type"] + (["n_pairs", "n_up", "n_down", "n_zero", "median_case", "median_control",
                             "median_difference"] if paired else
                            ["n_case", "n_control", "median_case", "median_control"])
    cols += ["hl_shift", "rank_biserial", "statistic", "statistic_value", "p_value", "p_floor", "q_bh", "method"]
    lines = ["\t".join(cols)]
    for r in report["results"]:
        lines.append("\t".join(_fmt(r.get(c)) for c in cols))
    return "\n".join(lines) + "\n"


def values_tsv(report):
    types = sorted({t for row in report["samples"] for t in row["values"]})
    lines = ["\t".join(["sample", "group", "pair", "total_cells"] + types)]
    for row in report["samples"]:
        lines.append("\t".join([row["sample"], row["group"], row.get("pair", ""), str(row["total_cells"])]
                               + [_fmt(row["values"][t]) for t in types]))
    return "\n".join(lines) + "\n"


def build_parser():
    p = argparse.ArgumentParser(
        description="Exact small-sample comparison of per-sample cell-type proportions between two groups "
                    "(Wilcoxon rank-sum or signed-rank), with rank-biserial and Hodges-Lehmann effect sizes "
                    "and Benjamini-Hochberg adjustment across cell types.",
        epilog="Example: composition_test.py cells.tsv --group-col response --case R --control NR",
    )
    p.add_argument("input", help="TSV or CSV: one row per cell, or one row per sample x cell type with --count-col")
    p.add_argument("--sample-col", default="sample_id", help="sample column, the unit of replication (default: sample_id)")
    p.add_argument("--celltype-col", default="cell_type", help="cell-type column (default: cell_type)")
    p.add_argument("--group-col", required=True, help="column holding the group label of each sample")
    p.add_argument("--case", required=True, help="group whose values come first in every difference")
    p.add_argument("--control", required=True, help="reference group")
    p.add_argument("--count-col", help="column with the number of cells (rows are then sample x cell type)")
    p.add_argument("--pair-col", help="paired design: column identifying the pair, usually the patient")
    p.add_argument("--patient-col", help="unpaired design: check that each patient contributes one sample")
    p.add_argument("--within", help="comma-separated cell types forming the denominator (parent population)")
    p.add_argument("--transform", choices=["proportion", "alr"], default="proportion",
                   help="proportion (default) or alr: log2 ratio to --reference with a 0.5 offset")
    p.add_argument("--reference", help="reference cell type for --transform alr")
    p.add_argument("--method", choices=["auto", "exact", "normal"], default="auto",
                   help="auto (default): exact up to --exact-max samples or non-zero pairs, normal above")
    p.add_argument("--exact-max", type=int, default=DEFAULT_EXACT_MAX,
                   help=f"largest design run exactly under --method auto (default: {DEFAULT_EXACT_MAX})")
    p.add_argument("--min-sample-cells", type=int, default=0,
                   help="exclude samples with fewer cells than this in the denominator (default: 0)")
    p.add_argument("--exclude-samples", help="comma-separated samples to leave out (sensitivity analysis)")
    p.add_argument("--out", help="write the JSON report here instead of stdout")
    p.add_argument("--tsv-out", help="also write the per-cell-type results as TSV")
    p.add_argument("--proportions-out", help="also write the per-sample values (proportions or log-ratios) as TSV")
    return p


def _split(text):
    return [t.strip() for t in text.split(",") if t.strip()] if text else []


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.exact_max < 1:
            raise InputError("--exact-max must be at least 1")
        if args.reference and args.transform != "alr":
            raise InputError("--reference is only used with --transform alr")
        fields, rows = read_rows(args.input)
        counts, labels = collect(rows, fields, args.sample_col, args.celltype_col, args.group_col,
                                 args.count_col, args.pair_col, args.patient_col)
        report = analyse(counts, labels, args.case, args.control, paired=bool(args.pair_col),
                         use_patient=bool(args.patient_col), within=_split(args.within),
                         transform=args.transform, reference=args.reference, method=args.method,
                         exact_max=args.exact_max, min_sample_cells=args.min_sample_cells,
                         exclude=tuple(_split(args.exclude_samples)))
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    report["design"]["unit_of_replication"] = args.sample_col
    text = json.dumps(_clean(report), indent=2, ensure_ascii=False) + "\n"
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)
    if args.tsv_out:
        with open(args.tsv_out, "w", encoding="utf-8") as fh:
            fh.write(results_tsv(report))
    if args.proportions_out:
        with open(args.proportions_out, "w", encoding="utf-8") as fh:
            fh.write(values_tsv(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
