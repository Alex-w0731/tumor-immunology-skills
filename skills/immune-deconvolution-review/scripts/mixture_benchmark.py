#!/usr/bin/env python3
"""Compare estimated cell-type fractions with known fractions (mixture benchmark).

Use on mixtures of known composition (in silico pseudo-bulk, sorted-cell or cell-line
mixtures) or on samples with an orthogonal measurement (flow cytometry, IHC counts). Both
inputs are samples x cell types tables (first column: sample ID). Only samples and cell
types present in both tables are compared; extra columns (for example coef_sum or rmse from
nnls_deconvolve.py --format tsv) are ignored and listed. 'NA' or an empty cell is missing
and is dropped pairwise; an infinite value is an error.

Metrics:
  Per cell type, across samples (can the method compare samples for this cell type?):
    n, Pearson r, Spearman rho (average ranks for ties), RMSE = sqrt(mean((est - true)^2)),
    bias = mean(est - true), OLS slope and intercept of est on true, and the mean estimate
    in samples where the truth is 0 (background prediction, the quantity studied in
    [Sturm2019]).
  Per sample, across cell types (is the composition within a sample right?):
    Pearson r, RMSE, and the sums of true and estimated values over the compared cell types.
  Overall: Pearson r and RMSE pooled over all pairs, median per-sample r.

Correlation needs at least --min-n pairs (default 3) and is null for a constant vector.
--merge sums cell types into a group in both tables before comparison (for collinear
types). --renormalize rescales rows to sum to 1 over the compared cell types, for example
to compare fractions relative to immune cells with a truth relative to all cells; rows with
a missing value are rescaled over fewer cell types and are flagged in the warnings.

Usage examples:
  python3 mixture_benchmark.py --estimated fractions.tsv --truth truth.tsv > benchmark.json
  python3 mixture_benchmark.py --estimated fractions.tsv --truth truth.tsv \\
      --merge "Macrophage=Macrophage_M1+Macrophage_M2" --renormalize truth --format tsv

Exit code 2 on invalid input, with a message naming the file, sample or column.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys

MISSING = {"", "NA", "NAN", "NULL", "NONE"}


class InputError(ValueError):
    """Invalid input; main() reports it and exits with code 2."""


def read_table(path: str, label: str):
    """Read a samples x columns table. Returns (sample ids, columns, {sample: {col: value}})."""
    delim = "," if path.lower().endswith(".csv") else "\t"
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            rows = [r for r in csv.reader(fh, delimiter=delim) if r and any(c.strip() for c in r)]
    except OSError as exc:
        raise InputError(f"{label}: cannot read {path}: {exc.strerror}") from None
    if len(rows) < 2:
        raise InputError(f"{label}: {path} needs a header row and at least one data row")
    header = [c.strip() for c in rows[0]]
    cols = header[1:]
    if not cols or any(not c for c in cols):
        raise InputError(f"{label}: header needs non-empty column names after the sample column")
    dups = sorted({c for c in cols if cols.count(c) > 1})
    if dups:
        raise InputError(f"{label}: duplicated column names: {', '.join(dups)}")
    ids, data = [], {}
    for n, r in enumerate(rows[1:], start=2):
        if len(r) != len(header):
            raise InputError(f"{label}: line {n} has {len(r)} fields but the header has {len(header)}")
        sid = r[0].strip()
        if not sid or sid in data:
            raise InputError(f"{label}: line {n} has an empty or repeated sample ID '{sid}'")
        rec = {}
        for col, cell in zip(cols, r[1:]):
            cell = cell.strip()
            if cell.upper() in MISSING:
                rec[col] = None
                continue
            try:
                v = float(cell)
            except ValueError:
                raise InputError(f"{label}: non-numeric value '{cell}' for sample '{sid}', column '{col}'") from None
            if not math.isfinite(v):
                raise InputError(f"{label}: infinite value for sample '{sid}', column '{col}'")
            rec[col] = v
        ids.append(sid)
        data[sid] = rec
    return ids, cols, data


def parse_merge(specs):
    """Parse 'Group=A+B' strings into [(group, [A, B])]."""
    out = []
    for spec in specs or []:
        name, sep, rhs = spec.partition("=")
        parts = [p.strip() for p in rhs.split("+") if p.strip()]
        if not sep or not name.strip() or len(parts) < 2:
            raise InputError(f"--merge '{spec}' must look like Group=TypeA+TypeB")
        out.append((name.strip(), parts))
    return out


def apply_merge(cols, data, merges, label):
    """Replace component columns by their sum. A table that already has the group keeps it."""
    cols = list(cols)
    for name, parts in merges:
        present = [p for p in parts if p in cols]
        if not present:
            continue
        if len(present) != len(parts):
            absent = [p for p in parts if p not in cols]
            raise InputError(f"{label}: --merge {name} needs {', '.join(absent)} as well")
        if name in cols:
            raise InputError(f"{label}: --merge {name} would overwrite an existing column")
        for rec in data.values():
            vals = [rec[p] for p in parts]
            rec[name] = None if any(v is None for v in vals) else sum(vals)
            for p in parts:
                del rec[p]
        idx = min(cols.index(p) for p in parts)
        cols = [c for c in cols if c not in parts]
        cols.insert(min(idx, len(cols)), name)
    return cols


def ranks(values):
    """Ranks starting at 1, ties receive the average rank."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            out[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return out


def pearson(x, y, min_n: int = 3):
    """Pearson correlation; None for fewer than min_n pairs or a constant vector."""
    n = len(x)
    if n < max(2, min_n) or max(x) == min(x) or max(y) == min(y):
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy)))


def spearman(x, y, min_n: int = 3):
    if len(x) < max(2, min_n):
        return None
    return pearson(ranks(x), ranks(y), min_n)


def rmse(est, true):
    return math.sqrt(sum((e - t) ** 2 for e, t in zip(est, true)) / len(est)) if est else None


def ols(x, y):
    """Slope and intercept of y on x; (None, None) when x is constant or too short."""
    n = len(x)
    if n < 2 or max(x) == min(x):
        return None, None
    mx, my = sum(x) / n, sum(y) / n
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)
    return slope, my - slope * mx


def _round(v, d: int = 4):
    return None if v is None else round(v, d)


def renormalize(data, samples, cell_types):
    """Rescale each row to sum to 1 over cell_types. Returns the samples with a missing value,
    whose rows are rescaled over fewer cell types and so use a different denominator."""
    incomplete = []
    for s in samples:
        rec = data[s]
        vals = [rec.get(c) for c in cell_types]
        if any(v is None for v in vals):
            incomplete.append(s)
        total = sum(v for v in vals if v is not None)
        if total > 0:
            for c in cell_types:
                if rec.get(c) is not None:
                    rec[c] = rec[c] / total
    return incomplete


def benchmark(est_ids, est_cols, est, true_ids, true_cols, true, merges=(), renorm="none", min_n=3):
    """Compute the metrics described in the module docstring. Returns a result dict."""
    est_cols = apply_merge(est_cols, est, merges, "estimated")
    true_cols = apply_merge(true_cols, true, merges, "truth")
    cell_types = [c for c in true_cols if c in est_cols]
    samples = [s for s in true_ids if s in est]
    if not cell_types:
        raise InputError("no cell-type column is shared by the two tables (check names, or use --merge)")
    if not samples:
        raise InputError("no sample ID is shared by the two tables")
    warnings = []
    for label, table, on in (("truth", true, ("truth", "both")), ("estimated", est, ("estimated", "both"))):
        if renorm in on:
            incomplete = renormalize(table, samples, cell_types)
            if incomplete:
                warnings.append(f"{label}: {len(incomplete)} sample(s) with missing values were renormalized over "
                                f"fewer cell types (for example {incomplete[0]}), so their denominator differs")

    per_type, all_e, all_t, n_missing = [], [], [], 0
    for c in cell_types:
        pairs = [(est[s][c], true[s][c]) for s in samples if est[s].get(c) is not None and true[s].get(c) is not None]
        n_missing += len(samples) - len(pairs)
        e = [p[0] for p in pairs]
        t = [p[1] for p in pairs]
        all_e += e
        all_t += t
        slope, intercept = ols(t, e)
        zero = [ev for ev, tv in pairs if tv == 0]
        per_type.append({
            "cell_type": c, "n": len(pairs),
            "pearson_r": _round(pearson(t, e, min_n)), "spearman_rho": _round(spearman(t, e, min_n)),
            "rmse": _round(rmse(e, t), 5), "bias": _round(sum(ev - tv for ev, tv in pairs) / len(pairs), 5) if pairs else None,
            "mean_true": _round(sum(t) / len(t), 5) if t else None, "mean_est": _round(sum(e) / len(e), 5) if e else None,
            "slope": _round(slope), "intercept": _round(intercept, 5),
            "n_true_zero": len(zero), "mean_est_when_true_zero": _round(sum(zero) / len(zero), 5) if zero else None,
        })

    per_sample = []
    for s in samples:
        pairs = [(est[s][c], true[s][c]) for c in cell_types if est[s].get(c) is not None and true[s].get(c) is not None]
        e = [p[0] for p in pairs]
        t = [p[1] for p in pairs]
        per_sample.append({
            "sample": s, "n": len(pairs), "pearson_r": _round(pearson(t, e, min_n)), "rmse": _round(rmse(e, t), 5),
            "sum_true": _round(sum(t), 4), "sum_est": _round(sum(e), 4),
        })

    sample_r = [p["pearson_r"] for p in per_sample if p["pearson_r"] is not None]
    sums_t = [p["sum_true"] for p in per_sample if p["n"]]
    sums_e = [p["sum_est"] for p in per_sample if p["n"]]
    if sums_t and sums_e:
        med_t, med_e = statistics.median(sums_t), statistics.median(sums_e)
        if (abs(med_e - 1) < 0.02 and med_t < 0.9) or (abs(med_t - 1) < 0.02 and med_e < 0.9):
            warnings.append(
                f"median row sums differ (truth {med_t:.3f}, estimate {med_e:.3f}) over the compared cell types: "
                "one table is relative to these cell types and the other to all cells, so bias and RMSE mix "
                "accuracy with a change of denominator; compare with --renormalize truth or use absolute estimates")
    for rec in per_type:
        if rec["pearson_r"] is None:
            warnings.append(f"{rec['cell_type']}: correlation undefined (fewer than {min_n} pairs or no variation)")
        elif rec["pearson_r"] >= 0.8 and rec["mean_true"] and rec["bias"] is not None and abs(rec["bias"]) > 0.5 * rec["mean_true"]:
            warnings.append(f"{rec['cell_type']}: high correlation (r = {rec['pearson_r']}) but bias of "
                            f"{rec['bias']} against a mean of {rec['mean_true']}: ranks samples well, wrong scale")
    return {
        "n_samples": len(samples), "cell_types": cell_types,
        "unmatched_samples": {"estimated_only": [s for s in est_ids if s not in true],
                              "truth_only": [s for s in true_ids if s not in est]},
        "ignored_columns": {"estimated_only": [c for c in est_cols if c not in cell_types],
                            "truth_only": [c for c in true_cols if c not in cell_types]},
        "n_missing_pairs": n_missing, "renormalize": renorm, "merged": [m[0] for m in merges],
        "per_cell_type": per_type, "per_sample": per_sample,
        "overall": {"n_pairs": len(all_e), "pearson_r": _round(pearson(all_t, all_e, min_n)),
                    "rmse": _round(rmse(all_e, all_t), 5),
                    "median_sample_r": _round(statistics.median(sample_r)) if sample_r else None},
        "warnings": warnings,
    }


def to_tsv(result) -> str:
    keys = ["cell_type", "n", "pearson_r", "spearman_rho", "rmse", "bias", "mean_true", "mean_est",
            "slope", "intercept", "n_true_zero", "mean_est_when_true_zero"]
    lines = ["\t".join(keys)]
    for rec in result["per_cell_type"]:
        lines.append("\t".join("NA" if rec[k] is None else str(rec[k]) for k in keys))
    return "\n".join(lines) + "\n"


def build_parser():
    p = argparse.ArgumentParser(
        description="Benchmark estimated cell-type fractions against known fractions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Both tables: rows = samples, columns = cell types, first column = sample ID.")
    p.add_argument("--estimated", required=True, help="estimated fractions or scores (samples x cell types)")
    p.add_argument("--truth", required=True, help="known fractions (samples x cell types)")
    p.add_argument("--merge", action="append", metavar="GROUP=A+B",
                   help="sum cell types into a group in both tables (repeatable)")
    p.add_argument("--renormalize", choices=["none", "truth", "estimated", "both"], default="none",
                   help="rescale rows to sum to 1 over the compared cell types before comparing")
    p.add_argument("--min-n", type=int, default=3, help="minimum pairs for a correlation (default 3)")
    p.add_argument("--format", choices=["json", "tsv"], default="json",
                   help="json: full result; tsv: the per-cell-type table only")
    p.add_argument("--out", help="write the result to this file instead of stdout")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.min_n < 2:
            raise InputError("--min-n must be at least 2")
        merges = parse_merge(args.merge)
        e_ids, e_cols, e = read_table(args.estimated, "estimated")
        t_ids, t_cols, t = read_table(args.truth, "truth")
        result = benchmark(e_ids, e_cols, e, t_ids, t_cols, t, merges, args.renormalize, args.min_n)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    text = json.dumps(result, indent=2) + "\n" if args.format == "json" else to_tsv(result)
    try:
        if args.out:
            with open(args.out, "w", encoding="utf-8") as fh:
                fh.write(text)
        else:
            sys.stdout.write(text)
    except OSError as exc:
        print(f"error: cannot write output: {exc}", file=sys.stderr)
        return 2
    for w in result["warnings"]:
        print(f"warning: {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
