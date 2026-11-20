#!/usr/bin/env python3
"""Reference-based deconvolution of bulk expression by non-negative least squares (NNLS).

For every mixture sample m (one value per gene) the script finds coefficients f >= 0 that
minimise ||S f - m||_2, where S is the signature matrix (genes x cell types). The solver is
the Lawson-Hanson active-set algorithm for NNLS [Lawson1995]; each unconstrained
sub-problem on the current passive set P is solved from the normal equations
(S_P' S_P) z = S_P' m by Cholesky factorisation. Fractions are f / sum(f): they describe the
mix of the cell types in the signature only, never the share of all cells in the sample.
sum(f) is reported as coef_sum.

Preprocessing, in this order:
  1. Gene identifiers are matched case-insensitively after trimming. At least
     2 x (number of cell types) shared genes are required (the tmelens convention).
     Matches that differ only in letter case are counted and flagged: human symbols are
     upper case and mouse symbols sentence case (CD8A, Cd8a), and a case-insensitive match
     across species is not an ortholog mapping.
  2. A matrix declared log2(x + 1) is linearised with 2^x - 1 (--mixture-scale log2,
     --signature-scale log2). With 'auto' (the default) a matrix that looks log-scaled stops
     the run, because regression deconvolution assumes linear-scale data
     [Newman2015; AvilaCobos2020]. The test is the NCBI GEO2R rule, as in tmelens: values
     look linear if the 99th percentile exceeds 100, or if max - min > 50 and the 25th
     percentile is positive.
  3. Signature rows that are zero for every cell type are dropped.
  4. Row scaling (default on, as in tmelens): each gene of S and m is divided by the gene's
     mean across the signature columns, so that highly expressed genes do not dominate.

Diagnostics (on the system actually solved, i.e. after row scaling):
  * 2-norm condition number of S, sqrt(lambda_max / lambda_min) of S'S (cyclic Jacobi);
    larger means less stable estimates [Abbas2009; Newman2015]. Reported as null above
    1e6 (numerically singular). The unscaled value is given for comparison with papers.
  * Pearson correlation between signature columns; pairs at or above --collinear-r.
  * Marker coverage: per cell type, the signature genes whose highest value is in that
    cell type, and how many were found in the mixture.
  * Per sample: Pearson r and RMSE between fitted and observed values, relative residual
    ||m - S f|| / ||m||, coef_sum and the number of non-zero coefficients.
  * Genes with the largest mean squared residual across samples.

Input: tab-separated (comma-separated if the name ends in .csv). First column: gene
identifier. Header: cell types (signature) or sample IDs (mixture). Numeric values only.
Output: JSON (default) or TSV with one row per sample. Exit code 2 on invalid input.

Usage examples:
  python3 nnls_deconvolve.py --signature signature.tsv --mixture tpm.tsv > deconvolution.json
  python3 nnls_deconvolve.py --signature signature.tsv --mixture log2_tpm.tsv \\
      --mixture-scale log2 --format tsv --out fractions.tsv --residuals-out residuals.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

SINGULAR_COND = 1e6


class InputError(ValueError):
    """Invalid input; main() reports it and exits with code 2."""


# --------------------------------------------------------------------------- input

def read_matrix(path: str, label: str):
    """Read a table whose first column is a row identifier. Returns (ids, columns, rows)."""
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
        raise InputError(f"{label}: header needs non-empty column names after the identifier column")
    dups = sorted({c for c in cols if cols.count(c) > 1})
    if dups:
        raise InputError(f"{label}: duplicated column names: {', '.join(dups)}")
    ids, values, seen = [], [], {}
    for n, r in enumerate(rows[1:], start=2):
        if len(r) != len(header):
            raise InputError(f"{label}: line {n} has {len(r)} fields but the header has {len(header)}")
        gid = r[0].strip()
        if not gid:
            raise InputError(f"{label}: line {n} has an empty identifier")
        key = gid.upper()
        if key in seen:
            raise InputError(f"{label}: identifier '{gid}' appears twice (lines {seen[key]} and {n}); collapse duplicates first")
        seen[key] = n
        vals = []
        for col, cell in zip(cols, r[1:]):
            try:
                v = float(cell)
            except ValueError:
                raise InputError(f"{label}: non-numeric value '{cell}' for '{gid}', column '{col}'") from None
            if not math.isfinite(v):
                raise InputError(f"{label}: missing or infinite value for '{gid}', column '{col}'")
            vals.append(v)
        ids.append(gid)
        values.append(vals)
    return ids, cols, values


def quantile(sorted_vals, q):
    """Quantile by linear interpolation between order statistics (R's default, type 7)."""
    h = (len(sorted_vals) - 1) * q
    lo = int(math.floor(h))
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (h - lo) * (sorted_vals[hi] - sorted_vals[lo])


def looks_linear(values) -> bool:
    """NCBI GEO2R rule (as in tmelens): True when the values do not look log-transformed."""
    flat = sorted(v for row in values for v in row)
    return quantile(flat, 0.99) > 100 or (flat[-1] - flat[0] > 50 and quantile(flat, 0.25) > 0)


def to_linear(values, scale: str, label: str):
    """Return (linear values, number of entries clipped at zero, scale used)."""
    if scale == "auto":
        if not looks_linear(values):
            flat = sorted(v for row in values for v in row)
            raise InputError(
                f"{label}: values look log-transformed (99th percentile {quantile(flat, 0.99):.3g}, "
                f"maximum {flat[-1]:.3g}). Regression deconvolution needs linear-scale data. Re-run with "
                f"--{label}-scale log2 to linearise log2(x + 1) values, or --{label}-scale linear if they are linear.")
        scale = "linear"
    if scale == "linear":
        for row in values:
            if any(v < 0 for v in row):
                n_neg = sum(1 for r in values for v in r if v < 0)
                raise InputError(f"{label}: {n_neg} negative value(s) on a linear scale; z-scores, log ratios "
                                 "and batch-corrected values cannot be deconvolved")
        return [list(row) for row in values], 0, "linear"
    if scale != "log2":
        raise InputError(f"{label}: unknown scale '{scale}'")
    out, clipped = [], 0
    for row in values:
        new = []
        for v in row:
            if v > 1000:
                raise InputError(f"{label}: value {v:g} is far too large for log2 data; is the matrix already linear?")
            x = 2.0 ** v - 1.0
            if x < 0:
                x, clipped = 0.0, clipped + 1
            new.append(x)
        out.append(new)
    return out, clipped, "log2"


# --------------------------------------------------------------------------- linear algebra

def dot(a, b) -> float:
    return sum(x * y for x, y in zip(a, b))


def gram(rows, k: int):
    """Cross-product matrix A'A for A given as a list of rows with k columns."""
    g = [[0.0] * k for _ in range(k)]
    for row in rows:
        for a in range(k):
            ra = row[a]
            if ra:
                ga = g[a]
                for b in range(a, k):
                    ga[b] += ra * row[b]
    for a in range(k):
        for b in range(a):
            g[a][b] = g[b][a]
    return g


def cholesky_solve(m, rhs):
    """Solve m z = rhs for symmetric positive-definite m; None when m is not (numerically) PD."""
    n = len(m)
    if n == 0:
        return []
    top = max(m[i][i] for i in range(n))
    if top <= 0:
        return None
    low = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = m[i][j] - sum(low[i][p] * low[j][p] for p in range(j))
            if i == j:
                if s <= 1e-13 * top:
                    return None
                low[i][i] = math.sqrt(s)
            else:
                low[i][j] = s / low[j][j]
    y = [0.0] * n
    for i in range(n):
        y[i] = (rhs[i] - sum(low[i][p] * y[p] for p in range(i))) / low[i][i]
    z = [0.0] * n
    for i in reversed(range(n)):
        z[i] = (y[i] - sum(low[p][i] * z[p] for p in range(i + 1, n))) / low[i][i]
    return z


def _solve_passive(g, h, passive):
    idx = [j for j, on in enumerate(passive) if on]
    sub = cholesky_solve([[g[a][b] for b in idx] for a in idx], [h[a] for a in idx])
    if sub is None:
        return None
    z = [0.0] * len(h)
    for j, v in zip(idx, sub):
        z[j] = v
    return z


def nnls_gram(g, h, tol=None, max_iter=None):
    """Lawson-Hanson active-set NNLS from the cross products g = A'A and h = A'b.

    Returns (x, iterations). tol is the threshold on the dual vector w = h - g x
    (default 1e-10 x max|h|). A column whose addition makes the passive sub-problem
    singular, or whose new coefficient is not positive, is skipped for that step (the
    safeguard of the original algorithm), so exactly collinear columns cannot loop forever.
    """
    n = len(h)
    if tol is None:
        tol = 1e-10 * max([abs(v) for v in h] + [0.0])
    if max_iter is None:
        max_iter = 30 * n + 50
    x = [0.0] * n
    passive = [False] * n
    w = list(h)
    excluded = set()
    it = 0
    while True:
        cand = [j for j in range(n) if not passive[j] and j not in excluded and w[j] > tol]
        if not cand:
            break
        t = max(cand, key=lambda j: w[j])
        passive[t] = True
        z = _solve_passive(g, h, passive)
        if z is None or z[t] <= 0.0:
            passive[t] = False
            excluded.add(t)
            continue
        while True:
            neg = [j for j in range(n) if passive[j] and z[j] <= 0.0]
            if not neg:
                break
            it += 1
            if it > max_iter:
                raise RuntimeError("NNLS did not converge")

            def ratio(j):
                d = x[j] - z[j]
                return x[j] / d if d > 0 else 0.0

            q = min(neg, key=ratio)
            alpha = ratio(q)
            x = [x[j] + alpha * (z[j] - x[j]) if passive[j] else 0.0 for j in range(n)]
            x[q] = 0.0
            for j in range(n):
                if passive[j] and x[j] <= 0.0:
                    passive[j] = False
                    x[j] = 0.0
            z = _solve_passive(g, h, passive)
            if z is None:  # cannot happen: a principal submatrix of a PD matrix is PD
                raise RuntimeError("singular NNLS sub-problem")
        x = z
        w = [h[j] - dot(g[j], x) for j in range(n)]
        excluded = set()
        it += 1
        if it > max_iter:
            raise RuntimeError("NNLS did not converge")
    return x, it


def nnls(a_rows, b):
    """NNLS for a matrix given as rows; returns (x, residual 2-norm)."""
    k = len(a_rows[0])
    g = gram(a_rows, k)
    h = [sum(row[j] * bi for row, bi in zip(a_rows, b)) for j in range(k)]
    x, _ = nnls_gram(g, h)
    rnorm = math.sqrt(sum((bi - dot(row, x)) ** 2 for row, bi in zip(a_rows, b)))
    return x, rnorm


def symmetric_eigenvalues(m, max_sweeps: int = 100):
    """Eigenvalues (ascending) of a symmetric matrix by the cyclic Jacobi method."""
    n = len(m)
    a = [list(row) for row in m]
    for _ in range(max_sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        total = sum(a[i][j] ** 2 for i in range(n) for j in range(n))
        if off <= 1e-30 * total or total == 0:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if apq == 0.0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * apq)
                t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
    return sorted(a[i][i] for i in range(n))


def condition_number(rows, k: int):
    """2-norm condition number of a matrix (rows x k); None when above SINGULAR_COND."""
    eig = symmetric_eigenvalues(gram(rows, k))
    if not eig or eig[-1] <= 0:
        return None
    if eig[0] <= eig[-1] / SINGULAR_COND ** 2:
        return None
    return math.sqrt(eig[-1] / eig[0])


def pearson(x, y):
    """Pearson correlation; None for fewer than two points or a constant vector."""
    n = len(x)
    if n < 2 or max(x) == min(x) or max(y) == min(y):
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy)))


def _round(v, d: int = 6):
    return None if v is None else round(v, d)


# --------------------------------------------------------------------------- deconvolution

def deconvolve(sig_ids, cell_types, sig_values, mix_ids, samples, mix_values,
               row_scale: bool = True, collinear_r: float = 0.9,
               min_genes_factor: float = 2.0, top_residuals: int = 10):
    """Deconvolve linear-scale mixtures. Returns (result dict, shared genes, residual rows)."""
    k = len(cell_types)
    mix_index = {gid.upper(): i for i, gid in enumerate(mix_ids)}
    coverage = {ct: {"n_marker_genes": 0, "n_found": 0} for ct in cell_types}
    shared, missing, zero_rows = [], [], 0
    for i, gid in enumerate(sig_ids):
        row = sig_values[i]
        found = gid.upper() in mix_index
        if max(row) > 0:
            top = cell_types[row.index(max(row))]
            coverage[top]["n_marker_genes"] += 1
            coverage[top]["n_found"] += int(found)
        if not found:
            missing.append(gid)
        elif max(row) <= 0:
            zero_rows += 1
        else:
            shared.append(i)
    need = int(math.ceil(min_genes_factor * k))
    if len(shared) < need:
        raise InputError(f"only {len(shared)} informative signature genes are present in the mixture; "
                         f"need at least {need} for {k} cell types (check identifiers and platform)")
    genes = [sig_ids[i] for i in shared]
    case_mismatch = [(sig_ids[i], mix_ids[mix_index[sig_ids[i].upper()]]) for i in shared
                     if sig_ids[i] != mix_ids[mix_index[sig_ids[i].upper()]]]
    s_raw = [list(sig_values[i]) for i in shared]
    m_rows = [list(mix_values[mix_index[sig_ids[i].upper()]]) for i in shared]
    s_fit = s_raw
    if row_scale:
        s_fit, scaled_m = [], []
        for srow, mrow in zip(s_raw, m_rows):
            scale = sum(srow) / k
            s_fit.append([v / scale for v in srow])
            scaled_m.append([v / scale for v in mrow])
        m_rows = scaled_m

    g = gram(s_fit, k)
    cond_fit = condition_number(s_fit, k)
    cond_raw = condition_number(s_raw, k)
    pairs = []
    for a in range(k):
        for b in range(a + 1, k):
            r = pearson([row[a] for row in s_fit], [row[b] for row in s_fit])
            if r is not None and r >= collinear_r:
                pairs.append({"cell_type_a": cell_types[a], "cell_type_b": cell_types[b], "r": _round(r, 4)})

    per_sample, residual_rows = [], [[0.0] * len(samples) for _ in genes]
    warnings = []
    for j, name in enumerate(samples):
        m = [row[j] for row in m_rows]
        h = [sum(row[c] * mv for row, mv in zip(s_fit, m)) for c in range(k)]
        try:
            x, _ = nnls_gram(g, h)
        except RuntimeError as exc:
            raise InputError(f"sample {name}: {exc}; check for extreme values or rescale the input") from None
        fitted = [dot(row, x) for row in s_fit]
        resid = [mv - fv for mv, fv in zip(m, fitted)]
        for gi, rv in enumerate(resid):
            residual_rows[gi][j] = rv
        total = sum(x)
        norm_m = math.sqrt(sum(v * v for v in m))
        rec = {
            "sample": name,
            "fractions": {ct: _round(v / total) for ct, v in zip(cell_types, x)} if total > 0 else None,
            "coefficients": {ct: _round(v) for ct, v in zip(cell_types, x)},
            "coef_sum": _round(total),
            "n_nonzero": sum(1 for v in x if v > 0),
            "pearson_r": _round(pearson(fitted, m), 4),
            "rmse": _round(math.sqrt(sum(v * v for v in resid) / len(resid))),
            "rel_residual": _round(math.sqrt(sum(v * v for v in resid)) / norm_m, 4) if norm_m > 0 else None,
        }
        if total <= 0:
            warnings.append(f"sample {name}: every coefficient is zero, so no fractions can be formed")
        per_sample.append(rec)

    n_s = len(samples)
    gene_stats = sorted(
        ({"gene": gname, "mean_sq_residual": sum(v * v for v in row) / n_s, "mean_residual": sum(row) / n_s}
         for gname, row in zip(genes, residual_rows)),
        key=lambda d: (-d["mean_sq_residual"], d["gene"]))[:max(0, top_residuals)]
    for d in gene_stats:
        d["mean_sq_residual"] = _round(d["mean_sq_residual"])
        d["mean_residual"] = _round(d["mean_residual"])

    if case_mismatch:
        a, b = case_mismatch[0]
        warnings.append(f"{len(case_mismatch)} of {len(genes)} genes matched only after ignoring letter case "
                        f"(for example '{a}' and '{b}'); check that both matrices use the same species and "
                        "identifier type")
    if cond_fit is None:
        warnings.append("the signature matrix is numerically singular: some cell types cannot be separated")
    for p in pairs:
        warnings.append(f"{p['cell_type_a']} and {p['cell_type_b']} have correlated signatures (r = {p['r']}); "
                        "their split is unstable, report their sum or merge them")
    for ct, cov in coverage.items():
        if cov["n_marker_genes"] and cov["n_found"] < 0.5 * cov["n_marker_genes"]:
            warnings.append(f"{ct}: only {cov['n_found']} of {cov['n_marker_genes']} marker genes found in the mixture")
        if cov["n_marker_genes"] == 0:
            warnings.append(f"{ct}: no signature gene has its highest value in this cell type")

    result = {
        "method": "NNLS, Lawson-Hanson active set",
        "settings": {"row_scale": row_scale, "collinear_r": collinear_r, "min_genes_factor": min_genes_factor},
        "cell_types": list(cell_types),
        "n_signature_genes": len(sig_ids),
        "n_genes_used": len(genes),
        "n_missing_genes": len(missing),
        "n_case_mismatched_genes": len(case_mismatch),
        "missing_genes": missing[:50],
        "n_zero_signature_rows_dropped": zero_rows,
        "marker_coverage": coverage,
        "condition_number_fit": _round(cond_fit, 3),
        "condition_number_unscaled": _round(cond_raw, 3),
        "collinear_pairs": pairs,
        "samples": per_sample,
        "top_residual_genes": gene_stats,
        "notes": [
            "fractions = coefficients / coef_sum: shares of the cell types in the signature, not of all cells",
            "pearson_r, rmse and residuals are on the fitted scale (after row scaling when it is on)",
        ],
        "warnings": warnings,
    }
    return result, genes, residual_rows


# --------------------------------------------------------------------------- output

def _fmt(v) -> str:
    return "NA" if v is None else (f"{v:.6g}" if isinstance(v, float) else str(v))


def to_tsv(result) -> str:
    cts = result["cell_types"]
    lines = ["\t".join(["sample"] + cts + ["coef_sum", "n_nonzero", "pearson_r", "rmse", "rel_residual"])]
    for rec in result["samples"]:
        fr = rec["fractions"] or {}
        vals = [fr.get(ct) for ct in cts] + [rec["coef_sum"], rec["n_nonzero"], rec["pearson_r"], rec["rmse"], rec["rel_residual"]]
        lines.append("\t".join([rec["sample"]] + [_fmt(v) for v in vals]))
    return "\n".join(lines) + "\n"


def build_parser():
    p = argparse.ArgumentParser(
        description="NNLS deconvolution of bulk expression with a linear-scale signature matrix, "
                    "with collinearity, coverage and goodness-of-fit diagnostics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Fractions are relative to the cell types in the signature. See the module docstring for formulas.")
    p.add_argument("--signature", required=True, help="genes x cell types table (first column: gene)")
    p.add_argument("--mixture", required=True, help="genes x samples table (first column: gene)")
    p.add_argument("--mixture-scale", choices=["auto", "linear", "log2"], default="auto",
                   help="'log2' means log2(x + 1) and is linearised; 'auto' stops if values look log-scaled")
    p.add_argument("--signature-scale", choices=["auto", "linear", "log2"], default="auto",
                   help="as --mixture-scale, for the signature matrix")
    p.add_argument("--no-row-scale", action="store_true", help="do not divide genes by their signature mean")
    p.add_argument("--collinear-r", type=float, default=0.9,
                   help="flag signature column pairs with Pearson r at or above this (default 0.9, a heuristic)")
    p.add_argument("--min-genes-factor", type=float, default=2.0,
                   help="require at least this many shared genes per cell type (default 2)")
    p.add_argument("--top-residual-genes", type=int, default=10, help="genes to list by mean squared residual")
    p.add_argument("--format", choices=["json", "tsv"], default="json")
    p.add_argument("--out", help="write the result to this file instead of stdout")
    p.add_argument("--residuals-out", help="write genes x samples residuals (observed - fitted, fitted scale) as TSV")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if not 0 < args.collinear_r <= 1:
            raise InputError("--collinear-r must be in (0, 1]")
        if args.min_genes_factor < 1:
            raise InputError("--min-genes-factor must be at least 1")
        sig_ids, cell_types, sig_vals = read_matrix(args.signature, "signature")
        mix_ids, samples, mix_vals = read_matrix(args.mixture, "mixture")
        sig_lin, sig_clip, sig_scale = to_linear(sig_vals, args.signature_scale, "signature")
        mix_lin, mix_clip, mix_scale = to_linear(mix_vals, args.mixture_scale, "mixture")
        result, genes, resid = deconvolve(sig_ids, cell_types, sig_lin, mix_ids, samples, mix_lin,
                                          row_scale=not args.no_row_scale, collinear_r=args.collinear_r,
                                          min_genes_factor=args.min_genes_factor,
                                          top_residuals=args.top_residual_genes)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    result["settings"].update({"signature_scale": sig_scale, "mixture_scale": mix_scale})
    for label, n in (("signature", sig_clip), ("mixture", mix_clip)):
        if n:
            result["warnings"].append(f"{label}: {n} value(s) below zero after 2^x - 1 were set to 0")
    text = json.dumps(result, indent=2) + "\n" if args.format == "json" else to_tsv(result)
    try:
        if args.out:
            with open(args.out, "w", encoding="utf-8") as fh:
                fh.write(text)
        else:
            sys.stdout.write(text)
        if args.residuals_out:
            with open(args.residuals_out, "w", encoding="utf-8") as fh:
                fh.write("\t".join(["gene"] + samples) + "\n")
                for gname, row in zip(genes, resid):
                    fh.write("\t".join([gname] + [f"{v:.6g}" for v in row]) + "\n")
    except OSError as exc:
        print(f"error: cannot write output: {exc}", file=sys.stderr)
        return 2
    for w in result["warnings"]:
        print(f"warning: {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
