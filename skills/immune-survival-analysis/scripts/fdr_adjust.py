#!/usr/bin/env python3
"""Benjamini-Hochberg and Holm adjusted p-values for a column of p-values.

Standard library only. Reads a TSV (or .csv) with one row per test, keeps every column and
appends p_bh and p_holm, in the input row order (or sorted by p with --sort).

Usage examples:
  python3 fdr_adjust.py --input cox_per_feature.tsv --column p > cox_per_feature_adj.tsv
  python3 fdr_adjust.py --input top_hits.tsv --column pvalue --n-tests 64 --out top_hits_adj.tsv

Formulas (the same as R p.adjust with methods 'BH' and 'holm'; source keys resolve in
references/sources.md). Sort the m non-missing p-values, p(1) <= ... <= p(m). N is the number
of tests that were performed: m by default, or --n-tests when the file holds only some of them
(for example the rows that were kept after a filter); N must be >= m.
  BH    q(i) = min over k >= i of min(1, N / k * p(k))           [Benjamini1995]
  Holm  h(i) = max over k <= i of min(1, (N - k + 1) * p(k))     [Holm1979]
BH controls the false discovery rate (the expected share of false findings among the tests
called significant); Holm controls the family-wise error rate (the chance of at least one
false finding). Holm-adjusted values are never smaller than BH-adjusted ones. Missing
p-values (empty, NA, NaN) stay missing and do not count towards m.

Output: the table as TSV on stdout (or --out), plus a one-line summary on stderr.
Exit code 2 on invalid input.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys

MISSING = {"", "na", "nan", "n/a", "null", "none", "."}
NEW_COLUMNS = ("p_bh", "p_holm")


class InputError(ValueError):
    """Raised for problems in the input table or the options; main() turns it into exit code 2."""


def _tests_performed(m: int, n_tests) -> int:
    if n_tests is None:
        return m
    if int(n_tests) != n_tests or n_tests < 1:
        raise InputError(f"--n-tests must be a positive whole number, got {n_tests}")
    if n_tests < m:
        raise InputError(f"--n-tests ({n_tests}) is smaller than the number of p-values in the file ({m})")
    return int(n_tests)


def _ascending(pvalues) -> list:
    """Indices of the non-missing p-values, sorted by p (stable, so ties keep input order)."""
    return sorted((i for i, p in enumerate(pvalues) if p is not None), key=lambda i: pvalues[i])


def adjust_bh(pvalues, n_tests=None) -> list:
    """Benjamini-Hochberg step-up adjusted p-values; None stays None."""
    order = _ascending(pvalues)
    n = _tests_performed(len(order), n_tests)
    out = [None] * len(pvalues)
    running = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        running = min(running, n / rank * pvalues[i])
        out[i] = min(1.0, running)
    return out


def adjust_holm(pvalues, n_tests=None) -> list:
    """Holm step-down adjusted p-values; None stays None."""
    order = _ascending(pvalues)
    n = _tests_performed(len(order), n_tests)
    out = [None] * len(pvalues)
    running = 0.0
    for rank in range(1, len(order) + 1):
        i = order[rank - 1]
        running = max(running, (n - rank + 1) * pvalues[i])
        out[i] = min(1.0, running)
    return out


def parse_pvalues(rows, column: str) -> list:
    """Return the p-values of one column as floats in [0, 1], or None when missing."""
    values = []
    for number, row in enumerate(rows, start=2):  # line 1 is the header
        text = row.get(column)
        if text is None or text.strip().lower() in MISSING:
            values.append(None)
            continue
        try:
            p = float(text)
        except ValueError:
            raise InputError(f"line {number}: '{text}' in column '{column}' is not a number") from None
        if not math.isfinite(p) or p < 0.0 or p > 1.0:
            raise InputError(f"line {number}: p-value {text} in column '{column}' is outside [0, 1]")
        values.append(p)
    return values


def read_table(path: str, delimiter=None) -> tuple:
    if delimiter is None:
        delimiter = "," if path.lower().endswith(".csv") else "\t"
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            header = list(reader.fieldnames or [])
            rows = [dict(row) for row in reader]
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc.strerror or exc}") from exc
    if not header:
        raise InputError(f"{path} is empty: a header line with column names is required")
    return header, rows


def _fmt(value) -> str:
    return "NA" if value is None else format(value, ".6g")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Append Benjamini-Hochberg (p_bh) and Holm (p_holm) adjusted p-values to a table.",
        epilog="Adjust over every test that was run (all features x endpoints x models), not only the rows "
               "you plan to show: pass --n-tests when the file is a subset.")
    parser.add_argument("--input", required=True, help="table with one row per test (TSV, or CSV if it ends in .csv)")
    parser.add_argument("--column", required=True, help="column holding the raw p-values")
    parser.add_argument("--n-tests", type=int, help="number of tests performed, if larger than the rows in the file")
    parser.add_argument("--alpha", type=float, default=0.05, help="threshold used only for the summary line (default 0.05)")
    parser.add_argument("--sort", action="store_true", help="sort the output by raw p-value (missing last)")
    parser.add_argument("--delimiter", help="field delimiter if not tab / .csv comma")
    parser.add_argument("--out", help="write the TSV here instead of stdout")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if not 0.0 < args.alpha < 1.0:
            raise InputError(f"--alpha must be between 0 and 1, got {args.alpha}")
        header, rows = read_table(args.input, args.delimiter)
        if args.column not in header:
            raise InputError(f"column '{args.column}' not found; available columns: {', '.join(header)}")
        clash = [c for c in NEW_COLUMNS if c in header]
        if clash:
            raise InputError(f"the input already has column(s) {', '.join(clash)}; rename them first")
        pvalues = parse_pvalues(rows, args.column)
        bh = adjust_bh(pvalues, args.n_tests)
        holm = adjust_holm(pvalues, args.n_tests)
        index = list(range(len(rows)))
        if args.sort:
            index.sort(key=lambda i: (pvalues[i] is None, pvalues[i] if pvalues[i] is not None else 0.0))
        try:
            handle = open(args.out, "w", newline="", encoding="utf-8") if args.out else sys.stdout
        except OSError as exc:
            raise InputError(f"cannot write {args.out}: {exc.strerror or exc}") from exc
        try:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header + list(NEW_COLUMNS))
            for i in index:
                writer.writerow([rows[i].get(c, "") for c in header] + [_fmt(bh[i]), _fmt(holm[i])])
        finally:
            if args.out:
                handle.close()
        m = sum(p is not None for p in pvalues)
        n = _tests_performed(m, args.n_tests)
        k_bh = sum(1 for q in bh if q is not None and q <= args.alpha)
        k_holm = sum(1 for q in holm if q is not None and q <= args.alpha)
        print(f"{m} p-values (adjusted for N = {n} tests, {len(rows) - m} missing); "
              f"BH <= {args.alpha}: {k_bh}; Holm <= {args.alpha}: {k_holm}", file=sys.stderr)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
