#!/usr/bin/env python3
"""Sensitivity, specificity, PPV and NPV of a response biomarker at ONE prespecified threshold.

Standard library only; deterministic. Research use only.

What is computed
  * A patient is test-positive when score >= threshold (--direction higher, the default) or
    score <= threshold (--direction lower). Label 1 = responder (the condition of interest).
  * sensitivity = TP/(TP+FN), specificity = TN/(TN+FP), PPV = TP/(TP+FP), NPV = TN/(TN+FN),
    each with an exact Clopper-Pearson interval [Clopper1934]: for x successes in n trials
    the limits are the alpha/2 quantile of Beta(x, n-x+1) and the 1-alpha/2 quantile of
    Beta(x+1, n-x), with 0 when x = 0 and 1 when x = n. The beta quantile is found by
    bisection on the regularized incomplete beta function (continued fraction).
  * PPV and NPV depend on the response rate of the cohort. --prevalence recomputes them for
    another response rate by Bayes' theorem (point estimates only).
  * --cohort repeats the table per cohort.
  * warnings (heuristic flag): a sensitivity or specificity interval wider than 0.30; labels
    excluded because they are in neither --positive nor --negative.

The script takes a threshold fixed in advance (a published cut-off or one learned in a
separate training cohort). It deliberately has no search for an "optimal" cut-off: a
cut-off chosen on the same patients overstates sensitivity and specificity [Leeflang2008].

Usage
  python3 threshold_metrics.py cohort.tsv --score tmb --label responder --threshold 10
  python3 threshold_metrics.py cohort.tsv --score tide --direction lower --threshold 0 \
      --label bor --positive CR,PR --negative SD,PD --cohort study --prevalence 0.2
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

MISSING = {"", "na", "nan", "n/a", "null", "none", "."}


class InputError(Exception):
    """Unusable input; main() prints the message and exits with code 2."""


# ----------------------------------------------------------------------------- statistics

def _betacf(a, b, x):
    """Continued fraction for the incomplete beta function (modified Lentz method)."""
    tiny, eps = 1e-300, 1e-15
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc(a, b, x):
    """Regularized incomplete beta function I_x(a, b) for a, b > 0."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                 + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def beta_quantile(q, a, b):
    """x with I_x(a, b) = q, by bisection (I_x is increasing in x)."""
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if betainc(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-15:
            break
    return (lo + hi) / 2


def clopper_pearson(x, n, alpha=0.05):
    """Exact two-sided (1 - alpha) interval for a binomial proportion [Clopper1934]."""
    if n == 0:
        return None, None
    if not 0 <= x <= n:
        raise ValueError("need 0 <= x <= n")
    lower = 0.0 if x == 0 else beta_quantile(alpha / 2, x, n - x + 1)
    upper = 1.0 if x == n else beta_quantile(1 - alpha / 2, x + 1, n - x)
    return lower, upper


def proportion(x, n, alpha):
    lo, hi = clopper_pearson(x, n, alpha)
    return {"x": x, "n": n, "estimate": x / n if n else None, "ci_low": lo, "ci_high": hi}


def confusion(records, threshold, direction="higher"):
    """Counts (tp, fp, fn, tn); test-positive means score >= threshold (or <= for 'lower')."""
    tp = fp = fn = tn = 0
    for r in records:
        positive = r["s"] >= threshold if direction == "higher" else r["s"] <= threshold
        if r["y"] == 1:
            tp, fn = (tp + 1, fn) if positive else (tp, fn + 1)
        else:
            fp, tn = (fp + 1, tn) if positive else (fp, tn + 1)
    return tp, fp, fn, tn


def metrics(tp, fp, fn, tn, alpha=0.05):
    n = tp + fp + fn + tn
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "sensitivity": proportion(tp, tp + fn, alpha),
            "specificity": proportion(tn, tn + fp, alpha),
            "ppv": proportion(tp, tp + fp, alpha),
            "npv": proportion(tn, tn + fn, alpha),
            "response_rate": proportion(tp + fn, n, alpha),
            "test_positive_rate": proportion(tp + fp, n, alpha)}


def predictive_values(sens, spec, prevalence):
    """PPV and NPV at another response rate (Bayes' theorem)."""
    ppv_den = sens * prevalence + (1 - spec) * (1 - prevalence)
    npv_den = spec * (1 - prevalence) + (1 - sens) * prevalence
    return (sens * prevalence / ppv_den if ppv_den else None,
            spec * (1 - prevalence) / npv_den if npv_den else None)


# ----------------------------------------------------------------------------- input

def _missing(text):
    return text is None or text.strip().casefold() in MISSING


def _categories(text):
    return {v.strip().casefold() for v in text.split(",") if v.strip()} if text else set()


def read_records(path, score, label, positive=None, negative=None, cohort=None):
    """Read a TSV/CSV and return (records {y, s, c}, exclusion counts, rows read)."""
    delim = "," if path.lower().endswith(".csv") else "\t"
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            if not reader.fieldnames:
                raise InputError(f"{path}: no header line")
            header = [h.strip() for h in reader.fieldnames]
            rows = [{(k or "").strip(): v for k, v in row.items()} for row in reader]
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc}") from None
    for col in (score, label, cohort):
        if col and col not in header:
            raise InputError(f"column '{col}' not found; columns are: {', '.join(header)}")
    pos_set, neg_set = _categories(positive), _categories(negative)
    if bool(pos_set) != bool(neg_set):
        raise InputError("give both --positive and --negative, or neither (then labels must be 1/0)")
    if pos_set & neg_set:
        raise InputError(f"label values in both --positive and --negative: {sorted(pos_set & neg_set)}")
    records, excluded, missing_score = [], {}, 0
    for i, row in enumerate(rows):
        line = i + 2
        raw = (row.get(label) or "").strip()
        if _missing(raw):
            excluded["(missing)"] = excluded.get("(missing)", 0) + 1
            continue
        if pos_set:
            key = raw.casefold()
            if key not in pos_set and key not in neg_set:
                excluded[raw] = excluded.get(raw, 0) + 1
                continue
            y = 1 if key in pos_set else 0
        else:
            try:
                value = float(raw)
            except ValueError:
                value = None
            if value not in (0.0, 1.0):
                raise InputError(f"column '{label}' line {line}: '{raw}' is not 1/0; "
                                 "map categories with --positive and --negative")
            y = int(value)
        if _missing(row.get(score)):
            missing_score += 1
            continue
        try:
            s = float(row[score])
        except ValueError:
            raise InputError(f"column '{score}' line {line}: '{row[score]}' is not a number") from None
        if not math.isfinite(s):
            raise InputError(f"column '{score}' line {line}: non-finite value")
        c = None
        if cohort:
            if _missing(row.get(cohort)):
                raise InputError(f"column '{cohort}' line {line}: empty cohort value")
            c = row[cohort].strip()
        records.append({"y": y, "s": s, "c": c})
    return records, {"excluded_by_label": excluded, "excluded_missing_score": missing_score}, len(rows)


# ----------------------------------------------------------------------------- command line

def _round(obj):
    if isinstance(obj, float):
        return round(obj, 6)
    if isinstance(obj, dict):
        return {k: _round(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_round(v) for v in obj]
    return obj


def build_parser():
    p = argparse.ArgumentParser(
        description="Sensitivity, specificity, PPV and NPV with exact Clopper-Pearson intervals at one "
                    "prespecified threshold. No cut-off search. Research use only.",
        epilog="Test-positive: score >= threshold (or <= with --direction lower). Label 1 = responder. "
               "Exit code 2 on unusable input.")
    p.add_argument("table", help="TSV (or .csv) file, one row per patient")
    p.add_argument("--score", required=True, help="column holding the marker")
    p.add_argument("--label", required=True, help="column holding the response label (1/0 or categories)")
    p.add_argument("--threshold", required=True, type=float, help="prespecified cut-off")
    p.add_argument("--direction", choices=("higher", "lower"), default="higher",
                   help="higher: score >= threshold is test-positive (default); lower: score <= threshold")
    p.add_argument("--positive", help="comma-separated label values that count as responders, e.g. CR,PR")
    p.add_argument("--negative", help="comma-separated label values that count as non-responders, e.g. SD,PD")
    p.add_argument("--cohort", help="column naming the cohort; adds one table per cohort")
    p.add_argument("--prevalence", type=float, help="response rate at which to recompute PPV and NPV")
    p.add_argument("--alpha", type=float, default=0.05, help="1 - confidence level (default 0.05)")
    p.add_argument("--out", help="write the JSON here instead of standard output")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")  # legacy Windows code pages cannot print every character
    try:
        if not 0 < args.alpha < 1:
            raise InputError("--alpha must lie strictly between 0 and 1")
        if not math.isfinite(args.threshold):
            raise InputError(f"--threshold must be a finite number, not {args.threshold}")
        if args.prevalence is not None and not 0 < args.prevalence < 1:
            raise InputError("--prevalence must lie strictly between 0 and 1")
        records, counts, n_rows = read_records(args.table, args.score, args.label, args.positive,
                                               args.negative, args.cohort)
        if not any(r["y"] == 1 for r in records) or not any(r["y"] == 0 for r in records):
            raise InputError("need at least one responder and one non-responder after exclusions")
        overall = metrics(*confusion(records, args.threshold, args.direction), alpha=args.alpha)
        rule = ("score >= " if args.direction == "higher" else "score <= ") + f"{args.threshold:g}"
        result = {"input": {"score": args.score, "label": args.label, "threshold": args.threshold,
                            "direction": args.direction, "test_positive_rule": rule,
                            "positive": args.positive, "negative": args.negative,
                            "cohort": args.cohort, "alpha": args.alpha},
                  "rows": dict(read=n_rows, used=len(records), **counts),
                  "overall": overall, "warnings": []}
        if args.cohort:
            result["overall"]["note"] = "pools all patients and ignores cohort"
            result["by_cohort"] = []
            for c in sorted({r["c"] for r in records}):
                row = {"cohort": c}
                row.update(metrics(*confusion([r for r in records if r["c"] == c], args.threshold,
                                              args.direction), alpha=args.alpha))
                result["by_cohort"].append(row)
        if args.prevalence is not None:
            ppv, npv = predictive_values(overall["sensitivity"]["estimate"],
                                         overall["specificity"]["estimate"], args.prevalence)
            result["at_prevalence"] = {"prevalence": args.prevalence, "ppv": ppv, "npv": npv,
                                       "note": "point estimates from this cohort's sensitivity and specificity"}
        for name in ("sensitivity", "specificity"):
            m = overall[name]
            if m["ci_high"] - m["ci_low"] > 0.3:
                result["warnings"].append(f"{name} {100 * (1 - args.alpha):g}% CI is {m['ci_high'] - m['ci_low']:.2f} wide "
                                          f"({m['x']}/{m['n']}): too few patients for a precise estimate "
                                          "(heuristic flag above 0.30)")
        real = [k for k in counts["excluded_by_label"] if k != "(missing)"]
        if real:
            result["warnings"].append("excluded by label: " + ", ".join(
                f"{k} = {counts['excluded_by_label'][k]}" for k in real)
                + "; rerun with them added to --negative as a sensitivity analysis")
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    text = json.dumps(_round(result), indent=2, ensure_ascii=False)
    if args.out:
        try:
            with open(args.out, "w", encoding="utf-8") as fh:
                fh.write(text + "\n")
        except OSError as exc:
            print(f"error: cannot write {args.out}: {exc}", file=sys.stderr)
            return 2
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
