#!/usr/bin/env python3
"""Kaplan-Meier estimates, median survival and the log-rank test from a patient table.

Standard library only. One row per patient; a time column, a 0/1 event column and either a
group column or a continuous score that is split at a prespecified point.

Usage examples:
  python3 km_logrank.py --input cohort.tsv --time os_months --event os_event --group subtype --id patient_id
  python3 km_logrank.py --input cohort.tsv --time pfi_months --event pfi_event \
      --score cd8_score --split median --times 12,24,60 --curve-out km_curve.tsv

With --id, an ID that occurs on more than one row stops the run (exit code 2): one sample per
patient is required, because several samples of one patient are pseudo-replicates.

What is computed (source keys resolve in references/sources.md):
  Kaplan-Meier   S(t) = prod over event times t_k <= t of (1 - d_k / n_k)        [Kaplan1958]
  Variance       Greenwood: var(S) = S^2 * sum d_k / (n_k * (n_k - d_k)); se = sqrt(var)
  Interval       'log' (default; also the default of survfit in the R survival package):
                     S * exp(-/+ z * sqrt(sum)), upper limit capped at 1
                 'log-log': exp(-exp(log(-log S) +/- z * sqrt(sum) / |log S|))
                 'plain':   S -/+ z * se, truncated to [0, 1]
  Median         smallest observed time with S(t) <= 0.5; its interval is the smallest time
                 at which the lower / upper confidence band is <= 0.5 (null if never)
  Follow-up      reverse Kaplan-Meier median: the same estimator with the event indicator
                 flipped, so that censoring is the event (a common convention)
  Log-rank       for G groups, O_g - E_g with E_g = sum_k d_k * n_gk / n_k and the
                 hypergeometric covariance
                     V_gh = sum_k d_k (n_k - d_k) / (n_k - 1) * (I[g=h] n_gk/n_k - n_gk n_hk / n_k^2);
                 chi-square = (O-E)' V^-1 (O-E) on G - 1 degrees of freedom
                 (R survdiff with rho = 0)                              [Mantel1966] [Peto1972]

Splitting a score: only prespecified rules are offered ('median', 'tertile', or a cut value
fixed before looking at the outcome). There is no option that searches for the cut-point with
the smallest p-value, because that search inflates the false-positive rate [Altman1994]; run
cutpoint_inflation.py to see by how much.

Output: JSON on stdout (or --out). Exit code 2 on invalid input.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import Counter

MISSING = {"", "na", "nan", "n/a", "null", "none", "."}
MEDIAN_TOL = 1e-12


class InputError(ValueError):
    """Raised for problems in the input table or the options; main() turns it into exit code 2."""


# --------------------------------------------------------------------------- statistics

def chi2_sf(x: float, df: int) -> float:
    """Upper tail probability of a chi-square variable with an integer number of degrees of freedom.

    Closed forms: even df = 2m gives exp(-x/2) * sum_{i<m} (x/2)^i / i!; odd df = 2m + 1 gives
    erfc(sqrt(x/2)) + exp(-x/2) * sum_{i=1..m} (x/2)^(i - 1/2) / Gamma(i + 1/2).
    """
    if df < 1 or int(df) != df:
        raise ValueError("df must be a positive integer")
    if x <= 0:
        return 1.0
    half = x / 2.0
    if df % 2 == 0:
        term, total = 1.0, 1.0
        for i in range(1, df // 2):
            term *= half / i
            total += term
        return min(1.0, math.exp(-half) * total)
    total = 0.0
    for i in range(1, (df - 1) // 2 + 1):
        total += math.exp((i - 0.5) * math.log(half) - math.lgamma(i + 0.5))
    return min(1.0, math.erfc(math.sqrt(half)) + math.exp(-half) * total)


def z_quantile(conf_level: float) -> float:
    if not 0.0 < conf_level < 1.0:
        raise InputError(f"--conf-level must be between 0 and 1, got {conf_level}")
    return statistics.NormalDist().inv_cdf(0.5 + conf_level / 2.0)


def kaplan_meier(times, events, conf_level: float = 0.95, conf_type: str = "log") -> list:
    """Kaplan-Meier table: one dict per distinct observed time.

    Keys: time, n_risk, n_event, n_censor, survival, se (Greenwood standard error of S),
    ci_lower, ci_upper. Limits are None where they are undefined (S = 0).
    """
    if conf_type not in ("log", "log-log", "plain"):
        raise InputError("--conf-type must be one of: log, log-log, plain")
    if len(times) != len(events):
        raise ValueError("times and events must have the same length")
    z = z_quantile(conf_level)
    n_event = Counter()
    n_censor = Counter()
    for t, e in zip(times, events):
        if e:
            n_event[t] += 1
        else:
            n_censor[t] += 1
    at_risk = len(times)
    surv = 1.0
    cum = 0.0  # sum of d / (n (n - d)): the variance of log S
    table = []
    for t in sorted(set(n_event) | set(n_censor)):
        d, c = n_event[t], n_censor[t]
        if d:
            surv *= 1.0 - d / at_risk
            cum = math.inf if d == at_risk else cum + d / (at_risk * (at_risk - d))
        if d == at_risk and d:
            surv = 0.0
        if surv <= 0.0 or math.isinf(cum):
            se = lower = upper = None
        else:
            root = math.sqrt(cum)
            se = surv * root
            if conf_type == "log":
                lower, upper = surv * math.exp(-z * root), min(1.0, surv * math.exp(z * root))
            elif conf_type == "plain":
                lower, upper = max(0.0, surv - z * se), min(1.0, surv + z * se)
            elif surv >= 1.0:
                lower, upper = 1.0, 1.0
            else:
                theta = math.log(-math.log(surv))
                half = z * root / abs(math.log(surv))
                lower, upper = math.exp(-math.exp(theta + half)), math.exp(-math.exp(theta - half))
        table.append({"time": t, "n_risk": at_risk, "n_event": d, "n_censor": c,
                      "survival": surv, "se": se, "ci_lower": lower, "ci_upper": upper})
        at_risk -= d + c
    return table


def _first_time_at_or_below_half(table, key):
    for row in table:
        value = row[key]
        if value is None and key == "ci_lower" and row["survival"] <= 0.0:
            value = 0.0  # S(t) = 0: the lower band has reached zero as well
        if value is not None and value <= 0.5 + MEDIAN_TOL:
            return row["time"]
    return None


def median_survival(table) -> dict:
    """Median survival time and its confidence limits from a Kaplan-Meier table.

    The median is the smallest time with S(t) <= 0.5 (None when the curve stays above 0.5).
    The limits are the smallest times at which the lower and the upper confidence band are
    <= 0.5, which is how R survival prints them. R reports the midpoint of the flat segment
    when S(t) is exactly 0.5 on an interval; this function reports the start of that segment.
    """
    return {"median": _first_time_at_or_below_half(table, "survival"),
            "ci_lower": _first_time_at_or_below_half(table, "ci_lower"),
            "ci_upper": _first_time_at_or_below_half(table, "ci_upper")}


def survival_at(table, t: float) -> dict:
    """Step-function value of a Kaplan-Meier table at time t, with the number still at risk."""
    row_at = None
    for row in table:
        if row["time"] <= t:
            row_at = row
        else:
            break
    n_risk = 0
    for row in table:
        if row["time"] >= t:
            n_risk = row["n_risk"]
            break
    last = table[-1]["time"] if table else None
    out = {"time": t, "n_risk": n_risk, "beyond_follow_up": bool(table) and t > last}
    if row_at is None:
        out.update({"survival": 1.0, "ci_lower": 1.0, "ci_upper": 1.0})
    else:
        out.update({k: row_at[k] for k in ("survival", "ci_lower", "ci_upper")})
    return out


def reverse_km_median(times, events):
    """Median follow-up by the reverse Kaplan-Meier method (censoring treated as the event)."""
    flipped = [0 if e else 1 for e in events]
    return median_survival(kaplan_meier(times, flipped))["median"]


def _solve(matrix, rhs):
    """Solve a small linear system by Gaussian elimination with partial pivoting."""
    n = len(rhs)
    a = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    scale = max((abs(v) for row in matrix for v in row), default=0.0)
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) <= 1e-12 * max(scale, 1e-300):
            raise ZeroDivisionError("singular matrix")
        a[col], a[pivot] = a[pivot], a[col]
        for r in range(col + 1, n):
            factor = a[r][col] / a[col][col]
            for k in range(col, n + 1):
                a[r][k] -= factor * a[col][k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (a[i][n] - sum(a[i][k] * x[k] for k in range(i + 1, n))) / a[i][i]
    return x


def logrank_test(times, events, groups, order=None) -> dict:
    """Log-rank test for two or more groups (hypergeometric variance, as R survdiff, rho = 0).

    Returns chi2, df, p and a per-group list with n, observed and expected events. chi2, df and
    p are None, with a 'reason', when the test is undefined (one group, or no events).
    A group with zero expected events carries no information and does not count towards df.
    """
    if not len(times) == len(events) == len(groups):
        raise ValueError("times, events and groups must have the same length")
    levels = list(order) if order is not None else sorted(set(groups), key=str)
    index = {g: i for i, g in enumerate(levels)}
    k = len(levels)
    n_group = [0] * k
    observed = [0.0] * k
    for g, e in zip(groups, events):
        n_group[index[g]] += 1
        observed[index[g]] += 1 if e else 0
    expected = [0.0] * k
    cov = [[0.0] * k for _ in range(k)]
    at_risk = list(n_group)
    by_time = {}
    for t, e, g in zip(times, events, groups):
        slot = by_time.setdefault(t, ([0] * k, [0] * k))  # (events, leaving) per group
        slot[1][index[g]] += 1
        if e:
            slot[0][index[g]] += 1
    for t in sorted(by_time):
        deaths, leaving = by_time[t]
        d, n = sum(deaths), sum(at_risk)
        if d:
            for g in range(k):
                expected[g] += d * at_risk[g] / n
            if n > 1:
                w = d * (n - d) / (n - 1)
                for g in range(k):
                    pg = at_risk[g] / n
                    for h in range(k):
                        cov[g][h] += w * ((pg if g == h else 0.0) - pg * at_risk[h] / n)
        for g in range(k):
            at_risk[g] -= leaving[g]
    table = [{"group": levels[g], "n": n_group[g], "observed": observed[g], "expected": expected[g]}
             for g in range(k)]
    result = {"chi2": None, "df": None, "p": None, "groups": table}
    informative = [g for g in range(k) if expected[g] > 0]
    if k < 2:
        result["reason"] = "only one group: the log-rank test needs at least two"
        return result
    if sum(observed) == 0:
        result["reason"] = "no events: the log-rank test is undefined"
        return result
    if len(informative) < 2:
        result["reason"] = "fewer than two groups had patients at risk at an event time"
        return result
    keep = informative[:-1]
    diff = [observed[g] - expected[g] for g in keep]
    sub = [[cov[g][h] for h in keep] for g in keep]
    try:
        solved = _solve(sub, diff)
    except ZeroDivisionError:
        result["reason"] = "the variance matrix is singular (a group has no variation in risk-set share)"
        return result
    chi2 = max(0.0, sum(d * s for d, s in zip(diff, solved)))
    df = len(informative) - 1
    result.update({"chi2": chi2, "df": df, "p": chi2_sf(chi2, df)})
    return result


# --------------------------------------------------------------------------- splitting a score

def quantile_type7(sorted_values, prob: float) -> float:
    """Sample quantile by linear interpolation (the default 'type 7' of R and NumPy)."""
    if not sorted_values:
        raise ValueError("no values")
    h = (len(sorted_values) - 1) * prob
    lo = int(math.floor(h))
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo])


def split_score(values, method: str, cut=None) -> tuple:
    """Assign prespecified groups to a continuous score.

    'median': 'high' if value > median, else 'low' (ties at the median go to 'low').
    'tertile': 'low' if value <= Q(1/3), 'mid' if value <= Q(2/3), else 'high' (type 7 quantiles).
    'cut': 'high' if value > cut, else 'low'; the cut must have been fixed in advance.
    Returns (labels, order, info).
    """
    ordered = sorted(values)
    if method == "median":
        point = statistics.median(ordered)
        labels = ["high" if v > point else "low" for v in values]
        return labels, ["low", "high"], {"method": "median", "cut_values": [point],
                                         "rule": "high if score > median, otherwise low"}
    if method == "tertile":
        q1, q2 = quantile_type7(ordered, 1 / 3), quantile_type7(ordered, 2 / 3)
        labels = ["low" if v <= q1 else "mid" if v <= q2 else "high" for v in values]
        return labels, ["low", "mid", "high"], {"method": "tertile", "cut_values": [q1, q2],
                                                "rule": "low if score <= Q(1/3), mid if <= Q(2/3), otherwise high"}
    if method == "cut":
        if cut is None:
            raise InputError("--cut needs a value")
        labels = ["high" if v > cut else "low" for v in values]
        return labels, ["low", "high"], {"method": "fixed cut-point supplied by the user", "cut_values": [cut],
                                         "rule": "high if score > cut, otherwise low"}
    raise InputError("--split must be 'median' or 'tertile'")


# --------------------------------------------------------------------------- input

def read_table(path: str, delimiter=None) -> tuple:
    """Read a delimited text table. Returns (header, rows as dicts). Tab unless the file ends in .csv."""
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


def _is_missing(text) -> bool:
    return text is None or text.strip().lower() in MISSING


def parse_records(header, rows, time_col: str, event_col: str, group_col=None, score_col=None, id_col=None) -> tuple:
    """Validate the table and return (records, n_dropped). Each record is (time, event, group_or_score).

    With id_col, an ID that occurs on more than one complete row raises InputError: several
    samples of one patient are pseudo-replicates and must be reduced to one before the analysis.
    """
    needed = [time_col, event_col] + [c for c in (group_col, score_col, id_col) if c]
    for col in needed:
        if col not in header:
            raise InputError(f"column '{col}' not found; available columns: {', '.join(header)}")
    records, dropped, ids = [], 0, Counter()
    for number, row in enumerate(rows, start=2):  # line 1 is the header
        cells = [row.get(c) for c in needed]
        if any(_is_missing(c) for c in cells):
            dropped += 1
            continue
        try:
            t = float(row[time_col])
        except ValueError:
            raise InputError(f"line {number}: time '{row[time_col]}' in column '{time_col}' is not a number") from None
        if not math.isfinite(t) or t < 0:
            raise InputError(f"line {number}: time {row[time_col]} in column '{time_col}' must be finite and >= 0")
        raw = row[event_col].strip()
        if raw in ("1", "1.0"):
            e = 1
        elif raw in ("0", "0.0"):
            e = 0
        else:
            raise InputError(f"line {number}: event '{raw}' in column '{event_col}' must be 1 (event) or 0 "
                             "(censored); recode labels such as '1:DECEASED' or 'Dead' first")
        if score_col:
            try:
                third = float(row[score_col])
            except ValueError:
                raise InputError(f"line {number}: score '{row[score_col]}' in column '{score_col}' is not a number") from None
            if not math.isfinite(third):
                raise InputError(f"line {number}: score in column '{score_col}' must be finite")
        elif group_col:
            third = row[group_col].strip()
        else:
            third = "all"
        records.append((t, e, third))
        if id_col:
            ids[row[id_col].strip()] += 1
    repeated = sorted(k for k, v in ids.items() if v > 1)
    if repeated:
        shown = ", ".join(f"'{k}' ({ids[k]} rows)" for k in repeated[:5])
        raise InputError(f"{len(repeated)} ID(s) in column '{id_col}' occur on more than one row, for example {shown}; "
                         "keep one sample per patient first (in TCGA the patient is the first 12 characters of the barcode)")
    return records, dropped


# --------------------------------------------------------------------------- report

def analyse(records, order=None, conf_level: float = 0.95, conf_type: str = "log", report_times=()) -> dict:
    """Kaplan-Meier summary per group, the log-rank test and the reverse-KM median follow-up."""
    if not records:
        raise InputError("no complete rows: every row has a missing time, event or group")
    times = [r[0] for r in records]
    events = [r[1] for r in records]
    groups = [r[2] for r in records]
    levels = list(order) if order else sorted(set(groups), key=str)
    levels = [g for g in levels if g in set(groups)]
    warnings = []
    summaries, curves = [], {}
    for g in levels:
        gt = [t for t, grp in zip(times, groups) if grp == g]
        ge = [e for e, grp in zip(events, groups) if grp == g]
        table = kaplan_meier(gt, ge, conf_level, conf_type)
        curves[g] = table
        med = median_survival(table)
        last_ten = [row["time"] for row in table if row["n_risk"] >= 10]
        item = {"group": g, "n": len(gt), "events": sum(ge), "censored": len(gt) - sum(ge),
                "median_survival": med["median"], "median_ci": [med["ci_lower"], med["ci_upper"]],
                "max_time": max(gt), "last_time_with_10_at_risk": last_ten[-1] if last_ten else None}
        if report_times:
            item["survival_at"] = [survival_at(table, t) for t in report_times]
        summaries.append(item)
        if sum(ge) < 10:
            warnings.append(f"group '{g}' has {sum(ge)} events (fewer than 10): its curve and any "
                            "comparison with it are imprecise")
    test = logrank_test(times, events, groups, order=levels)
    zero_times = sum(1 for t in times if t == 0)
    if zero_times:
        warnings.append(f"{zero_times} row(s) have time 0: check the time origin and whether they should be kept")
    if len(levels) > 1 and min(s["n"] for s in summaries) < 0.1 * len(times):
        warnings.append("the smallest group holds less than 10% of the patients")
    return {
        "n": len(times), "events": sum(events),
        "follow_up": {"median_reverse_km": reverse_km_median(times, events), "max_observed_time": max(times)},
        "conf_level": conf_level, "conf_type": conf_type,
        "groups": summaries, "logrank": test, "warnings": warnings, "_curves": curves,
    }


def write_curves(path: str, curves: dict, order) -> None:
    columns = ["group", "time", "n_risk", "n_event", "n_censor", "survival", "se", "ci_lower", "ci_upper"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for g in order:
            for row in curves[g]:
                writer.writerow([g] + ["NA" if row[c] is None else (repr(row[c]) if isinstance(row[c], float) else row[c])
                                       for c in columns[1:]])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Kaplan-Meier estimates with Greenwood confidence intervals, median survival per group, "
                    "reverse Kaplan-Meier median follow-up and the log-rank test.",
        epilog="Input: a TSV (or .csv) with one row per patient. The event column must be coded 1 = event, "
               "0 = censored. Rows with a missing time, event or group are dropped and counted in the output. "
               "Scores are split only at prespecified points; no cut-point search is offered.")
    parser.add_argument("--input", required=True, help="patient table (tab-separated, or comma-separated if it ends in .csv)")
    parser.add_argument("--time", required=True, help="column with follow-up time (any unit, >= 0)")
    parser.add_argument("--event", required=True, help="column with the event indicator: 1 = event, 0 = censored")
    parser.add_argument("--id", help="patient ID column; stops with exit code 2 if an ID occurs on more than one row")
    parser.add_argument("--group", help="column with group labels (omit for a single curve)")
    parser.add_argument("--score", help="continuous column to split into groups; needs --split or --cut")
    parser.add_argument("--split", choices=["median", "tertile"], help="prespecified split of --score")
    parser.add_argument("--cut", type=float, help="split --score at this value, fixed before looking at outcome (high if > cut)")
    parser.add_argument("--times", help="comma-separated times at which to report survival and number at risk, e.g. 12,24,60")
    parser.add_argument("--conf-level", type=float, default=0.95, help="confidence level (default 0.95)")
    parser.add_argument("--conf-type", choices=["log", "log-log", "plain"], default="log",
                        help="confidence interval transformation (default log, as R survfit)")
    parser.add_argument("--delimiter", help="field delimiter if not tab / .csv comma")
    parser.add_argument("--curve-out", help="also write the full step-function table (TSV) to this file, for plotting")
    parser.add_argument("--out", help="write the JSON report to this file instead of stdout")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.group and args.score:
            raise InputError("use either --group or --score, not both")
        if args.score and (args.split is None) == (args.cut is None):
            raise InputError("--score needs exactly one of --split or --cut")
        if not args.score and (args.split or args.cut is not None):
            raise InputError("--split and --cut need --score")
        report_times = []
        if args.times:
            try:
                report_times = [float(x) for x in args.times.split(",") if x.strip()]
            except ValueError:
                raise InputError(f"--times must be comma-separated numbers, got '{args.times}'") from None
            if any(not math.isfinite(x) or x < 0 for x in report_times):
                raise InputError(f"--times must be finite and >= 0, got '{args.times}'")
        header, rows = read_table(args.input, args.delimiter)
        records, dropped = parse_records(header, rows, args.time, args.event, args.group, args.score, args.id)
        order, grouping = None, {"column": args.group, "method": "column" if args.group else "none"}
        if args.score:
            if not records:
                raise InputError("no complete rows: every row has a missing time, event or score")
            labels, order, info = split_score([r[2] for r in records], "cut" if args.cut is not None else args.split, args.cut)
            records = [(t, e, lab) for (t, e, _), lab in zip(records, labels)]
            grouping = {"column": args.score, **info}
        report = analyse(records, order, args.conf_level, args.conf_type, report_times)
        curves = report.pop("_curves")
        report = {"input": {"file": args.input, "rows": len(rows), "rows_dropped_missing": dropped,
                            "id_column_checked": args.id},
                  "endpoint": {"time": args.time, "event": args.event}, "grouping": grouping, **report}
        if dropped:
            report["warnings"].append(f"{dropped} row(s) dropped for a missing value in a required column")
        text = json.dumps(report, indent=2)
        try:
            if args.curve_out:
                write_curves(args.curve_out, curves, [g["group"] for g in report["groups"]])
            if args.out:
                with open(args.out, "w", encoding="utf-8") as handle:
                    handle.write(text + "\n")
        except OSError as exc:
            raise InputError(f"cannot write {exc.filename}: {exc.strerror or exc}") from exc
        if not args.out:
            print(text)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
