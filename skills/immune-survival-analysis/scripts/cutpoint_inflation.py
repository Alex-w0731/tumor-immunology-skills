#!/usr/bin/env python3
"""How often does the 'optimal' cut-point of a marker reach p < alpha when the marker is noise?

Standard library only; seeded, so the same command gives the same numbers.

Usage examples:
  python3 cutpoint_inflation.py --n 150 --censoring 0.4 --n-sim 1000 --seed 1
  python3 cutpoint_inflation.py --input cohort.tsv --time os_months --event os_event \
      --score cd8_score --n-perm 1000 --seed 1

Simulation mode (default). Each simulated cohort has --n patients whose marker is independent
of survival by construction: marker ~ N(0, 1), event time ~ Exponential(1), censoring time ~
Exponential(c / (1 - c)), so the expected censored fraction is c = --censoring. For every
cohort the script computes three tests and reports how often each gives p < alpha. Because the
null is true, that share is the false-positive rate; its Monte Carlo standard error is
sqrt(r (1 - r) / n_sim).
  best cut-point   the two-group log-rank test at every cut that leaves at least
                   --min-group-fraction of the patients in each group, keeping the smallest p
                   (the 'optimal cut-point' approach)                         [Altman1994]
  median split     one prespecified log-rank test, high = marker > median
  continuous       the log-rank (Cox partial-likelihood score) test of the untransformed
                   marker at beta = 0, 1 degree of freedom

Permutation mode (--input). For a real cohort the script finds the best cut-point of --score and
its naive log-rank p, then shuffles the score across patients --n-perm times (seeded) and
reports the search-adjusted p = (1 + #{permutations whose best chi-square >= observed best
chi-square}) / (n_perm + 1). This max-statistic permutation test corrects only for the search
over cut-points within the stated range; it does not validate the cut in new patients. A
closed-form alternative is the maximally selected rank statistic                [Lausen1992]

Statistics. Every log-rank chi-square uses the hypergeometric variance of R survdiff (rho = 0)
[Mantel1966] [Peto1972]; for any single cut it equals km_logrank.py. All cuts are evaluated in
O(n log n) through the identities (A = the 'high' group, t_j = distinct event times with d_j
events among n_j at risk, w_j = d_j (n_j - d_j) / (n_j - 1))
  O_A - E_A = sum over i in A of (delta_i - H(t_i)),   H(t) = sum over t_j <= t of d_j / n_j
  V_A = sum over i in A of c(t_i)  -  sum over i, l in A of g(min(t_i, t_l))
  c(t) = sum over t_j <= t of w_j / n_j,   g(t) = sum over t_j <= t of w_j / n_j^2
with the double sum updated by two Fenwick trees as patients join A in decreasing marker order.

Output: JSON on stdout (or --out). Exit code 2 on invalid input.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys

MISSING = {"", "na", "nan", "n/a", "null", "none", "."}


class InputError(ValueError):
    """Raised for problems in the input table or the options; main() turns it into exit code 2."""


def chi2_1df_sf(x) -> float:
    """Upper tail of a chi-square variable with 1 degree of freedom; None (undefined test) gives 1."""
    if x is None or x <= 0.0:
        return 1.0
    return math.erfc(math.sqrt(x / 2.0))


class _Fenwick:
    """Binary indexed tree: point add, prefix sum over positions < i."""

    def __init__(self, size: int) -> None:
        self.tree = [0.0] * (size + 1)

    def add(self, i: int, value: float) -> None:
        i += 1
        while i < len(self.tree):
            self.tree[i] += value
            i += i & -i

    def prefix(self, i: int) -> float:
        total = 0.0
        while i > 0:
            total += self.tree[i]
            i -= i & -i
        return total


def prepare(times, events) -> dict:
    """Per-patient terms of the log-rank decomposition (see the module docstring)."""
    if len(times) != len(events):
        raise ValueError("times and events must have the same length")
    distinct = sorted(set(times))
    rank_of = {t: r for r, t in enumerate(distinct)}
    deaths = [0] * len(distinct)
    leaving = [0] * len(distinct)
    for t, e in zip(times, events):
        leaving[rank_of[t]] += 1
        deaths[rank_of[t]] += 1 if e else 0
    hazard, c_cum, g_cum = [], [], []
    at_risk, h, c, g = len(times), 0.0, 0.0, 0.0
    for r in range(len(distinct)):
        d = deaths[r]
        if d:
            h += d / at_risk
            if at_risk > 1:
                w = d * (at_risk - d) / (at_risk - 1)
                c += w / at_risk
                g += w / (at_risk * at_risk)
        hazard.append(h)
        c_cum.append(c)
        g_cum.append(g)
        at_risk -= leaving[r]
    ranks = [rank_of[t] for t in times]
    return {"a": [(1.0 if e else 0.0) - hazard[r] for e, r in zip(events, ranks)],
            "c": [c_cum[r] for r in ranks], "g": [g_cum[r] for r in ranks],
            "rank": ranks, "size": len(distinct), "events": sum(1 for e in events if e)}


def scan_cuts(marker, prep) -> list:
    """Two-group log-rank chi-square for every cut 'high if marker > cut'.

    Returns one dict per cut that separates distinct marker values: cut, n_high, n_low,
    o_minus_e_high (observed minus expected events in the high group), chi2 (None when the
    variance is zero, for example when there are no events).
    """
    n = len(marker)
    a, c, g, rank = prep["a"], prep["c"], prep["g"], prep["rank"]
    order = sorted(range(n), key=lambda i: -marker[i])
    count_tree, g_tree = _Fenwick(prep["size"]), _Fenwick(prep["size"])
    sum_a = sum_c = pairs = 0.0
    out = []
    for m in range(1, n):
        i = order[m - 1]
        r = rank[i]
        earlier = count_tree.prefix(r)  # members of A with a strictly earlier time
        pairs += g[i] + 2.0 * (g[i] * ((m - 1) - earlier) + g_tree.prefix(r))
        count_tree.add(r, 1.0)
        g_tree.add(r, g[i])
        sum_a += a[i]
        sum_c += c[i]
        nxt = order[m]
        if marker[nxt] == marker[i]:
            continue  # tied marker values cannot be split by a cut
        var = sum_c - pairs
        chi2 = sum_a * sum_a / var if var > 1e-12 * max(1.0, sum_c) else None
        out.append({"cut": marker[nxt], "n_high": m, "n_low": n - m, "o_minus_e_high": sum_a, "chi2": chi2})
    return out


def min_group_size(n: int, fraction: float) -> int:
    return max(1, math.ceil(fraction * n - 1e-9))


def summarise_cuts(marker, prep, min_group: int) -> dict:
    """Best eligible cut (largest chi-square), the median split and the number of eligible cuts."""
    cuts = scan_cuts(marker, prep)
    eligible = [k for k in cuts if k["n_high"] >= min_group and k["n_low"] >= min_group]
    best = None
    for k in eligible:
        if k["chi2"] is not None and (best is None or k["chi2"] > best["chi2"]):
            best = k
    median = statistics.median(marker)
    n_above = sum(1 for v in marker if v > median)
    split = next((k for k in cuts if k["n_high"] == n_above), None)
    return {"best": best, "median_split": split, "median": median, "n_eligible": len(eligible)}


def score_test(times, events, x) -> dict:
    """Log-rank (Cox score) test of a continuous covariate at beta = 0, 1 degree of freedom.

    U = sum over event times of (sum of x over the d_j events - d_j * mean of x in the risk set);
    V = sum over event times of w_j * (variance of x in the risk set, divisor n_j). With a 0/1
    covariate this is exactly the two-group log-rank test. U > 0 means higher x, higher hazard.
    """
    n = len(times)
    if n == 0:
        return {"u": 0.0, "var": 0.0, "chi2": None, "p": 1.0}
    centre = sum(x) / n  # U and V do not change when x is shifted; centring protects precision
    xc = [v - centre for v in x]
    idx = sorted(range(n), key=lambda i: times[i], reverse=True)
    s1 = s2 = u = var = 0.0
    n_risk = k = 0
    while k < n:
        t = times[idx[k]]
        event_sum, d = 0.0, 0
        while k < n and times[idx[k]] == t:
            i = idx[k]
            s1 += xc[i]
            s2 += xc[i] * xc[i]
            n_risk += 1
            if events[i]:
                event_sum += xc[i]
                d += 1
            k += 1
        if d:
            mean = s1 / n_risk
            u += event_sum - d * mean
            if n_risk > 1:
                var += d * (n_risk - d) / (n_risk - 1) * max(0.0, s2 / n_risk - mean * mean)
    chi2 = u * u / var if var > 1e-12 else None
    return {"u": u, "var": var, "chi2": chi2, "p": chi2_1df_sf(chi2)}


def _rate(hits: int, total: int) -> dict:
    r = hits / total
    return {"rate": r, "mc_se": math.sqrt(r * (1.0 - r) / total)}


def simulate(n: int, n_sim: int, censoring: float, min_fraction: float, alpha: float, seed: int) -> dict:
    """Null simulation: false-positive rates of the best cut, the median split and the continuous test."""
    rng = random.Random(seed)
    censor_rate = censoring / (1.0 - censoring) if censoring > 0 else 0.0
    min_group = min_group_size(n, min_fraction)
    hits = {"best": 0, "median": 0, "continuous": 0}
    min_p, events, n_cuts = [], [], []
    for _ in range(n_sim):
        marker = [rng.gauss(0.0, 1.0) for _ in range(n)]
        event_time = [rng.expovariate(1.0) for _ in range(n)]
        censor_time = [rng.expovariate(censor_rate) if censor_rate > 0 else math.inf for _ in range(n)]
        times = [min(t, c) for t, c in zip(event_time, censor_time)]
        status = [1 if t <= c else 0 for t, c in zip(event_time, censor_time)]
        prep = prepare(times, status)
        res = summarise_cuts(marker, prep, min_group)
        p_best = chi2_1df_sf(res["best"]["chi2"]) if res["best"] else 1.0
        p_median = chi2_1df_sf(res["median_split"]["chi2"]) if res["median_split"] else 1.0
        p_cont = score_test(times, status, marker)["p"]
        hits["best"] += p_best < alpha
        hits["median"] += p_median < alpha
        hits["continuous"] += p_cont < alpha
        min_p.append(p_best)
        events.append(prep["events"])
        n_cuts.append(res["n_eligible"])
    min_p.sort()
    return {
        "mode": "simulation",
        "settings": {"n": n, "n_sim": n_sim, "censoring_target": censoring, "min_group_fraction": min_fraction,
                     "min_group_size": min_group, "alpha": alpha, "seed": seed},
        "mean_events_per_cohort": statistics.fmean(events),
        "median_candidate_cuts": statistics.median(n_cuts),
        "false_positive_rate": {"best_cut_point": _rate(hits["best"], n_sim),
                                "median_split": _rate(hits["median"], n_sim),
                                "continuous_score_test": _rate(hits["continuous"], n_sim)},
        "best_cut_p_quantiles": {"q10": min_p[int(0.10 * (n_sim - 1))], "q25": min_p[int(0.25 * (n_sim - 1))],
                                 "q50": min_p[int(0.50 * (n_sim - 1))]},
        "note": "Synthetic null data: the marker is independent of survival by construction, so every "
                "p < alpha is a false positive. These rates describe the method, not any real cohort.",
    }


def permutation_test(times, events, marker, n_perm: int, min_fraction: float, seed: int) -> dict:
    """Best cut-point of a real marker with its naive and its search-adjusted (permutation) p-value."""
    n = len(times)
    min_group = min_group_size(n, min_fraction)
    if 2 * min_group > n:
        raise InputError(f"{n} patients are too few for groups of at least {min_group}")
    prep = prepare(times, events)
    if prep["events"] == 0:
        raise InputError("no events: a log-rank test is undefined")
    observed = summarise_cuts(marker, prep, min_group)
    best = observed["best"]
    if best is None:
        raise InputError("no eligible cut-point: too many tied scores for the minimum group size")
    rng = random.Random(seed)
    shuffled = list(marker)
    exceed = 0
    for _ in range(n_perm):
        rng.shuffle(shuffled)
        perm_best = summarise_cuts(shuffled, prep, min_group)["best"]
        if perm_best is not None and perm_best["chi2"] >= best["chi2"] * (1.0 - 1e-9):
            exceed += 1
    split = observed["median_split"]
    cont = score_test(times, events, marker)
    warnings = []
    if prep["events"] < 20:
        warnings.append(f"only {prep['events']} events: every test here is imprecise")
    return {
        "mode": "permutation",
        "n": n, "events": prep["events"],
        "settings": {"n_perm": n_perm, "min_group_fraction": min_fraction, "min_group_size": min_group, "seed": seed},
        "best_cut": {"cut": best["cut"], "rule": "high if score > cut", "n_high": best["n_high"], "n_low": best["n_low"],
                     "chi2": best["chi2"], "p_naive": chi2_1df_sf(best["chi2"]),
                     "high_group_observed_minus_expected": best["o_minus_e_high"],
                     "eligible_cuts_searched": observed["n_eligible"]},
        "search_adjusted_p": (1 + exceed) / (n_perm + 1),
        "median_split": None if split is None else {
            "cut": observed["median"], "n_high": split["n_high"], "n_low": split["n_low"], "chi2": split["chi2"],
            "p": chi2_1df_sf(split["chi2"]), "high_group_observed_minus_expected": split["o_minus_e_high"]},
        "continuous_score_test": {"chi2": cont["chi2"], "p": cont["p"],
                                  "direction": "higher score, higher hazard" if cont["u"] > 0 else "higher score, lower hazard"},
        "warnings": warnings,
    }


def read_cohort(path: str, time_col: str, event_col: str, score_col: str, delimiter=None) -> tuple:
    """Read (times, events, scores, rows_dropped_missing) from a delimited patient table."""
    if delimiter is None:
        delimiter = "," if path.lower().endswith(".csv") else "\t"
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            header = list(reader.fieldnames or [])
            rows = [dict(row) for row in reader]
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc.strerror or exc}") from exc
    for col in (time_col, event_col, score_col):
        if col not in header:
            raise InputError(f"column '{col}' not found; available columns: {', '.join(header) or '(none)'}")
    times, events, scores, dropped = [], [], [], 0
    for number, row in enumerate(rows, start=2):
        cells = [row.get(c) for c in (time_col, event_col, score_col)]
        if any(c is None or c.strip().lower() in MISSING for c in cells):
            dropped += 1
            continue
        try:
            t = float(cells[0])
        except ValueError:
            raise InputError(f"line {number}: time '{cells[0]}' in column '{time_col}' is not a number") from None
        try:
            s = float(cells[2])
        except ValueError:
            raise InputError(f"line {number}: score '{cells[2]}' in column '{score_col}' is not a number") from None
        if not (math.isfinite(t) and t >= 0):
            raise InputError(f"line {number}: time {cells[0]} in column '{time_col}' must be finite and >= 0")
        if not math.isfinite(s):
            raise InputError(f"line {number}: score in column '{score_col}' must be finite")
        if cells[1].strip() not in ("0", "1", "0.0", "1.0"):
            raise InputError(f"line {number}: event '{cells[1]}' in column '{event_col}' must be 1 (event) or 0 (censored)")
        times.append(t)
        events.append(1 if cells[1].strip() in ("1", "1.0") else 0)
        scores.append(s)
    if not times:
        raise InputError("no complete rows: every row has a missing time, event or score")
    return times, events, scores, dropped


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="False-positive rate of the 'optimal' survival cut-point under the null (simulation), or the "
                    "search-adjusted permutation p-value of the best cut-point in a real cohort (--input).",
        epilog="All randomness comes from --seed. Simulated numbers describe the method, never a real cohort.")
    parser.add_argument("--n", type=int, default=150, help="patients per simulated cohort (default 150)")
    parser.add_argument("--censoring", type=float, default=0.4, help="expected censored fraction, 0 to 0.9 (default 0.4)")
    parser.add_argument("--n-sim", type=int, default=1000, help="simulated cohorts (default 1000)")
    parser.add_argument("--alpha", type=float, default=0.05, help="significance level (default 0.05)")
    parser.add_argument("--min-group-fraction", type=float, default=0.1,
                        help="smallest share of patients allowed in either group when searching cuts (default 0.1)")
    parser.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    parser.add_argument("--input", help="patient table for permutation mode (TSV, or CSV if it ends in .csv)")
    parser.add_argument("--time", help="permutation mode: follow-up time column")
    parser.add_argument("--event", help="permutation mode: event column, 1 = event, 0 = censored")
    parser.add_argument("--score", help="permutation mode: continuous marker column")
    parser.add_argument("--n-perm", type=int, default=1000, help="permutation mode: number of permutations (default 1000)")
    parser.add_argument("--delimiter", help="field delimiter if not tab / .csv comma")
    parser.add_argument("--out", help="write the JSON report here instead of stdout")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if not 0.0 < args.min_group_fraction < 0.5:
            raise InputError(f"--min-group-fraction must be between 0 and 0.5, got {args.min_group_fraction}")
        if args.input:
            if not (args.time and args.event and args.score):
                raise InputError("permutation mode needs --time, --event and --score")
            if not 1 <= args.n_perm <= 1_000_000:
                raise InputError(f"--n-perm must be between 1 and 1000000, got {args.n_perm}")
            times, events, scores, dropped = read_cohort(args.input, args.time, args.event, args.score, args.delimiter)
            report = permutation_test(times, events, scores, args.n_perm, args.min_group_fraction, args.seed)
            report["input"] = {"file": args.input, "time": args.time, "event": args.event, "score": args.score,
                               "rows_dropped_missing": dropped}
        else:
            if args.n < 10:
                raise InputError(f"--n must be at least 10, got {args.n}")
            if not 0.0 <= args.censoring <= 0.9:
                raise InputError(f"--censoring must be between 0 and 0.9, got {args.censoring}")
            if not 0.0 < args.alpha < 1.0:
                raise InputError(f"--alpha must be between 0 and 1, got {args.alpha}")
            if not 1 <= args.n_sim <= 1_000_000:
                raise InputError(f"--n-sim must be between 1 and 1000000, got {args.n_sim}")
            if 2 * min_group_size(args.n, args.min_group_fraction) > args.n:
                raise InputError("--min-group-fraction leaves no room for a cut at this --n")
            report = simulate(args.n, args.n_sim, args.censoring, args.min_group_fraction, args.alpha, args.seed)
        text = json.dumps(report, indent=2)
        if args.out:
            try:
                with open(args.out, "w", encoding="utf-8") as handle:
                    handle.write(text + "\n")
            except OSError as exc:
                raise InputError(f"cannot write {args.out}: {exc.strerror or exc}") from exc
        else:
            print(text)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
