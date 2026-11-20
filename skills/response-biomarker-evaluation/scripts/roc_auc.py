#!/usr/bin/env python3
"""ROC AUC of a response biomarker, with uncertainty, per-cohort and paired comparisons.

Standard library only; deterministic for a given --seed. Research use only.

What is computed
  * AUC = U / (n1 * n0). U is the Mann-Whitney statistic: over every (responder,
    non-responder) pair it adds 1 when the responder's score lies on the response side and
    1/2 when the two scores tie [Hanley1982]. Label 1 = responder, as in tmelens
    `biomarker_auc`. --direction says which side was prespecified: higher (default) or
    lower scores in responders. The AUC is never flipped after looking at the data.
  * ci_low / ci_high: percentile bootstrap [Efron1986]. Responders and non-responders are
    resampled separately (within each cohort for the within-cohort AUC), so the class sizes
    stay fixed, as in tmelens. Quantiles use linear interpolation (numpy default, R type 7).
  * p_perm: two-sided permutation test of AUC = 0.5; labels are shuffled (within cohort for
    the within-cohort AUC) and p = (b + 1) / (m + 1), where b counts the m permutations whose
    |AUC - 0.5| is at least the observed one [Phipson2010].
  * se_delong: nonparametric standard error of the AUC from placement values [DeLong1988].
  * --cohort: AUC per cohort, the pooled AUC (all patients, cohort ignored) and the
    within-cohort AUC = sum_k U_k / sum_k n1_k * n0_k, which compares each responder only
    with non-responders of the same cohort.
  * --score2: paired comparison of two markers on the same patients: AUC difference with a
    paired stratified bootstrap CI and the DeLong z-test for correlated AUCs [DeLong1988].
    With --cohort it is also given per cohort and within cohort: the within-cohort
    difference is sum_k w_k (AUC1_k - AUC2_k) with w_k = n1_k n0_k / sum_j n1_j n0_j, its
    variance sum_k w_k^2 Var_k from the per-cohort paired DeLong variances (cohorts are
    independent), and its bootstrap resamples within cohort and class. This stratified
    combination is this skill's own extension of [DeLong1988], not a published test.
  * --plan: Hanley-McNeil standard error for an anticipated AUC and class sizes
    [Hanley1982], to show the precision a cohort can deliver before it is analysed.
  * warnings (heuristic flags): fewer than 10 patients in the smaller class overall or in a
    cohort; complete separation (AUC 0 or 1), where the bootstrap interval and se_delong
    collapse to a point; labelled patients without a score, by class; heavy ties; a pooled
    AUC more than 0.05 from the within-cohort AUC.

Input: TSV (comma-separated if the file name ends in .csv), one row per patient. The label
column holds 1/0, or categories mapped with --positive and --negative; rows matching
neither list are excluded and counted, never silently recoded. Output: JSON.

Usage
  python3 roc_auc.py cohort.tsv --score ifng6 --label responder --seed 1
  python3 roc_auc.py cohort.tsv --score tide --direction lower --label bor \
      --positive CR,PR --negative SD,PD --cohort study --score2 cd274 --seed 7 --out auc.json
  python3 roc_auc.py --plan --auc 0.70 --n-pos 15 --n-neg 45
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import random
import sys
from statistics import NormalDist

MISSING = {"", "na", "nan", "n/a", "null", "none", "."}
MIN_PER_CLASS = 3          # as in tmelens.biomarker_auc
SMALL_CLASS = 10           # heuristic: below this, the AUC is called exploratory


class InputError(Exception):
    """Unusable input; main() prints the message and exits with code 2."""


# ----------------------------------------------------------------------------- statistics

def mann_whitney_u(pos, neg):
    """U of `pos` over `neg`: every pair with pos > neg adds 1, every tie adds 1/2."""
    if not pos or not neg:
        raise ValueError("both groups need at least one value")
    srt = sorted(neg)
    u = 0.0
    for x in pos:
        lo = bisect.bisect_left(srt, x)
        hi = bisect.bisect_right(srt, x)
        u += lo + 0.5 * (hi - lo)
    return u


def auc(pos, neg):
    """Area under the empirical ROC curve, U / (n1 * n0) [Hanley1982]."""
    return mann_whitney_u(pos, neg) / (len(pos) * len(neg))


def midranks(values):
    """Ranks 1..n; tied values share the mean of the positions they occupy."""
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def within_auc(strata):
    """sum_k U_k / sum_k n1_k n0_k over (pos, neg) strata; one stratum gives the plain AUC."""
    u = sum(mann_whitney_u(p, q) for p, q in strata)
    return u / sum(len(p) * len(q) for p, q in strata)


def placements(pos, neg):
    """DeLong structural components: V10 per responder, V01 per non-responder [DeLong1988]."""
    n1, n0 = len(pos), len(neg)
    sneg, spos = sorted(neg), sorted(pos)
    v10 = []
    for x in pos:
        lo, hi = bisect.bisect_left(sneg, x), bisect.bisect_right(sneg, x)
        v10.append((lo + 0.5 * (hi - lo)) / n0)
    v01 = []
    for y in neg:
        lo, hi = bisect.bisect_left(spos, y), bisect.bisect_right(spos, y)
        v01.append((n1 - hi + 0.5 * (hi - lo)) / n1)
    return v10, v01


def _cov(a, b):
    n = len(a)
    if n < 2:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (n - 1)


def delong_se(pos, neg):
    """DeLong standard error of one AUC: sqrt(S10 / n1 + S01 / n0)."""
    v10, v01 = placements(pos, neg)
    var = _cov(v10, v10) / len(pos) + _cov(v01, v01) / len(neg)
    return math.sqrt(var) if var >= 0 else float("nan")


def _delong_parts(pos_pairs, neg_pairs):
    """AUC1, AUC2, Var(AUC1), Var(AUC2), Cov(AUC1, AUC2) from placement values [DeLong1988]."""
    a10, a01 = placements([a for a, _ in pos_pairs], [a for a, _ in neg_pairs])
    b10, b01 = placements([b for _, b in pos_pairs], [b for _, b in neg_pairs])
    n1, n0 = len(pos_pairs), len(neg_pairs)
    var1 = _cov(a10, a10) / n1 + _cov(a01, a01) / n0
    var2 = _cov(b10, b10) / n1 + _cov(b01, b01) / n0
    cov12 = _cov(a10, b10) / n1 + _cov(a01, b01) / n0
    return sum(a10) / n1, sum(b10) / n1, var1, var2, cov12


def _z_test(out, var_diff):
    """Add se_diff, z and the two-sided normal p to `out` when the variance is positive."""
    if var_diff > 1e-12:
        se = math.sqrt(var_diff)
        z = out["diff"] / se
        out.update(se_diff=se, z=z, p_delong=math.erfc(abs(z) / math.sqrt(2)))
    return out


def delong_paired(pos_pairs, neg_pairs):
    """DeLong test for two correlated AUCs; each pair holds (marker 1, marker 2) of one patient."""
    auc1, auc2, var1, var2, cov12 = _delong_parts(pos_pairs, neg_pairs)
    out = {"auc1": auc1, "auc2": auc2, "diff": auc1 - auc2,
           "se_diff": None, "z": None, "p_delong": None, "corr_auc": None}
    if var1 > 0 and var2 > 0:
        out["corr_auc"] = cov12 / math.sqrt(var1 * var2)
    return _z_test(out, var1 + var2 - 2 * cov12)


def delong_paired_within(strata_pairs):
    """Within-cohort AUC difference and its stratified z-test over (pos_pairs, neg_pairs) strata.

    diff = sum_k w_k (AUC1_k - AUC2_k), w_k = n1_k n0_k / sum_j n1_j n0_j (the same weights as
    within_auc); Var = sum_k w_k^2 Var_k(paired DeLong), cohorts being independent. This
    skill's extension of [DeLong1988]; each stratum needs at least 2 patients per class.
    """
    total = sum(len(p) * len(q) for p, q in strata_pairs)
    auc1 = auc2 = var = 0.0
    for p, q in strata_pairs:
        w = len(p) * len(q) / total
        a1, a2, v1, v2, c12 = _delong_parts(p, q)
        auc1, auc2, var = auc1 + w * a1, auc2 + w * a2, var + w * w * (v1 + v2 - 2 * c12)
    out = {"auc1": auc1, "auc2": auc2, "diff": auc1 - auc2, "se_diff": None, "z": None, "p_delong": None}
    return _z_test(out, var)


def hanley_mcneil_se(a, n_pos, n_neg):
    """SE of an AUC under the Hanley-McNeil approximation, Q1 = A/(2-A), Q2 = 2A^2/(1+A)."""
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n_pos - 1) * (q1 - a * a) + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg)
    return math.sqrt(max(var, 0.0))


def _quantile(sorted_vals, q):
    """Linear-interpolation quantile of an already sorted list (R type 7)."""
    h = (len(sorted_vals) - 1) * q
    lo = math.floor(h)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (h - lo) * (sorted_vals[hi] - sorted_vals[lo])


def _resample(rng, values):
    n = len(values)
    return [values[min(int(rng.random() * n), n - 1)] for _ in range(n)]


def _sample_indices(rng, n, k):
    """k distinct indices out of n (partial Fisher-Yates driven by rng.random only)."""
    pool = list(range(n))
    for i in range(k):
        j = i + min(int(rng.random() * (n - i)), n - i - 1)
        pool[i], pool[j] = pool[j], pool[i]
    return pool[:k]


def bootstrap_ci(strata, n_boot, alpha, rng):
    """Percentile CI of within_auc, resampling responders and non-responders within each stratum."""
    stats = sorted(within_auc([(_resample(rng, p), _resample(rng, q)) for p, q in strata])
                   for _ in range(n_boot))
    return _quantile(stats, alpha / 2), _quantile(stats, 1 - alpha / 2)


def bootstrap_diff_ci(strata_pairs, n_boot, alpha, rng):
    """Percentile CI of the (within-cohort) AUC1 - AUC2, resampling whole patients within each
    class of each (pos_pairs, neg_pairs) stratum; one stratum gives the plain paired bootstrap."""
    diffs = []
    for _ in range(n_boot):
        rs = [(_resample(rng, p), _resample(rng, q)) for p, q in strata_pairs]
        diffs.append(within_auc([([a for a, _ in p], [a for a, _ in q]) for p, q in rs])
                     - within_auc([([b for _, b in p], [b for _, b in q]) for p, q in rs]))
    diffs.sort()
    return _quantile(diffs, alpha / 2), _quantile(diffs, 1 - alpha / 2)


def permutation_p(strata, n_perm, rng):
    """Two-sided permutation p of AUC = 0.5, labels shuffled within each stratum [Phipson2010].

    The statistic |2 * sum U - sum n1 n0| is a whole number, so comparisons are exact.
    """
    prepared, u2_obs, total = [], 0.0, 0
    for p, q in strata:
        ranks = midranks(list(p) + list(q))
        n1 = len(p)
        prepared.append((ranks, n1))
        u2_obs += 2 * sum(ranks[:n1]) - n1 * (n1 + 1)
        total += n1 * len(q)
    d_obs = abs(u2_obs - total)
    b = 0
    for _ in range(n_perm):
        u2 = 0.0
        for ranks, n1 in prepared:
            u2 += 2 * sum(ranks[i] for i in _sample_indices(rng, len(ranks), n1)) - n1 * (n1 + 1)
        if abs(u2 - total) >= d_obs - 1e-9:
            b += 1
    return (b + 1) / (n_perm + 1)


def summarise(strata, n_boot, n_perm, alpha, seed, tag):
    """AUC, bootstrap CI and permutation p for a list of (pos, neg) strata."""
    lo, hi = bootstrap_ci(strata, n_boot, alpha, random.Random(f"{seed}/{tag}/boot"))
    out = {"n_pos": sum(len(p) for p, _ in strata), "n_neg": sum(len(q) for _, q in strata),
           "auc": within_auc(strata), "ci_low": lo, "ci_high": hi,
           "p_perm": permutation_p(strata, n_perm, random.Random(f"{seed}/{tag}/perm"))}
    if len(strata) == 1:
        out["se_delong"] = delong_se(*strata[0])
    return out


def compare(pos_pairs, neg_pairs, n_boot, alpha, seed, tag):
    """Paired comparison of two markers: DeLong test plus paired bootstrap CI of the difference."""
    out = {"n_pos": len(pos_pairs), "n_neg": len(neg_pairs)}
    out.update(delong_paired(pos_pairs, neg_pairs))
    lo, hi = bootstrap_diff_ci([(pos_pairs, neg_pairs)], n_boot, alpha, random.Random(f"{seed}/{tag}/diff"))
    out.update(diff_ci_low=lo, diff_ci_high=hi)
    return out


def compare_within(strata_pairs, n_boot, alpha, seed, tag):
    """Within-cohort paired comparison: stratified DeLong z-test plus stratified paired bootstrap CI."""
    out = {"n_pos": sum(len(p) for p, _ in strata_pairs), "n_neg": sum(len(q) for _, q in strata_pairs),
           "cohorts_used": len(strata_pairs)}
    out.update(delong_paired_within(strata_pairs))
    lo, hi = bootstrap_diff_ci(strata_pairs, n_boot, alpha, random.Random(f"{seed}/{tag}/diff"))
    out.update(diff_ci_low=lo, diff_ci_high=hi,
               note="patients compared only within their own cohort; stratified DeLong is this skill's extension")
    return out


def _flags(summary, where, warnings):
    """Heuristic warnings for one AUC summary: small smaller class, complete separation."""
    small = min(summary["n_pos"], summary["n_neg"])
    if small < SMALL_CLASS:
        warnings.append(f"{where}: only {small} patients in the smaller class (heuristic flag below "
                        f"{SMALL_CLASS}); asymptotic intervals for one AUC under-cover in small samples and "
                        "no small-sample interval method is uniformly best [Obuchowski1998]: quote "
                        "se_delong beside the bootstrap interval and treat the AUC as exploratory")
    if summary["auc"] in (0.0, 1.0):
        warnings.append(f"{where}: complete separation (AUC = {summary['auc']:.0f}); the bootstrap interval "
                        "and se_delong collapse to a point and do not show the real uncertainty: report "
                        "the counts and treat the AUC as exploratory")


# ----------------------------------------------------------------------------- input

def _missing(text):
    return text is None or text.strip().casefold() in MISSING


def _float(text, column, line):
    try:
        value = float(text)
    except ValueError:
        raise InputError(f"column '{column}' line {line}: '{text}' is not a number") from None
    if not math.isfinite(value):
        raise InputError(f"column '{column}' line {line}: non-finite value '{text}'")
    return value


def _categories(text):
    return {v.strip().casefold() for v in text.split(",") if v.strip()} if text else set()


def read_table(path):
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
    return header, rows


def build_records(header, rows, score, label, positive=None, negative=None, cohort=None, score2=None):
    """Rows -> records {y, s, s2, c}; returns (records, exclusion counts)."""
    for col in (score, label, cohort, score2):
        if col and col not in header:
            raise InputError(f"column '{col}' not found; columns are: {', '.join(header)}")
    pos_set, neg_set = _categories(positive), _categories(negative)
    if bool(pos_set) != bool(neg_set):
        raise InputError("give both --positive and --negative, or neither (then labels must be 1/0)")
    if pos_set & neg_set:
        raise InputError(f"label values in both --positive and --negative: {sorted(pos_set & neg_set)}")
    records, excluded, missing = [], {}, {"responder": 0, "non_responder": 0}
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
                                 "map categories with --positive and --negative (for example CR,PR and SD,PD)")
            y = int(value)
        if _missing(row.get(score)):
            missing["responder" if y == 1 else "non_responder"] += 1
            continue
        rec = {"y": y, "s": _float(row[score].strip(), score, line), "s2": None, "c": None}
        if score2 and not _missing(row.get(score2)):
            rec["s2"] = _float(row[score2].strip(), score2, line)
        if cohort:
            if _missing(row.get(cohort)):
                raise InputError(f"column '{cohort}' line {line}: empty cohort value")
            rec["c"] = row[cohort].strip()
        records.append(rec)
    return records, {"excluded_by_label": excluded, "excluded_missing_score": sum(missing.values()),
                     "missing_score_by_class": missing}


# ----------------------------------------------------------------------------- analysis

def _sign(direction):
    return 1.0 if direction == "higher" else -1.0


def _groups(recs, sign):
    return ([sign * r["s"] for r in recs if r["y"] == 1], [sign * r["s"] for r in recs if r["y"] == 0])


def _pairs(recs, s1, s2):
    cc = [r for r in recs if r["s2"] is not None]
    return ([(s1 * r["s"], s2 * r["s2"]) for r in cc if r["y"] == 1],
            [(s1 * r["s"], s2 * r["s2"]) for r in cc if r["y"] == 0])


def analyse(records, direction="higher", direction2="higher", use_cohort=False, use_score2=False,
            n_boot=2000, n_perm=5000, seed=1, alpha=0.05):
    sign, sign2 = _sign(direction), _sign(direction2)
    pos, neg = _groups(records, sign)
    if len(pos) < MIN_PER_CLASS or len(neg) < MIN_PER_CLASS:
        raise InputError(f"need at least {MIN_PER_CLASS} responders and {MIN_PER_CLASS} non-responders "
                         f"after exclusions; found {len(pos)} and {len(neg)}")
    warnings = []
    overall = summarise([(pos, neg)], n_boot, n_perm, alpha, seed, "overall")
    side = "higher" if direction == "higher" else "lower"
    if overall["auc"] > 0.5:
        overall["observed_direction"] = f"{side} in responders (as prespecified)"
    elif overall["auc"] < 0.5:
        overall["observed_direction"] = f"opposite to the prespecified direction ({side} in responders)"
    else:
        overall["observed_direction"] = "no separation"
    result = {"overall": overall}
    _flags(overall, "pooled" if use_cohort else "overall", warnings)
    scores = [r["s"] for r in records]
    if len(set(scores)) < 0.5 * len(scores):
        warnings.append(f"{len(set(scores))} distinct scores among {len(scores)} patients: ties count 1/2, "
                        "and the ROC curve has few operating points")
    if use_cohort:
        overall["note"] = "pools all patients and ignores cohort"
        by_cohort, strata, no_pairs = [], [], []
        for c in sorted({r["c"] for r in records}):
            sub = [r for r in records if r["c"] == c]
            p, q = _groups(sub, sign)
            if p and q:
                strata.append((p, q))
            else:
                no_pairs.append(c)
            if len(p) >= MIN_PER_CLASS and len(q) >= MIN_PER_CLASS:
                row = {"cohort": c}
                row.update(summarise([(p, q)], n_boot, n_perm, alpha, seed, f"cohort:{c}"))
                _flags(row, f"cohort {c}", warnings)
            else:
                row = {"cohort": c, "n_pos": len(p), "n_neg": len(q), "auc": None,
                       "note": f"fewer than {MIN_PER_CLASS} responders or non-responders: not estimated"}
            by_cohort.append(row)
        result["by_cohort"] = by_cohort
        if strata:
            within = summarise(strata, n_boot, n_perm, alpha, seed, "within")
            within["cohorts_used"] = len(strata)
            within["note"] = "responders compared only with non-responders of the same cohort"
            result["within_cohort"] = within
            gap = overall["auc"] - within["auc"]
            if abs(gap) > 0.05:
                warnings.append(f"pooled AUC differs from the within-cohort AUC by {gap:+.3f}: cohort "
                                "differences in score level or response rate distort the pooled value")
        if no_pairs:
            warnings.append(f"cohorts without both classes contribute no pairs: {', '.join(no_pairs)}")
    if use_score2:
        pp, nn = _pairs(records, sign, sign2)
        if len(pp) < MIN_PER_CLASS or len(nn) < MIN_PER_CLASS:
            raise InputError(f"--score2: need at least {MIN_PER_CLASS} responders and non-responders with "
                             f"both markers; found {len(pp)} and {len(nn)}")
        comp = {"overall": compare(pp, nn, n_boot, alpha, seed, "compare")}
        if len(pp) + len(nn) < len(records):
            warnings.append(f"--score2 comparison uses the {len(pp) + len(nn)} patients with both markers")
        if use_cohort:
            comp["overall"]["note"] = "pools all patients and ignores cohort"
            rows, strata_pairs = [], []
            for c in sorted({r["c"] for r in records}):
                p2, n2 = _pairs([r for r in records if r["c"] == c], sign, sign2)
                if len(p2) >= MIN_PER_CLASS and len(n2) >= MIN_PER_CLASS:
                    row = {"cohort": c}
                    row.update(compare(p2, n2, n_boot, alpha, seed, f"compare:{c}"))
                    strata_pairs.append((p2, n2))
                else:
                    row = {"cohort": c, "n_pos": len(p2), "n_neg": len(n2), "diff": None,
                           "note": f"fewer than {MIN_PER_CLASS} responders or non-responders: not estimated"}
                rows.append(row)
            comp["by_cohort"] = rows
            if strata_pairs:
                comp["within_cohort"] = compare_within(strata_pairs, n_boot, alpha, seed, "compare:within")
                if len(strata_pairs) < len(rows):
                    warnings.append("within-cohort comparison uses only the cohorts with at least "
                                    f"{MIN_PER_CLASS} responders and non-responders with both markers")
        result["comparison"] = comp
    result["warnings"] = warnings
    return result


def plan(a, n_pos, n_neg, alpha=0.05):
    if a is None or n_pos is None or n_neg is None:
        raise InputError("--plan needs --auc, --n-pos and --n-neg")
    if not 0 < a < 1:
        raise InputError("--auc must lie strictly between 0 and 1")
    if n_pos < 1 or n_neg < 1:
        raise InputError("--n-pos and --n-neg must be positive")
    se = hanley_mcneil_se(a, n_pos, n_neg)
    z = NormalDist().inv_cdf(1 - alpha / 2)
    return {"auc": a, "n_pos": n_pos, "n_neg": n_neg, "se_hanley_mcneil": se,
            "approx_ci_low": max(0.0, a - z * se), "approx_ci_high": min(1.0, a + z * se),
            "note": "Hanley-McNeil approximation for planning; report the data-based CI, not this one"}


def _round(obj):
    if isinstance(obj, float):
        return None if math.isnan(obj) else round(obj, 6)
    if isinstance(obj, dict):
        return {k: _round(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_round(v) for v in obj]
    return obj


def build_parser():
    p = argparse.ArgumentParser(
        description="ROC AUC of a response biomarker: Mann-Whitney AUC, stratified percentile bootstrap "
                    "CI, permutation p-value, DeLong SE; optional per-cohort, within-cohort and paired "
                    "two-marker analyses. Research use only.",
        epilog="Label 1 = responder. Rows whose label is in neither --positive nor --negative are excluded "
               "and counted in the output. Exit code 2 on unusable input.")
    p.add_argument("table", nargs="?", help="TSV (or .csv) file, one row per patient")
    p.add_argument("--score", help="column holding the marker")
    p.add_argument("--label", help="column holding the response label (1/0 or categories)")
    p.add_argument("--positive", help="comma-separated label values that count as responders, e.g. CR,PR")
    p.add_argument("--negative", help="comma-separated label values that count as non-responders, e.g. SD,PD")
    p.add_argument("--direction", choices=("higher", "lower"), default="higher",
                   help="prespecified side of the marker that indicates response (default: higher)")
    p.add_argument("--cohort", help="column naming the cohort; adds per-cohort and within-cohort AUCs")
    p.add_argument("--score2", help="second marker on the same patients; adds a paired comparison")
    p.add_argument("--direction2", choices=("higher", "lower"), default="higher",
                   help="prespecified direction of --score2 (default: higher)")
    p.add_argument("--n-boot", type=int, default=2000, help="bootstrap resamples (default 2000)")
    p.add_argument("--n-perm", type=int, default=5000, help="label permutations (default 5000)")
    p.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    p.add_argument("--alpha", type=float, default=0.05, help="1 - confidence level (default 0.05)")
    p.add_argument("--out", help="write the JSON here instead of standard output")
    g = p.add_argument_group("planning mode (no data needed)")
    g.add_argument("--plan", action="store_true", help="Hanley-McNeil SE for --auc, --n-pos, --n-neg")
    g.add_argument("--auc", type=float, help="anticipated AUC")
    g.add_argument("--n-pos", type=int, help="anticipated number of responders")
    g.add_argument("--n-neg", type=int, help="anticipated number of non-responders")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")  # legacy Windows code pages cannot print every character
    try:
        if not 0 < args.alpha < 1:
            raise InputError("--alpha must lie strictly between 0 and 1")
        if args.plan:
            result = plan(args.auc, args.n_pos, args.n_neg, args.alpha)
        else:
            if not (args.table and args.score and args.label):
                raise InputError("give a table, --score and --label (or use --plan)")
            if args.n_boot < 1 or args.n_perm < 1:
                raise InputError("--n-boot and --n-perm must be positive")
            header, rows = read_table(args.table)
            records, counts = build_records(header, rows, args.score, args.label, args.positive,
                                            args.negative, args.cohort, args.score2)
            result = {"input": {"score": args.score, "label": args.label, "direction": args.direction,
                                "positive": args.positive, "negative": args.negative,
                                "cohort": args.cohort, "score2": args.score2,
                                "direction2": args.direction2 if args.score2 else None,
                                "n_boot": args.n_boot, "n_perm": args.n_perm, "seed": args.seed,
                                "alpha": args.alpha},
                      "rows": dict(read=len(rows), used=len(records), **counts)}
            result.update(analyse(records, args.direction, args.direction2, bool(args.cohort),
                                  bool(args.score2), args.n_boot, args.n_perm, args.seed, args.alpha))
            real = [k for k in counts["excluded_by_label"] if k != "(missing)"]
            if real:
                result["warnings"].append(
                    "excluded by label: " + ", ".join(f"{k} = {counts['excluded_by_label'][k]}" for k in real)
                    + "; rerun with these categories added to --negative as a sensitivity analysis")
            miss = counts["missing_score_by_class"]
            if counts["excluded_missing_score"]:
                result["warnings"].append(
                    f"labelled patients without a score: responders {miss['responder']}, "
                    f"non-responders {miss['non_responder']}; check that missingness is unrelated to response")
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
