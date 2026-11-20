"""Unit tests for skills/response-biomarker-evaluation/scripts/roc_auc.py (standard library only)."""
import importlib.util
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "response-biomarker-evaluation" / "scripts" / "roc_auc.py"
_spec = importlib.util.spec_from_file_location("roc_auc", SCRIPT)
ra = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ra)


def brute_u(pos, neg):
    return sum(1.0 if x > y else 0.5 if x == y else 0.0 for x in pos for y in neg)


class TestAUC(unittest.TestCase):
    def test_hand_example_with_ties(self):
        # pos = 3, 2, 2 ; neg = 1, 2.  Pairs: (3,1)=1 (3,2)=1 (2,1)=1 (2,2)=1/2 (2,1)=1 (2,2)=1/2
        # U = 5 of n1*n0 = 6 pairs, AUC = 5/6.
        self.assertEqual(ra.mann_whitney_u([3, 2, 2], [1, 2]), 5.0)
        self.assertAlmostEqual(ra.auc([3, 2, 2], [1, 2]), 5 / 6, places=12)

    def test_ordinal_rating_table(self):
        # Five-category ordinal score with heavy ties.
        # responders per category 1..5:      3,  2, 2, 11, 33   (n1 = 51)
        # non-responders per category 1..5: 33,  6, 6, 11,  2   (n0 = 58)
        # U = sum_r count_r * (#non-responders below r + 0.5 * #non-responders at r)
        #   = 3*(0 + 16.5) + 2*(33 + 3) + 2*(39 + 3) + 11*(45 + 5.5) + 33*(56 + 1)
        #   = 49.5 + 72 + 84 + 555.5 + 1881 = 2642 ;  AUC = 2642 / 2958 = 0.893171
        pos = [1] * 3 + [2] * 2 + [3] * 2 + [4] * 11 + [5] * 33
        neg = [1] * 33 + [2] * 6 + [3] * 6 + [4] * 11 + [5] * 2
        self.assertEqual(ra.mann_whitney_u(pos, neg), 2642.0)
        self.assertAlmostEqual(ra.auc(pos, neg), 2642 / 2958, places=12)

    def test_matches_brute_force_and_midranks(self):
        rng = random.Random(3)
        for _ in range(25):
            pos = [rng.randint(0, 6) for _ in range(rng.randint(1, 12))]
            neg = [rng.randint(0, 6) for _ in range(rng.randint(1, 12))]
            u = ra.mann_whitney_u(pos, neg)
            self.assertEqual(u, brute_u(pos, neg))
            ranks = ra.midranks(pos + neg)
            self.assertAlmostEqual(sum(ranks[:len(pos)]) - len(pos) * (len(pos) + 1) / 2, u, places=9)

    def test_extremes_and_all_tied(self):
        self.assertEqual(ra.auc([5, 6], [1, 2]), 1.0)
        self.assertEqual(ra.auc([1, 2], [5, 6]), 0.0)
        self.assertEqual(ra.auc([4, 4], [4, 4, 4]), 0.5)

    def test_empty_group_raises(self):
        with self.assertRaises(ValueError):
            ra.auc([], [1.0])

    def test_midranks(self):
        self.assertEqual(ra.midranks([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])

    def test_within_cohort_vs_pooled(self):
        # Within each cohort the marker is uninformative, but cohort A has higher scores and more
        # responders, so pooling manufactures discrimination.
        a = ([10, 12, 14], [11, 13])   # U = 0 + 1 + 2 = 3 of 6 pairs
        b = ([1, 3], [0, 2, 4])        # U = 1 + 2 = 3 of 6 pairs
        self.assertEqual(ra.within_auc([a, b]), 0.5)
        # pooled: 10 beats 0,2,4 (3); 12 beats 11,0,2,4 (4); 14 beats all (5); 1 beats 0 (1); 3 beats 0,2 (2)
        self.assertEqual(ra.auc(a[0] + b[0], a[1] + b[1]), 15 / 25)


class TestDeLong(unittest.TestCase):
    def test_hand_se(self):
        # pos = 3, 2, 2 ; neg = 1, 2
        # V10 = (1, 0.75, 0.75): S10 = (1/36 + 1/144 + 1/144) / 2 = 1/48
        # V01 = (1, 2/3):        S01 = 2 * (1/6)^2 / 1 = 1/18
        # Var = S10/3 + S01/2 = 1/144 + 4/144 = 5/144 ; SE = sqrt(5)/12 = 0.186339
        self.assertAlmostEqual(ra.delong_se([3, 2, 2], [1, 2]), math.sqrt(5) / 12, places=12)

    def test_paired_negated_marker(self):
        # Marker 2 = -marker 1: AUC2 = 1 - AUC1, diff = 2*AUC1 - 1 = 2/3, SE(diff) = 2*SE(AUC1),
        # z = (2/3) / (sqrt(5)/6) = 1.788854, and the two AUC estimates have correlation -1.
        pos, neg = [3, 2, 2], [1, 2]
        res = ra.delong_paired([(x, -x) for x in pos], [(y, -y) for y in neg])
        self.assertAlmostEqual(res["diff"], 2 / 3, places=12)
        self.assertAlmostEqual(res["se_diff"], math.sqrt(5) / 6, places=12)
        self.assertAlmostEqual(res["z"], 4 / math.sqrt(5), places=9)
        self.assertAlmostEqual(res["p_delong"], math.erfc(4 / math.sqrt(5) / math.sqrt(2)), places=12)
        self.assertAlmostEqual(res["corr_auc"], -1.0, places=9)

    def test_identical_markers_give_no_test(self):
        pos, neg = [3, 2, 2, 5], [1, 2, 0]
        res = ra.delong_paired([(x, x) for x in pos], [(y, y) for y in neg])
        self.assertEqual(res["diff"], 0.0)
        self.assertIsNone(res["z"])
        self.assertIsNone(res["p_delong"])

    def test_hand_se_second_example(self):
        # pos = 0.9, 0.5, 0.5, 0.3 ; neg = 0.5, 0.2, 0.1 (one tie). U = 3 + 2.5 + 2.5 + 2 = 10 of 12.
        # V10 = (1, 5/6, 5/6, 2/3): S10 = 2 * (1/6)^2 / 3 = 1/54
        # V01 = (1/2, 1, 1):        S01 = ((1/3)^2 + 2 * (1/6)^2) / 2 = 1/12
        # Var = S10/4 + S01/3 = 1/216 + 6/216 = 7/216 ; SE = 0.180021
        pos, neg = [0.9, 0.5, 0.5, 0.3], [0.5, 0.2, 0.1]
        self.assertAlmostEqual(ra.auc(pos, neg), 10 / 12, places=12)
        self.assertAlmostEqual(ra.delong_se(pos, neg), math.sqrt(7 / 216), places=12)

    def test_stratified_within_cohort_difference(self):
        # Two identical strata, each pos = 3, 2, 2 ; neg = 1, 2 with marker 2 = -marker 1.
        # Per stratum: diff = 2/3 and Var = (sqrt(5)/6)^2 = 5/36 (test_paired_negated_marker).
        # Equal weights (6 pairs each, w = 1/2): diff = 2/3, Var = 2 * (1/4) * 5/36 = 5/72,
        # z = (2/3) / sqrt(5/72) = 2.529822.
        stratum = ([(x, -x) for x in (3, 2, 2)], [(y, -y) for y in (1, 2)])
        res = ra.delong_paired_within([stratum, stratum])
        self.assertAlmostEqual(res["diff"], 2 / 3, places=12)
        self.assertAlmostEqual(res["se_diff"], math.sqrt(5 / 72), places=12)
        self.assertAlmostEqual(res["z"], (2 / 3) / math.sqrt(5 / 72), places=9)

    def test_stratified_weights_follow_pair_counts(self):
        # The within-cohort difference equals within_auc(marker 1) - within_auc(marker 2).
        rng = random.Random(11)
        strata = []
        for n1, n0 in ((4, 9), (6, 5), (3, 12)):
            strata.append(([(rng.gauss(0.5, 1), rng.gauss(0, 1)) for _ in range(n1)],
                           [(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n0)]))
        res = ra.delong_paired_within(strata)
        a1 = ra.within_auc([([a for a, _ in p], [a for a, _ in q]) for p, q in strata])
        a2 = ra.within_auc([([b for _, b in p], [b for _, b in q]) for p, q in strata])
        self.assertAlmostEqual(res["diff"], a1 - a2, places=12)


class TestHanleyMcNeil(unittest.TestCase):
    def test_null_variance(self):
        # At A = 0.5, Q1 = Q2 = 1/3 and Var = (1/4 + (n1-1)/12 + (n0-1)/12) / (n1 n0)
        # = (n1 + n0 + 1) / (12 n1 n0), the exact null variance of U / (n1 n0).
        # n1 = n0 = 10: 21/1200, SE = 0.132288.
        self.assertAlmostEqual(ra.hanley_mcneil_se(0.5, 10, 10), math.sqrt(21 / 1200), places=12)

    def test_worked_value(self):
        # A = 0.893, n1 = 51, n0 = 58: Q1 = 0.893/1.107 = 0.806685, Q2 = 1.594898/1.893 = 0.842524,
        # A^2 = 0.797449; Var = (0.095551 + 50*0.009236 + 57*0.045075) / 2958 = 3.1266/2958
        # SE = 0.03251
        self.assertAlmostEqual(ra.hanley_mcneil_se(0.893, 51, 58), 0.03251, delta=0.0001)

    def test_plan_validation(self):
        with self.assertRaises(ra.InputError):
            ra.plan(1.2, 10, 10)
        with self.assertRaises(ra.InputError):
            ra.plan(0.7, None, 10)


class TestResampling(unittest.TestCase):
    def test_permutation_perfect_separation(self):
        # n1 = n0 = 5, perfect separation: exact two-sided p = 2 / C(10, 5) = 2/252 = 0.0079.
        p = ra.permutation_p([([6, 7, 8, 9, 10], [1, 2, 3, 4, 5])], 20000, random.Random(5))
        self.assertAlmostEqual(p, 2 / 252, delta=0.0025)

    def test_permutation_matches_full_enumeration_with_ties(self):
        # pos = 0.9, 0.5, 0.5, 0.3 ; neg = 0.5, 0.2, 0.1: observed U = 10 of 12, |U - 6| = 4.
        # Enumerate all C(7, 4) = 35 relabellings; the exact two-sided p is the share with |U - 6| >= 4.
        import itertools
        vals = [0.9, 0.5, 0.5, 0.3, 0.5, 0.2, 0.1]
        hits = 0
        for idx in itertools.combinations(range(7), 4):
            pos = [vals[i] for i in idx]
            neg = [vals[i] for i in range(7) if i not in idx]
            hits += abs(brute_u(pos, neg) - 6) >= 4
        exact = hits / 35                                   # 7/35 = 0.2
        p = ra.permutation_p([(vals[:4], vals[4:])], 20000, random.Random(8))
        self.assertAlmostEqual(exact, 0.2, places=12)
        self.assertAlmostEqual(p, exact, delta=4 * math.sqrt(0.2 * 0.8 / 20000))

    def test_permutation_all_tied_is_one(self):
        self.assertEqual(ra.permutation_p([([1, 1, 1], [1, 1, 1])], 200, random.Random(1)), 1.0)

    def test_permutation_never_zero(self):
        p = ra.permutation_p([([6, 7, 8, 9, 10], [1, 2, 3, 4, 5])], 9, random.Random(2))
        self.assertGreaterEqual(p, 0.1)    # (b + 1) / (m + 1) with m = 9

    def test_bootstrap_deterministic_and_brackets_estimate(self):
        pos = [2.1, 3.5, 1.2, 4.4, 2.9, 3.3]
        neg = [1.0, 2.2, 0.5, 1.9, 2.8, 1.1, 0.3]
        ci1 = ra.bootstrap_ci([(pos, neg)], 500, 0.05, random.Random("x"))
        ci2 = ra.bootstrap_ci([(pos, neg)], 500, 0.05, random.Random("x"))
        self.assertEqual(ci1, ci2)
        self.assertTrue(0 <= ci1[0] <= ra.auc(pos, neg) <= ci1[1] <= 1)

    def test_bootstrap_perfect_separation_is_degenerate(self):
        ci = ra.bootstrap_ci([([5, 6, 7], [1, 2, 3])], 200, 0.05, random.Random(1))
        self.assertEqual(ci, (1.0, 1.0))

    def test_quantile_matches_numpy_default(self):
        self.assertAlmostEqual(ra._quantile([1, 2, 3, 4], 0.5), 2.5)
        self.assertAlmostEqual(ra._quantile([1, 2, 3, 4], 0.25), 1.75)


class TestRecordsAndAnalysis(unittest.TestCase):
    HEADER = ["id", "bor", "score"]

    def rows(self, data):
        return [dict(zip(self.HEADER, r)) for r in data]

    def test_category_mapping_counts_exclusions(self):
        rows = self.rows([("a", "CR", "3"), ("b", "pr", "2"), ("c", "SD", "1"), ("d", "PD", "0"),
                          ("e", "NE", "5"), ("f", "", "4"), ("g", "PD", "NA")])
        recs, counts = ra.build_records(self.HEADER, rows, "score", "bor", "CR,PR", "SD,PD")
        self.assertEqual([r["y"] for r in recs], [1, 1, 0, 0])
        self.assertEqual(counts["excluded_by_label"], {"NE": 1, "(missing)": 1})
        self.assertEqual(counts["excluded_missing_score"], 1)
        self.assertEqual(counts["missing_score_by_class"], {"responder": 0, "non_responder": 1})

    def test_separation_and_small_cohort_flags(self):
        # Cohort A: 3 responders all above 4 non-responders (AUC 1); cohort B: 12 v 12, overlapping.
        recs = [{"y": 1, "s": s, "s2": None, "c": "A"} for s in (10.0, 11.0, 12.0)]
        recs += [{"y": 0, "s": s, "s2": None, "c": "A"} for s in (1.0, 2.0, 3.0, 4.0)]
        recs += [{"y": i % 2, "s": float(i), "s2": None, "c": "B"} for i in range(24)]
        res = ra.analyse(recs, use_cohort=True, n_boot=50, n_perm=50)
        joined = " | ".join(res["warnings"])
        self.assertIn("cohort A: only 3 patients", joined)
        self.assertIn("cohort A: complete separation", joined)
        self.assertNotIn("cohort B:", joined)

    def test_categorical_label_without_mapping_fails(self):
        with self.assertRaises(ra.InputError) as cm:
            ra.build_records(self.HEADER, self.rows([("a", "CR", "3")]), "score", "bor")
        self.assertIn("--positive", str(cm.exception))

    def test_bad_inputs_fail_loudly(self):
        with self.assertRaises(ra.InputError):
            ra.build_records(self.HEADER, self.rows([("a", "1", "high")]), "score", "bor")
        with self.assertRaises(ra.InputError):
            ra.build_records(self.HEADER, [], "tmb", "bor")
        with self.assertRaises(ra.InputError):
            ra.build_records(self.HEADER, [], "score", "bor", "CR,PR", None)

    def test_too_few_per_class(self):
        recs = [{"y": 1, "s": 1.0, "s2": None, "c": None}] * 2 + [{"y": 0, "s": 0.0, "s2": None, "c": None}] * 5
        with self.assertRaises(ra.InputError):
            ra.analyse(recs, n_boot=10, n_perm=10)

    def test_direction_lower_is_complement_without_ties(self):
        recs = [{"y": 1, "s": s, "s2": None, "c": None} for s in (5.0, 3.0, 4.5, 1.0)]
        recs += [{"y": 0, "s": s, "s2": None, "c": None} for s in (2.0, 0.5, 3.5, 1.5)]
        hi = ra.analyse(recs, "higher", n_boot=50, n_perm=50)["overall"]["auc"]
        lo = ra.analyse(recs, "lower", n_boot=50, n_perm=50)["overall"]["auc"]
        self.assertAlmostEqual(hi + lo, 1.0, places=12)

    def test_cohort_analysis_flags_pooling_gap(self):
        # Same data as test_within_cohort_vs_pooled: pooled AUC 0.6, within-cohort AUC 0.5.
        # Each cohort has fewer than 3 of one class, so per-cohort AUCs are not estimated.
        recs = []
        for c, pos, neg in (("A", [10, 12, 14], [11, 13]), ("B", [1, 3], [0, 2, 4])):
            recs += [{"y": 1, "s": float(s), "s2": None, "c": c} for s in pos]
            recs += [{"y": 0, "s": float(s), "s2": None, "c": c} for s in neg]
        res = ra.analyse(recs, use_cohort=True, n_boot=100, n_perm=100)
        self.assertEqual([r["cohort"] for r in res["by_cohort"]], ["A", "B"])
        self.assertTrue(all(r["auc"] is None for r in res["by_cohort"]))
        self.assertEqual(res["overall"]["auc"], 0.6)
        self.assertEqual(res["within_cohort"]["auc"], 0.5)
        self.assertTrue(any("within-cohort" in w for w in res["warnings"]))


class TestCommandLine(unittest.TestCase):
    DATA = """patient\tbor\tsig\tsig2\tstudy
p1\tCR\t2.5\t1.0\tA
p2\tPR\t1.9\t0.2\tA
p3\tPR\t1.2\t1.5\tA
p4\tSD\t0.8\t0.3\tA
p5\tPD\t0.4\t-0.5\tA
p6\tPD\t1.0\t0.9\tA
p7\tNE\t0.1\t0.0\tA
p8\tSD\t1.5\t1.1\tA
q1\tCR\t0.9\t0.4\tB
q2\tPR\t0.5\t0.8\tB
q3\tPR\t0.7\t-0.2\tB
q4\tPD\t-0.3\t0.1\tB
q5\tPD\t0.2\t0.5\tB
q6\tSD\t-0.1\t-0.4\tB
q7\tNE\t-0.5\t0.3\tB
q8\tPD\t0.6\t0.0\tB
"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "cohort.tsv")
        with open(self.path, "w", encoding="utf-8") as fh:
            fh.write(self.DATA)

    def tearDown(self):
        self.tmp.cleanup()

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)

    def test_full_run_is_deterministic(self):
        args = [self.path, "--score", "sig", "--label", "bor", "--positive", "CR,PR", "--negative", "SD,PD",
                "--cohort", "study", "--score2", "sig2", "--n-boot", "300", "--n-perm", "300", "--seed", "3"]
        first, second = self.run_cli(*args), self.run_cli(*args)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(first.stdout, second.stdout)
        out = json.loads(first.stdout)
        self.assertEqual(out["rows"]["excluded_by_label"], {"NE": 2})
        self.assertEqual((out["overall"]["n_pos"], out["overall"]["n_neg"]), (6, 8))
        self.assertEqual(len(out["by_cohort"]), 2)
        self.assertIn("within_cohort", out)
        self.assertIn("p_delong", out["comparison"]["overall"])
        self.assertTrue(any("NE = 2" in w for w in out["warnings"]))
        # Both cohorts have 3 responders and 4 non-responders, so both enter the stratified comparison,
        # whose difference is the difference of the two within-cohort AUCs.
        within = out["comparison"]["within_cohort"]
        self.assertEqual(within["cohorts_used"], 2)
        self.assertAlmostEqual(within["diff"], within["auc1"] - within["auc2"], places=5)
        self.assertAlmostEqual(within["auc1"], out["within_cohort"]["auc"], places=6)

    def test_unwritable_out_exit_2(self):
        proc = self.run_cli(self.path, "--score", "sig", "--label", "bor", "--positive", "CR,PR",
                            "--negative", "SD,PD", "--out", os.path.join(self.tmp.name, "missing", "x.json"))
        self.assertEqual(proc.returncode, 2)
        self.assertIn("cannot write", proc.stderr)

    def test_categorical_labels_without_mapping_exit_2(self):
        proc = self.run_cli(self.path, "--score", "sig", "--label", "bor")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--positive", proc.stderr)

    def test_plan_mode(self):
        proc = self.run_cli("--plan", "--auc", "0.5", "--n-pos", "10", "--n-neg", "10")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertAlmostEqual(json.loads(proc.stdout)["se_hanley_mcneil"], 0.132288, places=5)


if __name__ == "__main__":
    unittest.main()
