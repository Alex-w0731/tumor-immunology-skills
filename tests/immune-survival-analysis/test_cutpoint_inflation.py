"""Tests for skills/immune-survival-analysis/scripts/cutpoint_inflation.py."""
import importlib.util
import itertools
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "..", "skills", "immune-survival-analysis", "scripts")
SCRIPT = os.path.join(SCRIPTS, "cutpoint_inflation.py")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ci = _load("cutpoint_inflation")
km = _load("km_logrank")


class StatisticTests(unittest.TestCase):
    def test_hand_example_four_patients(self):
        # Same data as the hand-derived log-rank test in test_km_logrank: the high group fails at
        # 1 and 3, the low group at 2 and 4 -> chi2 = 8/13 for the cut that puts 2 patients high.
        times, events, marker = [1, 2, 3, 4], [1, 1, 1, 1], [0.9, 0.2, 0.8, 0.1]
        cuts = {k["n_high"]: k for k in ci.scan_cuts(marker, ci.prepare(times, events))}
        self.assertAlmostEqual(cuts[2]["chi2"], 8 / 13, places=12)
        self.assertAlmostEqual(cuts[2]["o_minus_e_high"], 2 - 4 / 3, places=12)
        self.assertEqual(cuts[2]["cut"], 0.2)

    def test_every_cut_equals_direct_logrank_with_ties(self):
        rng = random.Random(5)
        compared = 0
        for _ in range(12):
            n = rng.randint(8, 40)
            times = [round(rng.expovariate(0.3)) for _ in range(n)]       # ties and zeros
            events = [1 if rng.random() < 0.7 else 0 for _ in range(n)]
            marker = [round(rng.gauss(0, 1), 1) for _ in range(n)]         # tied markers
            for cut in ci.scan_cuts(marker, ci.prepare(times, events)):
                groups = ["high" if v > cut["cut"] else "low" for v in marker]
                ref = km.logrank_test(times, events, groups, order=["low", "high"])
                if ref["chi2"] is None:
                    continue
                compared += 1
                self.assertAlmostEqual(cut["chi2"], ref["chi2"], delta=1e-9 * max(1.0, ref["chi2"]))
        self.assertGreater(compared, 100)

    def test_score_test_continuous_hand_example(self):
        # x = 2, 1, 0 failing at 1, 2, 3.
        # t=1: risk set {2,1,0}, mean 1, U += 2 - 1 = 1, var = 5/3 - 1 = 2/3, w = 1
        # t=2: risk set {1,0}, mean 1/2, U += 1/2, var = 1/2 - 1/4 = 1/4
        # t=3: one patient left: U += 0, no variance.  chi2 = (3/2)^2 / (11/12) = 27/11
        res = ci.score_test([1, 2, 3], [1, 1, 1], [2.0, 1.0, 0.0])
        self.assertAlmostEqual(res["chi2"], 27 / 11, places=12)
        self.assertGreater(res["u"], 0)

    def test_score_test_binary_equals_logrank(self):
        times = [5, 8, 8, 12, 3, 9, 15, 2, 7, 7]
        events = [1, 0, 1, 1, 1, 0, 1, 1, 1, 0]
        x = [1, 0, 1, 0, 1, 0, 0, 1, 1, 0]
        ref = km.logrank_test(times, events, x)
        self.assertAlmostEqual(ci.score_test(times, events, [float(v) for v in x])["chi2"], ref["chi2"], places=10)

    def test_no_events_gives_undefined_tests(self):
        prep = ci.prepare([1, 2, 3, 4], [0, 0, 0, 0])
        self.assertTrue(all(k["chi2"] is None for k in ci.scan_cuts([1, 2, 3, 4], prep)))
        self.assertEqual(ci.score_test([1, 2, 3, 4], [0, 0, 0, 0], [1, 2, 3, 4])["p"], 1.0)

    def test_min_group_size(self):
        self.assertEqual(ci.min_group_size(150, 0.1), 15)
        self.assertEqual(ci.min_group_size(101, 0.1), 11)


class SimulationTests(unittest.TestCase):
    def test_deterministic_for_a_seed(self):
        one = ci.simulate(40, 30, 0.3, 0.1, 0.05, seed=3)
        two = ci.simulate(40, 30, 0.3, 0.1, 0.05, seed=3)
        three = ci.simulate(40, 30, 0.3, 0.1, 0.05, seed=4)
        self.assertEqual(one, two)
        self.assertNotEqual(one["best_cut_p_quantiles"], three["best_cut_p_quantiles"])

    def test_best_cut_inflates_while_prespecified_tests_hold_alpha(self):
        report = ci.simulate(100, 400, 0.3, 0.1, 0.05, seed=7)
        rates = {k: v["rate"] for k, v in report["false_positive_rate"].items()}
        self.assertGreater(rates["best_cut_point"], 0.2)
        self.assertGreater(rates["best_cut_point"], rates["median_split"] + 0.1)
        self.assertLess(rates["median_split"], 0.1)
        self.assertLess(rates["continuous_score_test"], 0.1)

    def test_permutation_mode(self):
        rng = random.Random(2)
        n = 80
        marker = [rng.gauss(0, 1) for _ in range(n)]
        strong = [rng.expovariate(math.exp(1.2 * m)) for m in marker]  # hazard rises steeply with the marker
        res = ci.permutation_test(strong, [1] * n, marker, 199, 0.1, seed=1)
        self.assertAlmostEqual(res["search_adjusted_p"], 1 / 200)
        self.assertEqual(res["continuous_score_test"]["direction"], "higher score, higher hazard")
        noise = [rng.expovariate(1.0) for _ in range(n)]
        res = ci.permutation_test(noise, [1] * n, marker, 199, 0.1, seed=1)
        self.assertGreaterEqual(res["search_adjusted_p"], res["best_cut"]["p_naive"])
        with self.assertRaises(ci.InputError):
            ci.permutation_test([1, 2, 3], [0, 0, 0], [1, 2, 3], 10, 0.1, seed=1)


class CommandLineTests(unittest.TestCase):
    def test_simulation_smoke(self):
        proc = subprocess.run([sys.executable, SCRIPT, "--n", "40", "--n-sim", "20", "--seed", "1"],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["settings"]["n_sim"], 20)

    def test_permutation_smoke_and_errors(self):
        rng = random.Random(9)
        rows = ["os_months\tos_event\tcd8"] + [f"{rng.expovariate(0.05):.2f}\t{int(rng.random() < 0.6)}\t{rng.gauss(0, 1):.3f}"
                                               for _ in range(60)] + ["12\t1\tNA"]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cohort.tsv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(rows) + "\n")
            base = [sys.executable, SCRIPT, "--input", path, "--time", "os_months", "--event", "os_event"]
            proc = subprocess.run(base + ["--score", "cd8", "--n-perm", "50"], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            report = json.loads(proc.stdout)
            self.assertEqual(report["input"]["rows_dropped_missing"], 1)
            self.assertIn("search_adjusted_p", report)
            proc = subprocess.run(base + ["--score", "cd4"], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("column 'cd4' not found", proc.stderr)


class ReviewAdditionTests(unittest.TestCase):
    """Added in the independent review: tied times in the score test, exact permutation p, errors."""

    def test_score_test_with_tied_times_by_hand(self):
        # times 1,1,2,3,3,4; events 1,1,0,1,1,0; x = .5, 2, 1, -1, 0, 3
        # t=1: all six at risk, mean 11/12, variance (divisor n) 15.25/6 - (11/12)^2 = 245/144;
        #      events x = .5 + 2 -> U += 2.5 - 2 * 11/12 = 2/3; w = 2*4/5 = 8/5 -> V += 8/5 * 245/144 = 49/18
        # t=3: at risk x = -1, 0, 3: mean 2/3, variance 10/3 - 4/9 = 26/9; events x = -1 + 0
        #      -> U += -1 - 4/3 = -7/3; w = 2*1/2 = 1 -> V += 26/9 = 52/18
        # U = -5/3, V = 101/18, chi2 = (25/9) / (101/18) = 50/101; U < 0: higher x, lower hazard
        res = ci.score_test([1, 1, 2, 3, 3, 4], [1, 1, 0, 1, 1, 0], [0.5, 2.0, 1.0, -1.0, 0.0, 3.0])
        self.assertAlmostEqual(res["u"], -5 / 3, places=12)
        self.assertAlmostEqual(res["var"], 101 / 18, places=12)
        self.assertAlmostEqual(res["chi2"], 50 / 101, places=12)

    def test_permutation_p_agrees_with_full_enumeration(self):
        # Exact permutation p from all 720 orderings of six marker values (minimum group size 1);
        # the seeded Monte Carlo p, (1 + exceedances) / (n_perm + 1), must agree within 4 standard errors.
        times, events = [3, 5, 6, 8, 11, 14], [1, 1, 0, 1, 1, 1]
        marker = [0.9, 0.7, 0.75, 0.3, 0.2, 0.1]
        prep = ci.prepare(times, events)
        observed = ci.summarise_cuts(marker, prep, 1)["best"]["chi2"]
        hits = total = 0
        for perm in itertools.permutations(marker):
            best = ci.summarise_cuts(list(perm), prep, 1)["best"]
            hits += best is not None and best["chi2"] >= observed * (1 - 1e-9)
            total += 1
        exact = hits / total
        res = ci.permutation_test(times, events, marker, 3000, 0.1, seed=4)
        self.assertLess(abs(res["search_adjusted_p"] - exact), 4 * math.sqrt(exact * (1 - exact) / 3000) + 1 / 3001)

    def test_bad_score_names_its_column_and_unwritable_out_exits_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "c.tsv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("t\te\tcd8\n5\t1\thigh\n")
            proc = subprocess.run([sys.executable, SCRIPT, "--input", path, "--time", "t", "--event", "e", "--score", "cd8"],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("score 'high' in column 'cd8'", proc.stderr)
            proc = subprocess.run([sys.executable, SCRIPT, "--n", "20", "--n-sim", "2", "--out", os.path.join(tmp, "missing", "x.json")],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("cannot write", proc.stderr)


if __name__ == "__main__":
    unittest.main()
