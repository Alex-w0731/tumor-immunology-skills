"""Unit tests for skills/response-biomarker-evaluation/scripts/threshold_metrics.py."""
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "response-biomarker-evaluation" / "scripts" / "threshold_metrics.py"
_spec = importlib.util.spec_from_file_location("threshold_metrics", SCRIPT)
tm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tm)


def binom_tail_ge(x, n, p):
    return sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(x, n + 1))


class TestClopperPearson(unittest.TestCase):
    def test_zero_successes_closed_form(self):
        # x = 0: the upper limit solves (1 - p)^n = alpha/2, so p = 1 - 0.025^(1/10) = 0.308497
        lo, hi = tm.clopper_pearson(0, 10)
        self.assertEqual(lo, 0.0)
        self.assertAlmostEqual(hi, 1 - 0.025 ** 0.1, places=10)

    def test_all_successes_closed_form(self):
        # x = n: the lower limit solves p^n = alpha/2, so p = 0.025^(1/10) = 0.691503
        lo, hi = tm.clopper_pearson(10, 10)
        self.assertAlmostEqual(lo, 0.025 ** 0.1, places=10)
        self.assertEqual(hi, 1.0)

    def test_defining_tail_equations(self):
        # The limits for x = 5 of n = 10 satisfy P(X >= 5 | p_L) = 0.025 and P(X <= 5 | p_U) = 0.025,
        # checked here directly with binomial probabilities. By symmetry p_L = 1 - p_U.
        # Expected values 0.187086 and 0.812914 (R binom.test(5, 10)); confirmed during review by an
        # independent bisection on the binomial tails that does not use this script's beta functions.
        lo, hi = tm.clopper_pearson(5, 10)
        self.assertAlmostEqual(binom_tail_ge(5, 10, lo), 0.025, places=9)
        self.assertAlmostEqual(1 - binom_tail_ge(6, 10, hi), 0.025, places=9)
        self.assertAlmostEqual(lo, 1 - hi, places=10)
        self.assertAlmostEqual(lo, 0.187086, places=6)
        self.assertAlmostEqual(hi, 0.812914, places=6)

    def test_other_alpha_and_large_n(self):
        lo, hi = tm.clopper_pearson(37, 120, alpha=0.1)
        self.assertAlmostEqual(binom_tail_ge(37, 120, lo), 0.05, places=8)
        self.assertAlmostEqual(1 - binom_tail_ge(38, 120, hi), 0.05, places=8)

    def test_betainc_identities(self):
        self.assertAlmostEqual(tm.betainc(1, 1, 0.3), 0.3, places=12)            # uniform
        self.assertAlmostEqual(tm.betainc(3, 1, 0.7), 0.7 ** 3, places=12)       # I_x(a, 1) = x^a
        self.assertAlmostEqual(tm.betainc(2.5, 4, 0.2) + tm.betainc(4, 2.5, 0.8), 1.0, places=12)

    def test_degenerate_inputs(self):
        self.assertEqual(tm.clopper_pearson(0, 0), (None, None))
        with self.assertRaises(ValueError):
            tm.clopper_pearson(4, 3)


class TestMetrics(unittest.TestCase):
    RECS = [{"y": 1, "s": s} for s in (5, 4, 3, 1)] + [{"y": 0, "s": s} for s in (4, 2, 1, 0, 0)]

    def test_hand_table(self):
        # threshold 3, score >= 3 positive. Responders 5,4,3 positive (TP 3), 1 negative (FN 1);
        # non-responders 4 positive (FP 1), 2,1,0,0 negative (TN 4).
        self.assertEqual(tm.confusion(self.RECS, 3), (3, 1, 1, 4))
        m = tm.metrics(3, 1, 1, 4)
        self.assertAlmostEqual(m["sensitivity"]["estimate"], 0.75)
        self.assertAlmostEqual(m["specificity"]["estimate"], 0.8)
        self.assertAlmostEqual(m["ppv"]["estimate"], 0.75)
        self.assertAlmostEqual(m["npv"]["estimate"], 0.8)
        self.assertAlmostEqual(m["response_rate"]["estimate"], 4 / 9)

    def test_boundary_is_inclusive_and_direction_lower(self):
        # score <= 3 positive: responders 3,1 (TP 2) and 5,4 (FN 2); non-responders 2,1,0,0 (FP 4), 4 (TN 1)
        self.assertEqual(tm.confusion(self.RECS, 3, "lower"), (2, 4, 2, 1))

    def test_no_test_positives(self):
        m = tm.metrics(0, 0, 3, 4)
        self.assertIsNone(m["ppv"]["estimate"])
        self.assertEqual(m["sensitivity"]["ci_low"], 0.0)

    def test_predictive_values_bayes(self):
        # sens 0.8, spec 0.9, prevalence 0.2: PPV = 0.16 / (0.16 + 0.08) = 2/3,
        # NPV = 0.72 / (0.72 + 0.04) = 0.947368
        ppv, npv = tm.predictive_values(0.8, 0.9, 0.2)
        self.assertAlmostEqual(ppv, 2 / 3, places=12)
        self.assertAlmostEqual(npv, 0.72 / 0.76, places=12)


class TestCommandLine(unittest.TestCase):
    DATA = "id\tbor\ttmb\tstudy\n" + "\n".join([
        "a\tCR\t14\tX", "b\tPR\t11\tX", "c\tPR\t6\tX", "d\tSD\t12\tX", "e\tPD\t3\tX", "f\tPD\t2\tX",
        "g\tNE\t20\tY", "h\tPR\t10\tY", "i\tSD\t4\tY", "j\tPD\t9\tY", "k\tPD\tNA\tY", "l\tCR\t8\tY"]) + "\n"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "tmb.tsv")
        with open(self.path, "w", encoding="utf-8") as fh:
            fh.write(self.DATA)

    def tearDown(self):
        self.tmp.cleanup()

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)

    def test_smoke(self):
        proc = self.run_cli(self.path, "--score", "tmb", "--label", "bor", "--threshold", "10",
                            "--positive", "CR,PR", "--negative", "SD,PD", "--cohort", "study", "--prevalence", "0.2")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out = json.loads(proc.stdout)
        # used: a b c d e f h i j l (g excluded as NE, k has no score)
        # responders a14 b11 c6 h10 l8 -> TP 3 (14, 11, 10), FN 2 ; non-responders d12 e3 f2 i4 j9 -> FP 1, TN 4
        o = out["overall"]
        self.assertEqual((o["tp"], o["fp"], o["fn"], o["tn"]), (3, 1, 2, 4))
        self.assertEqual(out["rows"]["excluded_by_label"], {"NE": 1})
        self.assertEqual(out["rows"]["excluded_missing_score"], 1)
        self.assertEqual([c["cohort"] for c in out["by_cohort"]], ["X", "Y"])
        self.assertIn("at_prevalence", out)
        self.assertTrue(any("too few patients" in w for w in out["warnings"]))

    def test_eval_case_values(self):
        # 18/20 and 17/20 (the threshold-and-leakage-audit eval): limits from the binomial tail equations.
        for x, n in ((18, 20), (17, 20), (18, 21), (17, 19)):
            lo, hi = tm.clopper_pearson(x, n)
            self.assertAlmostEqual(binom_tail_ge(x, n, lo), 0.025, places=9)
            self.assertAlmostEqual(1 - binom_tail_ge(x + 1, n, hi), 0.025, places=9)
        lo, hi = tm.clopper_pearson(18, 20)
        self.assertAlmostEqual(lo, 0.683017, places=6)
        self.assertAlmostEqual(hi, 0.987651, places=6)

    def test_non_finite_threshold_exit_2(self):
        proc = self.run_cli(self.path, "--score", "tmb", "--label", "bor", "--threshold", "nan",
                            "--positive", "CR,PR", "--negative", "SD,PD")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--threshold", proc.stderr)

    def test_missing_column_exit_2(self):
        proc = self.run_cli(self.path, "--score", "TMB_typo", "--label", "bor", "--threshold", "10",
                            "--positive", "CR,PR", "--negative", "SD,PD")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("TMB_typo", proc.stderr)


if __name__ == "__main__":
    unittest.main()
