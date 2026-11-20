"""Tests for skills/immune-survival-analysis/scripts/fdr_adjust.py."""
import importlib.util
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "..", "skills", "immune-survival-analysis", "scripts", "fdr_adjust.py")
spec = importlib.util.spec_from_file_location("fdr_adjust", SCRIPT)
fdr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fdr)


class AdjustTests(unittest.TestCase):
    def assertListAlmost(self, got, want):
        self.assertEqual(len(got), len(want))
        for g, w in zip(got, want):
            if w is None:
                self.assertIsNone(g)
            else:
                self.assertAlmostEqual(g, w, places=12)

    def test_evenly_spaced_hand_example(self):
        # p = .01 .02 .03 .04 .05, m = 5.
        # BH: 5/1*.01 = 5/2*.02 = 5/3*.03 = 5/4*.04 = 5/5*.05 = .05 -> all .05
        # Holm: 5*.01 = .05; 4*.02 = .08; 3*.03 = .09; 2*.04 = .08 -> running max .09; 1*.05 -> .09
        p = [0.01, 0.02, 0.03, 0.04, 0.05]
        self.assertListAlmost(fdr.adjust_bh(p), [0.05] * 5)
        self.assertListAlmost(fdr.adjust_holm(p), [0.05, 0.08, 0.09, 0.09, 0.09])

    def test_ties_and_missing_hand_example(self):
        # non-missing, sorted: .01 (k=1), .04 (k=2), .04 (k=3), .20 (k=4); m = 4
        # BH:   4*.01 = .04; 2*.04 = .08; 4/3*.04 = .053333; 1*.20 = .20
        #       running min from the top: .20, .053333, .053333, .04
        # Holm: 4*.01 = .04; 3*.04 = .12; 2*.04 = .08 -> .12; 1*.20 = .20
        p = [0.04, None, 0.01, 0.04, 0.20]
        self.assertListAlmost(fdr.adjust_bh(p), [0.04 * 4 / 3, None, 0.04, 0.04 * 4 / 3, 0.20])
        self.assertListAlmost(fdr.adjust_holm(p), [0.12, None, 0.04, 0.12, 0.20])

    def test_more_tests_performed_than_rows(self):
        # same p-values, N = 10 tests performed of which 4 are in the file
        # BH:   10*.01 = .1; 5*.04 = .2; 10/3*.04 = .133333; 2.5*.20 = .5 -> .5, .133333, .133333, .1
        # Holm: 10*.01 = .1; 9*.04 = .36; 8*.04 = .32 -> .36; 7*.20 = 1.4 -> capped at 1
        p = [0.04, None, 0.01, 0.04, 0.20]
        self.assertListAlmost(fdr.adjust_bh(p, n_tests=10), [0.4 / 3, None, 0.1, 0.4 / 3, 0.5])
        self.assertListAlmost(fdr.adjust_holm(p, n_tests=10), [0.36, None, 0.1, 0.36, 1.0])

    def test_edge_cases(self):
        self.assertListAlmost(fdr.adjust_bh([0.3]), [0.3])
        self.assertListAlmost(fdr.adjust_holm([0.3]), [0.3])
        self.assertListAlmost(fdr.adjust_bh([None, None]), [None, None])
        self.assertListAlmost(fdr.adjust_bh([]), [])
        self.assertListAlmost(fdr.adjust_bh([0.0, 1.0]), [0.0, 1.0])
        with self.assertRaises(fdr.InputError):
            fdr.adjust_bh([0.1, 0.2, 0.3], n_tests=2)

    def test_bh_never_exceeds_holm_and_both_keep_order(self):
        rng = random.Random(11)
        p = [rng.random() ** 3 for _ in range(200)]
        bh, holm = fdr.adjust_bh(p), fdr.adjust_holm(p)
        order = sorted(range(len(p)), key=lambda i: p[i])
        for a, b in zip(order, order[1:]):
            self.assertLessEqual(bh[a], bh[b] + 1e-15)
            self.assertLessEqual(holm[a], holm[b] + 1e-15)
        for q, h, raw in zip(bh, holm, p):
            self.assertLessEqual(raw, q + 1e-15)
            self.assertLessEqual(q, h + 1e-15)

    def test_parse_rejects_bad_values(self):
        with self.assertRaises(fdr.InputError):
            fdr.parse_pvalues([{"p": "0.2"}, {"p": "1.5"}], "p")
        with self.assertRaises(fdr.InputError):
            fdr.parse_pvalues([{"p": "abc"}], "p")
        self.assertEqual(fdr.parse_pvalues([{"p": "NA"}, {"p": ""}, {"p": "1e-3"}], "p"), [None, None, 0.001])


class CommandLineTests(unittest.TestCase):
    def run_cli(self, text, *args):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "tests.tsv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            return subprocess.run([sys.executable, SCRIPT, "--input", path, *args], capture_output=True, text=True)

    def test_smoke(self):
        proc = self.run_cli("feature\tp\nCD8\t0.04\nNK\tNA\nTreg\t0.01\nB\t0.04\nM2\t0.20\n", "--column", "p", "--sort")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = [line.split("\t") for line in proc.stdout.strip().split("\n")]
        self.assertEqual(lines[0], ["feature", "p", "p_bh", "p_holm"])
        self.assertEqual(lines[1], ["Treg", "0.01", "0.04", "0.04"])
        self.assertEqual(lines[-1], ["NK", "NA", "NA", "NA"])
        self.assertIn("4 p-values", proc.stderr)

    def test_bad_input_exits_2(self):
        proc = self.run_cli("feature\tp\nCD8\t1.2\n", "--column", "p")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("outside [0, 1]", proc.stderr)
        proc = self.run_cli("feature\tp\nCD8\t0.2\n", "--column", "pval")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("not found", proc.stderr)


class ReviewAdditionTests(unittest.TestCase):
    """Added in the independent review."""

    def test_fifteen_pvalues_by_hand(self):
        # BH step-up at q = 0.05: the largest k with p(k) <= k/15 * 0.05 is k = 4 (0.0095 <= 0.01333;
        # 0.0201 > 0.01667 and no later p(k) is below its bound) -> 4 rejections.
        # Holm step-down: 0.0001 <= .05/15, 0.0004 <= .05/14, 0.0019 <= .05/13, 0.0095 > .05/12 -> 3.
        # For p = 0.0095 (k = 4): BH = min over k >= 4 of 15 p(k) / k = 15 * 0.0095 / 4 = 0.035625;
        # Holm = max(15 * 0.0001, 14 * 0.0004, 13 * 0.0019, 12 * 0.0095) = 0.114.
        p = [0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459, 0.3240, 0.4262, 0.5719,
             0.6528, 0.7590, 1.0]
        bh, holm = fdr.adjust_bh(p), fdr.adjust_holm(p)
        self.assertEqual(sum(q <= 0.05 for q in bh), 4)
        self.assertEqual(sum(h <= 0.05 for h in holm), 3)
        self.assertAlmostEqual(bh[3], 0.035625, places=12)
        self.assertAlmostEqual(holm[3], 0.114, places=12)

    def test_unwritable_out_exits_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.tsv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("feature\tp\nCD8\t0.2\n")
            proc = subprocess.run([sys.executable, SCRIPT, "--input", path, "--column", "p", "--out",
                                   os.path.join(tmp, "missing", "o.tsv")], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("cannot write", proc.stderr)


if __name__ == "__main__":
    unittest.main()
