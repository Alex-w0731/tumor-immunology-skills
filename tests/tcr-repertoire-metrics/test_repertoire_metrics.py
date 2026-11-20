"""Unit tests for skills/tcr-repertoire-metrics/scripts/repertoire_metrics.py."""
import collections
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

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SCRIPT = os.path.join(ROOT, "skills", "tcr-repertoire-metrics", "scripts", "repertoire_metrics.py")

_spec = importlib.util.spec_from_file_location("repertoire_metrics", SCRIPT)
rm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rm)


def write(folder, name, text):
    path = os.path.join(folder, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return path


def run_cli(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True, timeout=120)


class GoldenValues(unittest.TestCase):
    def test_hand_derived_example(self):
        # counts 4, 2, 1, 1: N = 8, S = 4, p = 1/2, 1/4, 1/8, 1/8
        # H = (1/2) ln 2 + (1/4)(2 ln 2) + 2 (1/8)(3 ln 2) = 1.75 ln 2
        # J = 1.75 ln 2 / ln 4 = 0.875, so clonality = 0.125
        # lambda = 1/4 + 1/16 + 2/64 = 0.34375; 1/lambda = 32/11; exp(H) = 2 ** 1.75
        # Gini with x ascending 1, 1, 2, 4: (-3*1 - 1*1 + 1*2 + 3*4) / (4 * 8) = 10/32 = 0.3125
        # f1 = 2, f2 = 1: bias-corrected Chao1 = 4 + 2*1 / (2*2) = 4.5
        m = rm.diversity_metrics([4, 2, 1, 1], top=(1, 2))
        self.assertEqual(m["depth"], 8)
        self.assertEqual(m["richness"], 4)
        self.assertAlmostEqual(m["shannon"], 1.75 * math.log(2), places=12)
        self.assertAlmostEqual(m["pielou_evenness"], 0.875, places=12)
        self.assertAlmostEqual(m["clonality"], 0.125, places=12)
        self.assertAlmostEqual(m["simpson"], 0.34375, places=12)
        self.assertAlmostEqual(m["inverse_simpson"], 32 / 11, places=12)
        self.assertAlmostEqual(m["gini_simpson"], 0.65625, places=12)
        self.assertAlmostEqual(m["exp_shannon"], 2 ** 1.75, places=12)
        self.assertAlmostEqual(m["gini"], 0.3125, places=12)
        self.assertAlmostEqual(m["top1_fraction"], 0.5)
        self.assertAlmostEqual(m["top2_fraction"], 0.75)
        self.assertEqual((m["singletons"], m["doubletons"]), (2, 1))
        self.assertAlmostEqual(m["chao1"], 4.5)
        self.assertEqual(m["expanded_clonotypes"], 2)
        self.assertAlmostEqual(m["expanded_fraction"], 0.75)

    def test_expected_richness_hand_derived(self):
        # Hurlbert: E[S_2] = sum_i [1 - C(8 - c_i, 2) / C(8, 2)], C(8, 2) = 28
        # c = 4: 1 - 6/28; c = 2: 1 - 15/28; c = 1 (twice): 1 - 21/28  ->  49/28 = 1.75
        # check: 2 - P(both draws from the same clonotype) = 2 - (6 + 1) / 28 = 1.75
        self.assertAlmostEqual(rm.expected_richness([4, 2, 1, 1], 2), 1.75, places=12)
        self.assertAlmostEqual(rm.expected_richness([4, 2, 1, 1], 8), 4.0, places=12)
        self.assertAlmostEqual(rm.expected_richness([4, 2, 1, 1], 1), 1.0, places=12)

    def test_gini_matches_mean_absolute_difference(self):
        values = [7, 1, 3, 3, 12, 5]
        n, mean = len(values), sum(values) / len(values)
        mad = sum(abs(x - y) for x in values for y in values) / (n * n)
        self.assertAlmostEqual(rm.gini(values), mad / (2 * mean), places=12)


class EdgeCases(unittest.TestCase):
    def test_single_clonotype(self):
        m = rm.diversity_metrics([10])
        self.assertEqual(m["richness"], 1)
        self.assertIsNone(m["pielou_evenness"])
        self.assertIsNone(m["clonality"])
        self.assertEqual(m["simpson"], 1.0)
        self.assertEqual(m["gini"], 0.0)
        self.assertEqual(m["chao1"], 1.0)

    def test_zeros_are_ignored(self):
        self.assertEqual(rm.diversity_metrics([0, 4, 0, 2, 1, 1]), rm.diversity_metrics([4, 2, 1, 1]))

    def test_even_distribution_and_ties(self):
        m = rm.diversity_metrics([3, 3, 3], top=(1, 10))
        self.assertAlmostEqual(m["pielou_evenness"], 1.0, places=12)
        self.assertAlmostEqual(m["clonality"], 0.0, places=12)
        self.assertEqual(m["gini"], 0.0)
        for s in range(2, 40):   # evenness never leaves [0, 1] through rounding
            self.assertLessEqual(rm.diversity_metrics([7] * s)["pielou_evenness"], 1.0)
        self.assertAlmostEqual(m["top1_fraction"], 1 / 3)
        self.assertEqual(m["top10_fraction"], 1.0)

    def test_frequencies_give_no_count_based_metrics(self):
        m = rm.diversity_metrics([0.5, 0.25, 0.25])
        self.assertIsNone(m["chao1"])
        self.assertIsNone(m["singletons"])
        self.assertAlmostEqual(m["clonality"], 1 - 1.5 * math.log(2) / math.log(3), places=12)

    def test_empty_input_raises(self):
        with self.assertRaises(rm.InputError):
            rm.diversity_metrics([0, 0])
        with self.assertRaises(rm.InputError):
            rm.diversity_metrics([])


class Downsampling(unittest.TestCase):
    COUNTS = [40, 25, 10, 5, 5, 3, 2, 2] + [1] * 30

    def test_subsample_properties(self):
        rng = random.Random("1:S")
        drawn = rm.subsample(self.COUNTS, 50, rng)
        self.assertEqual(sum(drawn), 50)
        self.assertTrue(all(0 <= d <= c for d, c in zip(drawn, self.COUNTS)))
        self.assertEqual(rm.subsample(self.COUNTS, sum(self.COUNTS), rng), self.COUNTS)
        with self.assertRaises(rm.InputError):
            rm.subsample(self.COUNTS, sum(self.COUNTS) + 1, rng)

    def test_seed_reproducible(self):
        a = rm.downsampled_metrics(self.COUNTS, 40, 10, random.Random("7:x"))
        b = rm.downsampled_metrics(self.COUNTS, 40, 10, random.Random("7:x"))
        c = rm.downsampled_metrics(self.COUNTS, 40, 10, random.Random("8:x"))
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

    def test_mean_richness_matches_exact_rarefaction(self):
        result = rm.downsampled_metrics(self.COUNTS, 30, 400, random.Random("3:s"))
        exact = rm.expected_richness(self.COUNTS, 30)
        self.assertAlmostEqual(result["expected_richness"], exact)
        self.assertLess(abs(result["metrics"]["richness"]["mean"] - exact), 0.6)

    def test_floyd_sample_is_uniform(self):
        rng = random.Random(11)
        hits = [0] * 10
        for _ in range(20000):
            chosen = rm.floyd_sample(10, 3, rng)
            self.assertEqual(len(chosen), 3)
            for i in chosen:
                hits[i] += 1
        for h in hits:  # expected 6000 each; SD about 65
            self.assertLess(abs(h - 6000), 300)

    def test_quantile_and_summary(self):
        self.assertEqual(rm.quantile([1.0, 2.0, 3.0, 4.0, 5.0], 0.5), 3.0)
        self.assertAlmostEqual(rm.quantile([1.0, 2.0], 0.975), 1.975)
        s = rm.summarize([1.0, None, 3.0])
        self.assertEqual((s["mean"], s["n_undefined"]), (2.0, 1))


class LoadingAndCommandLine(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_mixed_chains_refused_then_filtered(self):
        path = write(self.dir, "mixed.tsv",
                     "sample\tlocus\tjunction_aa\tduplicate_count\tproductive\n"
                     "S1\tTRB\tCASSA\t5\tT\nS1\tTRA\tCAVR\t3\tT\nS1\tTRB\tCASSB\t2\tF\n"
                     "S1\tTRB\tCASSA\t1\tT\nS1\tTRB\tCASSC\t0\tT\n")
        proc = run_cli("--input", path, "--clone-col", "junction_aa", "--count-col", "duplicate_count")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("several receptor chains", proc.stderr)
        samples, qc, _ = rm.load_repertoires([path], ["junction_aa"], count_col="duplicate_count", chain="TRB")
        self.assertEqual(samples, {"S1": {"CASSA": 6}})
        q = qc["S1"]
        self.assertEqual((q["dropped_other_chain"], q["dropped_nonproductive"], q["dropped_zero_count"]), (1, 1, 1))
        self.assertEqual(q["rows_collapsed"], 1)
        self.assertEqual(q["chain"], "TRB")

    def test_cell_mode_counts_distinct_cells(self):
        path = write(self.dir, "cells.csv",
                     "barcode,sample,chain,cdr3,umis\n"
                     "c1,pre,TRB,CASSX,10\nc1,pre,TRB,CASSX,4\nc2,pre,TRB,CASSX,3\n"
                     "c3,pre,TRB,CASSY,8\nc3,pre,TRB,CASSZ,2\n")
        samples, qc, _ = rm.load_repertoires([path], ["cdr3"], cell_col="barcode")
        self.assertEqual(samples["pre"], {"CASSX": 2, "CASSY": 1, "CASSZ": 1})
        self.assertEqual((qc["pre"]["cells"], qc["pre"]["cells_in_several_clonotypes"]), (3, 1))

    def test_one_file_per_sample_and_order_independent_seeding(self):
        a = write(self.dir, "tumourA.tsv", "cdr3\tcount\n" + "".join(f"C{i}\t{1 + i % 4}\n" for i in range(60)))
        b = write(self.dir, "tumourB.tsv", "cdr3\tcount\n" + "".join(f"D{i}\t{1 + i % 2}\n" for i in range(30)))
        both = rm.analyse(*rm.load_repertoires([a, b], ["cdr3"], count_col="count"), downsample="40", iterations=5, seed=4)
        alone = rm.analyse(*rm.load_repertoires([a], ["cdr3"], count_col="count"), downsample="40", iterations=5, seed=4)
        self.assertEqual([s["sample"] for s in both["samples"]], ["tumourA", "tumourB"])
        self.assertEqual(both["samples"][0]["downsampled"], alone["samples"][0]["downsampled"])

    def test_below_target_excluded_and_min_depth(self):
        path = write(self.dir, "two.tsv", "sample\tcdr3\tcount\nA\tX\t10\nA\tY\t10\nB\tX\t3\nB\tZ\t2\n")
        loaded = rm.load_repertoires([path], ["cdr3"], count_col="count")
        res = rm.analyse(*loaded, downsample="10", iterations=3)
        self.assertIn("skipped", res["samples"][1]["downsampled"])
        self.assertTrue(any("below the downsampling depth" in n for n in res["notes"]))
        self.assertEqual(rm.analyse(*loaded, downsample="min", iterations=3)["settings"]["downsample_depth"], 5)
        self.assertTrue(any("differ" in n for n in rm.analyse(*loaded)["notes"]))

    def test_bad_inputs_exit_2(self):
        path = write(self.dir, "bad.tsv", "sample\tcdr3\tcount\nA\tX\t-1\n")
        proc = run_cli("--input", path, "--clone-col", "cdr3", "--count-col", "count")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("counts must be finite and >= 0", proc.stderr)
        proc = run_cli("--input", path, "--clone-col", "junction_aa", "--count-col", "count")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("junction_aa", proc.stderr)
        freq = write(self.dir, "freq.tsv", "sample\tcdr3\tfrequency\nA\tX\t0.6\nA\tY\t0.4\n")
        proc = run_cli("--input", freq, "--clone-col", "cdr3", "--count-col", "frequency", "--downsample", "min")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("integer counts", proc.stderr)

    def test_cli_json_and_tsv(self):
        path = write(self.dir, "clones.tsv", "sample\tcdr3\tcount\nS\tA\t4\nS\tB\t2\nS\tC\t1\nS\tD\t1\n")
        proc = run_cli("--input", path, "--clone-col", "cdr3", "--count-col", "count", "--top", "1,2")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out = json.loads(proc.stdout)
        self.assertAlmostEqual(out["samples"][0]["observed"]["clonality"], 0.125)
        proc = run_cli("--input", path, "--clone-col", "cdr3", "--count-col", "count",
                       "--downsample", "6", "--iterations", "5", "--format", "tsv")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        header, row = proc.stdout.strip().split("\n")
        self.assertIn("ds_clonality_mean", header.split("\t"))
        self.assertEqual(row.split("\t")[0], "S")


class ReviewAdditions(unittest.TestCase):
    """Cases added in the independent review (ties, exact draw law, file naming, empty files)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_ties_golden_values(self):
        # counts 3, 3, 2, 1, 1 (ties at the top and among singletons): N = 10, S = 5
        # H = -(2 * 0.3 ln 0.3 + 0.2 ln 0.2 + 2 * 0.1 ln 0.1) = 1.5047882836; J = H / ln 5
        # lambda = 0.09 + 0.09 + 0.04 + 0.01 + 0.01 = 0.24; Gini with 1, 1, 2, 3, 3:
        # (-4 - 2 + 0 + 6 + 12) / (5 * 10) = 0.24; Chao1 = 5 + 2 * 1 / (2 * 2) = 5.5
        # top-1 is 0.3 whichever of the tied clonotypes is ranked first
        m = rm.diversity_metrics([1, 3, 2, 3, 1], top=(1, 2))
        h = -(2 * 0.3 * math.log(0.3) + 0.2 * math.log(0.2) + 2 * 0.1 * math.log(0.1))
        self.assertAlmostEqual(m["shannon"], h, places=12)
        self.assertAlmostEqual(m["clonality"], 1 - h / math.log(5), places=12)
        self.assertAlmostEqual(m["inverse_simpson"], 1 / 0.24, places=12)
        self.assertAlmostEqual(m["gini"], 0.24, places=12)
        self.assertAlmostEqual(m["chao1"], 5.5)
        self.assertAlmostEqual(m["top1_fraction"], 0.3)
        self.assertAlmostEqual(m["top2_fraction"], 0.6)
        self.assertAlmostEqual(m["expanded_fraction"], 0.8)

    def test_expected_richness_equals_brute_force_enumeration(self):
        # all C(10, 3) = 120 subsets of the 10 units of counts 3, 3, 2, 1, 1:
        # mean number of distinct clonotypes = 306 / 120 = 2.55 (Hurlbert: 170/120 + 64/120 + 72/120)
        counts = [3, 3, 2, 1, 1]
        units = [i for i, c in enumerate(counts) for _ in range(c)]
        subsets = list(itertools.combinations(range(len(units)), 3))
        brute = sum(len({units[u] for u in s}) for s in subsets) / len(subsets)
        self.assertAlmostEqual(brute, 2.55, places=12)
        self.assertAlmostEqual(rm.expected_richness(counts, 3), brute, places=12)

    def test_subsample_follows_the_hypergeometric_law(self):
        # drawing 2 of the 4 units of counts 2, 1, 1: P(2,0,0) = P(0,1,1) = 1/6, P(1,1,0) = P(1,0,1) = 2/6
        rng = random.Random("law")
        n = 30000
        seen = collections.Counter(tuple(rm.subsample([2, 1, 1], 2, rng)) for _ in range(n))
        exact = {(2, 0, 0): 1 / 6, (1, 1, 0): 2 / 6, (1, 0, 1): 2 / 6, (0, 1, 1): 1 / 6}
        self.assertEqual(set(seen), set(exact))
        for outcome, prob in exact.items():   # SD of a frequency is at most 0.003 here
            self.assertLess(abs(seen[outcome] / n - prob), 0.012)

    def test_sample_names_from_files_keep_inner_dots(self):
        self.assertEqual(rm.sample_name_from_path("/x/P1.pre.tsv"), "P1.pre")
        self.assertEqual(rm.sample_name_from_path("P1.post.csv.gz"), "P1.post")
        self.assertEqual(rm.sample_name_from_path("run7.clones"), "run7.clones")
        pre = write(self.dir, "P1.pre.tsv", "cdr3\tcount\nX\t3\nY\t1\n")
        post = write(self.dir, "P1.post.tsv", "cdr3\tcount\nX\t1\nZ\t5\n")
        samples, _, _ = rm.load_repertoires([pre, post], ["cdr3"], count_col="count")
        self.assertEqual(samples, {"P1.pre": {"X": 3, "Y": 1}, "P1.post": {"X": 1, "Z": 5}})

    def test_same_sample_from_two_files_is_pooled_with_a_note(self):
        os.makedirs(os.path.join(self.dir, "lane1"))
        os.makedirs(os.path.join(self.dir, "lane2"))
        a = write(self.dir, os.path.join("lane1", "S.tsv"), "cdr3\tcount\nX\t3\n")
        b = write(self.dir, os.path.join("lane2", "S.tsv"), "cdr3\tcount\nX\t1\nY\t2\n")
        samples, _, notes = rm.load_repertoires([a, b], ["cdr3"], count_col="count")
        self.assertEqual(samples, {"S": {"X": 4, "Y": 2}})
        self.assertTrue(any("pools rows from 2 files" in n for n in notes))

    def test_header_only_file_exits_2(self):
        path = write(self.dir, "empty.tsv", "cdr3\tcount\n")
        proc = run_cli("--input", path, "--clone-col", "cdr3", "--count-col", "count", "--format", "tsv")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("no data rows", proc.stderr)

    def test_shallow_sample_note_only_with_three_samples(self):
        self.assertIsNone(rm.shallow_note({"A": 100, "B": 5}))
        self.assertIsNone(rm.shallow_note({"A": 100, "B": 90, "C": 30}))
        self.assertIn("'C'", rm.shallow_note({"A": 100, "B": 90, "C": 10}))
        path = write(self.dir, "three.tsv", "sample\tcdr3\tcount\nA\tx\t50\nA\ty\t50\nB\tx\t40\nB\tz\t40\nC\tx\t2\nC\tq\t1\n")
        res = rm.analyse(*rm.load_repertoires([path], ["cdr3"], count_col="count"), downsample="min", iterations=2)
        self.assertTrue(any("fifth of the median" in n for n in res["notes"]))


if __name__ == "__main__":
    unittest.main()
