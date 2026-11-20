"""Unit tests for skills/single-cell-patient-level-stats/scripts/composition_test.py."""
import importlib.util
import itertools
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "single-cell-patient-level-stats" / "scripts" / "composition_test.py"
_spec = importlib.util.spec_from_file_location("composition_test_under_test", SCRIPT)
ct = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ct)


def independent_midranks(values):
    # rank_i = #(values below v_i) + (#(values equal to v_i) + 1) / 2, computed without sorting
    return [sum(1 for w in values if w < v) + (sum(1 for w in values if w == v) + 1) / 2.0 for v in values]


def brute_ranksum_p(x, y):
    pooled = list(x) + list(y)
    ranks = independent_midranks(pooled)
    n1, big_n = len(x), len(pooled)
    mean = n1 * (big_n + 1) / 2.0
    observed = abs(sum(ranks[:n1]) - mean)
    hits = total = 0
    for combo in itertools.combinations(range(big_n), n1):
        total += 1
        if abs(sum(ranks[i] for i in combo) - mean) >= observed - 1e-9:
            hits += 1
    return hits / total


def brute_signedrank_p(diffs):
    nonzero = [d for d in diffs if d != 0]
    ranks = independent_midranks([abs(d) for d in nonzero])
    mean = sum(ranks) / 2.0
    observed = abs(sum(r for r, d in zip(ranks, nonzero) if d > 0) - mean)
    hits = total = 0
    for signs in itertools.product((0, 1), repeat=len(nonzero)):
        total += 1
        if abs(sum(r for r, s in zip(ranks, signs) if s) - mean) >= observed - 1e-9:
            hits += 1
    return hits / total


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


class TestRankStatistics(unittest.TestCase):
    def test_midranks_with_ties(self):
        self.assertEqual(ct.midranks([3, 1, 3, 2]), [3.5, 1.0, 3.5, 2.0])

    def test_ranksum_exact_golden_no_ties(self):
        # 3 vs 3, complete separation: of the C(6,3) = 20 equally likely splits only {1,2,3}
        # and {4,5,6} are this extreme, so p = 2/20 = 0.1, which is also the floor.
        p, floor = ct.ranksum_exact([1, 2, 3], [4, 5, 6])
        self.assertAlmostEqual(p, 0.1)
        self.assertAlmostEqual(floor, 0.1)

    def test_ranksum_floor_is_two_over_binomial(self):
        # Smallest two-sided exact p without ties = 2 / C(n1 + n2, n1).
        for n1, n2 in ((2, 2), (3, 4), (4, 4), (5, 5)):
            p, floor = ct.ranksum_exact(list(range(n1)), list(range(10, 10 + n2)))
            self.assertAlmostEqual(floor, 2 / math.comb(n1 + n2, n1))
            self.assertAlmostEqual(p, floor)

    def test_ranksum_exact_golden_with_ties(self):
        # x = [0, 0, 1], y = [0, 2, 3]. Midranks: the three zeros share rank 2; 1 -> 4, 2 -> 5,
        # 3 -> 6. Doubled scores {4,4,4,8,10,12}; T2(x) = 16; E[T2] = 3 * 7 = 21, |dev| = 5.
        # The 20 subsets of size 3 give sums 12(1) 16(3) 18(3) 20(3) 22(3) 24(3) 26(3) 30(1);
        # |sum - 21| >= 5 for 12, 16, 26, 30: (1 + 3 + 3 + 1) / 20 = 0.4. Floor: sums 12, 30 = 2/20.
        p, floor = ct.ranksum_exact([0, 0, 1], [0, 2, 3])
        self.assertAlmostEqual(p, 0.4)
        self.assertAlmostEqual(floor, 0.1)

    def test_ranksum_exact_matches_bruteforce(self):
        x = [0.1, 0.2, 0.2, 0.0, 0.5]
        y = [0.2, 0.3, 0.0, 0.4, 0.4, 0.05]
        p, _ = ct.ranksum_exact(x, y)
        self.assertAlmostEqual(p, brute_ranksum_p(x, y))
        p_swapped, _ = ct.ranksum_exact(y, x)
        self.assertAlmostEqual(p, p_swapped)

    def test_ranksum_normal_golden(self):
        # 1..5 vs 6..10: T = 15, E[T] = 27.5, Var = 5*5*11/12 = 22.9167;
        # z = (12.5 - 0.5) / 4.7871 = 2.5067; p = erfc(z / sqrt 2) = 0.01219.
        p = ct.ranksum_normal([1, 2, 3, 4, 5], [6, 7, 8, 9, 10])
        self.assertAlmostEqual(p, math.erfc(12 / math.sqrt(275 / 12) / math.sqrt(2)))
        self.assertAlmostEqual(p, 0.01219, places=4)

    def test_normal_close_to_exact_for_moderate_n(self):
        x = [i * 0.01 for i in range(20)]
        y = [0.05 + i * 0.01 for i in range(20)]
        exact = ct.compare_unpaired(x, y, method="exact")["p_value"]
        normal = ct.compare_unpaired(x, y, method="normal")["p_value"]
        self.assertLess(abs(exact - normal), 0.01)

    def test_all_values_tied(self):
        res = ct.compare_unpaired([0.0, 0.0], [0.0, 0.0, 0.0])
        self.assertEqual(res["p_value"], 1.0)
        self.assertEqual(res["rank_biserial"], 0.0)
        self.assertEqual(ct.ranksum_normal([1, 1], [1, 1]), 1.0)

    def test_unpaired_effect_sizes_golden(self):
        # x = [0,0,1] vs y = [0,2,3]: x > y in 1 pair, ties in 2 pairs, so U_x = 1 + 0.5 * 2 = 2;
        # rank-biserial = 2 * 2 / 9 - 1 = -5/9. The 9 differences x - y sorted:
        # -3,-3,-2,-2,-2,-1,0,0,1 -> Hodges-Lehmann shift = median = -2.
        res = ct.compare_unpaired([0, 0, 1], [0, 2, 3])
        self.assertAlmostEqual(res["statistic_value"], 2.0)
        self.assertAlmostEqual(res["rank_biserial"], -5 / 9)
        self.assertAlmostEqual(res["hl_shift"], -2.0)
        self.assertEqual(res["method"], "exact")

    def test_signedrank_golden(self):
        # Differences 0.1, 0.2, -0.05, 0.3, 0 (floating noise removed by rounding). Zero dropped:
        # ranks of |d| = 2, 3, 1, 4; W+ = 9, W- = 1; rank-biserial = (9 - 1) / 10 = 0.8.
        # n = 4: W+ in {0, 1, 9, 10} has probability 4/16 -> p = 0.25; floor = 2/16 = 0.125.
        # Walsh averages of all five differences: 15 values, median (8th) = 0.1.
        res = ct.compare_paired([0.3, 0.4, 0.15, 0.5, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2])
        self.assertEqual((res["n_up"], res["n_down"], res["n_zero"]), (3, 1, 1))
        self.assertAlmostEqual(res["statistic_value"], 9.0)
        self.assertAlmostEqual(res["p_value"], 0.25)
        self.assertAlmostEqual(res["p_floor"], 0.125)
        self.assertAlmostEqual(res["rank_biserial"], 0.8)
        self.assertAlmostEqual(res["hl_shift"], 0.1)

    def test_signedrank_five_pairs_cannot_reach_005(self):
        # 5 pairs all in one direction: p = 2 / 2^5 = 0.0625, the floor for 5 pairs.
        res = ct.compare_paired([2, 3, 4, 5, 6], [1, 1, 1, 1, 1])
        self.assertAlmostEqual(res["p_value"], 0.0625)
        self.assertAlmostEqual(res["p_floor"], 0.0625)

    def test_signedrank_matches_bruteforce_with_ties(self):
        diffs = [1, -1, 2, 2, -3, 4, 1, 0]
        p, _ = ct.signedrank_exact([d for d in diffs if d != 0])
        self.assertAlmostEqual(p, brute_signedrank_p(diffs))

    def test_paired_all_zero_differences(self):
        res = ct.compare_paired([0.1, 0.2], [0.1, 0.2])
        self.assertEqual(res["p_value"], 1.0)
        self.assertEqual(res["n_zero"], 2)

    def test_bh_golden(self):
        # Sorted p: 0.01, 0.03, 0.04, 0.20 (m = 4) -> raw m*p/rank: 0.04, 0.06, 0.0533, 0.20;
        # step-up cumulative minimum from the top: 0.04, 0.0533, 0.0533, 0.20.
        adj = ct.bh_adjust([0.01, 0.04, 0.03, 0.20])
        for got, want in zip(adj, [0.04, 0.16 / 3, 0.16 / 3, 0.20]):
            self.assertAlmostEqual(got, want)
        self.assertEqual(ct.bh_adjust([]), [])


class TestCompositionDesign(unittest.TestCase):
    def per_cell_rows(self):
        # Synthetic: A1, A2 in group A; B1, B2 in group B. B1 has no 'Bcell' rows at all.
        layout = {"A1": ("A", {"T": 6, "Bcell": 2, "Mye": 2}), "A2": ("A", {"T": 5, "Bcell": 1, "Mye": 4}),
                  "B1": ("B", {"T": 2, "Mye": 8}), "B2": ("B", {"T": 3, "Bcell": 1, "Mye": 6})}
        lines = ["sample_id\tgroup\tcell_type"]
        for s, (g, cts) in layout.items():
            for c, n in cts.items():
                lines += [f"{s}\t{g}\t{c}"] * n
        return "\n".join(lines) + "\n"

    def run_analysis(self, text, **kwargs):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "in.tsv")
            write(path, text)
            fields, rows = ct.read_rows(path)
        count_col = kwargs.pop("count_col", None)
        pair_col = kwargs.pop("pair_col", None)
        patient_col = kwargs.pop("patient_col", None)
        counts, labels = ct.collect(rows, fields, "sample_id", "cell_type", "group", count_col, pair_col, patient_col)
        return ct.analyse(counts, labels, kwargs.pop("case", "A"), kwargs.pop("control", "B"),
                          paired=bool(pair_col), use_patient=bool(patient_col), **kwargs)

    def test_absent_cell_type_is_zero_not_missing(self):
        report = self.run_analysis(self.per_cell_rows())
        bcell = next(r for r in report["results"] if r["cell_type"] == "Bcell")
        self.assertEqual((bcell["n_case"], bcell["n_control"]), (2, 2))
        b1 = next(s for s in report["samples"] if s["sample"] == "B1")
        self.assertEqual(b1["values"]["Bcell"], 0.0)
        self.assertAlmostEqual(b1["values"]["T"], 0.2)
        self.assertTrue(any("cannot reach" in n or "no cell type can reach" in n for n in report["notes"]))

    def test_within_changes_denominator(self):
        report = self.run_analysis(self.per_cell_rows(), within=["T", "Bcell"])
        a1 = next(s for s in report["samples"] if s["sample"] == "A1")
        self.assertAlmostEqual(a1["values"]["T"], 6 / 8)
        self.assertNotIn("Mye", a1["values"])

    def test_alr_transform(self):
        report = self.run_analysis(self.per_cell_rows(), transform="alr", reference="Mye")
        a1 = next(s for s in report["samples"] if s["sample"] == "A1")
        self.assertAlmostEqual(a1["values"]["T"], math.log2(6.5 / 2.5))
        self.assertNotIn("Mye", [r["cell_type"] for r in report["results"]])

    def paired_counts(self):
        # Synthetic long counts: P1-P3 complete, P4 has only the 'B' sample.
        rows = ["sample_id\tpatient\tgroup\tcell_type\tn_cells"]
        data = {("P1a", "P1", "A"): (30, 70), ("P1b", "P1", "B"): (10, 90), ("P2a", "P2", "A"): (25, 75),
                ("P2b", "P2", "B"): (20, 80), ("P3a", "P3", "A"): (40, 60), ("P3b", "P3", "B"): (15, 85),
                ("P4b", "P4", "B"): (5, 95)}
        for (s, p, g), (t, m) in data.items():
            rows += [f"{s}\t{p}\t{g}\tT\t{t}", f"{s}\t{p}\t{g}\tMye\t{m}"]
        return "\n".join(rows) + "\n"

    def test_paired_drops_incomplete_pair(self):
        report = self.run_analysis(self.paired_counts(), count_col="n_cells", pair_col="patient")
        self.assertEqual(report["design"]["n_pairs"], 3)
        self.assertEqual(report["design"]["incomplete_pairs"], ["P4"])
        t = next(r for r in report["results"] if r["cell_type"] == "T")
        self.assertEqual((t["n_up"], t["n_down"]), (3, 0))
        self.assertAlmostEqual(t["p_value"], 0.25)  # 3 pairs: floor 2/2^3

    def test_paired_duplicate_sample_in_group_is_error(self):
        text = self.paired_counts() + "P1c\tP1\tA\tT\t5\nP1c\tP1\tA\tMye\t5\n"
        with self.assertRaises(ct.InputError):
            self.run_analysis(text, count_col="n_cells", pair_col="patient")

    def test_unpaired_patient_in_both_groups_is_error(self):
        with self.assertRaisesRegex(ct.InputError, "paired design"):
            self.run_analysis(self.paired_counts(), count_col="n_cells", patient_col="patient")

    def test_duplicate_count_rows_and_bad_counts(self):
        text = "sample_id\tgroup\tcell_type\tn_cells\nS1\tA\tT\t3\nS1\tA\tT\t4\nS2\tB\tT\t1\n"
        with self.assertRaisesRegex(ct.InputError, "more than once"):
            self.run_analysis(text, count_col="n_cells")
        text = "sample_id\tgroup\tcell_type\tn_cells\nS1\tA\tT\t3.5\nS2\tB\tT\t1\n"
        with self.assertRaisesRegex(ct.InputError, "whole number"):
            self.run_analysis(text, count_col="n_cells")

    def test_inconsistent_group_and_missing_group(self):
        text = "sample_id\tgroup\tcell_type\nS1\tA\tT\nS1\tB\tT\nS2\tB\tT\n"
        with self.assertRaisesRegex(ct.InputError, "constant within a sample"):
            self.run_analysis(text)
        with self.assertRaisesRegex(ct.InputError, "not found"):
            self.run_analysis(self.per_cell_rows(), case="C")


class TestIndependentGoldens(unittest.TestCase):
    """Reviewer goldens: hand derivations, cross-checked by exact-fraction enumeration in a
    separate script that shares no code with composition_test.py."""

    def test_ranksum_ties_four_vs_three(self):
        # Case [0.1, 0.3, 0.3, 0] vs control [0.3, 0, 0.05]. Midranks: 0, 0 -> 1.5; 0.05 -> 3;
        # 0.1 -> 4; the three 0.3 -> 6. Control rank sum 6 + 1.5 + 3 = 10.5, E = 3 * 8 / 2 = 12,
        # |dev| = 1.5. Of the 35 control 3-subsets of {1.5, 1.5, 3, 4, 6, 6, 6}, only sums 11.5 (6
        # subsets) and 13 (3 subsets) lie closer to 12, so p = 26/35. Floor: sums 6 and 18 only,
        # 2/35. Pairs: 6 favour the case, 3 the control, 3 tied -> r = (6 - 3) / 12 = 0.25.
        # Differences sorted: -0.3 -0.2 -0.05 0 0 0 0.05 0.1 0.25 0.25 0.3 0.3 -> HL = 0.025.
        res = ct.compare_unpaired([0.1, 0.3, 0.3, 0.0], [0.3, 0.0, 0.05])
        self.assertAlmostEqual(res["p_value"], 26 / 35)
        self.assertAlmostEqual(res["p_floor"], 2 / 35)
        self.assertAlmostEqual(res["rank_biserial"], 0.25)
        self.assertAlmostEqual(res["hl_shift"], 0.025)

    def test_ranksum_two_vs_four_no_ties(self):
        # Complete separation: p = floor = 2 / C(6, 2) = 2/15; all 8 pairs favour the case (r = 1);
        # differences 1, 2, 3, 4, 5, 6, 7, 8 -> HL = 4.5.
        res = ct.compare_unpaired([5, 9], [1, 2, 3, 4])
        self.assertAlmostEqual(res["p_value"], 2 / 15)
        self.assertAlmostEqual(res["p_floor"], 2 / 15)
        self.assertAlmostEqual(res["rank_biserial"], 1.0)
        self.assertAlmostEqual(res["hl_shift"], 4.5)

    def test_signedrank_ties_and_zero(self):
        # d = 2, -2, 2, 1, 0, 3; the zero leaves the test. |d| = 2, 2, 2, 1, 3 -> midranks 3, 3, 3, 1, 5.
        # W+ = 3 + 3 + 1 + 5 = 12, W- = 3, E = 7.5, |dev| = 4.5. Over the 32 sign patterns of scores
        # {1, 3, 3, 3, 5}, sums <= 3 come from {}, {1} and each single 3 (5 patterns), and by symmetry
        # 5 patterns give sums >= 12: p = 10/32. Floor: sums 0 and 15, 2/32. r = (12 - 3) / 15 = 0.6.
        # Walsh averages of all six differences (zero included): 21 values, 11th smallest = 1.0.
        res = ct.compare_paired([2, -2, 2, 1, 0, 3], [0, 0, 0, 0, 0, 0])
        self.assertEqual((res["n_up"], res["n_down"], res["n_zero"]), (4, 1, 1))
        self.assertAlmostEqual(res["statistic_value"], 12.0)
        self.assertAlmostEqual(res["p_value"], 10 / 32)
        self.assertAlmostEqual(res["p_floor"], 2 / 32)
        self.assertAlmostEqual(res["rank_biserial"], 0.6)
        self.assertAlmostEqual(res["hl_shift"], 1.0)

    def test_signedrank_normal_golden(self):
        # d = 1..10, all positive: W+ = 55, E = 27.5, Var = 10 * 11 * 21 / 24 = 96.25;
        # z = (27.5 - 0.5) / sqrt(96.25) = 2.7521 -> p = erfc(z / sqrt 2) = 0.005922.
        self.assertAlmostEqual(ct.signedrank_normal(list(range(1, 11))), 0.005922, places=6)

    def test_bh_with_tied_pvalues(self):
        # Sorted p 0.01, 0.02, 0.04, 0.04, 0.5 (m = 5) -> m p / rank = 0.05, 0.05, 0.0667, 0.05, 0.5;
        # cumulative minimum from the top: 0.05, 0.05, 0.05, 0.05, 0.5.
        adj = ct.bh_adjust([0.02, 0.5, 0.01, 0.04, 0.04])
        for got, want in zip(adj, [0.05, 0.5, 0.05, 0.05, 0.05]):
            self.assertAlmostEqual(got, want)

    def test_field_count_errors_and_blank_lines(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "in.csv")
            write(path, "sample_id,group,cell_type,n_cells\nS1,A,T,3,EXTRA\nS2,B,T,1\n")
            with self.assertRaisesRegex(ct.InputError, "line 2: more fields"):
                ct.read_rows(path)
            write(path, "sample_id,group,cell_type,n_cells\nS1,A,T,3\nS2,B,T\n")
            with self.assertRaisesRegex(ct.InputError, "line 3: fewer fields"):
                ct.read_rows(path)
            write(path, "sample_id,group,cell_type,n_cells\nS1,A,T,3\n\nS2,B,T,1\n")
            self.assertEqual(len(ct.read_rows(path)[1]), 2)


class TestEvalInput(unittest.TestCase):
    def test_paired_eval_counts(self):
        # Synthetic eval input. CD8_T on minus pre, by hand from the counts: P1 420/1905 - 120/1600 =
        # 0.1455, P2 380/2100 - 200/2000 = 0.0810, P3 160/1075 - 90/1000 = 0.0588, P4 330/2820 -
        # 300/3000 = 0.0170, P5 150/915 - 60/800 = 0.0889. All five up: p = floor = 2/32; median and
        # Hodges-Lehmann shift 0.08095. Myeloid falls in 5 of 5 (p = 0.0625) although its counts
        # barely move. P6 has no on-treatment sample.
        path = ROOT / "skills" / "single-cell-patient-level-stats" / "evals" / "paired-composition-on-treatment" / "files" / "counts.tsv"
        fields, rows = ct.read_rows(str(path))
        counts, labels = ct.collect(rows, fields, "sample_id", "cell_type", "timepoint", "n_cells", "patient_id")
        report = ct.analyse(counts, labels, "on", "pre", paired=True)
        self.assertEqual(report["design"]["incomplete_pairs"], ["P6"])
        cd8 = next(r for r in report["results"] if r["cell_type"] == "CD8_T")
        self.assertEqual((cd8["n_pairs"], cd8["n_up"]), (5, 5))
        self.assertAlmostEqual(cd8["p_value"], 0.0625)
        self.assertAlmostEqual(cd8["p_floor"], 0.0625)
        self.assertAlmostEqual(cd8["median_difference"], 0.0809524, places=6)
        self.assertAlmostEqual(cd8["hl_shift"], 0.0809524, places=6)
        mye = next(r for r in report["results"] if r["cell_type"] == "Myeloid")
        self.assertEqual(mye["n_down"], 5)


class TestCommandLine(unittest.TestCase):
    def test_cli_smoke_and_exit_code(self):
        rows = ["sample_id,group,cell_type,n_cells"]
        for i, (t, m) in enumerate([(30, 70), (35, 65), (28, 72), (10, 90), (12, 88), (15, 85)]):
            g = "R" if i < 3 else "NR"
            rows += [f"S{i},{g},CD8_T,{t}", f"S{i},{g},Myeloid,{m}"]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "counts.csv")
            tsv = os.path.join(d, "res.tsv")
            write(path, "\n".join(rows) + "\n")
            proc = subprocess.run([sys.executable, str(SCRIPT), path, "--count-col", "n_cells", "--group-col",
                                   "group", "--case", "R", "--control", "NR", "--tsv-out", tsv],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            report = json.loads(proc.stdout)
            cd8 = next(r for r in report["results"] if r["cell_type"] == "CD8_T")
            self.assertAlmostEqual(cd8["p_value"], 0.1)  # 3 vs 3 complete separation
            self.assertAlmostEqual(cd8["rank_biserial"], 1.0)
            with open(tsv, encoding="utf-8") as fh:
                self.assertTrue(fh.readline().startswith("cell_type\tn_case"))
            bad = subprocess.run([sys.executable, str(SCRIPT), path, "--group-col", "arm", "--case", "R",
                                  "--control", "NR"], capture_output=True, text=True)
            self.assertEqual(bad.returncode, 2)
            self.assertIn("arm", bad.stderr)
            empty = os.path.join(d, "empty.tsv")
            write(empty, "")
            proc = subprocess.run([sys.executable, str(SCRIPT), empty, "--group-col", "g", "--case", "a",
                                   "--control", "b"], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)


if __name__ == "__main__":
    unittest.main()
