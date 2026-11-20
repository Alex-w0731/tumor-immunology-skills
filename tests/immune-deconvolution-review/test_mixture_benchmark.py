"""Unit tests for skills/immune-deconvolution-review/scripts/mixture_benchmark.py."""
import copy
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "immune-deconvolution-review" / "scripts" / "mixture_benchmark.py"
spec = importlib.util.spec_from_file_location("mixture_benchmark", SCRIPT)
mb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mb)


def tables(est_rows, true_rows, cols):
    """Build (ids, cols, data) triples from {sample: [values]} dictionaries."""
    def build(rows):
        return list(rows), list(cols), {s: dict(zip(cols, v)) for s, v in rows.items()}
    return build(est_rows), build(true_rows)


def run(est_rows, true_rows, cols, **kw):
    (ei, ec, e), (ti, tc, t) = tables(est_rows, true_rows, cols)
    return mb.benchmark(ei, ec, e, ti, tc, t, **kw)


class TestStatistics(unittest.TestCase):
    def test_spearman_with_ties_golden(self):
        # x = (1, 2, 2, 3) -> ranks (1, 2.5, 2.5, 4); y = (1, 3, 2, 4) -> ranks (1, 3, 2, 4).
        # Centred: (-1.5, 0, 0, 1.5) and (-1.5, 0.5, -0.5, 1.5); sxy = 4.5, sxx = 4.5, syy = 5,
        # rho = 4.5 / sqrt(22.5) = 0.9486833.
        self.assertEqual(mb.ranks([1, 2, 2, 3]), [1.0, 2.5, 2.5, 4.0])
        self.assertAlmostEqual(mb.spearman([1, 2, 2, 3], [1, 3, 2, 4]), 4.5 / math.sqrt(22.5), places=12)

    def test_edge_cases(self):
        self.assertIsNone(mb.pearson([0.1, 0.2], [0.1, 0.3]))          # below min_n = 3
        self.assertIsNone(mb.pearson([0.0, 0.0, 0.0], [0.1, 0.2, 0.3]))  # constant truth
        self.assertEqual(mb.ols([1.0, 1.0], [2.0, 3.0]), (None, None))
        self.assertIsNone(mb.rmse([], []))


class TestBenchmark(unittest.TestCase):
    def test_constant_shift_golden(self):
        # est = true + 0.1 for T: r = 1, bias = 0.1, RMSE = 0.1, slope = 1, intercept = 0.1.
        res = run({"s1": [0.2], "s2": [0.3], "s3": [0.4]}, {"s1": [0.1], "s2": [0.2], "s3": [0.3]}, ["T"])
        rec = res["per_cell_type"][0]
        self.assertEqual(rec["pearson_r"], 1.0)
        self.assertAlmostEqual(rec["bias"], 0.1, places=5)
        self.assertAlmostEqual(rec["rmse"], 0.1, places=5)
        self.assertAlmostEqual(rec["slope"], 1.0, places=4)
        self.assertAlmostEqual(rec["intercept"], 0.1, places=5)

    def test_per_sample_golden(self):
        # true (0.2, 0.3, 0.5), est (0.25, 0.25, 0.5): sxy = 1/24, sxx = 7/150, syy = 1/24,
        # r = sqrt(150/168) = 0.944911; RMSE = sqrt((0.05^2 + 0.05^2 + 0) / 3) = 0.040825.
        res = run({"s1": [0.25, 0.25, 0.5]}, {"s1": [0.2, 0.3, 0.5]}, ["A", "B", "C"])
        rec = res["per_sample"][0]
        self.assertAlmostEqual(rec["pearson_r"], 0.9449, places=4)
        self.assertAlmostEqual(rec["rmse"], 0.04082, places=5)
        self.assertIsNone(res["per_cell_type"][0]["pearson_r"])        # one sample only
        self.assertTrue(any("correlation undefined" in w for w in res["warnings"]))

    def test_background_prediction(self):
        res = run({"a": [0.05], "b": [0.03], "c": [0.25], "d": [0.35]},
                  {"a": [0.0], "b": [0.0], "c": [0.2], "d": [0.4]}, ["X"])
        rec = res["per_cell_type"][0]
        self.assertEqual(rec["n_true_zero"], 2)
        self.assertAlmostEqual(rec["mean_est_when_true_zero"], 0.04, places=5)

    def test_denominator_mismatch_and_renormalize(self):
        # Truth covers only half of each sample (the rest is tumour); estimates sum to 1.
        true = {"s1": [0.1, 0.4], "s2": [0.2, 0.3], "s3": [0.3, 0.2]}
        est = {s: [v[0] / 0.5, v[1] / 0.5] for s, v in true.items()}
        res = run(est, copy.deepcopy(true), ["T", "M"])
        self.assertTrue(any("row sums differ" in w for w in res["warnings"]))
        self.assertTrue(any("wrong scale" in w for w in res["warnings"]))
        self.assertGreater(res["per_cell_type"][0]["bias"], 0.09)
        fixed = run(est, copy.deepcopy(true), ["T", "M"], renorm="truth")
        self.assertAlmostEqual(fixed["per_cell_type"][0]["bias"], 0.0, places=5)
        self.assertFalse(any("row sums differ" in w for w in fixed["warnings"]))

    def test_merge_and_unmatched(self):
        (ei, ec, e), _ = tables({"s1": [0.1, 0.2, 0.7], "s2": [0.2, 0.2, 0.6], "s3": [0.3, 0.1, 0.6], "x": [0, 0, 1]},
                                {}, ["M1", "M2", "T"])
        ti, tc = ["s1", "s2", "s3"], ["Macrophage", "T"]
        t = {"s1": {"Macrophage": 0.3, "T": 0.7}, "s2": {"Macrophage": 0.4, "T": 0.6}, "s3": {"Macrophage": 0.4, "T": 0.6}}
        res = mb.benchmark(ei, ec, e, ti, tc, t, merges=mb.parse_merge(["Macrophage=M1+M2"]))
        self.assertEqual(res["cell_types"], ["Macrophage", "T"])
        self.assertAlmostEqual(res["per_cell_type"][0]["rmse"], 0.0, places=5)
        self.assertEqual(res["unmatched_samples"]["estimated_only"], ["x"])
        with self.assertRaisesRegex(mb.InputError, "needs M3"):
            mb.apply_merge(["M1", "M2"], {"s": {"M1": 1, "M2": 1}}, mb.parse_merge(["Mac=M1+M2+M3"]), "est")
        with self.assertRaisesRegex(mb.InputError, "must look like"):
            mb.parse_merge(["Macrophage"])

    def test_missing_values_and_no_overlap(self):
        res = run({"a": [None], "b": [0.2], "c": [0.3], "d": [0.4]}, {"a": [0.1], "b": [0.2], "c": [0.3], "d": [0.5]}, ["T"])
        self.assertEqual(res["per_cell_type"][0]["n"], 3)
        self.assertEqual(res["n_missing_pairs"], 1)
        (ei, ec, e), (ti, tc, t) = tables({"a": [0.1]}, {"a": [0.1]}, ["T"])
        with self.assertRaisesRegex(mb.InputError, "no cell-type column"):
            mb.benchmark(ei, ["B"], {"a": {"B": 0.1}}, ti, tc, t)


class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.est = d / "est.tsv"
        self.est.write_text("sample\tT\tB\tcoef_sum\ns1\t0.2\t0.8\t1.3\ns2\t0.3\t0.7\t1.1\ns3\t0.4\t0.6\t0.9\n")
        self.true = d / "truth.tsv"
        self.true.write_text("sample\tT\tB\ns1\t0.1\t0.9\ns2\t0.2\t0.8\ns3\t0.3\t0.7\n")

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, *extra):
        return subprocess.run([sys.executable, str(SCRIPT), "--estimated", str(self.est), "--truth", str(self.true), *extra],
                              capture_output=True, text=True, timeout=60)

    def test_json_and_tsv(self):
        proc = self.cli()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out = json.loads(proc.stdout)
        self.assertEqual(out["ignored_columns"]["estimated_only"], ["coef_sum"])
        self.assertAlmostEqual(out["per_cell_type"][0]["bias"], 0.1, places=5)
        proc = self.cli("--format", "tsv")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(proc.stdout.startswith("cell_type\tn\tpearson_r"))

    def test_bad_input_exits_2_and_help(self):
        self.true.write_text("sample\tT\ns1\tabc\n")
        proc = self.cli()
        self.assertEqual(proc.returncode, 2)
        self.assertIn("'abc'", proc.stderr)
        proc = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True, timeout=60)
        self.assertEqual(proc.returncode, 0)
        self.assertIn("usage", proc.stdout.lower())

    def test_other_errors_exit_2(self):
        proc = self.cli("--out", str(Path(self.tmp.name) / "missing_dir" / "b.json"))
        self.assertEqual(proc.returncode, 2)
        self.assertIn("cannot write output", proc.stderr)
        proc = self.cli("--min-n", "1")
        self.assertEqual(proc.returncode, 2)
        self.true.write_text("sample\tT\tB\ns1\t0.1\t0.9\ns1\t0.2\t0.8\n")
        proc = self.cli()
        self.assertEqual(proc.returncode, 2)
        self.assertIn("repeated sample ID", proc.stderr)
        self.true.write_text("sample\tT\tB\ns1\tinf\t0.9\n")
        proc = self.cli()
        self.assertEqual(proc.returncode, 2)
        self.assertIn("infinite", proc.stderr)


class TestReviewerGoldens(unittest.TestCase):
    """Golden values derived by hand during review, independent of this module's code."""

    def test_spearman_three_way_tie(self):
        # x = (1,1,2,3,3,3) -> ranks (1.5,1.5,3,5,5,5); y = (2,1,4,3,5,6) -> ranks equal to y.
        # Centred: dx = (-2,-2,-0.5,1.5,1.5,1.5), dy = (-1.5,-2.5,0.5,-0.5,1.5,2.5);
        # sxy = 13, sxx = 15, syy = 17.5, so rho = 13 / sqrt(262.5) = 0.8023774.
        self.assertEqual(mb.ranks([1, 1, 2, 3, 3, 3]), [1.5, 1.5, 3.0, 5.0, 5.0, 5.0])
        self.assertAlmostEqual(mb.spearman([1, 1, 2, 3, 3, 3], [2, 1, 4, 3, 5, 6]), 13 / math.sqrt(262.5), places=12)

    def test_renormalize_flags_rows_with_missing_values(self):
        true = {"s1": [0.1, 0.3, None], "s2": [0.2, 0.2, 0.1], "s3": [0.3, 0.1, 0.1]}
        est = {"s1": [0.25, 0.75, 0.0], "s2": [0.4, 0.4, 0.2], "s3": [0.6, 0.2, 0.2]}
        res = run(est, true, ["A", "B", "C"], renorm="truth")
        self.assertTrue(any("renormalized over fewer cell types" in w for w in res["warnings"]))
        # s2: truth (0.2, 0.2, 0.1) / 0.5 = (0.4, 0.4, 0.2), identical to the estimate.
        self.assertAlmostEqual(res["per_sample"][1]["rmse"], 0.0, places=9)

    def test_synthetic_eval_case(self):
        # Pins what the benchmark-pseudobulk-zh graders rely on. Recomputed in review with plain
        # formulas: B cell bias +0.1030 and slope 1.384 after renormalizing the truth; Treg mean
        # estimate where the truth is 0 (PB03, PB06, PB09, PB11) = 0.1386 / 4 = 0.03465.
        d = ROOT / "skills" / "immune-deconvolution-review" / "evals" / "benchmark-pseudobulk-zh" / "files"
        ei, ec, e = mb.read_table(str(d / "estimated_fractions_synthetic.tsv"), "estimated")
        ti, tc, t = mb.read_table(str(d / "truth_fractions_synthetic.tsv"), "truth")
        raw = mb.benchmark(ei, ec, copy.deepcopy(e), ti, tc, copy.deepcopy(t))
        self.assertTrue(any("row sums differ" in w for w in raw["warnings"]))
        res = mb.benchmark(ei, ec, e, ti, tc, t, renorm="truth")
        by = {r["cell_type"]: r for r in res["per_cell_type"]}
        self.assertAlmostEqual(by["B_cell"]["bias"], 0.1030, places=3)
        self.assertAlmostEqual(by["B_cell"]["slope"], 1.384, places=3)
        self.assertGreater(by["B_cell"]["pearson_r"], 0.99)
        self.assertEqual(by["Treg"]["n_true_zero"], 4)
        self.assertAlmostEqual(by["Treg"]["mean_est_when_true_zero"], 0.03465, places=5)
        self.assertFalse(any("row sums differ" in w for w in res["warnings"]))


if __name__ == "__main__":
    unittest.main()
