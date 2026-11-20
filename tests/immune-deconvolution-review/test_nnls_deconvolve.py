"""Unit tests for skills/immune-deconvolution-review/scripts/nnls_deconvolve.py."""
import importlib.util
import itertools
import json
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "immune-deconvolution-review" / "scripts" / "nnls_deconvolve.py"
spec = importlib.util.spec_from_file_location("nnls_deconvolve", SCRIPT)
nd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nd)

# Golden system used below. Signature (linear scale), genes x [T, B]:
#   G1 10 0 | G2 0 10 | G3 5 5 | G4 2 8
# Mixture = 0.3 * T + 0.7 * B = [3, 7, 5, 6.2]; the exact solution is f = (0.3, 0.7).
SIG_IDS = ["G1", "G2", "G3", "G4"]
SIG = [[10.0, 0.0], [0.0, 10.0], [5.0, 5.0], [2.0, 8.0]]
MIX = [[3.0], [7.0], [5.0], [6.2]]


def brute_force_nnls(a_rows, b):
    """Best feasible unconstrained fit over all column subsets (equals the NNLS optimum)."""
    k = len(a_rows[0])
    best = (math.sqrt(sum(v * v for v in b)), [0.0] * k)
    for size in range(1, k + 1):
        for subset in itertools.combinations(range(k), size):
            sub = [[row[j] for j in subset] for row in a_rows]
            g = nd.gram(sub, size)
            h = [sum(r[j] * bi for r, bi in zip(sub, b)) for j in range(size)]
            z = nd.cholesky_solve(g, h)
            if z is None or any(v <= 0 for v in z):
                continue
            x = [0.0] * k
            for j, v in zip(subset, z):
                x[j] = v
            res = math.sqrt(sum((bi - nd.dot(row, x)) ** 2 for row, bi in zip(a_rows, b)))
            if res < best[0] - 1e-12:
                best = (res, x)
    return best


class TestNNLSSolver(unittest.TestCase):
    def test_hand_derived_examples(self):
        # A = [[1,0],[1,0],[0,1]], b = [2,1,1]: x1 is the mean of 2 and 1 = 1.5, x2 = 1;
        # residual (0.5, -0.5, 0) has norm sqrt(0.5).
        x, r = nd.nnls([[1, 0], [1, 0], [0, 1]], [2, 1, 1])
        self.assertAlmostEqual(x[0], 1.5, places=12)
        self.assertAlmostEqual(x[1], 1.0, places=12)
        self.assertAlmostEqual(r, math.sqrt(0.5), places=12)
        # b = [-1,-1,-1]: A'b = (-2, -1) < 0, so x = 0 and the residual norm is sqrt(3).
        x, r = nd.nnls([[1, 0], [1, 0], [0, 1]], [-1, -1, -1])
        self.assertEqual(x, [0.0, 0.0])
        self.assertAlmostEqual(r, math.sqrt(3), places=12)

    def test_active_constraint(self):
        # Fit b = c1 + c2 * t at t = 1, 2, 3 with b = (3, 2, 1). Unconstrained: c = (4, -1).
        # With c2 >= 0 the optimum is c2 = 0, c1 = mean(b) = 2 (KKT: w = (0, -2)).
        x, r = nd.nnls([[1, 1], [1, 2], [1, 3]], [3, 2, 1])
        self.assertAlmostEqual(x[0], 2.0, places=12)
        self.assertEqual(x[1], 0.0)
        self.assertAlmostEqual(r, math.sqrt(2), places=12)

    def test_matches_brute_force_and_kkt(self):
        rng = random.Random(20261120)
        for _ in range(40):
            a = [[rng.uniform(-1, 3) for _ in range(4)] for _ in range(7)]
            b = [rng.uniform(-2, 6) for _ in range(7)]
            x, r = nd.nnls(a, b)
            best_r, best_x = brute_force_nnls(a, b)
            self.assertAlmostEqual(r, best_r, places=9)
            for u, v in zip(x, best_x):
                self.assertAlmostEqual(u, v, places=7)
            resid = [bi - nd.dot(row, x) for row, bi in zip(a, b)]
            w = [sum(row[j] * ri for row, ri in zip(a, resid)) for j in range(4)]
            for xj, wj in zip(x, w):
                self.assertGreaterEqual(xj, 0.0)
                self.assertLessEqual(wj, 1e-8)
                if xj > 0:
                    self.assertAlmostEqual(wj, 0.0, places=8)

    def test_exactly_collinear_columns_terminate(self):
        a = [[1, 1, 0], [2, 2, 1], [0, 0, 3], [1, 1, 1]]
        b = [0.4 * r[0] + 0.6 * r[2] for r in a]
        x, r = nd.nnls(a, b)
        self.assertAlmostEqual(x[0] + x[1], 0.4, places=10)
        self.assertAlmostEqual(x[2], 0.6, places=10)
        self.assertAlmostEqual(r, 0.0, places=10)


class TestLinearAlgebraAndScale(unittest.TestCase):
    def test_condition_number_golden(self):
        # [[1,1],[0,1]]'[[1,1],[0,1]] = [[1,1],[1,2]], eigenvalues (3 -/+ sqrt 5)/2, so the
        # 2-norm condition number is sqrt(lmax/lmin) = (3 + sqrt 5)/2 = 2.6180339887...
        eig = nd.symmetric_eigenvalues([[1.0, 1.0], [1.0, 2.0]])
        self.assertAlmostEqual(eig[0], (3 - math.sqrt(5)) / 2, places=12)
        self.assertAlmostEqual(eig[1], (3 + math.sqrt(5)) / 2, places=12)
        self.assertAlmostEqual(nd.condition_number([[1.0, 1.0], [0.0, 1.0]], 2), (3 + math.sqrt(5)) / 2, places=10)
        self.assertAlmostEqual(nd.condition_number([[1.0, 0.0], [0.0, 2.0]], 2), 2.0, places=12)
        self.assertIsNone(nd.condition_number([[1.0, 1.0], [2.0, 2.0]], 2))

    def test_geo2r_rule(self):
        self.assertTrue(nd.looks_linear([[0.0, 150.0]]))                    # 99th pct > 100
        self.assertTrue(nd.looks_linear([[float(v) for v in range(1, 61)]]))  # range > 50, q25 > 0
        self.assertFalse(nd.looks_linear([[0.0] * 10 + [60.0]]))           # q25 = 0
        self.assertFalse(nd.looks_linear([[2.5, 7.1, 12.0, 15.9]]))         # typical log2 values

    def test_scale_handling(self):
        log_mix = [[math.log2(v[0] + 1.0)] for v in MIX]
        lin, clipped, used = nd.to_linear(log_mix, "log2", "mixture")
        self.assertEqual((clipped, used), (0, "log2"))
        for got, want in zip(lin, MIX):
            self.assertAlmostEqual(got[0], want[0], places=12)
        with self.assertRaisesRegex(nd.InputError, "look log-transformed"):
            nd.to_linear(log_mix, "auto", "mixture")
        with self.assertRaisesRegex(nd.InputError, "negative"):
            nd.to_linear([[1.0, -0.5]], "linear", "mixture")
        lin, clipped, _ = nd.to_linear([[-1.0, 0.0]], "log2", "mixture")
        self.assertEqual((lin, clipped), ([[0.0, 0.0]], 1))
        with self.assertRaisesRegex(nd.InputError, "too large"):
            nd.to_linear([[5000.0]], "log2", "mixture")

    def test_pearson_edge_cases(self):
        self.assertIsNone(nd.pearson([1.0], [2.0]))
        self.assertIsNone(nd.pearson([0.1, 0.1, 0.1], [1.0, 2.0, 3.0]))
        # sxy = 4.5, sxx = 2, syy = 366/36, so r = 4.5 / sqrt(61/3) = 0.9979487...
        self.assertAlmostEqual(nd.pearson([1.0, 2.0, 3.0], [2.0, 4.0, 6.5]), 0.997949, places=6)


class TestDeconvolve(unittest.TestCase):
    def run_golden(self, row_scale):
        res, genes, resid = nd.deconvolve(SIG_IDS, ["T", "B"], SIG, SIG_IDS, ["S1"], MIX, row_scale=row_scale)
        rec = res["samples"][0]
        self.assertAlmostEqual(rec["fractions"]["T"], 0.3, places=6)
        self.assertAlmostEqual(rec["fractions"]["B"], 0.7, places=6)
        self.assertAlmostEqual(rec["coef_sum"], 1.0, places=6)
        self.assertAlmostEqual(rec["rmse"], 0.0, places=6)
        self.assertEqual(rec["pearson_r"], 1.0)
        self.assertEqual(genes, SIG_IDS)
        return res

    def test_exact_recovery_with_and_without_row_scaling(self):
        self.run_golden(True)
        res = self.run_golden(False)
        self.assertEqual(res["n_genes_used"], 4)
        self.assertEqual(res["collinear_pairs"], [])

    def test_coefficients_keep_the_scale(self):
        mix = [[2 * v[0]] for v in MIX]  # twice the material: same fractions, coef_sum 2
        res, _, _ = nd.deconvolve(SIG_IDS, ["T", "B"], SIG, SIG_IDS, ["S1"], mix)
        self.assertAlmostEqual(res["samples"][0]["coef_sum"], 2.0, places=6)
        self.assertAlmostEqual(res["samples"][0]["fractions"]["T"], 0.3, places=6)

    def test_duplicate_signature_columns_are_flagged(self):
        sig = [[r[0], r[0], r[1]] for r in SIG]
        res, _, _ = nd.deconvolve(SIG_IDS, ["T", "T2", "B"], sig, SIG_IDS, ["S1"], MIX, min_genes_factor=1)
        rec = res["samples"][0]
        self.assertAlmostEqual(rec["fractions"]["T"] + rec["fractions"]["T2"], 0.3, places=6)
        self.assertIsNone(res["condition_number_fit"])
        self.assertEqual(res["collinear_pairs"][0]["r"], 1.0)
        self.assertTrue(any("numerically singular" in w for w in res["warnings"]))

    def test_marker_coverage_missing_and_zero_rows(self):
        sig_ids = SIG_IDS + ["G5", "G6"]
        sig = SIG + [[1.0, 9.0], [0.0, 0.0]]
        mix_ids = SIG_IDS + ["G6"]
        res, _, _ = nd.deconvolve(sig_ids, ["T", "B"], sig, mix_ids, ["S1"], MIX + [[0.0]])
        # G1 -> T, G3 (tie) -> T (first maximum), G2, G4, G5 -> B; G5 is missing; G6 is all zero.
        self.assertEqual(res["marker_coverage"]["T"], {"n_marker_genes": 2, "n_found": 2})
        self.assertEqual(res["marker_coverage"]["B"], {"n_marker_genes": 3, "n_found": 2})
        self.assertEqual(res["missing_genes"], ["G5"])
        self.assertEqual(res["n_zero_signature_rows_dropped"], 1)
        self.assertEqual(res["n_genes_used"], 4)

    def test_unexplained_gene_tops_the_residual_list(self):
        sig_ids = SIG_IDS + ["TUMOUR1"]
        sig = SIG + [[0.1, 0.1]]
        mix = MIX + [[50.0]]
        res, _, _ = nd.deconvolve(sig_ids, ["T", "B"], sig, sig_ids, ["S1"], mix, row_scale=False)
        top = res["top_residual_genes"][0]
        self.assertEqual(top["gene"], "TUMOUR1")
        self.assertGreater(top["mean_residual"], 0)

    def test_errors_and_zero_sample(self):
        with self.assertRaisesRegex(nd.InputError, "need at least 4"):
            nd.deconvolve(SIG_IDS, ["T", "B"], SIG, SIG_IDS[:3], ["S1"], MIX[:3])
        res, _, _ = nd.deconvolve(SIG_IDS, ["T", "B"], SIG, SIG_IDS, ["S1", "Z"], [r + [0.0] for r in MIX])
        self.assertIsNone(res["samples"][1]["fractions"])
        self.assertEqual(res["samples"][1]["coef_sum"], 0.0)
        self.assertTrue(any("every coefficient is zero" in w for w in res["warnings"]))


class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.sig = d / "sig.tsv"
        self.sig.write_text("gene\tT\tB\n" + "".join(f"{g}\t{r[0]}\t{r[1]}\n" for g, r in zip(SIG_IDS, SIG)))
        self.mix = d / "mix.tsv"
        self.mix.write_text("gene\tS1\n" + "".join(f"{g}\t{r[0]}\n" for g, r in zip(SIG_IDS, MIX)))
        self.dir = d

    def tearDown(self):
        self.tmp.cleanup()

    def run_cli(self, *extra):
        cmd = [sys.executable, str(SCRIPT), "--signature", str(self.sig), "--mixture", str(self.mix), *extra]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60)

    def test_json_tsv_and_residuals(self):
        proc = self.run_cli("--signature-scale", "linear", "--mixture-scale", "linear")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out = json.loads(proc.stdout)
        self.assertAlmostEqual(out["samples"][0]["fractions"]["B"], 0.7, places=6)
        self.assertEqual(out["settings"]["mixture_scale"], "linear")
        resid = self.dir / "resid.tsv"
        proc = self.run_cli("--signature-scale", "linear", "--mixture-scale", "linear", "--format", "tsv",
                            "--residuals-out", str(resid))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        header, row = proc.stdout.strip().splitlines()
        self.assertEqual(header.split("\t")[:4], ["sample", "T", "B", "coef_sum"])
        self.assertAlmostEqual(float(row.split("\t")[1]), 0.3, places=6)
        self.assertEqual(resid.read_text().splitlines()[0], "gene\tS1")

    def test_auto_scale_refuses_small_values_and_help_works(self):
        proc = self.run_cli()
        self.assertEqual(proc.returncode, 2)
        self.assertIn("look log-transformed", proc.stderr)
        proc = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True, timeout=60)
        self.assertEqual(proc.returncode, 0)
        self.assertIn("usage", proc.stdout.lower())

    def test_bad_value_exits_2(self):
        self.mix.write_text("gene\tS1\nG1\tNA\nG2\t1\nG3\t1\nG4\t1\n")
        proc = self.run_cli("--signature-scale", "linear", "--mixture-scale", "linear")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("'G1'", proc.stderr)

    def test_other_errors_exit_2(self):
        proc = self.run_cli("--signature-scale", "linear", "--mixture-scale", "linear",
                            "--out", str(self.dir / "no_such_dir" / "out.json"))
        self.assertEqual(proc.returncode, 2)
        self.assertIn("cannot write output", proc.stderr)
        proc = self.run_cli("--collinear-r", "1.5")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--collinear-r", proc.stderr)
        self.mix.write_text("gene\tS1\nG1\t3\ng1\t3\nG2\t7\nG3\t5\nG4\t6.2\n")
        proc = self.run_cli("--signature-scale", "linear", "--mixture-scale", "linear")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("appears twice", proc.stderr)


class TestReviewerGoldens(unittest.TestCase):
    """Golden values derived by hand during review, independent of the solver code."""

    def test_active_constraint_on_second_column(self):
        # A = [[1,0],[0,1],[1,1]], b = (3,-1,2). Unconstrained: A'A = [[2,1],[1,2]], A'b = (5,1),
        # x = (1/3)[[2,-1],[-1,2]](5,1) = (3,-1), infeasible. With x2 = 0: x1 = 5/2 = 2.5.
        # KKT: w2 = a2'(b - 2.5 a1) = (0,1,1).(0.5,-1,-0.5) = -1.5 <= 0, so (2.5, 0) is optimal;
        # the residual (0.5,-1,-0.5) has norm sqrt(1.5).
        x, r = nd.nnls([[1, 0], [0, 1], [1, 1]], [3, -1, 2])
        self.assertAlmostEqual(x[0], 2.5, places=12)
        self.assertEqual(x[1], 0.0)
        self.assertAlmostEqual(r, math.sqrt(1.5), places=12)

    def test_condition_number_of_a_tall_matrix(self):
        # S = [[2,0],[0,1],[1,1]]: S'S = [[5,1],[1,2]], eigenvalues (7 +/- sqrt 13)/2, so the
        # 2-norm condition number is sqrt((7 + sqrt 13)/(7 - sqrt 13)) = 1.7675918...
        want = math.sqrt((7 + math.sqrt(13)) / (7 - math.sqrt(13)))
        self.assertAlmostEqual(nd.condition_number([[2.0, 0.0], [0.0, 1.0], [1.0, 1.0]], 2), want, places=10)

    def test_matches_that_differ_only_in_case_are_flagged(self):
        res, _, _ = nd.deconvolve(SIG_IDS, ["T", "B"], SIG, [g.lower() for g in SIG_IDS], ["S1"], MIX)
        self.assertEqual(res["n_case_mismatched_genes"], 4)
        self.assertTrue(any("ignoring letter case" in w for w in res["warnings"]))
        self.assertAlmostEqual(res["samples"][0]["fractions"]["T"], 0.3, places=6)
        res, _, _ = nd.deconvolve(SIG_IDS, ["T", "B"], SIG, SIG_IDS, ["S1"], MIX)
        self.assertEqual(res["n_case_mismatched_genes"], 0)


class TestSyntheticEvalCase(unittest.TestCase):
    """Pins what the nnls-macrophage-split graders rely on. The same fractions were recomputed
    in review by exhaustive search over column subsets with Gauss-Jordan solves (max difference
    5e-7), and the M1/M2 column correlation by a plain Pearson formula on the row-scaled matrix."""

    def test_macrophage_case(self):
        d = ROOT / "skills" / "immune-deconvolution-review" / "evals" / "nnls-macrophage-split" / "files"
        sig_ids, cts, sig = nd.read_matrix(str(d / "signature_synthetic.tsv"), "signature")
        mix_ids, samples, mix = nd.read_matrix(str(d / "expression_log2tpm_synthetic.tsv"), "mixture")
        with self.assertRaisesRegex(nd.InputError, "look log-transformed"):
            nd.to_linear(mix, "auto", "mixture")
        lin, _, _ = nd.to_linear(mix, "log2", "mixture")
        res, _, _ = nd.deconvolve(sig_ids, cts, sig, mix_ids, samples, lin)
        self.assertEqual([(p["cell_type_a"], p["cell_type_b"]) for p in res["collinear_pairs"]],
                         [("Macrophage_M1", "Macrophage_M2")])
        self.assertAlmostEqual(res["collinear_pairs"][0]["r"], 0.9027, places=4)
        self.assertEqual([g["gene"] for g in res["top_residual_genes"][:2]], ["IDO1", "CXCL9"])
        self.assertTrue(all(g["mean_residual"] > 0 for g in res["top_residual_genes"][:2]))
        by = {r["sample"]: r for r in res["samples"]}
        self.assertTrue(all(by[s]["coefficients"]["Macrophage_M2"] == 0.0 for s in samples))
        self.assertLess(max(by[s]["coef_sum"] for s in ("N1", "N2", "N3")),
                        min(by[s]["coef_sum"] for s in ("R1", "R2", "R3")))
        self.assertAlmostEqual(by["R1"]["coef_sum"], 0.635, places=3)
        self.assertAlmostEqual(by["R1"]["fractions"]["Macrophage_M1"], 0.3989, places=4)


if __name__ == "__main__":
    unittest.main()
