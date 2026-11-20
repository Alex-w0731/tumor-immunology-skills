"""Unit tests for skills/spatial-neighborhood-statistics/scripts/neighborhood_enrichment.py."""
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
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "spatial-neighborhood-statistics" / "scripts" / "neighborhood_enrichment.py"
_spec = importlib.util.spec_from_file_location("neighborhood_enrichment_under_test", SCRIPT)
ne = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ne)


def brute_knn(xs, ys, k):
    out = []
    for i in range(len(xs)):
        d = sorted(((xs[j] - xs[i]) ** 2 + (ys[j] - ys[i]) ** 2, j) for j in range(len(xs)) if j != i)
        out.append([j for _, j in d[:k]])
    return out


def brute_radius(xs, ys, radius):
    # same rule as the script: distance <= radius with a 1e-9 relative tolerance for rounding
    return [[j for j in range(len(xs)) if j != i and math.hypot(xs[j] - xs[i], ys[j] - ys[i]) <= radius * (1 + 1e-9)]
            for i in range(len(xs))]


def directed_count(labels, nbrs, focal, a, b):
    return sum(1 for i in range(len(labels)) if focal[i] and labels[i] == a for j in nbrs[i] if labels[j] == b)


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def table(rows, header=("x", "y", "cell_type", "sample", "patient")):
    return "\t".join(header) + "\n" + "".join("\t".join(str(v) for v in r) + "\n" for r in rows)


def compartment_tissue(seed=3, nx=20, ny=10):
    # synthetic: left half mostly Tumor, right half mostly Immune; labels random within each half
    rng = random.Random(seed)
    xs, ys, types, comp = [], [], [], []
    for gx in range(nx):
        for gy in range(ny):
            xs.append(gx * 10 + rng.uniform(-2, 2))
            ys.append(gy * 10 + rng.uniform(-2, 2))
            left = gx < nx // 2
            comp.append("left" if left else "right")
            p_tumor = 0.9 if left else 0.1
            types.append("Tumor" if rng.random() < p_tumor else "Immune")
    return xs, ys, types, comp


class GoldenPath(unittest.TestCase):
    # Four cells on a line, 1 unit apart, labelled A B A B; radius 1 links neighbours on the line.
    # Undirected edges (0,1), (1,2), (2,3) -> 6 directed edges, every cell focal, so pairs are unordered.
    # Observed C(A,B): cell 0 has one B neighbour, cell 2 has two -> 3.
    # Expected under label permutation: 6 edges x n_A n_B / (N (N - 1)) = 6 x 2 x 2 / 12 = 2;
    #   C(A,A) and C(B,B): 6 x 2 x 1 / 12 = 1 each.
    # Exact null over the 6 distinct labellings (A positions {0,1},{0,2},{0,3},{1,2},{1,3},{2,3}):
    #   C(A,B) = 1, 3, 2, 2, 3, 1 -> mean 2 and P(C >= 3) = 2/6 = 1/3.
    xs, ys, types = [0.0, 1.0, 2.0, 3.0], [0.0] * 4, ["A", "B", "A", "B"]

    def test_counts_expectation_and_tail(self):
        res = ne.analyze_sample(self.xs, self.ys, self.types, graph="radius", radius=1.0, min_cells=1,
                                n_perm=20000, seed=7)
        rows = {(r["type_a"], r["type_b"]): r for r in res["rows"]}
        self.assertEqual(set(rows), {("A", "A"), ("A", "B"), ("B", "B")})
        ab = rows[("A", "B")]
        self.assertEqual(ab["observed"], 3)
        self.assertAlmostEqual(ab["expected"], 2.0, places=12)
        self.assertAlmostEqual(rows[("A", "A")]["expected"], 1.0, places=12)
        self.assertEqual(rows[("A", "A")]["observed"], 0)
        self.assertAlmostEqual(ab["p_enrich"], 1 / 3, delta=0.015)
        self.assertAlmostEqual(ab["log2_ratio"], math.log2(3.5 / 2.5), places=12)
        self.assertFalse(res["info"]["pairs_ordered"])


class GoldenGrid(unittest.TestCase):
    # 3 x 3 unit lattice (index = 3 * row + column), labels   A B A / B B B / A B C.
    # Radius 1.0 sits exactly on the lattice spacing: the 12 horizontal and vertical edges count
    # (24 directed), the diagonals (1.41) do not. Observed C(A,B): each corner A has 2 B
    # neighbours (corner 0: cells 1 and 3; corner 2: 1 and 5; corner 6: 3 and 7) -> 6.
    # Expected: 24 x n_A n_B / (N (N - 1)) = 24 x 3 x 5 / 72 = 5; C(A,A): 24 x 3 x 2 / 72 = 2.
    # The exact tails are enumerated below over the 9! / (3! 5! 1!) = 504 distinct labellings.
    xs = [float(i % 3) for i in range(9)]
    ys = [float(i // 3) for i in range(9)]
    labels = ["A", "B", "A", "B", "B", "B", "A", "B", "C"]

    def test_counts_expectation_and_exact_tails(self):
        edges = [(i, j) for i in range(9) for j in range(9)
                 if i != j and abs(self.xs[i] - self.xs[j]) + abs(self.ys[i] - self.ys[j]) == 1]
        self.assertEqual(len(edges), 24)
        labellings = set(itertools.permutations(self.labels))
        self.assertEqual(len(labellings), 504)
        counts = [sum(1 for i, j in edges if lab[i] == "A" and lab[j] == "B") for lab in labellings]
        tail_ge = sum(c >= 6 for c in counts) / 504
        tail_le = sum(c <= 6 for c in counts) / 504
        self.assertAlmostEqual(sum(counts) / 504, 5.0, places=12)
        res = ne.analyze_sample(self.xs, self.ys, self.labels, graph="radius", radius=1.0, min_cells=1,
                                n_perm=20000, seed=5, pairs=[("A", "B"), ("A", "A")])
        rows = {(r["type_a"], r["type_b"]): r for r in res["rows"]}
        ab = rows[("A", "B")]
        self.assertEqual(ab["observed"], 6)
        self.assertAlmostEqual(ab["expected"], 5.0, places=12)
        self.assertAlmostEqual(rows[("A", "A")]["expected"], 2.0, places=12)
        self.assertAlmostEqual(ab["log2_ratio"], math.log2(6.5 / 5.5), places=12)
        self.assertAlmostEqual(ab["p_enrich"], tail_ge, delta=0.012)     # about 4 Monte Carlo SE
        self.assertAlmostEqual(ab["p_deplete"], tail_le, delta=0.012)

    def test_knn_ties_go_to_lower_index(self):
        # the centre cell 4 has four neighbours at distance 1 (cells 1, 3, 5, 7): k = 2 keeps 1 and 3
        nbrs = ne.knn_lists(self.xs, self.ys, 2)
        self.assertEqual(nbrs[4], [1, 3])
        self.assertEqual(nbrs[8], [5, 7])


class ScaleDependence(unittest.TestCase):
    def test_z_grows_with_tissue_size_log2_ratio_does_not(self):
        # two far-apart copies of one tissue in one sample: observed and expected double, the SD
        # grows by about sqrt(2), so z grows by about sqrt(2) while log2(O/E) barely moves
        rng = random.Random(12)
        xs = [rng.uniform(0, 200) for _ in range(150)]
        ys = [rng.uniform(0, 200) for _ in range(150)]
        types = ["A" if x < 100 and rng.random() < 0.6 else rng.choice(["A", "B", "C"]) for x in xs]
        one = ne.analyze_sample(xs, ys, types, graph="radius", radius=30, pairs=[("A", "A")],
                                n_perm=2000, seed=1)["rows"][0]
        two = ne.analyze_sample(xs + [x + 10000 for x in xs], ys + ys, types + types, graph="radius",
                                radius=30, pairs=[("A", "A")], n_perm=2000, seed=1)["rows"][0]
        self.assertEqual(two["observed"], 2 * one["observed"])
        self.assertGreater(one["z"], 2)
        self.assertTrue(1.2 < two["z"] / one["z"] < 1.65, (one["z"], two["z"]))
        self.assertLess(abs(two["log2_ratio"] - one["log2_ratio"]), 0.05)


class ExactExpectation(unittest.TestCase):
    def test_stratified_expectation_matches_enumeration(self):
        # 7 cells in two strata (4 + 3), a focal subset, directed kNN graph; enumerate every
        # within-stratum permutation (4! x 3! = 144) and average the directed counts.
        rng = random.Random(11)
        xs = [rng.uniform(0, 10) for _ in range(7)]
        ys = [rng.uniform(0, 10) for _ in range(7)]
        labels = [0, 1, 1, 2, 0, 2, 1]
        strata = ["s", "s", "s", "s", "t", "t", "t"]
        focal = [True, True, False, True, True, False, True]
        nbrs, _, _ = ne.build_graph(xs, ys, "knn", k=2, knn_mode="directed")
        n_types = 3
        codes = [a * n_types + b for a in range(3) for b in range(3)]
        stats = ne.permutation_counts(labels, nbrs, focal, strata, codes, n_types, 1, random.Random(1))
        idx_s, idx_t = [0, 1, 2, 3], [4, 5, 6]
        totals = dict.fromkeys(codes, 0)
        n_perm = 0
        for ps in itertools.permutations([labels[i] for i in idx_s]):
            for pt in itertools.permutations([labels[i] for i in idx_t]):
                lab = list(labels)
                for i, v in zip(idx_s, ps):
                    lab[i] = v
                for i, v in zip(idx_t, pt):
                    lab[i] = v
                n_perm += 1
                for c in codes:
                    totals[c] += directed_count(lab, nbrs, focal, *divmod(c, n_types))
        self.assertEqual(n_perm, 144)
        for c in codes:
            self.assertAlmostEqual(stats[c][1], totals[c] / n_perm, places=10, msg=f"code {c}")
            self.assertEqual(stats[c][0], directed_count(labels, nbrs, focal, *divmod(c, n_types)))

    def test_permutation_mean_close_to_exact(self):
        xs, ys, types, _ = compartment_tissue()
        res = ne.analyze_sample(xs, ys, types, graph="knn", k=6, n_perm=400, seed=2)
        labels = sorted(set(types))
        index = {t: i for i, t in enumerate(labels)}
        nbrs, _, _ = ne.build_graph(xs, ys, "knn", 6, "union")
        codes = [index["Immune"] * 2 + index["Tumor"]]
        stats = ne.permutation_counts([index[t] for t in types], nbrs, [True] * len(xs), [""] * len(xs),
                                      codes, 2, 400, random.Random(5))
        obs, exp, sd, _, _, perm_mean = stats[codes[0]]
        self.assertLess(abs(perm_mean - exp), 4 * sd / math.sqrt(400))
        self.assertTrue(res["rows"])


class Graphs(unittest.TestCase):
    def test_knn_matches_brute_force_with_ties_and_random_points(self):
        grid_x = [float(i % 7) for i in range(49)]
        grid_y = [float(i // 7) for i in range(49)]
        for k in (1, 3, 6, 10):
            self.assertEqual(ne.knn_lists(grid_x, grid_y, k), brute_knn(grid_x, grid_y, k))
        rng = random.Random(4)
        xs = [rng.uniform(0, 100) for _ in range(150)]
        ys = [rng.uniform(0, 30) for _ in range(150)]
        for k in (1, 6, 15):
            self.assertEqual(ne.knn_lists(xs, ys, k), brute_knn(xs, ys, k))
        clustered_x = [rng.gauss(0, 1) for _ in range(60)] + [rng.gauss(500, 1) for _ in range(60)]
        clustered_y = [rng.gauss(0, 1) for _ in range(120)]
        self.assertEqual(ne.knn_lists(clustered_x, clustered_y, 5), brute_knn(clustered_x, clustered_y, 5))

    def test_knn_degenerate_inputs(self):
        self.assertEqual(ne.knn_lists([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 6), [[1, 2], [0, 2], [0, 1]])
        self.assertEqual(ne.knn_lists([1.0], [1.0], 3), [[]])
        line = [float(i) for i in range(20)]
        self.assertEqual(ne.knn_lists(line, [0.0] * 20, 2), brute_knn(line, [0.0] * 20, 2))

    def test_radius_includes_exact_boundary_after_scaling(self):
        # integer grid scaled by 0.1: neighbours exactly 0.1 apart must be found despite rounding
        xs = [(i % 6) * 0.1 for i in range(36)]
        ys = [(i // 6) * 0.1 for i in range(36)]
        got = ne.radius_lists(xs, ys, 0.1)
        self.assertEqual(got, brute_radius(xs, ys, 0.1))
        self.assertEqual(len(got[14]), 4)

    def test_union_is_symmetric_directed_is_not_and_cap_prunes(self):
        rng = random.Random(8)
        xs = [rng.uniform(0, 50) for _ in range(40)] + [500.0]
        ys = [rng.uniform(0, 50) for _ in range(40)] + [500.0]
        union, sym, _ = ne.build_graph(xs, ys, "knn", k=3, knn_mode="union")
        self.assertTrue(sym)
        for i, lst in enumerate(union):
            for j in lst:
                self.assertIn(i, union[j])
        directed, sym_d, _ = ne.build_graph(xs, ys, "knn", k=3, knn_mode="directed")
        self.assertFalse(sym_d)
        self.assertTrue(all(len(lst) == 3 for lst in directed))
        capped, _, pruned = ne.build_graph(xs, ys, "knn", k=3, knn_mode="union", max_dist=100.0)
        self.assertEqual(capped[40], [])
        self.assertGreater(pruned, 0)

    def test_focal_by_buffer(self):
        xs = [0.0, 10.0, 0.0, 10.0, 5.0]
        ys = [0.0, 0.0, 10.0, 10.0, 5.0]
        self.assertEqual(ne.focal_by_buffer(xs, ys, 3.0, "bbox"), [False] * 4 + [True])
        self.assertEqual(ne.focal_by_buffer(xs, ys, 3.0, "hull"), [False] * 4 + [True])
        self.assertEqual(ne.focal_by_buffer(xs, ys, 6.0, "hull"), [False] * 5)


class Statistics(unittest.TestCase):
    def test_bh_and_by_golden(self):
        # sorted p 0.01, 0.03, 0.04, 0.2 (m = 4): BH = min over higher ranks of p m / rank
        #   0.2 -> 0.2; 0.04 -> 0.0533; 0.03 -> min(0.06, 0.0533); 0.01 -> 0.04
        q = ne.adjust_p([0.01, 0.04, 0.03, 0.2], "bh")
        for got, want in zip(q, [0.04, 0.16 / 3, 0.16 / 3, 0.2]):
            self.assertAlmostEqual(got, want, places=12)
        c = 1 + 1 / 2 + 1 / 3 + 1 / 4
        q_by = ne.adjust_p([0.01, 0.04, 0.03, 0.2], "by")
        self.assertAlmostEqual(q_by[0], 0.04 * c, places=12)
        self.assertEqual(ne.adjust_p([], "bh"), [])

    def test_sign_test_golden_and_floor(self):
        # 6 of 6 patients in one direction: 2 x 0.5^6 = 0.03125; 5 of 5: 2 x 0.5^5 = 0.0625
        self.assertAlmostEqual(ne.sign_test(6, 0), 0.03125, places=12)
        self.assertAlmostEqual(ne.sign_test(0, 5), 0.0625, places=12)
        self.assertEqual(ne.sign_test(3, 3), 1.0)
        self.assertEqual(ne.sign_test(0, 0), 1.0)
        # 8 vs 2: 2 x (1 + 10 + 45) / 1024
        self.assertAlmostEqual(ne.sign_test(8, 2), 2 * 56 / 1024, places=12)
        self.assertAlmostEqual(ne.sign_test_floor(5), 0.0625, places=12)
        self.assertEqual(ne.sign_test_floor(0), 1.0)

    def test_summary_uses_one_effect_per_patient(self):
        rows = []
        effects = {("s1", "P1"): 1.0, ("s2", "P1"): 3.0, ("s3", "P2"): -1.0, ("s4", "P3"): 0.5}
        expected = {"s1": 2.0, "s2": 8.0, "s3": 4.0, "s4": 30.0}
        for (s, p), v in effects.items():
            rows.append({"sample": s, "patient": p, "type_a": "A", "type_b": "B", "log2_ratio": v,
                         "expected": expected[s], "call": "ns"})
        out = ne.summarize_across(rows)
        self.assertEqual(len(out), 1)
        o = out[0]
        self.assertEqual(o["n_samples"], 4)
        self.assertEqual(o["n_units"], 3)
        self.assertEqual(o["median_expected"], 6.0)      # median of 2, 4, 8, 30 over samples
        self.assertEqual(o["unit_effects"], {"P1": 2.0, "P2": -1.0, "P3": 0.5})
        self.assertEqual((o["n_pos"], o["n_neg"]), (2, 1))
        self.assertEqual(o["median_log2_ratio"], 0.5)
        self.assertEqual(o["sign_test_p"], 1.0)


class Behaviour(unittest.TestCase):
    def test_compartments_create_segregation_that_strata_remove(self):
        xs, ys, types, comp = compartment_tissue()
        plain = ne.analyze_sample(xs, ys, types, graph="knn", k=6, n_perm=300, seed=1)
        row = [r for r in plain["rows"] if (r["type_a"], r["type_b"]) == ("Immune", "Tumor")][0]
        self.assertEqual(row["call"], "depleted")
        self.assertLess(row["z"], -4)
        strat = ne.analyze_sample(xs, ys, types, graph="knn", k=6, n_perm=300, seed=1, strata=comp)
        srow = [r for r in strat["rows"] if (r["type_a"], r["type_b"]) == ("Immune", "Tumor")][0]
        self.assertEqual(srow["call"], "ns")
        self.assertLess(abs(srow["log2_ratio"]), 0.3)
        self.assertGreater(abs(row["log2_ratio"]), 1.0)

    def test_random_labels_give_no_calls(self):
        rng = random.Random(21)
        xs = [rng.uniform(0, 300) for _ in range(300)]
        ys = [rng.uniform(0, 300) for _ in range(300)]
        types = [rng.choice(["T", "B", "M"]) for _ in range(300)]
        res = ne.analyze_sample(xs, ys, types, graph="radius", radius=25, n_perm=300, seed=3)
        self.assertEqual(len(res["rows"]), 6)
        self.assertTrue(all(r["call"] == "ns" for r in res["rows"]))

    def test_edge_cases(self):
        tiny = ne.analyze_sample([0.0, 1.0], [0.0, 0.0], ["A", "B"])
        self.assertTrue(tiny["info"]["skipped"])
        one_type = ne.analyze_sample([float(i) for i in range(12)], [0.0] * 12, ["A"] * 12, n_perm=10)
        self.assertEqual(len(one_type["rows"]), 1)       # A:A only, labels cannot move
        self.assertEqual(one_type["rows"][0]["p_two_sided"], 1.0)
        few = ne.analyze_sample([float(i) for i in range(12)], [0.0] * 12, ["A"] * 6 + ["B"] * 6,
                                n_perm=10, min_cells=7)
        self.assertEqual(few["rows"], [])
        self.assertTrue(any("no testable pair" in n for n in few["notes"]))
        with self.assertRaises(ne.InputError):
            ne.analyze_sample([0.0] * 5, [0.0] * 5, ["A", "B", "A", "B", "A"], max_dist_factor=1.5)

    def test_force_ordered_and_seed_determinism(self):
        xs, ys, types, _ = compartment_tissue(nx=8, ny=8)
        a = ne.analyze_sample(xs, ys, types, n_perm=50, seed=9, force_ordered=True)
        pairs = {(r["type_a"], r["type_b"]) for r in a["rows"]}
        self.assertIn(("Tumor", "Immune"), pairs)
        self.assertIn(("Immune", "Tumor"), pairs)
        b = ne.analyze_sample(xs, ys, types, n_perm=50, seed=9, force_ordered=True)
        self.assertEqual(a["rows"], b["rows"])


class CommandLine(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True, timeout=120)

    def test_smoke_two_patients_and_outputs(self):
        xs, ys, types, _ = compartment_tissue(nx=10, ny=6)
        rows = []
        for s, p in (("s1", "P1"), ("s2", "P1"), ("s3", "P2")):
            rows += [(round(x * 2, 3), round(y * 2, 3), t, s, p) for x, y, t in zip(xs, ys, types)]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cells.tsv")
            write(path, table(rows))
            tsv, summ = os.path.join(tmp, "pairs.tsv"), os.path.join(tmp, "summary.tsv")
            proc = self.run_cli(path, "--sample-col", "sample", "--patient-col", "patient", "--um-per-unit",
                                "0.5", "--k", "6", "--max-dist-factor", "1.8", "--n-perm", "99", "--seed", "4",
                                "--tsv-out", tsv, "--summary-tsv", summ)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            out = json.loads(proc.stdout)
            self.assertEqual(out["distance_unit"], "um")
            self.assertEqual(len(out["samples"]), 3)
            summary = {(r["type_a"], r["type_b"]): r for r in out["across_samples"]}
            self.assertEqual(summary[("Immune", "Tumor")]["n_units"], 2)
            self.assertEqual(summary[("Immune", "Tumor")]["n_samples"], 3)
            self.assertTrue(any("smallest attainable" in n for n in out["notes"]))
            with open(tsv, encoding="utf-8") as fh:
                header = fh.readline().rstrip("\n").split("\t")
            self.assertEqual(header, ne.PAIR_COLUMNS)
            self.assertTrue(os.path.getsize(summ) > 0)
            # micrometre scaling: median nearest-neighbour distance is reported in um (half the pixels)
            med = out["samples"][0]["graph"]["median_nn_distance"]
            self.assertTrue(0 < med < 12)

    def test_focal_column_makes_all_samples_ordered(self):
        rows = []
        rng = random.Random(2)
        for s in ("a", "b"):
            for i in range(40):
                focal = "1" if (s == "a" or i % 2) else "0"     # sample a: every cell focal
                rows.append((rng.uniform(0, 100), rng.uniform(0, 100), rng.choice(["X", "Y"]), s, focal))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cells.tsv")
            write(path, table(rows, ("x", "y", "cell_type", "sample", "focal")))
            proc = self.run_cli(path, "--sample-col", "sample", "--focal-col", "focal", "--n-perm", "20")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            out = json.loads(proc.stdout)
            self.assertTrue(all(s["pairs_ordered"] for s in out["samples"]))

    def test_errors_exit_2_with_message(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = os.path.join(tmp, "empty.tsv")
            write(empty, "x\ty\tcell_type\n")
            proc = self.run_cli(empty)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("no data rows", proc.stderr)
            bad = os.path.join(tmp, "bad.tsv")
            write(bad, table([(1, "a", "T", "s1", "P1")]))
            proc = self.run_cli(bad)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("coordinates must be numbers", proc.stderr)
            proc = self.run_cli(bad, "--type-col", "phenotype")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("'phenotype'", proc.stderr)
            two = os.path.join(tmp, "two.tsv")
            write(two, table([(1, 1, "T", "s1", "P1"), (2, 2, "B", "s1", "P2"), (3, 3, "T", "s1", "P1")]))
            proc = self.run_cli(two, "--sample-col", "sample", "--patient-col", "patient")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("maps to 2 patients", proc.stderr)
            proc = self.run_cli(two, "--pairs", "T:NK")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("NK", proc.stderr)
            proc = self.run_cli(two, "--graph", "radius")
            self.assertEqual(proc.returncode, 2)
            # a radius typed without --graph radius must not silently run the default kNN graph
            proc = self.run_cli(two, "--radius", "25")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--graph radius", proc.stderr)
            proc = self.run_cli(two, "--graph", "radius", "--radius", "25", "--k", "6")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("--k", proc.stderr)
            proc = self.run_cli(two, "--min-cells", "1", "--n-perm", "5",
                                "--out", os.path.join(tmp, "no_such_dir", "out.json"))
            self.assertEqual(proc.returncode, 2)
            self.assertIn("cannot write", proc.stderr)


if __name__ == "__main__":
    unittest.main()
