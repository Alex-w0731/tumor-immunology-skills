"""Unit tests for skills/tcr-repertoire-metrics/scripts/repertoire_overlap.py."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SCRIPT = os.path.join(ROOT, "skills", "tcr-repertoire-metrics", "scripts", "repertoire_overlap.py")

_spec = importlib.util.spec_from_file_location("repertoire_overlap", SCRIPT)
ro = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ro)


def write(folder, name, text):
    path = os.path.join(folder, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return path


class OverlapValues(unittest.TestCase):
    def test_hand_derived_example(self):
        # A = {x: 2, y: 1, z: 1}, B = {y: 2, z: 1, w: 1}; N_A = N_B = 4
        # shared {y, z}, union {x, y, z, w}: Jaccard = 2/4 = 0.5
        # sum p q = (1/4)(2/4) + (1/4)(1/4) = 0.1875; sum p^2 = sum q^2 = 1/4 + 1/16 + 1/16 = 0.375
        # Morisita-Horn = 2 * 0.1875 / 0.75 = 0.5
        # shared fraction of A = (1 + 1)/4 = 0.5; of B = (2 + 1)/4 = 0.75
        m = ro.overlap_metrics({"x": 2, "y": 1, "z": 1}, {"y": 2, "z": 1, "w": 1})
        self.assertEqual(m["shared_clonotypes"], 2)
        self.assertAlmostEqual(m["jaccard"], 0.5)
        self.assertAlmostEqual(m["morisita_horn"], 0.5)
        self.assertAlmostEqual(m["shared_fraction_a"], 0.5)
        self.assertAlmostEqual(m["shared_fraction_b"], 0.75)

    def test_abundance_weighting(self):
        # one shared clone holding 90 % of each sample plus 10 private singletons per sample
        # Jaccard = 1/21; sum p q = 0.81, sum p^2 = sum q^2 = 0.81 + 10 * 0.0001 = 0.811
        # Morisita-Horn = 1.62 / 1.622
        a = {"big": 90, **{f"a{i}": 1 for i in range(10)}}
        b = {"big": 90, **{f"b{i}": 1 for i in range(10)}}
        m = ro.overlap_metrics(a, b)
        self.assertAlmostEqual(m["jaccard"], 1 / 21)
        self.assertAlmostEqual(m["morisita_horn"], 1.62 / 1.622)

    def test_identity_scaling_disjoint_symmetry(self):
        a = {"x": 5, "y": 3, "z": 1}
        self.assertEqual((ro.overlap_metrics(a, a)["jaccard"], ro.overlap_metrics(a, a)["morisita_horn"]), (1.0, 1.0))
        scaled = ro.overlap_metrics(a, {k: 10 * v for k, v in a.items()})
        self.assertAlmostEqual(scaled["morisita_horn"], 1.0)
        disjoint = ro.overlap_metrics(a, {"q": 4})
        self.assertEqual((disjoint["jaccard"], disjoint["morisita_horn"]), (0.0, 0.0))
        b = {"x": 1, "w": 7, "y": 2}
        self.assertAlmostEqual(ro.overlap_metrics(a, b)["morisita_horn"], ro.overlap_metrics(b, a)["morisita_horn"])

    def test_empty_sample_raises(self):
        with self.assertRaises(ro.InputError):
            ro.overlap_metrics({"x": 0}, {"x": 1})


class OverlapAnalysis(unittest.TestCase):
    SAMPLES = {"tumour": {"A": 50, "B": 20, "C": 5, "D": 1, "E": 1},
               "blood": {"A": 2, "B": 1, "F": 100, "G": 30, **{f"s{i}": 1 for i in range(67)}}}
    QC = {"tumour": {"chain": "TRB"}, "blood": {"chain": "TRB"}}

    def test_requires_two_samples_and_valid_pairs(self):
        with self.assertRaises(ro.InputError):
            ro.analyse_overlap({"tumour": self.SAMPLES["tumour"]}, self.QC, [])
        with self.assertRaises(ro.InputError):
            ro.analyse_overlap(self.SAMPLES, self.QC, [], pairs="tumour:normal")

    def test_downsampling_reproducible_and_depth_note(self):
        res = ro.analyse_overlap(self.SAMPLES, self.QC, [])
        self.assertTrue(any("depth" in n for n in res["notes"]))
        a = ro.analyse_overlap(self.SAMPLES, self.QC, [], downsample="min", iterations=8, seed=5)
        b = ro.analyse_overlap(self.SAMPLES, self.QC, [], downsample="min", iterations=8, seed=5)
        self.assertEqual(a, b)
        ds = a["pairs"][0]["downsampled"]
        self.assertEqual(ds["depth"], 77)
        self.assertLessEqual(ds["metrics"]["jaccard"]["p97_5"], 1.0)
        self.assertEqual(ds["metrics"]["richness_a"]["mean"], 5.0)   # tumour is drawn at full depth

    def test_min_depth_uses_only_the_requested_pairs(self):
        # sample C (depth 3) is not in the requested pair, so it must not set the common depth;
        # a repeated pair is reported once
        samples = {"A": {"x": 50, "y": 50}, "B": {"x": 40, "z": 40}, "C": {"x": 2, "q": 1}}
        qc = {n: {} for n in samples}
        res = ro.analyse_overlap(samples, qc, [], pairs="A:B,A:B", downsample="min", iterations=3)
        self.assertEqual(res["settings"]["downsample_depth"], 80)
        self.assertEqual(len(res["pairs"]), 1)
        self.assertFalse(any("'C'" in n for n in res["notes"]))

    def test_depth_moves_jaccard_in_either_direction(self):
        # worked example of references/formulas.md, each sample subsampled to 55 of its 110 units:
        # A, B share one clonotype of 100 units and hold 10 private singletons each:
        #   full depth J = 1/21 = 0.048; each private singleton survives with P = 55/110 = 0.5,
        #   so about 1 shared of 1 + 5 + 5 clonotypes: J rises to about 0.09
        # C, D share 10 singletons and hold one private clonotype of 100 units each:
        #   full depth J = 10/12 = 0.83; a shared singleton is seen in both with P = 0.25:
        #   about 2.5 shared of 2 + 7.5 clonotypes: J falls to about 0.26
        samples = {"A": {"s": 100, **{f"a{i}": 1 for i in range(10)}},
                   "B": {"s": 100, **{f"b{i}": 1 for i in range(10)}},
                   "C": {"a": 100, **{f"t{i}": 1 for i in range(10)}},
                   "D": {"b": 100, **{f"t{i}": 1 for i in range(10)}}}
        res = ro.analyse_overlap(samples, {n: {} for n in samples}, [], pairs="A:B,C:D",
                                 downsample="55", iterations=400, seed=1)
        ab, cd = res["pairs"]
        self.assertAlmostEqual(ab["observed"]["jaccard"], 1 / 21)
        self.assertAlmostEqual(cd["observed"]["jaccard"], 10 / 12)
        self.assertGreater(ab["downsampled"]["metrics"]["jaccard"]["mean"], 0.08)
        self.assertLess(cd["downsampled"]["metrics"]["jaccard"]["mean"], 0.32)
        # the shared fraction of A stays at 100/110 because the shared clonotype is always drawn
        self.assertAlmostEqual(ab["downsampled"]["metrics"]["shared_fraction_a"]["mean"], 100 / 110, places=2)

    def test_second_hand_example(self):
        # A = {x: 3, y: 1}, B = {x: 1, z: 1}: J = 1/3; sum p q = 0.75 * 0.5 = 0.375,
        # sum p^2 = 0.5625 + 0.0625 = 0.625, sum q^2 = 0.5: MH = 0.75 / 1.125 = 2/3;
        # shared fractions 3/4 (A) and 1/2 (B)
        m = ro.overlap_metrics({"x": 3, "y": 1}, {"x": 1, "z": 1})
        self.assertAlmostEqual(m["jaccard"], 1 / 3)
        self.assertAlmostEqual(m["morisita_horn"], 2 / 3)
        self.assertEqual((m["shared_fraction_a"], m["shared_fraction_b"]), (0.75, 0.5))

    def test_full_depth_identical_samples(self):
        same = {"s1": {"x": 3, "y": 2}, "s2": {"x": 3, "y": 2}}
        res = ro.analyse_overlap(same, {"s1": {}, "s2": {}}, [], downsample="min", iterations=3)
        metrics = res["pairs"][0]["downsampled"]["metrics"]
        self.assertEqual((metrics["jaccard"]["mean"], metrics["morisita_horn"]["mean"]), (1.0, 1.0))


class OverlapCommandLine(unittest.TestCase):
    def test_cli_two_files(self):
        with tempfile.TemporaryDirectory() as folder:
            a = write(folder, "pre.tsv", "cdr3_aa\tcount\nX\t2\nY\t1\nZ\t1\n")
            b = write(folder, "post.tsv", "cdr3_aa\tcount\nY\t2\nZ\t1\nW\t1\n")
            proc = subprocess.run([sys.executable, SCRIPT, "--input", a, b, "--clone-col", "cdr3_aa",
                                   "--count-col", "count"], capture_output=True, text=True, timeout=120)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            pair = json.loads(proc.stdout)["pairs"][0]
            self.assertEqual((pair["sample_a"], pair["sample_b"]), ("pre", "post"))
            self.assertAlmostEqual(pair["observed"]["morisita_horn"], 0.5)
            proc = subprocess.run([sys.executable, SCRIPT, "--input", a, "--clone-col", "cdr3_aa",
                                   "--count-col", "count"], capture_output=True, text=True, timeout=120)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("at least two samples", proc.stderr)
            proc = subprocess.run([sys.executable, SCRIPT, "--input", a, b, "--clone-col", "cdr3_aa",
                                   "--count-col", "count", "--downsample", "min", "--iterations", "4",
                                   "--format", "tsv"], capture_output=True, text=True, timeout=120)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("ds_morisita_horn_mean", proc.stdout.splitlines()[0].split("\t"))


if __name__ == "__main__":
    unittest.main()
