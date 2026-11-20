"""Unit tests for skills/immune-signature-curation/scripts/signature_overlap.py (standard library only)."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "immune-signature-curation" / "scripts" / "signature_overlap.py"
_spec = importlib.util.spec_from_file_location("signature_overlap", SCRIPT)
so = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(so)

# Two curated tmelens panels (src/tmelens/data/signatures.json, read November 2026).
CHECKPOINTS = ["PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1"]
EXHAUSTION = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1", "CXCL13", "LAYN"]
CYTOTOXIC = ["CTSW", "GNLY", "GZMA", "GZMB", "GZMH", "KLRB1", "KLRD1", "KLRK1", "NKG7", "PRF1"]
CYT = ["GZMA", "PRF1"]
PANEL = ("Immune_Checkpoints\tcurated\t" + "\t".join(CHECKPOINTS) + "\n"
         "T_cell_Exhaustion\tcurated\t" + "\t".join(EXHAUSTION) + "\n"
         "Danaher_Cytotoxic_cells\t10.1186/s40425-017-0215-8\t" + "\t".join(CYTOTOXIC) + "\n"
         "CYT\t10.1016/j.cell.2014.12.033\t" + "\t".join(CYT) + "\n")


def prep(text, **kw):
    return so.prepare(so.parse_gmt(text)[0], **kw)


class TestFormulas(unittest.TestCase):
    def test_golden_tmelens_panels(self):
        # Hand derivation: shared = {PDCD1, CTLA4, LAG3, HAVCR2, TIGIT} = 5 genes;
        # union = 8 + 9 - 5 = 12, so J = 5/12 = 0.4167; OC = 5 / min(8, 9) = 5/8 = 0.625.
        st = so.overlap_stats(set(CHECKPOINTS), set(EXHAUSTION))
        self.assertEqual(st["n_shared"], 5)
        self.assertAlmostEqual(st["jaccard"], 5 / 12)
        self.assertAlmostEqual(st["overlap_coefficient"], 0.625)
        self.assertEqual(so.classify(st), "overlapping")

    def test_nested_and_identical(self):
        # CYT (GZMA, PRF1) lies inside the 10-gene cytotoxic panel: J = 2/10, OC = 2/2 = 1.
        st = so.overlap_stats(set(CYTOTOXIC), set(CYT))
        self.assertAlmostEqual(st["jaccard"], 0.2)
        self.assertEqual(st["overlap_coefficient"], 1.0)
        self.assertEqual(so.classify(st), "nested")
        self.assertEqual(so.classify(so.overlap_stats({"A", "B"}, {"B", "A"})), "identical")

    def test_threshold_edges(self):
        # {A,B,C,D} vs {A,B,C,E}: shared 3, union 5, J = 0.6 >= 0.5 -> near_duplicate.
        st = so.overlap_stats(set("ABCD"), set("ABCE"))
        self.assertAlmostEqual(st["jaccard"], 0.6)
        self.assertEqual(so.classify(st), "near_duplicate")
        self.assertEqual(so.classify(st, jaccard_min=0.7, overlap_min=0.8), "shared")
        # Exactly at the threshold counts (thresholds are inclusive).
        st2 = so.overlap_stats(set("ABC"), set("BCD"))     # J = 2/4 = 0.5, OC = 2/3
        self.assertEqual(so.classify(st2), "near_duplicate")
        self.assertEqual(so.classify(so.overlap_stats({"A"}, {"B"})), "distinct")

    def test_ayers_lists_from_eval_panel(self):
        # Gene lists as in the tmelens library (IFNG_Ayers6, TIS_Ayers18). Hand derivation:
        # TIS vs IFNG share {STAT1, IDO1, CXCL9} = 3; union 18 + 6 - 3 = 21; J = 3/21 = 1/7; OC = 3/6
        # = 0.5, exactly the default cut-off, so 'overlapping' (inclusive; EnrichmentMap's preset
        # 'overlap > 0.5' would not link this pair). TIS vs CHECKPOINTS share {CD274, PDCD1LG2, LAG3,
        # TIGIT, IDO1} = 5; union 18 + 8 - 5 = 21; J = 5/21; OC = 5/8.
        ifng = {"IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"}
        tis = {"CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
               "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"}
        st = so.overlap_stats(tis, ifng)
        self.assertEqual((st["n_shared"], st["shared_genes"]), (3, ["CXCL9", "IDO1", "STAT1"]))
        self.assertAlmostEqual(st["jaccard"], 1 / 7)
        self.assertEqual(st["overlap_coefficient"], 0.5)
        self.assertEqual(so.classify(st), "overlapping")
        st = so.overlap_stats(tis, set(CHECKPOINTS))
        self.assertAlmostEqual(st["jaccard"], 5 / 21)
        self.assertAlmostEqual(st["overlap_coefficient"], 5 / 8)

    def test_empty_sets(self):
        st = so.overlap_stats(set(), {"A"})
        self.assertIsNone(st["jaccard"])
        self.assertIsNone(st["overlap_coefficient"])
        self.assertEqual(so.classify(st), "empty")


class TestCompare(unittest.TestCase):
    def test_within_file(self):
        r = so.compare(prep(PANEL))
        flags = {(p["set_a"], p["set_b"]): p["flag"] for p in r["pairs"]}
        self.assertEqual(flags[("Danaher_Cytotoxic_cells", "CYT")], "nested")
        self.assertEqual(flags[("Immune_Checkpoints", "T_cell_Exhaustion")], "overlapping")
        self.assertEqual(len(r["pairs"]), 2)                  # four distinct pairs are not listed
        self.assertEqual(r["summary"]["n_pairs_compared"], 6)
        self.assertEqual(r["pairs"][0]["flag"], "nested")      # sorted by flag severity
        best = {s["name"]: s["best_match"] for s in r["sets"]}
        self.assertEqual(best["CYT"]["set"], "Danaher_Cytotoxic_cells")
        self.assertEqual(r["clusters"], [])
        self.assertEqual(len(so.compare(prep(PANEL), all_pairs=True)["pairs"]), 6)

    def test_best_match_ties_and_duplicate_genes(self):
        # X shares {a, b} with Z and with Y: J = 2/6 = 1/3 and OC = 2/4 = 1/2 for both, a tie; the best
        # match is the name first in sort order (Y), whatever the file order.
        r = so.compare(prep("X\tx\ta\tb\tc\td\nZ\tx\ta\tb\tg\th\nY\tx\ta\tb\te\tf\n"))
        best = {s["name"]: s["best_match"] for s in r["sets"]}
        self.assertEqual(best["X"], {"set": "Y", "jaccard": 0.3333, "overlap_coefficient": 0.5})
        self.assertEqual(best["Z"]["set"], "X")
        # A gene listed twice counts once: {a, a, b} vs {a, c} -> shared 1, union 3, J = 1/3, OC = 1/2.
        r = so.compare(prep("D\tx\ta\ta\tb\nE\tx\ta\tc\n"))
        self.assertEqual((r["sets"][0]["n_input"], r["pairs"][0]["jaccard"], r["pairs"][0]["flag"]),
                         (2, 0.3333, "overlapping"))

    def test_clusters_of_near_duplicates(self):
        text = "A\tx\tG1\tG2\tG3\tG4\nB\tx\tG1\tG2\tG3\tG5\nC\tx\tG4\tG3\tG2\tG1\nD\tx\tG9\tG8\n"
        r = so.compare(prep(text))
        self.assertEqual(r["clusters"], [["A", "B", "C"]])

    def test_universe_restriction(self):
        # With only GZMA measured, CYT keeps one gene and stays nested in the cytotoxic panel.
        r = so.compare(prep(PANEL, universe={"GZMA", "PDCD1", "CTLA4"}))
        sets = {s["name"]: s for s in r["sets"]}
        self.assertEqual((sets["CYT"]["n_input"], sets["CYT"]["n_used"]), (2, 1))
        self.assertEqual(sets["CYT"]["dropped_by_universe"], ["PRF1"])
        flags = {(p["set_a"], p["set_b"]): p["flag"] for p in r["pairs"]}
        self.assertEqual(flags[("Immune_Checkpoints", "T_cell_Exhaustion")], "identical")

    def test_against_and_ignore_case(self):
        mine = "CYT_local\tuser\tgzma\tprf1\n"
        published = "CYT\t10.1016/j.cell.2014.12.033\tGZMA\tPRF1\n"
        a, b = so.parse_gmt(mine)[0], so.parse_gmt(published)[0]
        exact = so.compare(so.prepare(a), so.prepare(b))
        self.assertEqual(exact["pairs"], [])
        folded = so.compare(so.prepare(a, ignore_case=True), so.prepare(b, ignore_case=True))
        self.assertEqual(folded["pairs"][0]["flag"], "identical")
        self.assertEqual(folded["summary"]["n_sets_b"], 1)

    def test_single_set_and_parsing(self):
        r = so.compare(prep("Only\tx\tCD3E\tCD3D\n"))
        self.assertEqual(r["summary"]["n_pairs_compared"], 0)
        sets, skipped = so.parse_gmt("# note\nshort\tx\n\nOK\tx\t CD3E \t\tCD3D\n")
        self.assertEqual(sets, {"OK": ["CD3E", "CD3D"]})
        self.assertEqual(skipped, [1, 2])
        with self.assertRaises(so.InputError):
            so.parse_gmt("A\tx\tG1\nA\tx\tG2\n")
        with self.assertRaises(ValueError):
            so.compare(prep(PANEL), jaccard_min=0)

    def test_tsv(self):
        tsv = so.to_tsv(so.compare(prep(PANEL)))
        lines = tsv.strip().split("\n")
        self.assertEqual(lines[0].split("\t"), list(so.TSV_COLUMNS))
        self.assertEqual(lines[1].split("\t")[:8],
                         ["Danaher_Cytotoxic_cells", "CYT", "10", "2", "2", "0.2000", "1.0000", "nested"])


class TestCommandLine(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)

    def test_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            gmt = os.path.join(tmp, "panel.gmt")
            uni = os.path.join(tmp, "genes.txt")
            dup = os.path.join(tmp, "dup.gmt")
            Path(gmt).write_text(PANEL, encoding="utf-8")
            Path(uni).write_text("gene\nGZMA\nPRF1\nPDCD1\n", encoding="utf-8")
            Path(dup).write_text("A\tx\tG1\nA\tx\tG2\n", encoding="utf-8")
            p = self.run_cli(gmt)
            self.assertEqual(p.returncode, 0, p.stderr)
            rep = json.loads(p.stdout)
            self.assertEqual(rep["summary"]["by_flag"], {"nested": 1, "overlapping": 1})
            p = self.run_cli(gmt, "--universe", uni, "--format", "tsv")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(p.stdout.startswith("set_a\tset_b"))
            bad = self.run_cli(dup)
            self.assertEqual(bad.returncode, 2)
            self.assertIn("more than once", bad.stderr)
            self.assertEqual(self.run_cli(os.path.join(tmp, "none.gmt")).returncode, 2)
            unwritable = self.run_cli(gmt, "--out", os.path.join(tmp, "no_such_dir", "o.json"))
            self.assertEqual(unwritable.returncode, 2)
            self.assertIn("cannot write", unwritable.stderr)
            self.assertEqual(self.run_cli(gmt, "--overlap", "1.5").returncode, 2)
            self.assertIn("usage", self.run_cli("--help").stdout.lower())


if __name__ == "__main__":
    unittest.main()
