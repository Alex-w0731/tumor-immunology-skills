"""Unit tests for skills/tme-bulk-profiling/scripts/signature_coverage.py (synthetic data only)."""
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "tme-bulk-profiling" / "scripts" / "signature_coverage.py"
sys.dont_write_bytecode = True          # keep __pycache__ out of the skill folder
spec = importlib.util.spec_from_file_location("signature_coverage", SCRIPT)
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)

MATRIX = [
    "gene\tS1\tS2\tS3",
    "CD8A\t5.1\t6.2\t7.3",
    "cd8b\t4.0\t4.5\t5.0",
    "GZMA\t3.0\t3.0\t3.0",      # constant across samples
    "PRF1\tNA\tNA\tNA",         # entirely missing: not measured (tmelens drops it)
    "CD3E\t2.0\t2.5\tNA",
    "ACTB\t10\t11\t12",
]
GMT = [
    "TCELL\tfirst definition\tCD3D",
    "CD8\tdesc\tCD8A\tCD8B",
    "CYT\tdesc\tGZMA\tPRF1",
    "TCELL\tdesc\tCD3D\tCD3E\tCD3G\tCD3E",   # repeated name: this last definition wins
    "EMPTY\tdesc\t",
    "MOUSE\tdesc\tCd8a",
    "BROKEN line separated by spaces",
]


def by_name(rows):
    return {r["signature"]: r for r in rows}


def write(tmp, name, lines):
    path = Path(tmp) / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class TestPureFunctions(unittest.TestCase):
    def setUp(self):
        self.sets, self.gmt_warnings = sc.parse_gmt(GMT)
        records = [line.split("\t") for line in MATRIX]
        self.measured, self.constant, self.ids = sc.scan_matrix(records)

    def test_parse_gmt(self):
        self.assertEqual(list(self.sets), ["TCELL", "CD8", "CYT", "EMPTY", "MOUSE"])
        self.assertEqual(self.sets["TCELL"], ["CD3D", "CD3E", "CD3G", "CD3E"])
        self.assertEqual(self.sets["EMPTY"], [])
        codes = {w["code"] for w in self.gmt_warnings}
        self.assertEqual(codes, {"GMT_LINES_SKIPPED", "GMT_DUPLICATE_NAMES"})

    def test_scan_matrix(self):
        self.assertEqual(self.measured, {"CD8A", "CD8B", "GZMA", "CD3E", "ACTB"})
        self.assertEqual(self.constant, {"GZMA"})

    def test_coverage_golden_values(self):
        rows = by_name(sc.coverage_table(self.sets, self.measured, self.constant))
        # CD8: {CD8A, CD8B}, both measured (cd8b matches after upper-casing) -> 2/2 = 1.0
        self.assertEqual((rows["CD8"]["n_found"], rows["CD8"]["coverage"], rows["CD8"]["flagged"]), (2, 1.0, False))
        # CYT: GZMA measured but constant, PRF1 all missing -> 1/2 = 0.5, not below 0.5 (tmelens keeps it);
        # informative coverage (1 - 1) / 2 = 0
        cyt = rows["CYT"]
        self.assertEqual((cyt["coverage"], cyt["flagged"], cyt["n_constant"], cyt["informative_coverage"]),
                         (0.5, False, 1, 0.0))
        self.assertEqual(cyt["missing"], ["PRF1"])
        # TCELL: CD3E repeated counts once -> distinct {CD3D, CD3E, CD3G}; 1 found -> 1/3 = 0.3333, flagged
        t = rows["TCELL"]
        self.assertEqual((t["n_genes"], t["n_found"], t["coverage"], t["flagged"]), (3, 1, 0.3333, True))
        self.assertTrue(rows["EMPTY"]["flagged"])
        self.assertEqual(rows["EMPTY"]["coverage"], 0.0)

    def test_min_genes_rule(self):
        rows = by_name(sc.coverage_table(self.sets, self.measured, self.constant, 0.5, 2))
        self.assertTrue(rows["CYT"]["flagged"])      # one gene found, two required
        self.assertFalse(rows["CD8"]["flagged"])

    def test_warnings(self):
        rows = sc.coverage_table(self.sets, self.measured, self.constant)
        codes = {w["code"] for w in sc.coverage_warnings(self.sets, rows, self.ids, 0.5)}
        self.assertIn("GMT_LOWER_CASE_GENES", codes)
        self.assertIn("GMT_SET_WITHOUT_GENES", codes)
        self.assertIn("COVERAGE_PASSES_ONLY_WITH_CONSTANT_GENES", codes)
        self.assertNotIn("IDENTIFIER_MISMATCH", codes)

    def test_identifier_mismatch(self):
        records = [["id", "S1", "S2"], ["ENSG00000153563.16", "1", "2"], ["ENSG00000172116.23", "3", "4"]]
        measured, constant, ids = sc.scan_matrix(records)
        sets = {"CD8": ["CD8A", "CD8B"]}
        rows = sc.coverage_table(sets, measured, constant)
        self.assertEqual(rows[0]["n_found"], 0)
        codes = {w["code"] for w in sc.coverage_warnings(sets, rows, ids, 0.5)}
        self.assertIn("IDENTIFIER_MISMATCH", codes)
        self.assertEqual(sc.identifier_kind(["7157", "1956", "3845"]), "entrez_or_integer")

    def test_orientation_and_layout(self):
        samples_by_genes = [["sample", "CD8A", "GZMA"], ["S1", "1", "2"], ["S2", "3", "2"]]
        measured, constant, _ = sc.scan_matrix(samples_by_genes, "samples-by-genes")
        self.assertEqual((measured, constant), ({"CD8A", "GZMA"}, {"GZMA"}))
        r_layout = [["S1", "S2"], ["CD8A", "1", "2"], ["GZMA", "4", "4"]]   # header one field short
        measured, constant, _ = sc.scan_matrix(r_layout)
        self.assertEqual((measured, constant), ({"CD8A", "GZMA"}, {"GZMA"}))
        with self.assertRaises(sc.InputError):
            sc.scan_matrix([["gene", "S1", "S2"], ["CD8A", "1"]])
        with self.assertRaises(sc.InputError):
            sc.scan_matrix([["gene", "S1"]])

    def test_single_sample_is_never_constant(self):
        measured, constant, _ = sc.scan_matrix([["gene", "S1"], ["CD8A", "4"], ["GZMA", "4"]])
        self.assertEqual((measured, constant), ({"CD8A", "GZMA"}, set()))

    def test_tmelens_aliases(self):
        # tmelens 0.1.0 renames CTGF -> CCN2 only when CCN2 is absent (preprocess.prepare_expression)
        records = [["gene", "S1", "S2"], ["CTGF", "5", "5"], ["FAM101B", "1", "2"], ["RFLNB", "3", "4"],
                   ["ACTA2", "6", "7"]]
        measured, constant, _ = sc.scan_matrix(records)
        measured, constant, applied = sc.apply_aliases(measured, constant)
        self.assertEqual(applied, {"CTGF": "CCN2"})                 # RFLNB already present: no rename
        self.assertEqual(measured, {"CCN2", "FAM101B", "RFLNB", "ACTA2"})
        self.assertEqual(constant, {"CCN2"})                         # constant status follows the rename
        # Pan-F-TBRS-like set of 4 genes: ACTA2, CCN2, RFLNB, TAGLN -> 3 found / 4 = 0.75
        rows = sc.coverage_table({"PANF": ["ACTA2", "CCN2", "RFLNB", "TAGLN"]}, measured, constant)
        self.assertEqual((rows[0]["n_found"], rows[0]["coverage"], rows[0]["n_constant"]), (3, 0.75, 1))
        self.assertEqual(rows[0]["missing"], ["TAGLN"])

    def test_tied_columns_samples_by_genes(self):
        # columns as genes; GZMA tied in every sample, CD8A has one missing value but two numbers
        recs = [["sample", "CD8A", "GZMA", "PRF1"], ["S1", "1", "2", "NA"], ["S2", "NA", "2", ""],
                ["S3", "3", "2", "NA"]]
        measured, constant, ids = sc.scan_matrix(recs, "samples-by-genes")
        self.assertEqual((measured, constant, ids), ({"CD8A", "GZMA"}, {"GZMA"}, ["CD8A", "GZMA", "PRF1"]))


class TestCommandLine(unittest.TestCase):
    def test_json_and_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            gmt = write(tmp, "sets.gmt", GMT)
            mat = write(tmp, "expr.tsv", MATRIX)
            proc = subprocess.run([sys.executable, str(SCRIPT), str(gmt), str(mat)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            rep = json.loads(proc.stdout)
            self.assertEqual((rep["n_matrix_genes"], rep["n_signatures"], rep["n_flagged"]), (5, 5, 2))
            out = Path(tmp) / "cov.tsv"
            proc = subprocess.run([sys.executable, str(SCRIPT), str(gmt), str(mat), "--format", "tsv",
                                   "--min-coverage", "0.8", "--out", str(out)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            lines = out.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0].split("\t")[0], "signature")
            cyt = next(line.split("\t") for line in lines if line.startswith("CYT\t"))
            self.assertEqual(cyt[6], "yes")          # 0.5 < 0.8
            self.assertEqual(cyt[7], "PRF1")

    def test_gzip_comments_blank_lines_and_alias_in_report(self):
        import gzip
        with tempfile.TemporaryDirectory() as tmp:
            gmt = write(tmp, "sets.gmt", ["PANF\tdesc\tCCN2\tACTA2\tTAGLN"])
            mat = Path(tmp) / "expr.tsv.gz"
            with gzip.open(mat, "wt", encoding="utf-8") as fh:
                fh.write("# synthetic\n\ngene\tS1\tS2\nCTGF\t1\t2\n\nACTA2\t3\t4\n")
            proc = subprocess.run([sys.executable, str(SCRIPT), str(gmt), str(mat)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            rep = json.loads(proc.stdout)
            self.assertEqual(rep["renamed_symbols"], {"CTGF": "CCN2"})
            # CCN2 (renamed from CTGF) and ACTA2 found, TAGLN missing -> 2/3 = 0.6667, not flagged
            row = rep["signatures"][0]
            self.assertEqual((row["n_found"], row["coverage"], row["flagged"]), (2, 0.6667, False))

    def test_bad_inputs_exit_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            mat = write(tmp, "expr.tsv", MATRIX)
            empty_gmt = write(tmp, "empty.gmt", ["only one field"])
            for args in ([str(empty_gmt), str(mat)],
                         [str(Path(tmp) / "absent.gmt"), str(mat)],
                         [str(write(tmp, "ok.gmt", GMT)), str(mat), "--min-coverage", "1.5"]):
                proc = subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)
                self.assertEqual(proc.returncode, 2, args)
                self.assertIn("error:", proc.stderr)


if __name__ == "__main__":
    unittest.main()
