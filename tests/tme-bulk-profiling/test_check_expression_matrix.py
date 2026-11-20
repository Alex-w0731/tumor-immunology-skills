"""Unit tests for skills/tme-bulk-profiling/scripts/check_expression_matrix.py (synthetic data only)."""
import gzip
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "tme-bulk-profiling" / "scripts" / "check_expression_matrix.py"
sys.dont_write_bytecode = True          # keep __pycache__ out of the skill folder
spec = importlib.util.spec_from_file_location("check_expression_matrix", SCRIPT)
cem = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cem)

ANCHORS = ["ACTB", "GAPDH", "PTPRC", "CD8A", "CD3E", "GZMB"]


def records(header, rows):
    """Header list plus (label, values) rows -> list of string records."""
    out = [list(header)]
    for label, values in rows:
        out.append([label] + [str(v) for v in values])
    return out


def codes(report):
    return {i["code"]: i["severity"] for i in report["issues"]}


def log_rows(labels, n_samples=4):
    """Deterministic log2-scale values between 1 and 14 (not integers, no negatives)."""
    rows = []
    for i, label in enumerate(labels):
        rows.append((label, [round(1.25 + ((i * 3 + j * 5) % 11) * 1.1 + 0.1 * j, 3) for j in range(n_samples)]))
    return rows


class TestPureHelpers(unittest.TestCase):
    def test_quantile_golden(self):
        # numpy 'linear' (type 7): h = (n - 1) q; [1, 2, 3, 4] at q = 0.25 -> h = 0.75 -> 1 + 0.75 * (2 - 1) = 1.75
        self.assertAlmostEqual(cem.quantile([1, 2, 3, 4], 0.25), 1.75)
        # 101 values 0..100 at q = 0.99 -> h = 99 exactly -> 99
        self.assertAlmostEqual(cem.quantile(list(range(101)), 0.99), 99.0)
        self.assertEqual(cem.quantile([7.0], 0.5), 7.0)
        with self.assertRaises(ValueError):
            cem.quantile([], 0.5)

    def test_parse_value(self):
        self.assertEqual(cem.parse_value("3.5"), ("number", 3.5))
        self.assertEqual(cem.parse_value(" -2 "), ("number", -2.0))
        for token in ("NA", "", "nan", "#N/A", "null"):
            self.assertEqual(cem.parse_value(token)[0], "missing", token)
        self.assertEqual(cem.parse_value("-Inf")[0], "infinite")
        for token in ("1,5", "1_000", "high"):
            self.assertEqual(cem.parse_value(token)[0], "text", token)

    def test_classify_identifier(self):
        cases = {
            "ENSG00000141510.17": "ensembl_gene",
            "ENSMUSG00000059552": "ensembl_gene",
            "ENST00000269305": "ensembl_transcript",
            "NM_000546.6": "refseq",
            "201422_at": "probe",
            "ILMN_1651199": "probe",
            "1-Mar": "excel_date",
            "ACTB|60": "composite",
            "7157": "integer",
            "HLA-DRA": "other",
            "": "empty",
        }
        for label, expected in cases.items():
            self.assertEqual(cem.classify_identifier(label), expected, label)

    def test_geo2r_rule_golden(self):
        base = {"min": 0.0, "q25": 0.0, "median": 2.0, "q75": 10.0, "q99": 99.0, "max": 99.9}
        self.assertFalse(cem.geo2r_linear(base))             # q99 <= 100 and q25 not > 0
        self.assertTrue(cem.geo2r_linear({**base, "q25": 0.1}))  # range 99.9 > 50 with q25 > 0
        self.assertTrue(cem.geo2r_linear({**base, "q99": 150.0, "max": 200.0}))

    def test_infer_scale_rules(self):
        log_q = {"min": 0.0, "q25": 3.0, "median": 5.0, "q75": 8.0, "q99": 14.0, "max": 16.0}
        self.assertEqual(cem.infer_scale(log_q, False)["call"], "log")
        # integer counts below 50 with many zeros: GEO2R misses them, the count rule catches them
        counts_q = {"min": 0.0, "q25": 0.0, "median": 2.0, "q75": 6.0, "q99": 30.0, "max": 40.0}
        sc = cem.infer_scale(counts_q, True)
        self.assertEqual(sc["call"], "linear")
        self.assertFalse(sc["tmelens_auto_log_would_transform"])
        centered_q = {"min": -2.0, "q25": -0.6, "median": 0.0, "q75": 0.6, "q99": 2.0, "max": 2.5}
        self.assertEqual(cem.infer_scale(centered_q, False, "gene")["call"], "centered")
        flat = {"min": 1.0, "q25": 1.0, "median": 1.0, "q75": 1.0, "q99": 1.0, "max": 1.0}
        self.assertEqual(cem.infer_scale(flat, False)["call"], "undetermined")

    def test_find_duplicates(self):
        d = cem.find_duplicates(["CD8A", "CD8A", "GZMA", "gzma", "PRF1"])
        self.assertEqual(d["n_duplicated_ids"], 1)
        self.assertEqual(d["n_entries_involved"], 2)
        self.assertEqual(d["n_more_after_upper_casing"], 1)
        self.assertEqual(d["examples_after_upper_casing"], ["GZMA"])


class TestProfileMatrix(unittest.TestCase):
    def test_clean_log_matrix(self):
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], log_rows(ANCHORS)))
        self.assertEqual(rep["orientation"]["call"], "genes_x_samples")
        self.assertEqual(rep["orientation"]["confidence"], "high")
        self.assertEqual(rep["identifiers"]["type"], "gene_symbol")
        self.assertEqual(rep["scale"]["call"], "log")
        c = codes(rep)
        self.assertIn("SCALE_LOG", c)
        self.assertNotIn("error", c.values())
        message = next(i["message"] for i in rep["issues"] if i["code"] == "SCALE_LOG")
        self.assertIn("log=False", message)
        # tmelens 0.1.0 cli.main maps --log no to log=False, so the advice must be to pass it
        self.assertIn("pass --log no", message)

    def test_quantiles_with_ties_and_scale_golden(self):
        # 3 genes x 4 samples with heavy ties; the 12 values sorted: 1 1 2 2 2 2 2 2 2 3 60 60.
        # Type 7 with n = 12: q25 -> h = 11 * 0.25 = 2.75 -> 2 + 0.75 * (2 - 2) = 2;
        # median -> h = 5.5 -> 2; q99 -> h = 10.89 -> 60 + 0.89 * (60 - 60) = 60.
        # GEO2R: q99 60 <= 100, but range 59 > 50 with q25 2 > 0 -> linear. All values are integers
        # with max > 1, so the counts rule holds too.
        rows = [("ACTB", [1, 2, 60, 2]), ("GAPDH", [1, 3, 60, 2]), ("PTPRC", [2, 2, 2, 2])]
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], rows))
        q = rep["scale"]["quantiles"]
        self.assertEqual((q["q25"], q["median"], q["q99"], q["max"]), (2.0, 2.0, 60.0, 60.0))
        self.assertEqual(rep["scale"]["call"], "linear")
        self.assertTrue(rep["scale"]["geo2r_rule_linear"])        # range 59 > 50 with q25 2 > 0
        self.assertEqual(rep["constant"]["n_constant_genes"], 1)   # PTPRC
        # column sums 4, 7, 122, 6: not one million -> counts, totals differ 122 / 4 = 30.5-fold
        c = codes(rep)
        self.assertIn("COUNT_LIKE", c)
        self.assertEqual(rep["sample_totals"]["max_over_min"], 30.5)

    def test_integer_values_summing_to_one_million_are_rounded_tpm(self):
        # each column is a permutation of 600000, 300000, 99000, 1000 -> sums exactly 1e6
        col = [600000, 300000, 99000, 1000]
        rows = [(g, [col[(i + j) % 4] for j in range(3)]) for i, g in enumerate(ANCHORS[:4])]
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3"], rows))
        c = codes(rep)
        self.assertIn("SUM_TO_ONE_MILLION", c)
        self.assertNotIn("COUNT_LIKE", c)
        msg = next(i["message"] for i in rep["issues"] if i["code"] == "SUM_TO_ONE_MILLION")
        self.assertIn("rounded TPM or CPM", msg)

    def test_transposed_matrix_is_an_error(self):
        rows = log_rows(ANCHORS)
        header = ["sample"] + [g for g, _ in rows]
        transposed = [(f"S{j + 1}", [r[1][j] for r in rows]) for j in range(4)]
        rep = cem.profile_matrix(records(header, transposed))
        self.assertEqual(rep["orientation"]["call"], "samples_x_genes")
        self.assertEqual(rep["orientation"]["n_genes"], len(ANCHORS))
        self.assertEqual(codes(rep)["ORIENTATION_SAMPLES_BY_GENES"], "error")
        forced = cem.profile_matrix(records(header, transposed), orientation="genes-by-samples")
        self.assertEqual(forced["orientation"]["basis"], "stated by the user")

    def test_versioned_ensembl_ids(self):
        labels = [f"ENSG{n:011d}.{v}" for n, v in ((141510, 17), (111640, 15), (81237, 20), (153563, 16))]
        rep = cem.profile_matrix(records(["id", "S1", "S2", "S3", "S4"], log_rows(labels)))
        self.assertEqual(rep["identifiers"]["type"], "ensembl_gene_id")
        self.assertEqual(rep["identifiers"]["ensembl_versioned"], 4)
        msg = next(i["message"] for i in rep["issues"] if i["code"] == "IDENTIFIERS_ENSEMBL")
        self.assertIn("version suffix", msg)

    def test_duplicates_missing_and_constant(self):
        rows = log_rows(ANCHORS + ["CD8A", "gzmb", "PRF1", "KRT18"])
        rows[-2] = ("PRF1", ["NA", "", "NA", "NA"])        # entirely missing gene
        rows[-1] = ("KRT18", [5.0, "NA", 5.0, 5.0])        # constant gene with one gap
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], rows))
        dup = rep["duplicates"]["gene_ids"]
        self.assertEqual(dup["n_duplicated_ids"], 1)
        self.assertEqual(dup["n_more_after_upper_casing"], 1)
        v = rep["values"]
        self.assertEqual(v["n_missing"], 5)
        self.assertEqual(v["genes_all_missing"], 1)
        self.assertEqual(v["genes_with_missing"], 2)
        self.assertEqual(rep["constant"]["n_constant_genes"], 1)
        c = codes(rep)
        for code in ("DUPLICATE_GENE_IDS", "DUPLICATE_GENE_IDS_AFTER_UPPER_CASING", "MISSING_VALUES", "CONSTANT_GENES"):
            self.assertIn(code, c)

    def test_tpm_like_linear_matrix(self):
        # every sample sums to exactly 1e6: 400000.5 + 299999.5 + 200000.25 + 98999.75 + 1000 = 1000000
        col = [400000.5, 299999.5, 200000.25, 98999.75, 1000.0]
        rows = [(g, [col[(i + j) % 5] for j in range(3)]) for i, g in enumerate(ANCHORS[:5])]
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3"], rows))
        self.assertEqual(rep["scale"]["call"], "linear")
        self.assertTrue(rep["sample_totals"]["sum_to_one_million"])
        c = codes(rep)
        self.assertIn("SCALE_LINEAR", c)
        self.assertIn("SUM_TO_ONE_MILLION", c)

    def test_counts_and_negative_linear(self):
        rows = [(g, [0, 3 * (i + 1), 10 * (i + 1), 150 * (i + 1)]) for i, g in enumerate(ANCHORS)]
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], rows))
        self.assertEqual(rep["scale"]["call"], "linear")
        self.assertIn("COUNT_LIKE", codes(rep))
        rows[0] = (rows[0][0], [-5, 3, 10, 150])
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], rows))
        self.assertEqual(codes(rep)["NEGATIVE_ON_LINEAR_SCALE"], "error")

    def test_gene_centered_standardised_matrix(self):
        base = [-1.5, -0.5, 0.5, 1.5]
        sd = math.sqrt(sum(x * x for x in base) / 3)        # sample SD of base = sqrt(5/3)
        rows = []
        for i, g in enumerate(ANCHORS):
            perm = base[i % 4:] + base[:i % 4]
            rows.append((g, [round(x / sd, 6) for x in perm]))
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], rows))
        self.assertEqual(rep["scale"]["call"], "centered")
        self.assertEqual(rep["scale"]["centered_axis"], "gene")
        self.assertTrue(rep["scale"]["standardized"])
        self.assertEqual(codes(rep)["SCALE_CENTERED_PER_GENE"], "error")

    def test_mouse_case_and_excel_dates(self):
        labels = [g[0] + g[1:].lower() for g in ANCHORS] + ["1-Mar", "2-Sep"]
        rep = cem.profile_matrix(records(["gene", "S1", "S2", "S3", "S4"], log_rows(labels)))
        self.assertEqual(rep["identifiers"]["symbol_case"], "sentence")
        c = codes(rep)
        self.assertIn("SYMBOL_CASE_NOT_HUMAN", c)
        self.assertEqual(c["EXCEL_DATE_SYMBOLS"], "error")

    def test_r_layout_and_annotation_column(self):
        rows = [(g, [f"desc {g}"] + vals) for g, vals in log_rows(ANCHORS)]
        recs = [["Description", "S1", "S2", "S3", "S4"]] + [[g] + [str(v) for v in vals] for g, vals in rows]
        rep = cem.profile_matrix(recs)
        self.assertEqual(rep["layout"], "header_without_corner")
        self.assertEqual(rep["annotation_columns"], ["Description"])
        self.assertEqual(rep["shape"]["n_columns"], 4)
        self.assertIn("ANNOTATION_COLUMNS", codes(rep))

    def test_unreadable_inputs(self):
        with self.assertRaises(cem.MatrixError):
            cem.profile_matrix([])
        with self.assertRaises(cem.MatrixError):
            cem.profile_matrix([["gene", "S1", "S2"], ["ACTB", "1", "2"], ["GAPDH", "1"]])
        with self.assertRaises(cem.MatrixError):
            cem.profile_matrix([["gene", "S1"], ["ACTB", "x"], ["GAPDH", "y"]])


class TestCommandLine(unittest.TestCase):
    def test_cli_tsv_and_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            lines = ["gene\tS1\tS2\tS3"] + ["\t".join([g] + [str(v) for v in vals]) for g, vals in log_rows(ANCHORS, 3)]
            plain = Path(tmp) / "expr.tsv"
            plain.write_text("\n".join(lines) + "\n", encoding="utf-8")
            zipped = Path(tmp) / "expr.csv.gz"
            with gzip.open(zipped, "wt", encoding="utf-8") as fh:
                fh.write("\n".join(line.replace("\t", ",") for line in lines) + "\n")
            out = Path(tmp) / "report.json"
            proc = subprocess.run([sys.executable, str(SCRIPT), str(plain), "--out", str(out)],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            rep = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(rep["orientation"]["n_samples"], 3)
            proc = subprocess.run([sys.executable, str(SCRIPT), str(zipped)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(json.loads(proc.stdout)["scale"]["call"], "log")

    def test_cli_bad_input_exits_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.tsv"
            empty.write_text("", encoding="utf-8")
            proc = subprocess.run([sys.executable, str(SCRIPT), str(empty)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)
            self.assertIn("error:", proc.stderr)
            proc = subprocess.run([sys.executable, str(SCRIPT), str(Path(tmp) / "missing.tsv")],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)


if __name__ == "__main__":
    unittest.main()
