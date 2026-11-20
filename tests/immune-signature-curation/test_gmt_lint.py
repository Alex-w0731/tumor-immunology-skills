"""Unit tests for skills/immune-signature-curation/scripts/gmt_lint.py (standard library only)."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "immune-signature-curation" / "scripts" / "gmt_lint.py"
_spec = importlib.util.spec_from_file_location("gmt_lint", SCRIPT)
gl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gl)


def codes(report):
    return {i["code"] for i in report["issues"]}


def issue(report, code):
    return [i for i in report["issues"] if i["code"] == code]


class TestTokens(unittest.TestCase):
    def test_spreadsheet_forms(self):
        # Ziemann 2016: SEPT2 -> '2-Sep', MARCH1 -> '1-Mar', SEPT2 -> '2006/09/02',
        # RIKEN accession 2310009E13 -> '2.31E+13'. Abeysooriya 2021: 9.33E+22, five-digit serials.
        self.assertEqual(gl.classify_token("2-Sep")[0], "date")
        self.assertEqual(gl.classify_token("1-Mar")[0], "date")
        self.assertEqual(gl.classify_token("Mar-3")[0], "date")
        self.assertEqual(gl.classify_token("2006/09/02")[0], "date")
        self.assertEqual(gl.classify_token("2.31E+13")[0], "number")
        self.assertEqual(gl.classify_token("9.33E+22")[0], "number")
        self.assertEqual(gl.classify_token("43891")[0], "integer")
        # The unconverted RIKEN accession is not a spreadsheet number.
        self.assertEqual(gl.classify_token("2310009E13")[0], "symbol")

    def test_identifier_kinds(self):
        self.assertEqual(gl.classify_token("ENSG00000141510.17")[0], "ensembl_gene")
        self.assertEqual(gl.classify_token("NM_000546.6")[0], "refseq")
        self.assertEqual(gl.classify_token("1007_s_at")[0], "probe")
        self.assertEqual(gl.classify_token("HLA-DRA")[0], "symbol")
        self.assertEqual(gl.classify_token("CD8A/CD8B")[0], "multi")
        self.assertEqual(gl.classify_token("NA")[0], "missing")
        self.assertEqual(gl.classify_token(' "GZMB" ')[1:], ("GZMB", ["whitespace", "quoted"]))

    def test_case_style(self):
        self.assertEqual(gl.case_style("CD8A"), "upper")
        self.assertEqual(gl.case_style("C3orf18"), "upper")   # HGNC C#orf# exception
        self.assertEqual(gl.case_style("Cd8a"), "sentence")
        self.assertEqual(gl.case_style("cd8a"), "lower")
        self.assertEqual(gl.case_style("NKp46"), "mixed")

    def test_date_candidates_documented_only(self):
        self.assertEqual(gl.date_candidates("3-Sep"), [{"original": "SEPT3", "current": "SEPTIN3"}])
        self.assertEqual(gl.date_candidates("1-Mar"), [{"original": "MARCH1", "current": "MARCHF1"}])
        self.assertEqual(gl.date_candidates("Mar-8"), [{"original": "MARCH8", "current": "MARCHF8"}])
        # HGNC REST, November 2026: SEPT2 is the previous symbol of SEPTIN2 and an alias of SEPTIN6;
        # DEC1 is the previous symbol of DELEC1 and an alias of BHLHE40. Every target is listed.
        self.assertEqual(gl.date_candidates("2-Sep"), [{"original": "SEPT2", "current": "SEPTIN2 or SEPTIN6"}])
        self.assertEqual(gl.date_candidates("1-Dec"), [{"original": "DEC1", "current": "DELEC1 or BHLHE40"}])
        self.assertEqual(gl.date_candidates("20-Sep"), [])
        self.assertEqual(gl.date_candidates("CD8A"), [])


class TestSymbolMapping(unittest.TestCase):
    def test_builtin_tables(self):
        self.assertEqual(gl.map_symbol("CTGF"), ("retired", ["CCN2"]))
        self.assertEqual(gl.map_symbol("SEPT13"), ("retired", ["SEPTIN7P2"]))
        self.assertEqual(gl.map_symbol("MARCH11"), ("retired", ["MARCHF11"]))
        self.assertEqual(gl.map_symbol("PD-L1"), ("alias", ["CD274"]))
        self.assertEqual(gl.map_symbol("Tim-3"), ("alias", ["HAVCR2"]))
        # HGNC: PD-1 is an alias of PDCD1 only, PD1 of PDCD1, SNCA and SPATA2.
        self.assertEqual(gl.map_symbol("PD-1"), ("alias", ["PDCD1"]))
        self.assertEqual(gl.map_symbol("PD1")[0], "ambiguous_alias")
        self.assertEqual(gl.map_symbol("FCGR3")[0], "ambiguous_previous")
        self.assertEqual(gl.map_symbol("CD8A"), ("not_checked", ["CD8A"]))

    def test_hgnc_table(self):
        table = ("hgnc_id\tsymbol\tstatus\talias_symbol\tprev_symbol\n"
                 "HGNC:1\tCCN2\tApproved\tIGFBP8\tCTGF\n"
                 "HGNC:2\tCD274\tApproved\t\"PD-L1|PDL1|B7-H1\"\t\n"
                 "HGNC:3\tOLDGENE\tSymbol Withdrawn\t\t\n"
                 "HGNC:4\tSNCA\tApproved\tPD1\t\n"
                 "HGNC:5\tPDCD1\tApproved\tPD1|PD-1\t\n")
        h = gl.load_hgnc_table(table)
        self.assertEqual(gl.map_symbol("CTGF", h), ("retired", ["CCN2"]))
        self.assertEqual(gl.map_symbol("B7-H1", h), ("alias", ["CD274"]))
        self.assertEqual(gl.map_symbol("PD1", h), ("ambiguous_alias", ["PDCD1", "SNCA"]))
        self.assertEqual(gl.map_symbol("CD274", h), ("approved", ["CD274"]))
        self.assertEqual(gl.map_symbol("cd274", h), ("case_variant", ["CD274"]))
        self.assertEqual(gl.map_symbol("OLDGENE", h), ("unknown", []))

    def test_previous_symbol_that_is_also_an_alias_is_ambiguous(self):
        # HGNC REST, November 2026: DEC1 is the previous symbol of DELEC1 and an alias of BHLHE40, so a
        # 'retired' call (the old behaviour) would silently pick DELEC1. CD11B is a previous symbol and
        # an alias of ITGAM only, so it stays a plain rename.
        table = ("symbol\tstatus\talias_symbol\tprev_symbol\n"
                 "DELEC1\tApproved\tCTS9\tDEC1\n"
                 "BHLHE40\tApproved\tDEC1|SHARP2\tBHLHB2|STRA13\n"
                 "ITGAM\tApproved\tCD11B\tCD11B\n")
        h = gl.load_hgnc_table(table)
        self.assertEqual(gl.map_symbol("DEC1", h), ("ambiguous_previous", ["DELEC1", "BHLHE40"]))
        self.assertEqual(gl.map_symbol("CD11B", h), ("retired", ["ITGAM"]))
        self.assertEqual(gl.map_symbol("SHARP2", h), ("alias", ["BHLHE40"]))
        # The built-in table agrees with HGNC.
        self.assertEqual(gl.map_symbol("DEC1"), ("ambiguous_previous", ["DELEC1", "BHLHE40"]))
        self.assertEqual(gl.map_symbol("SEPT2"), ("ambiguous_previous", ["SEPTIN2", "SEPTIN6"]))
        self.assertEqual(gl.map_symbol("CD11B"), ("retired", ["ITGAM"]))
        r = gl.lint_gmt("T\tcurated\tDEC1\tCD3E\tCD8A\n", species="human")
        self.assertEqual(issue(r, "AMBIGUOUS_ALIAS")[0]["items"], ["DEC1 -> DELEC1 or BHLHE40"])
        self.assertNotIn("RETIRED_SYMBOL", codes(r))

    def test_hgnc_table_errors(self):
        with self.assertRaises(gl.InputError):
            gl.load_hgnc_table("")
        with self.assertRaises(gl.InputError):
            gl.load_hgnc_table("gene\tvalue\nA\t1\n")


class TestLint(unittest.TestCase):
    def test_spreadsheet_damage_is_an_error(self):
        text = "IMMUNE_X\tcurated\tCD3E\t1-Mar\t2-Sep\t2.31E+13\t43891\tCD8A\tGZMB\n"
        r = gl.lint_gmt(text)
        self.assertFalse(r["summary"]["ok"])
        self.assertTrue({"SPREADSHEET_DATE", "SPREADSHEET_NUMBER", "SPREADSHEET_SERIAL"} <= codes(r))
        items = issue(r, "SPREADSHEET_DATE")[0]["items"]
        self.assertIn("1-Mar: was MARCH1 (now MARCHF1)?", items)
        self.assertIn("2-Sep: was SEPT2 (now SEPTIN2 or SEPTIN6)?", items)
        # Golden serial: 2021-01-01 is serial 44197 in the 1900 date system, so 2020-01-01 is
        # 44197 - 366 = 43831; + 31 (January) + 29 (February 2020) = 43891 = 2020-03-01 = '1-Mar'.
        self.assertEqual(issue(r, "SPREADSHEET_SERIAL")[0]["items"],
                         ["43891 (2020-03-01): was MARCH1 (now MARCHF1)?"])
        self.assertEqual(gl.serial_note("43525"), "43525 (2019-03-01): was MARCH1 (now MARCHF1)?")
        # 43831 + 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 (January to August 2020) = 44075 = 2020-09-01.
        self.assertEqual(gl.serial_note("44075"), "44075 (2020-09-01): was SEPT1 (now SEPTIN1)?")
        self.assertEqual(gl.serial_note("45000"), "45000 (2023-03-15): origin unknown")

    def test_duplicate_set_names_and_genes(self):
        text = ("T_cells\tcurated\tCD3D\tCD3E\tCD3E\tCD3G\n"
                "T_cells\tcurated\tCD3D\tCD3E\n"
                "t_cells\tcurated\tCD2\tCD5\n")
        r = gl.lint_gmt(text)
        dup = issue(r, "DUPLICATE_SET_NAME")[0]
        self.assertEqual(dup["severity"], "error")
        self.assertEqual(dup["items"], ["1", "2"])
        self.assertEqual(issue(r, "DUPLICATE_GENE")[0]["items"], ["CD3E"])
        self.assertIn("SET_NAME_CASE_DUPLICATE", codes(r))

    def test_retired_alias_and_clash(self):
        text = "Stroma\t10.1038/nature25501\tCTGF\tCCN2\tACTA2\tPD-L1\tPD1\n"
        r = gl.lint_gmt(text, species="human")
        self.assertEqual(issue(r, "RETIRED_SYMBOL")[0]["items"], ["CTGF -> CCN2"])
        self.assertEqual(issue(r, "ALIAS_SYMBOL")[0]["items"], ["PD-L1 -> CD274"])
        self.assertEqual(issue(r, "AMBIGUOUS_ALIAS")[0]["items"], ["PD1 -> PDCD1 or SNCA or SPATA2"])
        self.assertEqual(issue(r, "DUPLICATE_AFTER_MAPPING")[0]["items"], ["CCN2"])
        self.assertEqual(r["sets"][0]["provenance"], "doi")

    def test_case_problems_in_human_file(self):
        # 3 of 34 symbols are not upper case: 91% upper, so the file is called human-style (>= 90%).
        text = ("A\tcurated\tCD8A\tCd8a\tGZMB\tPRF1\tNKG7\tGNLY\tCD3E\tCD3D\tCD2\tCD5\n"
                "B\tcurated\tcd274\tNKp46\tIFNG\tSTAT1\tCXCL9\tCXCL10\tIDO1\tHLA-DRA\tGBP1\tIRF1\n"
                "C\tcurated\tCD19\tMS4A1\tCD79A\tCD79B\tBLK\tFCRL2\tTCL1A\tSPIB\tPNOC\tTNFRSF17\n")
        r = gl.lint_gmt(text)
        self.assertEqual(r["identifiers"]["case_style"], "human-style")
        self.assertEqual(r["identifiers"]["species_assumed"], "human")
        self.assertEqual(issue(r, "CASE_DUPLICATE")[0]["items"], ["CD8A/Cd8a"])
        self.assertEqual(issue(r, "LOWERCASE_SYMBOLS")[0]["items"], ["cd274"])
        mism = {s for i in issue(r, "SPECIES_CASE_MISMATCH") for s in i["items"]}
        self.assertEqual(mism, {"Cd8a"})          # NKp46 is reported as an alias instead
        self.assertIn("NKp46 -> NCR1", issue(r, "ALIAS_SYMBOL")[0]["items"])

    def test_mouse_file(self):
        # H2-K1 does not vote; 10 of the other 11 symbols are sentence case (91%): mouse-style.
        text = ("Mouse_T\tcurated\tCd3e\tCd8a\tGzmb\tPrf1\tNkg7\tH2-K1\tCD274\tIfng\tStat1\tCxcl9"
                "\tCxcl10\tIl2ra\n")
        r = gl.lint_gmt(text)
        self.assertEqual(r["identifiers"]["species_assumed"], "mouse")
        self.assertEqual(r["settings"]["symbol_check"], "skipped_mouse")
        self.assertEqual(issue(r, "SPECIES_CASE_MISMATCH")[0]["items"], ["CD274"])   # H2-K1 is mouse MHC

    def test_mouse_mhc_class_ii(self):
        # MGI class II symbols (H2-Aa, H2-Ab1, H2-Eb1) are neither sentence nor upper case. Before the
        # fix they voted 'mixed', so this mouse antigen-presentation set was called species 'unknown'.
        r = gl.lint_gmt("Mouse_MHCII\tcurated\tCd74\tCiita\tH2-Aa\tH2-Ab1\tH2-Eb1\tCd86\n")
        self.assertEqual(r["identifiers"]["case_style"], "mouse-style")
        self.assertEqual(r["identifiers"]["species_assumed"], "mouse")
        self.assertNotIn("SPECIES_CASE_MISMATCH", codes(r))
        # In a human file every H2- symbol is flagged, including the all-upper-case H2-K1.
        r = gl.lint_gmt("Human_APC\tcurated\tCD74\tCIITA\tH2-K1\tH2-Ab1\tHLA-DRA\n", species="human")
        self.assertEqual(issue(r, "SPECIES_CASE_MISMATCH")[0]["items"], ["H2-Ab1", "H2-K1"])

    def test_sizes(self):
        big = "\t".join(f"G{i}" for i in range(1, 502))
        text = ("One\tcurated\tFOXP3\n"
                "Small\tcurated\tCD8A\tCD8B\n"
                f"Big\tcurated\t{big}\n")
        r = gl.lint_gmt(text)
        self.assertIn("SINGLE_GENE_SET", codes(r))
        self.assertEqual(issue(r, "SMALL_SET")[0]["severity"], "info")
        self.assertIn("LARGE_SET", codes(r))
        self.assertEqual(r["sizes"], {"min": 1, "median": 2, "max": 501, "n_single_gene": 1,
                                      "n_below_min": 1, "n_above_max": 1})

    def test_format_problems(self):
        text = ("﻿# exported from a spreadsheet\tx\ty\r\n"
                "\r\n"
                "Bad line with spaces CD3E CD3D\r\n"
                "NoGenes\tcurated\r\n"
                "Cytotoxic\tGZMA\tPRF1\tGZMB\tNKG7\r\n")
        r = gl.lint_gmt(text)
        c = codes(r)
        self.assertTrue({"BOM", "CRLF_LINE_ENDINGS", "BLANK_LINES", "COMMENT_LINE",
                         "WRONG_DELIMITER", "TOO_FEW_FIELDS", "DESCRIPTION_LOOKS_LIKE_GENE"} <= c)
        self.assertEqual(r["format"]["n_sets"], 1)

    def test_identical_sets_and_provenance(self):
        text = ("CYT\t10.1016/j.cell.2014.12.033\tGZMA\tPRF1\n"
                "Cytolytic\tpublished\tPRF1\tGZMA\n"
                "Mine\t\tGZMA\tGZMB\n")
        r = gl.lint_gmt(text)
        self.assertEqual(issue(r, "IDENTICAL_SETS")[0]["items"], ["CYT", "Cytolytic"])
        self.assertEqual(r["provenance"], {"doi": 1, "empty": 1, "published": 1})
        self.assertIn("PROVENANCE_WITHOUT_REFERENCE", codes(r))
        self.assertIn("NO_PROVENANCE", codes(r))

    def test_entrez_file_has_no_serial_errors(self):
        text = "Set\tcurated\t925\t926\t3002\t5551\t10563\t11005\n"
        r = gl.lint_gmt(text)
        self.assertEqual(r["identifiers"]["file_type"], "entrez")
        self.assertNotIn("SPREADSHEET_SERIAL", codes(r))

    def test_empty_input_and_bad_arguments(self):
        r = gl.lint_gmt("")
        self.assertEqual(issue(r, "NO_SETS")[0]["severity"], "error")
        self.assertEqual(r["sizes"], {})
        with self.assertRaises(ValueError):
            gl.lint_gmt("A\tx\tCD3E\n", min_size=20, max_size=5)
        with self.assertRaises(ValueError):
            gl.lint_gmt("A\tx\tCD3E\n", species="rat")


class TestCommandLine(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)

    def test_cli_strict_and_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = os.path.join(tmp, "good.gmt")
            bad = os.path.join(tmp, "bad.gmt")
            Path(good).write_text("CYT\t10.1016/j.cell.2014.12.033\tGZMA\tPRF1\n", encoding="utf-8")
            Path(bad).write_text("X\tcurated\t1-Mar\tCD3E\n", encoding="utf-8")
            p = self.run_cli(good)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(json.loads(p.stdout)["summary"]["ok"])
            self.assertEqual(self.run_cli(bad).returncode, 0)
            self.assertEqual(self.run_cli(bad, "--strict").returncode, 1)
            out = os.path.join(tmp, "r.json")
            self.assertEqual(self.run_cli(bad, "--out", out).returncode, 0)
            self.assertEqual(json.loads(Path(out).read_text())["file"], "bad.gmt")
            missing = self.run_cli(os.path.join(tmp, "nope.gmt"))
            self.assertEqual(missing.returncode, 2)
            self.assertIn("no such file", missing.stderr)
            unwritable = self.run_cli(good, "--out", os.path.join(tmp, "no_such_dir", "r.json"))
            self.assertEqual(unwritable.returncode, 2)
            self.assertIn("cannot write", unwritable.stderr)
            self.assertEqual(self.run_cli(good, "--min-size", "0").returncode, 2)

    def test_help(self):
        p = self.run_cli("--help")
        self.assertEqual(p.returncode, 0)
        self.assertIn("usage", p.stdout.lower())


if __name__ == "__main__":
    unittest.main()
