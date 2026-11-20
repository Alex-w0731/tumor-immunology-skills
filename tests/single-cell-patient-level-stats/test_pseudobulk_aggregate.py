"""Unit tests for skills/single-cell-patient-level-stats/scripts/pseudobulk_aggregate.py."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "skills" / "single-cell-patient-level-stats" / "scripts" / "pseudobulk_aggregate.py"
_spec = importlib.util.spec_from_file_location("pseudobulk_aggregate_under_test", SCRIPT)
pa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pa)

# Synthetic per-cell table. Hand sums: S1::T = (1+2, 0+5) = (3, 5) from 2 cells; S1::B = (0, 1)
# from 1 cell; S2::T = (3+0+4, 3+0+1) = (7, 4) from 3 cells.
WIDE = """cell_id\tsample_id\tcell_type\tcondition\tGENE_A\tGENE_B
c1\tS1\tT\tctrl\t1\t0
c2\tS1\tT\tctrl\t2\t5
c3\tS1\tB\tctrl\t0\t1
c4\tS2\tT\ttreat\t3\t3
c5\tS2\tT\ttreat\t0\t0
c6\tS2\tT\ttreat\t4\t1
"""

# The same cells in long format. Zero rows are omitted except one for c5, which would
# otherwise not be counted (a cell with no row is invisible in a long table).
LONG = """sample_id\tcell_type\tcell_id\tcondition\tgene\tcount
S1\tT\tc1\tctrl\tGENE_A\t1
S1\tT\tc2\tctrl\tGENE_A\t2
S1\tT\tc2\tctrl\tGENE_B\t5
S1\tB\tc3\tctrl\tGENE_B\t1
S2\tT\tc4\ttreat\tGENE_A\t3
S2\tT\tc4\ttreat\tGENE_B\t3
S2\tT\tc5\ttreat\tGENE_A\t0
S2\tT\tc6\ttreat\tGENE_A\t4
S2\tT\tc6\ttreat\tGENE_B\t1
"""


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def aggregate(text, fmt="wide", **kw):
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "in.tsv")
        write(path, text)
        fh, header, reader = pa.open_table(path)
        with fh:
            if fmt == "wide":
                return pa.aggregate_wide(header, reader, "sample_id", "cell_type", kw.get("cell_col", "cell_id"),
                                         kw.get("label_cols", ("condition",)), kw.get("exclude_cols", ()),
                                         kw.get("required", ("condition",)))
            return pa.aggregate_long(header, reader, "sample_id", "cell_type", kw.get("cell_col", "cell_id"),
                                     "gene", "count", ("condition",), ("condition",))


class TestAggregation(unittest.TestCase):
    def test_wide_golden_sums_and_cell_numbers(self):
        genes, sums, n_cells, labels, n_input = aggregate(WIDE)
        self.assertEqual(genes, ["GENE_A", "GENE_B"])
        self.assertEqual(sums[("S1", "T")], [3, 5])
        self.assertEqual(sums[("S1", "B")], [0, 1])
        self.assertEqual(sums[("S2", "T")], [7, 4])
        self.assertEqual(n_cells, {("S1", "T"): 2, ("S1", "B"): 1, ("S2", "T"): 3})
        self.assertEqual(n_input, 6)
        self.assertEqual(labels["S2"], {"condition": "treat"})

    def test_long_matches_wide(self):
        wide = aggregate(WIDE)
        long_ = aggregate(LONG, fmt="long")
        self.assertEqual(wide[0], long_[0])
        self.assertEqual(wide[1], long_[1])
        self.assertEqual(wide[2], long_[2])

    def test_min_cells_filter_and_summary(self):
        genes, sums, n_cells, labels, n_input = aggregate(WIDE)
        meta, kept, summary = pa.summarise(genes, sums, n_cells, labels, n_input, min_cells=2, min_units=1,
                                           group_col="condition", label_cols=("condition",))
        self.assertEqual(kept, ["S1::T", "S2::T"])
        self.assertEqual(summary["dropped"], [{"sample": "S1", "cell_type": "B", "n_cells": 1}])
        b_row = next(r for r in meta if r["pseudobulk_id"] == "S1::B")
        self.assertEqual((b_row["kept"], b_row["total_counts"]), ("no", 1))
        t_entry = next(e for e in summary["by_cell_type"] if e["cell_type"] == "T")
        self.assertEqual(t_entry["units_per_group"], {"ctrl": 1, "treat": 1})
        self.assertTrue(t_entry["testable"])
        b_entry = next(e for e in summary["by_cell_type"] if e["cell_type"] == "B")
        self.assertFalse(b_entry["testable"])

    def test_repeated_patient_counted_once(self):
        # Synthetic: patient P3 has two samples (H3 primary, H3m metastasis) in group HR.
        text = "sample_id\tcell_type\tpatient\tsubtype\tG1\n"
        for s, p, g in (("T1", "P1", "TN"), ("T2", "P2", "TN"), ("H3", "P3", "HR"), ("H3m", "P3", "HR")):
            text += f"{s}\tTreg\t{p}\t{g}\t2\n" * 3
        genes, sums, n_cells, labels, n_input = aggregate(text, cell_col=None, label_cols=("subtype", "patient"),
                                                          required=("subtype", "patient"))
        _, _, summary = pa.summarise(genes, sums, n_cells, labels, n_input, min_cells=3, min_units=2,
                                     group_col="subtype", patient_col="patient", label_cols=("subtype", "patient"))
        entry = summary["by_cell_type"][0]
        self.assertEqual(entry["units_per_group"], {"HR": 1, "TN": 2})
        self.assertEqual(entry["repeated_units"], [{"unit": "P3", "group": "HR", "samples": ["H3", "H3m"]}])
        self.assertFalse(entry["testable"])
        self.assertTrue(summary["notes"])

    def test_non_integer_count_is_rejected_with_column_name(self):
        bad = WIDE.replace("c2\tS1\tT\tctrl\t2\t5", "c2\tS1\tT\tctrl\t2\t1.37")
        with self.assertRaisesRegex(pa.InputError, "GENE_B"):
            aggregate(bad)
        with self.assertRaisesRegex(pa.InputError, "non-negative integer"):
            aggregate(WIDE.replace("\t3\t3\n", "\t-3\t3\n"))

    def test_inconsistent_label_and_duplicate_cell(self):
        with self.assertRaisesRegex(pa.InputError, "constant within a sample"):
            aggregate(WIDE.replace("c5\tS2\tT\ttreat", "c5\tS2\tT\tctrl"))
        with self.assertRaisesRegex(pa.InputError, "more than once"):
            aggregate(WIDE.replace("c6\tS2", "c5\tS2"))

    def test_long_format_needs_cell_column_and_unique_cell_gene(self):
        with self.assertRaisesRegex(pa.InputError, "--cell-col"):
            aggregate(LONG, fmt="long", cell_col=None)
        with self.assertRaisesRegex(pa.InputError, "more than one row"):
            aggregate(LONG + "S2\tT\tc6\ttreat\tGENE_B\t2\n", fmt="long")

    def test_integer_qc_columns_are_not_summed_as_genes(self):
        # nFeature_RNA holds integers, so the count check alone would let it through as a gene.
        text = "sample_id\tcell_type\tnFeature_RNA\tG1\nS1\tT\t1200\t3\nS1\tT\t900\t1\n"
        with self.assertRaisesRegex(pa.InputError, "nFeature_RNA.*--exclude-cols"):
            aggregate(text, cell_col=None, label_cols=(), required=())
        genes, sums, _, _, _ = aggregate(text, cell_col=None, label_cols=(), required=(),
                                         exclude_cols=("nFeature_RNA",))
        self.assertEqual((genes, sums[("S1", "T")]), (["G1"], [4]))
        self.assertEqual(pa.metadata_like(["CD8A", "total_counts_mt", "pct_counts_mt", "Leiden"]),
                         ["total_counts_mt", "pct_counts_mt", "Leiden"])

    def test_metadata_only_table(self):
        genes, sums, n_cells, _, _ = aggregate("sample_id\tcell_type\tcondition\nS1\tT\tc\nS1\tT\tc\n", cell_col=None)
        self.assertEqual(genes, [])
        self.assertEqual(n_cells[("S1", "T")], 2)


class TestCommandLine(unittest.TestCase):
    def test_cli_writes_matrix_meta_and_summary(self):
        with tempfile.TemporaryDirectory() as d:
            src, counts, meta = (os.path.join(d, n) for n in ("cells.tsv", "pb.tsv", "meta.tsv"))
            write(src, WIDE)
            proc = subprocess.run([sys.executable, str(SCRIPT), src, "--cell-col", "cell_id", "--group-col",
                                   "condition", "--min-cells", "2", "--min-units", "1", "--out", counts,
                                   "--meta-out", meta], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["pseudobulks"], {"total": 3, "kept": 2, "dropped": 1})
            with open(counts, encoding="utf-8") as fh:
                lines = fh.read().splitlines()
            self.assertEqual(lines, ["gene\tS1::T\tS2::T", "GENE_A\t3\t7", "GENE_B\t5\t4"])
            with open(meta, encoding="utf-8") as fh:
                self.assertEqual(fh.readline().strip().split("\t"),
                                 ["pseudobulk_id", "sample", "cell_type", "n_cells", "total_counts", "kept", "condition"])
            bad = subprocess.run([sys.executable, str(SCRIPT), src, "--sample-col", "donor"],
                                 capture_output=True, text=True)
            self.assertEqual(bad.returncode, 2)
            self.assertIn("donor", bad.stderr)
            write(src, "")
            empty = subprocess.run([sys.executable, str(SCRIPT), src], capture_output=True, text=True)
            self.assertEqual(empty.returncode, 2)

    def test_cli_on_treg_eval_input(self):
        # Synthetic eval input. Counted independently with awk: cells per sample H1 15, H2 4, H3 12,
        # H4 10, H4m 11, T1 14, T2 12, T3 6, T4 11; FOXP3 sums H1 23, H4m 31, T4 23. At 10 cells H2
        # and T3 drop; H4 (exactly 10) stays; H4 and H4m are one patient (P8), so 3 vs 3 patients.
        src = ROOT / "skills" / "single-cell-patient-level-stats" / "evals" / "treg-pseudobulk-subtype-zh" / "files" / "cells.tsv"
        with tempfile.TemporaryDirectory() as d:
            counts = os.path.join(d, "pb.tsv")
            proc = subprocess.run([sys.executable, str(SCRIPT), str(src), "--cell-col", "cell_id", "--group-col",
                                   "subtype", "--patient-col", "patient_id", "--covariates", "tissue",
                                   "--exclude-cols", "percent_mt,data_status", "--out", counts],
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            summary = json.loads(proc.stdout)
            with open(counts, encoding="utf-8") as fh:
                rows = {line.split("\t")[0]: line.rstrip("\n").split("\t")[1:] for line in fh}
        self.assertEqual(summary["dropped"], [{"sample": "H2", "cell_type": "Treg", "n_cells": 4},
                                              {"sample": "T3", "cell_type": "Treg", "n_cells": 6}])
        entry = summary["by_cell_type"][0]
        self.assertEqual(entry["units_per_group"], {"HR+": 3, "TNBC": 3})
        self.assertEqual(entry["repeated_units"], [{"unit": "P8", "group": "HR+", "samples": ["H4", "H4m"]}])
        self.assertFalse(entry["testable"])
        foxp3 = dict(zip(rows["gene"], rows["FOXP3"]))
        self.assertEqual((foxp3["H1::Treg"], foxp3["H4m::Treg"], foxp3["T4::Treg"]), ("23", "31", "23"))
        self.assertEqual(summary["first_genes"][:2], ["FOXP3", "IL2RA"])


if __name__ == "__main__":
    unittest.main()
