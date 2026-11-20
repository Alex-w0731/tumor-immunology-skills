#!/usr/bin/env python3
"""Check a bulk expression matrix before immune-signature scoring. Standard library only.

Reads a delimited text matrix (TSV or CSV, optionally gzip-compressed). The first row holds the
column identifiers and the first column holds the row identifiers; a header that is one field
shorter than the data rows (the R `write.table` layout) is accepted. Leading lines that start
with '#' are skipped. The JSON report covers:

  orientation   genes x samples or samples x genes
  identifiers   gene symbol, Ensembl, RefSeq, probe, integer or composite; symbol case;
                Ensembl version suffixes; symbols turned into dates by a spreadsheet
  scale         linear, log or centered (heuristic below)
  duplicates    repeated gene and sample identifiers, exact and after upper-casing
  values        missing, non-numeric, infinite, negative and zero values; constant genes
  issues        coded findings, each with a severity and the action it calls for

Usage examples:

  python3 check_expression_matrix.py expression.tsv
  python3 check_expression_matrix.py cohort.csv.gz --orientation genes-by-samples --out check.json

Orientation heuristic [heuristic]. An axis "holds genes" when at least half of its labels match
a gene-identifier pattern (Ensembl gene or transcript, RefSeq, array probe) or when at least 3
labels are anchor genes (48 widely measured genes such as ACTB, GAPDH, PTPRC and CD8A; approved
symbols at HGNC as of November 2026 [HGNC2026]; matching ignores case). If exactly one axis
holds genes the call has confidence "high". Otherwise the shape decides with confidence "low":
an axis at least 5 times longer than the other is taken as the gene axis.

Scale heuristic. Quantiles are linear-interpolation quantiles (the numpy default) of the finite
values, taken from every k-th row when the matrix holds more than --max-values values. The
rules are applied in this order:

  1. centered  negative values are present and at least 90% of the non-constant genes (or
               samples) have a mean within 0.05 standard deviations of zero       [heuristic]
  2. linear    the GEO2R rule holds: q99 > 100, or (max - min > 50 and q25 > 0). This is the
               rule that tmelens 0.1.0 applies with log="auto" and documents as the NCBI GEO2R
               rule [Wang2026]; a convention, not a validated classifier
  3. linear    max > 50; 2**50 is about 1e15, so no log2 expression value is that large
                                                                                  [heuristic]
  4. linear    every finite value is a non-negative integer and max > 1 (counts)  [heuristic]
  5. log       otherwise

Rule 2 misses linear data in which at least a quarter of the values are zero and the 99th
percentile is at most 100 (sparse or low-depth matrices, small panels). Rules 3 and 4 catch
most of these, and the report then says that tmelens log="auto" would leave the matrix
untransformed. Non-integer linear data whose values all lie below 50 are called "log": confirm
the scale from the data source whenever the call matters.

Exit status: 0 when a report was produced (read its "issues" list), 2 when the input could not
be read as a matrix.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from array import array
from pathlib import Path

GEO2R_Q99 = 100.0            # convention: tmelens.needs_log_transform, "GEO2R rule" [Wang2026]
GEO2R_RANGE = 50.0           # same convention
LOG_BOUND = 50.0             # heuristic: no log2 expression value exceeds 50
CENTERED_FRACTION = 0.9      # heuristic
CENTERED_TOLERANCE = 0.05    # heuristic: |mean| <= 0.05 SD
UNIT_SD_RANGE = (0.9, 1.1)   # heuristic: standardised features
PATTERN_MIN_FRACTION = 0.5   # heuristic
ANCHOR_MIN_HITS = 3          # heuristic
SHAPE_RATIO = 5              # heuristic
TYPE_MIN_FRACTION = 0.8      # heuristic
FEW_SAMPLES = 10             # heuristic
FEW_GENES = 5000             # heuristic
PER_MILLION = 1e6            # TPM and CPM sum to one million per sample by definition
PER_MILLION_TOLERANCE = 0.01
DEFAULT_MAX_VALUES = 2_000_000

MISSING_TOKENS = frozenset({"", "na", "n/a", "#n/a", "nan", "null", "none", "."})

# Approved HGNC symbols (checked at genenames.org, November 2026): widely measured
# housekeeping, tumour and immune genes, used only to recognise which axis holds genes.
ANCHOR_GENES = frozenset({
    "ACTB", "GAPDH", "B2M", "TUBB", "RPLP0", "RPL13A", "HPRT1", "TBP", "PGK1", "PPIA", "UBC",
    "YWHAZ", "SDHA", "TFRC", "GUSB", "HMBS", "EEF1A1", "RPS18", "POLR2A", "ALDOA",
    "TP53", "MYC", "EGFR", "KRAS", "PTEN", "VIM", "EPCAM", "KRT18", "MKI67", "STAT1", "COL1A1",
    "PTPRC", "CD3E", "CD4", "CD8A", "CD19", "CD68", "MS4A1", "GZMB", "PRF1", "FOXP3", "IFNG",
    "CD274", "PDCD1", "CTLA4", "HLA-A", "HLA-B", "HLA-DRA",
})

RE_ENSEMBL_GENE = re.compile(r"^ENS([A-Z]{0,6})G\d{11}(\.\d+)?(_PAR_Y)?$")
RE_ENSEMBL_TRANSCRIPT = re.compile(r"^ENS[A-Z]{0,6}T\d{11}(\.\d+)?(_PAR_Y)?$")
RE_REFSEQ = re.compile(r"^[NX][MR]_\d+(\.\d+)?$")
RE_PROBE = re.compile(r"^(AFFX-.+|.+_(at|st)|ILMN_\d+|A_\d{2}_P\d+|TC\d+\.hg\.\d+)$")
RE_INTEGER = re.compile(r"^\d+$")
RE_EXCEL_DATE = re.compile(
    r"^\d{1,2}-(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)$", re.IGNORECASE
)
GENE_PATTERN_CLASSES = ("ensembl_gene", "ensembl_transcript", "refseq", "probe")
DELIMITERS = {"tab": "\t", "comma": ",", "semicolon": ";"}
TALLY_KEYS = ("all_missing", "constant", "all_zero", "evaluable", "centered", "unit_sd", "negative")


class MatrixError(ValueError):
    """The input cannot be read as a labelled numeric matrix."""


# --------------------------------------------------------------------------- pure helpers

def quantile(sorted_values, q: float) -> float:
    """Linear-interpolation quantile (numpy default, 'type 7') of an ascending sequence."""
    n = len(sorted_values)
    if n == 0:
        raise ValueError("quantile of an empty sequence")
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must lie between 0 and 1")
    h = (n - 1) * q
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    return float(sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo]))


def parse_value(token: str):
    """Classify one cell: ('number', x), ('missing', None), ('infinite', None) or ('text', None)."""
    try:
        x = float(token)
    except ValueError:
        return ("missing", None) if token.strip().lower() in MISSING_TOKENS else ("text", None)
    if x != x:
        return "missing", None
    if x in (math.inf, -math.inf):
        return "infinite", None
    if "_" in token:                   # float() accepts '1_000'; a data file never means that
        return "text", None
    return "number", x


def classify_identifier(label: str) -> str:
    """Return the identifier class of one row or column label."""
    t = label.strip()
    if not t:
        return "empty"
    if RE_ENSEMBL_GENE.match(t):
        return "ensembl_gene"
    if RE_ENSEMBL_TRANSCRIPT.match(t):
        return "ensembl_transcript"
    if RE_REFSEQ.match(t):
        return "refseq"
    if RE_EXCEL_DATE.match(t):
        return "excel_date"
    if RE_PROBE.match(t):
        return "probe"
    if "|" in t:
        return "composite"
    if RE_INTEGER.match(t):
        return "integer"
    return "other"


def anchor_case(label: str):
    """Return 'upper', 'sentence' or 'other' if the label names an anchor gene, else None.

    Composite labels such as 'ACTB|60' are split on '|' and each part is tested.
    """
    for part in label.split("|"):
        p = part.strip()
        if p.upper() in ANCHOR_GENES:
            if p == p.upper():
                return "upper"
            if p[:1].isupper() and p[1:] == p[1:].lower():
                return "sentence"
            return "other"
    return None


def summarise_labels(labels, examples: int = 5) -> dict:
    """Count identifier classes, anchor genes and suspicious labels on one axis."""
    counts: dict = {}
    case_counts = {"upper": 0, "sentence": 0, "other": 0}
    excel, versioned, prefixes = [], 0, {}
    for label in labels:
        cls = classify_identifier(label)
        counts[cls] = counts.get(cls, 0) + 1
        if cls == "excel_date":
            excel.append(label.strip())
        elif cls == "ensembl_gene":
            m = RE_ENSEMBL_GENE.match(label.strip())
            prefix = "ENS" + m.group(1) + "G"
            prefixes[prefix] = prefixes.get(prefix, 0) + 1
            if m.group(2):
                versioned += 1
        case = anchor_case(label)
        if case:
            case_counts[case] += 1
    n = len(labels)
    pattern_n = sum(counts.get(c, 0) for c in GENE_PATTERN_CLASSES)
    anchors = sum(case_counts.values())
    return {
        "n": n,
        "class_counts": dict(sorted(counts.items())),
        "pattern_fraction": pattern_n / n if n else 0.0,
        "anchor_hits": anchors,
        "anchor_case": case_counts,
        "holds_genes": bool(n) and (pattern_n / n >= PATTERN_MIN_FRACTION or anchors >= ANCHOR_MIN_HITS),
        "excel_date": {"n": len(excel), "examples": excel[:examples]},
        "ensembl_versioned": versioned,
        "ensembl_prefix": max(sorted(prefixes), key=prefixes.get) if prefixes else None,
        "examples": [str(x).strip() for x in labels[:examples]],
    }


def call_orientation(rows: dict, cols: dict, forced: str = "auto") -> dict:
    """Decide which axis holds genes from two `summarise_labels` results."""
    n_rows, n_cols = rows["n"], cols["n"]
    if forced == "genes-by-samples":
        call, confidence, basis = "genes_x_samples", "high", "stated by the user"
    elif forced == "samples-by-genes":
        call, confidence, basis = "samples_x_genes", "high", "stated by the user"
    elif rows["holds_genes"] and not cols["holds_genes"]:
        call, confidence, basis = "genes_x_samples", "high", "identifiers"
    elif cols["holds_genes"] and not rows["holds_genes"]:
        call, confidence, basis = "samples_x_genes", "high", "identifiers"
    elif rows["holds_genes"] and cols["holds_genes"]:
        call, confidence, basis = "undetermined", "none", "both axes hold gene identifiers"
    elif n_rows >= SHAPE_RATIO * n_cols:
        call, confidence, basis = "genes_x_samples", "low", "shape"
    elif n_cols >= SHAPE_RATIO * n_rows:
        call, confidence, basis = "samples_x_genes", "low", "shape"
    else:
        call, confidence, basis = "undetermined", "none", "no gene identifiers recognised"
    if call == "samples_x_genes":
        n_genes, n_samples = n_cols, n_rows
    elif call == "genes_x_samples":
        n_genes, n_samples = n_rows, n_cols
    else:
        n_genes = n_samples = None
    return {
        "call": call,
        "confidence": confidence,
        "basis": basis,
        "n_genes": n_genes,
        "n_samples": n_samples,
        "row_labels": {"pattern_fraction": round(rows["pattern_fraction"], 4), "anchor_genes": rows["anchor_hits"]},
        "column_labels": {"pattern_fraction": round(cols["pattern_fraction"], 4), "anchor_genes": cols["anchor_hits"]},
    }


def identifier_type(summary: dict) -> str:
    """Name the identifier type of the gene axis from a `summarise_labels` result."""
    n = summary["n"]
    if n == 0:
        return "unknown"
    counts = summary["class_counts"]
    for cls, name in (
        ("ensembl_gene", "ensembl_gene_id"),
        ("ensembl_transcript", "ensembl_transcript_id"),
        ("refseq", "refseq_id"),
        ("probe", "probe_id"),
        ("integer", "integer_id"),
        ("composite", "composite_id"),
    ):
        if counts.get(cls, 0) / n >= TYPE_MIN_FRACTION:
            return name
    symbol_like = counts.get("other", 0) + counts.get("excel_date", 0)
    if symbol_like / n >= TYPE_MIN_FRACTION:
        return "gene_symbol" if summary["anchor_hits"] >= ANCHOR_MIN_HITS else "unrecognised"
    return "mixed"


def symbol_case(summary: dict):
    """'upper', 'sentence' or 'mixed' for the anchor genes found on an axis; None if none found."""
    case = summary["anchor_case"]
    total = sum(case.values())
    if total == 0:
        return None
    if case["upper"] / total >= TYPE_MIN_FRACTION:
        return "upper"
    if case["sentence"] / total >= TYPE_MIN_FRACTION:
        return "sentence"
    return "mixed"


def find_duplicates(labels, examples: int = 5) -> dict:
    """Identifiers that occur more than once after stripping, and those that only collide
    once upper-cased."""
    exact: dict = {}
    folded: dict = {}
    for label in labels:
        t = label.strip()
        exact[t] = exact.get(t, 0) + 1
        folded[t.upper()] = folded.get(t.upper(), 0) + 1
    dup = {k: v for k, v in exact.items() if v > 1 and k}
    exact_upper = {k.upper() for k in dup}
    case_only = sorted(k for k, v in folded.items() if v > 1 and k and k not in exact_upper)
    return {
        "n_duplicated_ids": len(dup),
        "n_entries_involved": sum(dup.values()),
        "examples": [{"id": k, "count": dup[k]} for k in sorted(dup)[:examples]],
        "n_more_after_upper_casing": len(case_only),
        "examples_after_upper_casing": case_only[:examples],
        "n_empty": exact.get("", 0),
    }


def geo2r_linear(q: dict) -> bool:
    """The GEO2R log2 auto-detection rule, as implemented by tmelens.needs_log_transform."""
    return bool(q["q99"] > GEO2R_Q99 or (q["max"] - q["min"] > GEO2R_RANGE and q["q25"] > 0))


def infer_scale(q: dict, integer_valued: bool, centered_axis=None) -> dict:
    """Apply the documented scale rules to a quantile summary (keys: min q25 median q75 q99 max).

    `integer_valued` says whether every finite value is an integer; `centered_axis` names the
    axis ('gene' or 'sample') whose means sit at zero, or None.
    """
    reasons = []
    rule = geo2r_linear(q)
    above = q["max"] > LOG_BOUND
    negative = q["min"] < 0
    counts = bool(integer_valued and not negative and q["max"] > 1)
    if q["max"] == q["min"]:
        call = "undetermined"
        reasons.append("all finite values are identical")
    elif negative and centered_axis:
        call = "centered"
        reasons.append(f"negative values, and {centered_axis} means sit at zero")
    elif rule or above or counts:
        call = "linear"
        if rule:
            reasons.append("GEO2R rule met (q99 > 100, or range > 50 with q25 > 0)")
        if above:
            reasons.append("maximum above 50, too large for a log2 expression value")
        if counts:
            reasons.append("all values are non-negative integers")
    else:
        call = "log"
        reasons.append("GEO2R rule not met, maximum at most 50, values not all integers")
    return {
        "call": call,
        "geo2r_rule_linear": rule,
        "max_above_log_bound": above,
        "integer_valued": counts,
        "tmelens_auto_log_would_transform": rule,
        "reasons": reasons,
    }


def feature_flags(n: int, low: float, high: float, mean: float, sd: float) -> dict:
    """Describe one gene or sample from its finite-value count, range, mean and sample SD."""
    constant = n >= 2 and low == high
    evaluable = n >= 2 and not constant
    return {
        "all_missing": n == 0,
        "constant": constant,
        "all_zero": constant and low == 0.0,
        "evaluable": evaluable,
        "centered": evaluable and abs(mean) <= CENTERED_TOLERANCE * sd,
        "unit_sd": evaluable and UNIT_SD_RANGE[0] <= sd <= UNIT_SD_RANGE[1],
        "negative": n > 0 and low < 0,
    }


def _is_centered(tally: dict) -> bool:
    return tally["evaluable"] > 0 and tally["centered"] >= CENTERED_FRACTION * tally["evaluable"]


def _is_standardized(tally: dict) -> bool:
    return tally["evaluable"] > 0 and tally["unit_sd"] >= CENTERED_FRACTION * tally["evaluable"]


class _Columns:
    """Running count, mean, sum of squared deviations (Welford) and range for every column."""

    def __init__(self, n: int) -> None:
        self.count = [0] * n
        self.mean = [0.0] * n
        self.m2 = [0.0] * n
        self.low = [math.inf] * n
        self.high = [-math.inf] * n
        self.missing = [0] * n
        self.text = [0] * n

    def add(self, j: int, x: float) -> None:
        c = self.count[j] + 1
        self.count[j] = c
        delta = x - self.mean[j]
        mean = self.mean[j] + delta / c
        self.mean[j] = mean
        self.m2[j] += delta * (x - mean)
        if x < self.low[j]:
            self.low[j] = x
        if x > self.high[j]:
            self.high[j] = x

    def flags(self, j: int) -> dict:
        n = self.count[j]
        sd = math.sqrt(self.m2[j] / (n - 1)) if n > 1 else 0.0
        return feature_flags(n, self.low[j], self.high[j], self.mean[j], sd)


def profile_matrix(records, orientation: str = "auto", max_values: int = DEFAULT_MAX_VALUES,
                   n_rows_hint=None, examples: int = 5) -> dict:
    """Profile a matrix given as an iterable of field lists (header record first).

    `n_rows_hint` (the number of data rows, if known) sets the row stride of the quantile
    sample. Raises MatrixError when the records do not form a labelled numeric matrix.
    """
    it = iter(records)
    try:
        header = [h.strip() for h in next(it)]
    except StopIteration:
        raise MatrixError("the file is empty") from None
    try:
        first = next(it)
    except StopIteration:
        raise MatrixError("the file has a header but no data rows") from None
    if len(header) < 2 and len(first) < 2:
        raise MatrixError("only one field per line; check the delimiter (--sep)")
    if len(first) == len(header) + 1:
        layout, col_labels = "header_without_corner", header
    elif len(first) == len(header):
        layout, col_labels = "header_with_corner", header[1:]
    else:
        raise MatrixError(
            f"data row 1 has {len(first)} fields but the header has {len(header)}; check the "
            "delimiter (--sep) and that the first row holds the column identifiers"
        )
    n_cols = len(col_labels)
    width = n_cols + 1
    stride = 1
    if n_rows_hint and n_rows_hint * n_cols > max_values:
        stride = int(math.ceil(n_rows_hint * n_cols / max_values))

    cols = _Columns(n_cols)
    row_labels, row_sums, row_missing = [], array("d"), array("l")
    row_tally = dict.fromkeys(TALLY_KEYS, 0)
    constant_rows, text_examples = [], {}
    sample = array("d")
    totals = {"missing": 0, "infinite": 0, "negative": 0, "zero": 0, "finite": 0, "non_integer": 0}
    state = {"min": math.inf, "max": -math.inf, "rows_with_values": 0}

    def consume(index: int, rec) -> None:
        if len(rec) != width:
            raise MatrixError(
                f"data row {index + 1} ('{str(rec[0])[:40] if rec else ''}') has {len(rec)} "
                f"fields; expected {width}"
            )
        label = rec[0].strip()
        values = []
        missing = 0
        for j in range(n_cols):
            kind, x = parse_value(rec[j + 1])
            if kind == "number":
                values.append(x)
                cols.add(j, x)
                if x < 0:
                    totals["negative"] += 1
                elif x == 0:
                    totals["zero"] += 1
                if x != math.floor(x):
                    totals["non_integer"] += 1
            elif kind == "missing":
                missing += 1
                cols.missing[j] += 1
            elif kind == "infinite":
                totals["infinite"] += 1
            else:
                cols.text[j] += 1
                if j in text_examples or len(text_examples) < 50:
                    seen = text_examples.setdefault(j, [])
                    if len(seen) < 2:
                        seen.append({"row": label, "column": col_labels[j], "value": rec[j + 1].strip()[:30]})
        n = len(values)
        totals["missing"] += missing
        totals["finite"] += n
        if n:
            low, high = min(values), max(values)
            mean = math.fsum(values) / n
            sd = math.sqrt(math.fsum((v - mean) ** 2 for v in values) / (n - 1)) if n > 1 else 0.0
            state["min"], state["max"] = min(state["min"], low), max(state["max"], high)
            if state["rows_with_values"] % stride == 0:
                sample.extend(values)
            state["rows_with_values"] += 1
        else:
            low = high = mean = sd = 0.0
        flags = feature_flags(n, low, high, mean, sd)
        for key in TALLY_KEYS:
            if flags[key]:
                row_tally[key] += 1
        if flags["constant"] and len(constant_rows) < examples:
            constant_rows.append(label)
        row_labels.append(label)
        row_sums.append(mean * n)
        row_missing.append(missing)

    consume(0, first)
    for i, rec in enumerate(it, start=1):
        consume(i, rec)

    n_rows = len(row_labels)
    annotation = [j for j in range(n_cols) if cols.text[j] == n_rows]
    keep = [j for j in range(n_cols) if cols.text[j] != n_rows]
    if not keep or totals["finite"] == 0:
        raise MatrixError("no numeric values found; check the delimiter (--sep) and the file layout")
    kept_labels = [col_labels[j] for j in keep]
    n_text = sum(cols.text[j] for j in keep)
    text_list = [e for j in sorted(text_examples) if j in set(keep) for e in text_examples[j]][:examples]

    col_flags = [cols.flags(j) for j in keep]
    col_tally = {k: sum(1 for f in col_flags if f[k]) for k in TALLY_KEYS}

    row_summary = summarise_labels(row_labels, examples)
    col_summary = summarise_labels(kept_labels, examples)
    orient = call_orientation(row_summary, col_summary, orientation)
    genes_in_rows = orient["call"] != "samples_x_genes"
    if genes_in_rows:
        gene_summary, gene_tally, sample_tally = row_summary, row_tally, col_tally
        gene_labels, sample_labels = row_labels, kept_labels
        gene_missing, sample_missing = list(row_missing), [cols.missing[j] for j in keep]
        constant_examples = constant_rows
        sums = [cols.mean[j] * cols.count[j] for j in keep if cols.count[j]]
    else:
        gene_summary, gene_tally, sample_tally = col_summary, col_tally, row_tally
        gene_labels, sample_labels = kept_labels, row_labels
        gene_missing, sample_missing = [cols.missing[j] for j in keep], list(row_missing)
        constant_examples = [kept_labels[i] for i, f in enumerate(col_flags) if f["constant"]][:examples]
        sums = [s for s, m in zip(row_sums, row_missing) if m < len(keep)]
    n_genes_axis, n_samples_axis = len(gene_labels), len(sample_labels)

    ordered = sorted(sample)
    q = {
        "min": state["min"],
        "q25": quantile(ordered, 0.25),
        "median": quantile(ordered, 0.5),
        "q75": quantile(ordered, 0.75),
        "q99": quantile(ordered, 0.99),
        "max": state["max"],
    }
    centered_axis = "gene" if _is_centered(gene_tally) else ("sample" if _is_centered(sample_tally) else None)
    scale = infer_scale(q, totals["non_integer"] == 0, centered_axis)
    if scale["call"] != "centered":
        centered_axis = None
    scale["centered_axis"] = centered_axis
    scale["standardized"] = bool(
        centered_axis and _is_standardized(gene_tally if centered_axis == "gene" else sample_tally)
    )
    scale["quantiles"] = {k: round(v, 6) for k, v in q.items()}
    scale["quantiles_from"] = "all values" if stride == 1 else f"every {stride}th row"

    sample_totals = None
    if orient["call"] != "undetermined" and sums:
        sums.sort()
        sample_totals = {
            "min": round(sums[0], 4),
            "median": round(quantile(sums, 0.5), 4),
            "max": round(sums[-1], 4),
            "max_over_min": round(sums[-1] / sums[0], 4) if sums[0] > 0 else None,
            "sum_to_one_million": all(abs(s - PER_MILLION) <= PER_MILLION_TOLERANCE * PER_MILLION for s in sums),
        }

    n_values = n_rows * len(keep)
    report = {
        "tool": "check_expression_matrix",
        "layout": layout,
        "shape": {"n_rows": n_rows, "n_columns": len(keep)},
        "annotation_columns": [col_labels[j] for j in annotation],
        "orientation": orient,
        "gene_axis": "rows" if genes_in_rows else "columns",
        "gene_axis_assumed": orient["call"] == "undetermined",
        "identifiers": {
            "type": identifier_type(gene_summary),
            "class_counts": gene_summary["class_counts"],
            "anchor_genes_found": gene_summary["anchor_hits"],
            "symbol_case": symbol_case(gene_summary),
            "ensembl_prefix": gene_summary["ensembl_prefix"],
            "ensembl_versioned": gene_summary["ensembl_versioned"],
            "excel_date_like": gene_summary["excel_date"],
            "examples": gene_summary["examples"],
        },
        "scale": scale,
        "sample_totals": sample_totals,
        "values": {
            "n_values": n_values,
            "n_missing": totals["missing"],
            "fraction_missing": round(totals["missing"] / n_values, 6),
            "n_non_numeric": n_text,
            "non_numeric_examples": text_list,
            "n_infinite": totals["infinite"],
            "n_negative": totals["negative"],
            "fraction_negative": round(totals["negative"] / totals["finite"], 6),
            "fraction_zero": round(totals["zero"] / totals["finite"], 6),
            "genes_with_missing": sum(1 for m in gene_missing if m),
            "genes_all_missing": gene_tally["all_missing"],
            "max_fraction_missing_per_gene": round(max(gene_missing) / n_samples_axis, 6),
            "samples_all_missing": sample_tally["all_missing"],
            "max_fraction_missing_per_sample": round(max(sample_missing) / n_genes_axis, 6),
            "genes_with_negative": gene_tally["negative"],
        },
        "constant": {
            "n_constant_genes": gene_tally["constant"],
            "n_all_zero_genes": gene_tally["all_zero"],
            "n_constant_samples": sample_tally["constant"],
            "examples": constant_examples,
        },
        "duplicates": {
            "gene_ids": find_duplicates(gene_labels, examples),
            "sample_ids": find_duplicates(sample_labels, examples),
        },
    }
    report["issues"] = build_issues(report)
    return report


def build_issues(report: dict) -> list:
    """Turn a profile into coded findings.

    Severity: 'error' (resolve before scoring), 'warning' (decide, and state the decision in
    the write-up), 'info'.
    """
    issues = []

    def add(code: str, severity: str, message: str) -> None:
        issues.append({"code": code, "severity": severity, "message": message})

    orient, ident, scale = report["orientation"], report["identifiers"], report["scale"]
    values, dups = report["values"], report["duplicates"]

    if orient["call"] == "samples_x_genes":
        add("ORIENTATION_SAMPLES_BY_GENES", "error",
            "Rows are samples and columns are genes. Transpose to genes x samples before scoring: "
            "tmelens expects genes in rows and would find no signature gene.")
    elif orient["call"] == "undetermined":
        add("ORIENTATION_UNDETERMINED", "error",
            f"Could not tell which axis holds genes ({orient['basis']}). The per-gene checks below "
            "assume genes in rows. Confirm the layout with the data source and rerun with --orientation.")
    elif orient["confidence"] == "low":
        add("ORIENTATION_BY_SHAPE_ONLY", "warning",
            "No gene identifiers were recognised, so the orientation was inferred from the matrix "
            "shape alone. Confirm it and rerun with --orientation.")

    if report["annotation_columns"]:
        add("ANNOTATION_COLUMNS", "warning",
            "Columns without any numeric value were left out of the checks and must be removed "
            f"before scoring: {', '.join(report['annotation_columns'][:5])}.")

    kind = ident["type"]
    if kind == "ensembl_gene_id":
        extra = " Strip the version suffix first." if ident["ensembl_versioned"] else ""
        species = ("" if ident["ensembl_prefix"] in (None, "ENSG")
                   else f" The prefix {ident['ensembl_prefix']} is not the human one (ENSG).")
        add("IDENTIFIERS_ENSEMBL", "error",
            "Gene identifiers are Ensembl gene IDs. Map them to HGNC symbols before scoring "
            f"symbol-based signatures; tmelens does not map identifiers.{extra}{species}")
    elif kind == "ensembl_transcript_id":
        add("IDENTIFIERS_TRANSCRIPTS", "error",
            "Identifiers are Ensembl transcript IDs. Aggregate to gene level and map to HGNC symbols first.")
    elif kind in ("refseq_id", "probe_id", "integer_id", "composite_id"):
        hint = {
            "refseq_id": "RefSeq accessions",
            "probe_id": "array probe IDs (map them with the platform annotation and choose how to collapse probes)",
            "integer_id": "integers (Entrez gene IDs or row numbers)",
            "composite_id": "composite labels such as SYMBOL|ID (split them and keep the symbol)",
        }[kind]
        add("IDENTIFIERS_NOT_SYMBOLS", "error",
            f"Gene identifiers are {hint}. Convert to HGNC symbols before scoring symbol-based signatures.")
    elif kind == "unrecognised":
        add("IDENTIFIERS_UNRECOGNISED", "warning",
            "Labels look like symbols but fewer than 3 anchor genes (ACTB, GAPDH, PTPRC, CD8A, ...) "
            "were found: a small panel, custom identifiers or another naming scheme. Check "
            "signature coverage before trusting any score.")
    elif kind == "mixed":
        add("IDENTIFIERS_MIXED", "warning",
            f"Identifier types are mixed ({ident['class_counts']}). Harmonise them to HGNC symbols.")
    if ident["symbol_case"] == "sentence":
        add("SYMBOL_CASE_NOT_HUMAN", "warning",
            "Recognised symbols are in sentence case (mouse or rat style). Upper-casing them is not "
            "an orthologue mapping: map to human orthologues before scoring human signatures.")
    if ident["excel_date_like"]["n"]:
        add("EXCEL_DATE_SYMBOLS", "error",
            f"{ident['excel_date_like']['n']} identifiers look like dates (for example "
            f"{ident['excel_date_like']['examples'][0]}): a spreadsheet converted gene symbols. "
            "Recover the identifiers from the original file; the conversion cannot be undone "
            "reliably (1-Mar may have been MARCH1 or MARC1, now MARCHF1 and MTARC1 at HGNC).")

    gene_dups, sample_dups = dups["gene_ids"], dups["sample_ids"]
    if gene_dups["n_empty"]:
        add("EMPTY_IDENTIFIERS", "error", f"{gene_dups['n_empty']} gene identifiers are empty.")
    if gene_dups["n_duplicated_ids"]:
        add("DUPLICATE_GENE_IDS", "warning",
            f"{gene_dups['n_duplicated_ids']} gene identifiers occur more than once "
            f"({gene_dups['n_entries_involved']} entries). Choose a collapsing rule and report it; "
            "tmelens averages duplicated symbols by default.")
    if gene_dups["n_more_after_upper_casing"]:
        add("DUPLICATE_GENE_IDS_AFTER_UPPER_CASING", "warning",
            f"{gene_dups['n_more_after_upper_casing']} more identifiers collide once upper-cased; "
            "tmelens upper-cases symbols and would then average them.")
    if sample_dups["n_duplicated_ids"]:
        add("DUPLICATE_SAMPLE_IDS", "error",
            f"{sample_dups['n_duplicated_ids']} sample identifiers are repeated. Resolve them "
            "before any join with clinical data.")

    if values["n_non_numeric"]:
        add("NON_NUMERIC_VALUES", "error",
            f"{values['n_non_numeric']} cells are neither numbers nor missing-value tokens (see "
            "non_numeric_examples). Fix the export: decimal commas, text or spreadsheet errors.")
    if values["n_infinite"]:
        add("INFINITE_VALUES", "error",
            f"{values['n_infinite']} infinite values, usually the log of zero upstream. Go back to "
            "the untransformed data or decide how to replace them.")
    if values["samples_all_missing"]:
        add("SAMPLES_ALL_MISSING", "error",
            f"{values['samples_all_missing']} samples have no numeric value at all.")
    if values["n_missing"]:
        add("MISSING_VALUES", "warning",
            f"{values['n_missing']} missing values in {values['genes_with_missing']} genes "
            f"({values['genes_all_missing']} entirely missing; at most "
            f"{values['max_fraction_missing_per_gene']:.0%} of a gene). tmelens drops all-missing genes "
            "and fills the rest with the gene median: decide whether genes with many gaps should be "
            "dropped instead, and report it.")

    if scale["call"] == "centered":
        how = "standardised to unit variance" if scale["standardized"] else "centred"
        if scale["centered_axis"] == "gene":
            add("SCALE_CENTERED_PER_GENE", "error",
                f"Each gene was {how} across samples. Within-sample ranks no longer reflect expression, "
                "so rank-based single-sample scores and mean scores are not interpretable, and the "
                "data cannot be log-transformed. Obtain the matrix as it was before centring; only a "
                "cohort-relative z-score reading remains valid, relative to whichever samples were "
                "used for the centring.")
        else:
            add("SCALE_CENTERED_PER_SAMPLE", "error",
                f"Each sample was {how} across genes. Values are no longer expression units and "
                "cannot be log-transformed; within-sample ranks are preserved. Obtain the matrix as "
                "it was before centring, or state that scores are computed on per-sample standardised values.")
    elif scale["call"] == "linear":
        if values["n_negative"]:
            add("NEGATIVE_ON_LINEAR_SCALE", "error",
                f"{values['n_negative']} negative values on a linear scale (minimum "
                f"{scale['quantiles']['min']}). log2(x + 1) is undefined below -1 and tmelens raises on "
                "negative input: find out how they arose before flooring or dropping anything.")
        if scale["tmelens_auto_log_would_transform"]:
            add("SCALE_LINEAR", "info",
                "Values are on a linear scale. Apply log2(x + 1) before z-score or mean scoring "
                "(tmelens log='auto' does so).")
        else:
            add("SCALE_AUTO_LOG_WOULD_NOT_TRIGGER", "warning",
                f"Values are on a linear scale ({'; '.join(scale['reasons'])}), but the GEO2R rule is "
                "not met, so tmelens log='auto' would leave them untransformed. Request the transform "
                "explicitly (--log yes, or log=True).")
        totals = report["sample_totals"] or {}
        if totals.get("sum_to_one_million"):
            rounded = " The values are integers: rounded TPM or CPM." if scale["integer_valued"] else ""
            add("SUM_TO_ONE_MILLION", "info",
                "Every sample sums to within 1% of one million: TPM or CPM, already normalised for "
                f"depth.{rounded}")
        elif scale["integer_valued"]:
            ratio = totals.get("max_over_min")
            spread = f" Sample totals differ by a factor of {ratio}." if ratio else ""
            add("COUNT_LIKE", "warning",
                "All values are non-negative integers: read counts, not normalised for sequencing "
                f"depth.{spread} Convert to CPM or TPM before z-score or mean scoring; tmelens "
                "log-transforms but does not normalise library size.")
    elif scale["call"] == "log":
        if values["n_negative"]:
            add("NEGATIVE_VALUES", "warning",
                f"{values['n_negative']} negative values (minimum {scale['quantiles']['min']}) without "
                "centring: expected for logs of values below 1 (log-CPM, log2(x + 0.001)) or for "
                "log-ratios. A minimum near -10 means log2(x + 0.001): zeros then sit far below every "
                "expressed value and dominate means and z-scores of lowly expressed genes, so "
                "re-express as log2(x + 1). Confirm the transform from the data source; do not log again.")
        add("SCALE_LOG", "info",
            "Values look log-transformed already. Do not log again: pass --log no to tmelens score "
            "or report (tmelens 0.1.0 maps it to log=False), or log=False to "
            "tmelens.prepare_expression, once the data source confirms the unit.")
    else:
        add("SCALE_UNDETERMINED", "warning", "All values are identical; the scale cannot be assessed.")

    const = report["constant"]
    if const["n_constant_genes"]:
        add("CONSTANT_GENES", "info",
            f"{const['n_constant_genes']} genes are constant across samples ({const['n_all_zero_genes']} "
            "all zero). A constant signature gene counts as covered but adds nothing to a z-score: "
            "check with signature_coverage.py.")
    if orient["n_samples"] is not None and orient["n_samples"] < FEW_SAMPLES:
        add("FEW_SAMPLES", "warning",
            f"Only {orient['n_samples']} samples. Gene-wise z-scores and cohort quantiles (tertile, "
            "median) are unstable at this size: prefer per-sample scores and report continuous values.")
    if orient["n_genes"] is not None and orient["n_genes"] < FEW_GENES:
        add("FEW_GENES", "warning",
            f"Only {orient['n_genes']} genes: a targeted panel or a filtered matrix. Rank-based "
            "single-sample scores rank a set against the other genes in the matrix, so they are not "
            "comparable with genome-wide scores.")
    return issues


# --------------------------------------------------------------------------- file reading

def open_text(path: Path):
    """Open a plain or gzip-compressed text file as UTF-8 (a byte-order mark is ignored)."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, "r", encoding="utf-8-sig", newline="")


def choose_delimiter(path: Path, first_line: str, sep: str = "auto") -> str:
    """Return the requested delimiter, or the most frequent of tab, comma and semicolon in
    the header line (file extension as the tie-breaker when none occurs)."""
    if sep != "auto":
        return DELIMITERS[sep]
    counts = {name: first_line.count(char) for name, char in DELIMITERS.items()}
    best = max(("tab", "comma", "semicolon"), key=lambda k: counts[k])
    if counts[best] == 0:
        return "," if ".csv" in [s.lower() for s in path.suffixes] else "\t"
    return DELIMITERS[best]


def iter_records(path: Path, sep: str = "auto"):
    """Yield field lists, skipping leading '#' lines, a GCT size line and blank lines."""
    with open_text(path) as fh:
        first = fh.readline()
        gct = first.startswith("#1.")
        while first.startswith("#"):
            first = fh.readline()
        if gct:
            first = fh.readline()
        if not first.strip():
            return
        delimiter = choose_delimiter(path, first, sep)
        for rec in csv.reader([first], delimiter=delimiter):
            yield rec
        for rec in csv.reader(fh, delimiter=delimiter):
            if rec and any(field.strip() for field in rec):
                yield rec


def count_lines(path: Path) -> int:
    """Number of lines in the file (used only to size the quantile sample)."""
    with open_text(path) as fh:
        return sum(1 for _ in fh)


def check_file(path, sep: str = "auto", orientation: str = "auto",
               max_values: int = DEFAULT_MAX_VALUES, examples: int = 5) -> dict:
    """Profile the matrix stored at `path` and return the report dictionary."""
    path = Path(path)
    if not path.is_file():
        raise MatrixError(f"file not found: {path}")
    try:
        hint = max(count_lines(path) - 1, 1)
        report = profile_matrix(iter_records(path, sep), orientation, max_values, hint, examples)
    except (UnicodeDecodeError, csv.Error, OSError, EOFError) as exc:
        raise MatrixError(f"could not read {path.name} as delimited UTF-8 text: {exc}") from None
    return {"file": path.name, **report}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check a bulk expression matrix (orientation, identifier type, linear or log "
                    "scale, duplicated identifiers, missing and negative values, constant genes) "
                    "and print a JSON report.",
        epilog="Exit status 0: report produced (read its 'issues' list). Exit status 2: unreadable input.",
    )
    parser.add_argument("matrix",
                        help="TSV or CSV matrix, optionally .gz; the first row and the first column hold identifiers")
    parser.add_argument("--sep", choices=["auto", "tab", "comma", "semicolon"], default="auto",
                        help="field delimiter (default: the most frequent of tab, comma and semicolon in the header)")
    parser.add_argument("--orientation", choices=["auto", "genes-by-samples", "samples-by-genes"], default="auto",
                        help="state the layout instead of inferring it")
    parser.add_argument("--max-values", type=int, default=DEFAULT_MAX_VALUES,
                        help="largest number of values used for the quantiles; larger matrices are "
                             "sampled by taking every k-th row (default: %(default)s)")
    parser.add_argument("--examples", type=int, default=5,
                        help="how many example identifiers to list per finding (default: %(default)s)")
    parser.add_argument("--out", help="write the JSON report to this file instead of standard output")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.max_values < 1 or args.examples < 0:
        print("error: --max-values must be at least 1 and --examples at least 0", file=sys.stderr)
        return 2
    try:
        report = check_file(args.matrix, args.sep, args.orientation, args.max_values, args.examples)
    except MatrixError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    text = json.dumps(report, indent=2, ensure_ascii=False)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
