#!/usr/bin/env python3
"""Per-signature gene coverage of a bulk expression matrix. Standard library only.

Reads gene sets from a GMT file (one set per line, tab-separated: name, description, genes) and
an expression matrix (TSV or CSV, optionally gzip-compressed; genes in the first column by
default), and reports for every signature how many of its genes are measured.

Definitions, chosen to reproduce tmelens 0.1.0 [Wang2026]:

  measured gene   the identifier occurs in the matrix with at least one numeric value
                  (tmelens drops all-missing genes before scoring). Matching ignores case and
                  surrounding spaces, because tmelens upper-cases matrix symbols.
  renamed symbols the three retired symbols that tmelens 0.1.0 renames in the matrix (CTGF to
                  CCN2, FAM101B to RFLNB, KIAA0125 to FAM30A) are renamed here too, under the
                  same condition: the old symbol is measured and the current one is not.
  n_genes        distinct genes of the set after upper-casing (a repeated gene counts once,
                  as in tmelens.signature_coverage)
  coverage        n_found / n_genes (0 for a set without genes)
  flagged         n_found < --min-genes or coverage < --min-coverage. With the defaults
                  (1 and 0.5) this is exactly the rule under which tmelens.score_signatures
                  returns NaN for the set. The 0.5 is a software default (convention), not a
                  validated threshold.
  n_constant      measured genes whose numeric values are identical in every sample. They count
                  as covered but add no between-sample information to z-score or mean scores
                  (tmelens sets the z-score of a zero-variance gene to 0).
  informative_coverage  (n_found - n_constant) / n_genes

Warnings cover identifier mismatches (Ensembl, Entrez or probe IDs on one side and symbols on the
other), GMT genes containing lower-case letters (tmelens upper-cases the matrix but not GMT
genes, so it would never match them), duplicated set names (tmelens keeps the last one), GMT
lines without genes, and sets whose coverage passes only because constant genes are counted.

Usage examples:

  python3 signature_coverage.py signatures.gmt expression.tsv
  python3 signature_coverage.py sets.gmt cohort.csv.gz --min-coverage 0.8 --format tsv --out coverage.tsv

Exit status: 0 when a report was produced (read "flagged" and "warnings"), 2 when an input
could not be read or an argument is invalid.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from pathlib import Path

DEFAULT_MIN_COVERAGE = 0.5   # tmelens score_signatures(min_coverage=0.5) [Wang2026]; convention
DEFAULT_MIN_GENES = 1        # tmelens score_signatures(min_genes=1) [Wang2026]
MISSING_TOKENS = frozenset({"", "na", "n/a", "#n/a", "nan", "null", "none", "."})
DELIMITERS = {"tab": "\t", "comma": ",", "semicolon": ";"}
ID_PATTERNS = (
    ("ensembl", re.compile(r"^ENS[A-Z]{0,6}[GT]\d{11}(\.\d+)?$", re.IGNORECASE)),
    ("entrez_or_integer", re.compile(r"^\d+$")),
    ("probe", re.compile(r"^(AFFX-.+|.+_(at|st)|ILMN_\d+|A_\d{2}_P\d+)$", re.IGNORECASE)),
)
MAJORITY = 0.8               # heuristic: share of identifiers that names an identifier type
# tmelens 0.1.0 prepare_expression(map_aliases=True) renames these retired symbols [Wang2026];
# the current symbols were checked at HGNC (November 2026) [HGNC2026].
TMELENS_ALIASES = {"CTGF": "CCN2", "FAM101B": "RFLNB", "KIAA0125": "FAM30A"}


class InputError(ValueError):
    """An input file or argument cannot be used."""


# --------------------------------------------------------------------------- pure functions

def normalise(gene: str) -> str:
    return gene.strip().upper()


def parse_gmt(lines) -> tuple:
    """Parse GMT lines into ({name: [genes as written]}, warnings).

    Lines with fewer than three tab-separated fields are skipped, as tmelens.read_gmt does.
    A repeated set name keeps the last definition, as a Python dict (and tmelens) does.
    """
    sets: dict = {}
    warnings = []
    skipped, duplicated = [], []
    for number, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n\r")
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            skipped.append(number)
            continue
        name = parts[0].strip()
        genes = [g.strip() for g in parts[2:] if g.strip()]
        if name in sets:
            duplicated.append(name)
        sets[name] = genes
    if skipped:
        warnings.append(_warning(
            "GMT_LINES_SKIPPED",
            f"{len(skipped)} GMT line(s) have fewer than three tab-separated fields and were skipped "
            f"(first: line {skipped[0]}). A space-separated file is the usual cause."))
    if duplicated:
        warnings.append(_warning(
            "GMT_DUPLICATE_NAMES",
            f"Set names defined more than once (the last definition is used, as in tmelens): "
            f"{', '.join(sorted(set(duplicated))[:5])}."))
    return sets, warnings


def identifier_kind(ids) -> str:
    """'ensembl', 'entrez_or_integer', 'probe' or 'symbol_or_other' by a MAJORITY share."""
    ids = [i for i in ids if i]
    if not ids:
        return "symbol_or_other"
    for kind, pattern in ID_PATTERNS:
        if sum(1 for i in ids if pattern.match(i)) / len(ids) >= MAJORITY:
            return kind
    return "symbol_or_other"


def coverage_table(sets: dict, measured: set, constant: set, min_coverage: float = DEFAULT_MIN_COVERAGE,
                   min_genes: int = DEFAULT_MIN_GENES) -> list:
    """One record per set. `measured` and `constant` hold upper-cased matrix identifiers."""
    rows = []
    for name, genes in sets.items():
        distinct = list(dict.fromkeys(normalise(g) for g in genes))
        found = [g for g in distinct if g in measured]
        flat = [g for g in found if g in constant]
        n = len(distinct)
        cov = len(found) / n if n else 0.0
        informative = (len(found) - len(flat)) / n if n else 0.0
        rows.append({
            "signature": name,
            "n_genes": n,
            "n_found": len(found),
            "coverage": round(cov, 4),
            "n_constant": len(flat),
            "informative_coverage": round(informative, 4),
            "flagged": len(found) < min_genes or cov < min_coverage or n == 0,
            "missing": [g for g in distinct if g not in measured],
            "constant": flat,
        })
    return rows


def coverage_warnings(sets: dict, rows: list, matrix_ids, min_coverage: float) -> list:
    """Warnings about identifiers, case and constant genes."""
    warnings = []
    gmt_genes = [g for genes in sets.values() for g in genes]
    matrix_kind = identifier_kind(list(matrix_ids))
    gmt_kind = identifier_kind([normalise(g) for g in gmt_genes])
    if matrix_kind != gmt_kind:
        warnings.append(_warning(
            "IDENTIFIER_MISMATCH",
            f"Matrix identifiers look like '{matrix_kind}' but GMT genes look like '{gmt_kind}'. "
            "Map both to the same identifier type (HGNC symbols for published signatures) before "
            "reading coverage: a mismatch makes every signature look absent."))
    lower = sorted({g for g in gmt_genes if g != g.upper()})
    if lower:
        warnings.append(_warning(
            "GMT_LOWER_CASE_GENES",
            f"{len(lower)} GMT gene(s) contain lower-case letters (for example {', '.join(lower[:3])}). "
            "They match here because matching ignores case, but tmelens upper-cases only the matrix, "
            "so it would miss them: upper-case the GMT first. Sentence case usually means mouse "
            "symbols, and upper-casing is not an orthologue mapping."))
    empty = [r["signature"] for r in rows if r["n_genes"] == 0]
    if empty:
        warnings.append(_warning("GMT_SET_WITHOUT_GENES", f"Sets without genes: {', '.join(empty[:5])}."))
    rescued = [r["signature"] for r in rows if not r["flagged"] and r["informative_coverage"] < min_coverage]
    if rescued:
        warnings.append(_warning(
            "COVERAGE_PASSES_ONLY_WITH_CONSTANT_GENES",
            f"These sets pass the coverage threshold only because constant genes are counted: "
            f"{', '.join(rescued[:5])}."))
    return warnings


def _warning(code: str, message: str) -> dict:
    return {"code": code, "message": message}


# --------------------------------------------------------------------------- matrix reading

def _is_number(token: str) -> bool:
    try:
        x = float(token)
    except ValueError:
        return False
    return x == x and not math.isinf(x) and "_" not in token


def scan_matrix(records, orientation: str = "genes-by-samples") -> tuple:
    """Return (measured ids, constant ids, all ids) from field lists (header record first).

    Identifiers are upper-cased. An identifier is measured when any of its values is numeric and
    constant when every row (or column) that carries it has two or more numeric values, all equal.
    Records are consumed one at a time, so a genome-wide matrix is never held in memory.
    """
    if orientation not in ("genes-by-samples", "samples-by-genes"):
        raise InputError(f"unknown orientation {orientation!r}")
    it = iter(records)
    try:
        header = [h.strip() for h in next(it)]
    except StopIteration:
        raise InputError("the matrix file is empty") from None

    stats: dict = {}   # id -> [n_numeric_any, all_vectors_constant]

    def merge(gene_id: str, n_numeric: int, flat: bool) -> None:
        entry = stats.setdefault(normalise(gene_id), [0, True])
        entry[0] += n_numeric
        entry[1] = entry[1] and flat

    width, col_labels, columns = None, [], []   # columns: [n, low, high] per column
    for i, rec in enumerate(it, start=1):
        if width is None:
            width = len(rec)
            if width == len(header) + 1:
                col_labels = header
            elif width == len(header):
                col_labels = header[1:]
            else:
                raise InputError(
                    f"data row 1 has {width} fields but the header has {len(header)}; check the delimiter")
            if not col_labels:
                raise InputError("the matrix has no value columns; check the delimiter (--sep)")
            columns = [[0, math.inf, -math.inf] for _ in col_labels]
        if len(rec) != width:
            raise InputError(f"data row {i} has {len(rec)} fields; expected {width}")
        if orientation == "genes-by-samples":
            nums = [float(v) for v in rec[1:] if _is_number(v.strip())]
            merge(rec[0], len(nums), len(nums) >= 2 and min(nums) == max(nums))
        else:
            for col, v in zip(columns, rec[1:]):
                if _is_number(v.strip()):
                    x = float(v)
                    col[0] += 1
                    col[1] = min(col[1], x)
                    col[2] = max(col[2], x)
    if width is None:
        raise InputError("the matrix has a header but no data rows")
    if orientation == "samples-by-genes":
        for label, (n, low, high) in zip(col_labels, columns):
            merge(label, n, n >= 2 and low == high)
    stats.pop("", None)
    measured = {g for g, (n, _) in stats.items() if n > 0}
    constant = {g for g, (n, flat) in stats.items() if n > 0 and flat}
    return measured, constant, sorted(stats)


def apply_aliases(measured: set, constant: set, aliases=None) -> tuple:
    """Rename retired symbols as tmelens.prepare_expression does: old -> new when the old symbol
    is measured and the new one is not. Returns (measured, constant, {old: new} applied)."""
    aliases = TMELENS_ALIASES if aliases is None else aliases
    measured, constant = set(measured), set(constant)
    applied = {old: new for old, new in aliases.items() if old in measured and new not in measured}
    for old, new in applied.items():
        measured.discard(old)
        measured.add(new)
        if old in constant:
            constant.discard(old)
            constant.add(new)
    return measured, constant, applied


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, "r", encoding="utf-8-sig", newline="")


def iter_matrix_records(path: Path, sep: str = "auto"):
    """Yield field lists from a delimited matrix file, one line at a time. '#' lines at the top
    (and a GCT size line) and blank lines are skipped."""
    with _open_text(path) as fh:
        first = fh.readline()
        gct = first.startswith("#1.")
        while first.startswith("#"):
            first = fh.readline()
        if gct:
            first = fh.readline()
        while first and not first.strip():
            first = fh.readline()
        if not first:
            raise InputError("the matrix file is empty")
        if sep == "auto":
            counts = {name: first.count(ch) for name, ch in DELIMITERS.items()}
            best = max(("tab", "comma", "semicolon"), key=lambda k: counts[k])
            delimiter = DELIMITERS[best] if counts[best] else ("," if ".csv" in path.name.lower() else "\t")
        else:
            delimiter = DELIMITERS[sep]
        yield from csv.reader([first], delimiter=delimiter)
        for rec in csv.reader(fh, delimiter=delimiter):
            if rec and any(field.strip() for field in rec):
                yield rec


def read_matrix(path: Path, sep: str = "auto", orientation: str = "genes-by-samples") -> tuple:
    """Scan a delimited matrix file (see `iter_matrix_records` and `scan_matrix`)."""
    if not path.is_file():
        raise InputError(f"file not found: {path}")
    try:
        return scan_matrix(iter_matrix_records(path, sep), orientation)
    except (UnicodeDecodeError, OSError, EOFError, csv.Error) as exc:
        raise InputError(f"could not read {path.name} as delimited UTF-8 text: {exc}") from None


def read_gmt(path: Path) -> tuple:
    if not path.is_file():
        raise InputError(f"file not found: {path}")
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            sets, warnings = parse_gmt(fh)
    except (UnicodeDecodeError, OSError) as exc:
        raise InputError(f"could not read {path.name} as UTF-8 text: {exc}") from None
    if not sets:
        raise InputError(f"no gene sets found in {path.name}; GMT lines are name<TAB>description<TAB>gene...")
    return sets, warnings


def build_report(gmt_path, matrix_path, min_coverage=DEFAULT_MIN_COVERAGE, min_genes=DEFAULT_MIN_GENES,
                 sep="auto", orientation="genes-by-samples") -> dict:
    if not 0.0 <= min_coverage <= 1.0:
        raise InputError("--min-coverage must lie between 0 and 1")
    if min_genes < 0:
        raise InputError("--min-genes must be at least 0")
    sets, warnings = read_gmt(Path(gmt_path))
    measured, constant, ids = read_matrix(Path(matrix_path), sep, orientation)
    measured, constant, renamed = apply_aliases(measured, constant)
    rows = coverage_table(sets, measured, constant, min_coverage, min_genes)
    warnings += coverage_warnings(sets, rows, ids, min_coverage)
    return {
        "tool": "signature_coverage",
        "gmt": Path(gmt_path).name,
        "matrix": Path(matrix_path).name,
        "orientation": orientation,
        "min_coverage": min_coverage,
        "min_genes": min_genes,
        "renamed_symbols": renamed,
        "n_matrix_genes": len(measured),
        "n_signatures": len(rows),
        "n_flagged": sum(1 for r in rows if r["flagged"]),
        "signatures": rows,
        "warnings": warnings,
    }


def to_tsv(report: dict) -> str:
    cols = ["signature", "n_genes", "n_found", "coverage", "n_constant", "informative_coverage", "flagged", "missing"]
    out = ["\t".join(cols)]
    for r in report["signatures"]:
        out.append("\t".join(
            ",".join(r[c]) if c == "missing" else ("yes" if r[c] is True else "no" if r[c] is False else str(r[c]))
            for c in cols))
    return "\n".join(out) + "\n"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Report per-signature gene coverage of an expression matrix and flag signatures "
                    "below a coverage threshold (defaults reproduce the tmelens NaN rule).",
        epilog="Exit status 0: report produced. Exit status 2: unreadable input or invalid argument.")
    p.add_argument("gmt", help="gene sets in GMT format (name<TAB>description<TAB>gene...)")
    p.add_argument("matrix", help="expression matrix, TSV or CSV, optionally .gz")
    p.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                   help="flag sets whose found fraction is below this (default: %(default)s, the tmelens default)")
    p.add_argument("--min-genes", type=int, default=DEFAULT_MIN_GENES,
                   help="flag sets with fewer found genes than this (default: %(default)s, the tmelens default)")
    p.add_argument("--orientation", choices=["genes-by-samples", "samples-by-genes"], default="genes-by-samples",
                   help="where the gene identifiers are: first column (default) or header row")
    p.add_argument("--sep", choices=["auto", "tab", "comma", "semicolon"], default="auto",
                   help="matrix delimiter (default: the most frequent of tab, comma and semicolon in the header)")
    p.add_argument("--format", choices=["json", "tsv"], default="json", help="output format (default: json)")
    p.add_argument("--out", help="write the report to this file instead of standard output")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = build_report(args.gmt, args.matrix, args.min_coverage, args.min_genes, args.sep, args.orientation)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    text = to_tsv(report) if args.format == "tsv" else json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
