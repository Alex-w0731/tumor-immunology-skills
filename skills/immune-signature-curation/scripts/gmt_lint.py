#!/usr/bin/env python3
"""Lint a GMT gene-set file before it is used for immune scoring or enrichment analysis.

GMT is the gene-set format of GSEA and MSigDB: a text file in which each line defines one gene set
by its name, a second field (a description or identifier, often a URL or DOI) and the genes, all
tab-separated; the gene identifiers must match those of the expression data [Reimand2019].

What is checked (issue codes appear in the JSON report):

  format       byte-order mark, CRLF line endings, blank lines, lines starting with '#', lines
               without a tab, lines with fewer than three fields
  set names    empty, surrounding whitespace, duplicated exactly or only by case (dictionary-based
               readers such as tmelens read_gmt keep the last copy [Wang2026]), identical genes
               under two names
  description  empty, a bare 'published' label, or shaped like a gene symbol (the description
               column is probably missing, so the first gene is read as the description and lost);
               the provenance kind is recorded: doi, pmid, url, published, curated, derived, user,
               empty, gene_like, text
  genes        empty fields, whitespace, quotes, several identifiers in one field, missing-value
               tokens, duplicates within a set (exact, by case, or after renaming), spreadsheet
               damage (dates, floating-point numbers, five-digit date serials) [Ziemann2016]
               [Abeysooriya2021], identifier types mixed in one file, case style (HGNC symbols are
               upper case except the 'orf' of C#orf# symbols [Bruford2020]; mouse symbols start
               with an upper-case letter followed by lower case [MGI2026]), previous symbols and
               common immune aliases from a small built-in table checked at HGNC in November 2026
               [HGNC2026], or every symbol against a full HGNC table given with --hgnc
  size         one-gene sets; sets below --min-size (default 10) or above --max-size (default
               500), the permissive ends of the bounds (fewer than 10-15 or more than 200-500
               genes) that Reimand et al. often recommend for excluding pathway gene sets when
               analysing human expression data [Reimand2019]; small marker panels are normal for
               mean or z-score scoring and are reported as info only

A symbol that HGNC links to more than one approved gene, counting previous symbols and aliases
together, is reported as ambiguous: DEC1 is the previous symbol of DELEC1 and an alias of BHLHE40
[HGNC2026].

Heuristics of this script (not published thresholds): a file is called human-style when at least
90% of its symbols are upper case and mouse-style when at least 90% are sentence case; mouse MHC
symbols (H2-K1, H2-Ab1) do not vote; a five-digit number in a file of symbols is read as a
spreadsheet date serial and decoded in the 1900 date system (day 0 = 1899-12-30), so that 43891
reads as 2020-03-01, the date a spreadsheet stores for '1-Mar' typed in 2020; a file saved in the
1904 date system decodes to wrong dates.

The script never rewrites the file: renaming a symbol or restoring a date-mangled gene needs the
original source.

Usage examples:
  python3 gmt_lint.py signatures.gmt --out lint.json
  python3 gmt_lint.py sets.gmt.gz --hgnc hgnc_complete_set.txt --species human --strict

Exit status: 0 report written; 1 report written, --strict given and the report has errors;
2 unreadable input, unwritable --out or bad arguments.
"""
from __future__ import annotations

import argparse
import datetime
import gzip
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

SEVERITY_RANK = {"error": 0, "warning": 1, "info": 2}
ID_KINDS = ("symbol", "integer", "ensembl_gene", "ensembl_other", "refseq", "probe")
MISSING_TOKENS = {"na", "n/a", "#n/a", "nan", "null", "none", "-", "?"}
NOT_GENE_WORDS = {"NA", "NONE", "NULL", "USER", "CURATED", "DERIVED", "PUBLISHED", "UNKNOWN",
                  "OTHER", "TCGA", "GEO", "GSEA", "MSIGDB", "KEGG", "GO"}

_MON = r"(?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]{0,6}"
RE_DATE_D_MON = re.compile(r"^(\d{1,2})-(" + _MON + r")(?:-\d{2,4})?$", re.IGNORECASE)
RE_DATE_MON_D = re.compile(r"^(" + _MON + r")-(\d{1,4})$", re.IGNORECASE)
RE_DATE_NUMERIC = re.compile(r"^\d{1,4}[/.\-]\d{1,2}[/.\-]\d{1,4}$")
# Spreadsheet floating-point display such as 2.31E+13 [Ziemann2016] or 9.33E+22 [Abeysooriya2021].
# An unconverted RIKEN clone accession such as 2310009E13 has no '.' and no exponent sign.
RE_SCIENTIFIC = re.compile(r"^[+-]?(?:\d+\.\d*|\.\d+)[eE][+-]?\d+$|^[+-]?\d+[eE][+-]\d+$")
RE_DECIMAL = re.compile(r"^[+-]?\d*\.\d+$")
RE_INTEGER = re.compile(r"^\d+$")
RE_ENSEMBL_GENE = re.compile(r"^ENS[A-Z]{0,6}G\d{11}(?:\.\d+)?$")
RE_ENSEMBL_OTHER = re.compile(r"^ENS[A-Z]{0,6}[TPE]\d{11}(?:\.\d+)?$")
RE_REFSEQ = re.compile(r"^[NX][MRP]_\d+(?:\.\d+)?$")
RE_PROBE = re.compile(r"^(?:AFFX-\S+|\d+(?:_[a-z])*_at|ILMN_\d+|A_\d{2}_P\d+)$")
RE_SYMBOL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.\-_@]*$")
RE_ORF = re.compile(r"^C(?:\d{1,2}|X|Y)orf\d+$")
# Mouse MHC symbols keep an upper-case letter after 'H2-' (H2-K1, H2-D1, H2-Aa, H2-Ab1, H2-Eb1),
# so they fit neither the sentence-case nor the upper-case pattern; no human symbol starts 'H2-'.
RE_MOUSE_MHC = re.compile(r"^H2-[A-Z][A-Za-z0-9]*$")
RE_DOI = re.compile(r"\b10\.\d{4,9}/\S+")
RE_PMID = re.compile(r"^PMID:?\s*\d+$", re.IGNORECASE)
RE_GENE_LIKE = re.compile(r"^[A-Z][A-Z0-9]{1,11}(?:-[A-Z0-9]{1,6})?$")

# Previous HGNC symbols -> current symbol(s), each checked at rest.genenames.org in November 2026
# [HGNC2026]. SEPT#, MARCH# and DEC1 were renamed because spreadsheets turn them into dates
# [Bruford2020]; CTGF, FAM101B and KIAA0125 are the three renames tmelens applies [Wang2026].
RETIRED_BUILTIN = {
    "CTGF": ("CCN2",), "FAM101B": ("RFLNB",), "KIAA0125": ("FAM30A",),
    "MARC1": ("MTARC1",), "MARC2": ("MTARC2",), "DEC1": ("DELEC1",), "SEPT13": ("SEPTIN7P2",),
    "CD8": ("CD8A",), "CD45": ("PTPRC",), "CD39": ("ENTPD1",), "IDO": ("IDO1",),
    "CD11B": ("ITGAM",), "CD11C": ("ITGAX",), "FCGR3": ("FCGR3A", "FCGR3B"),
}
RETIRED_BUILTIN.update({f"SEPT{i}": (f"SEPTIN{i}",) for i in range(1, 15) if i != 13})
RETIRED_BUILTIN.update({f"MARCH{i}": (f"MARCHF{i}",) for i in range(1, 12)})

# HGNC alias_symbol entries for common immune protein names, upper-cased, as HGNC lists them
# (November 2026 [HGNC2026]). Several targets = ambiguous: PD1 is an alias of PDCD1, SNCA and
# SPATA2, while PD-1 is an alias of PDCD1 only. A previous symbol of one gene that is an alias
# of another is ambiguous too: DEC1 (DELEC1 / BHLHE40), SEPT2 (SEPTIN2 / SEPTIN6).
ALIAS_BUILTIN = {
    "DEC1": ("BHLHE40",), "SEPT2": ("SEPTIN6",),
    "TIM3": ("HAVCR2",), "TIM-3": ("HAVCR2",), "PDL1": ("CD274",), "PD-L1": ("CD274",),
    "PD1": ("PDCD1", "SNCA", "SPATA2"), "PD-1": ("PDCD1",), "PD-L2": ("PDCD1LG2",),
    "B7H3": ("CD276",), "B7-H3": ("CD276",), "CTLA-4": ("CTLA4",), "VISTA": ("VSIR",),
    "CD16": ("FCGR3A", "FCGR3B"), "CD16A": ("FCGR3A",), "CD32": ("FCGR2A", "FCGR2B"),
    "CD56": ("NCAM1",), "NCAM": ("NCAM1",), "CD103": ("ITGAE",), "CD25": ("IL2RA", "ISG20"),
    "CD127": ("IL7R",), "CD31": ("PECAM1",), "CD73": ("NT5E",), "CD206": ("MRC1",),
    "GITR": ("TNFRSF18",), "OX40": ("TNFRSF4",), "CD137": ("TNFRSF9",), "4-1BB": ("TNFRSF9",),
    "NKP46": ("NCR1",), "NKG2D": ("KLRK1",), "TIM-1": ("HAVCR1",),
}

HGNC_COLUMNS = {
    "symbol": ("symbol", "approved_symbol"),
    "prev": ("prev_symbol", "previous_symbol", "previous_symbols"),
    "alias": ("alias_symbol", "alias_symbols"),
    "status": ("status",),
}


class InputError(ValueError):
    """The input cannot be read."""


# --------------------------------------------------------------------------- pure helpers

def classify_token(field: str):
    """Classify one gene field. Returns (kind, cleaned token, flags)."""
    flags = []
    tok = field
    if tok.strip() != tok:
        flags.append("whitespace")
        tok = tok.strip()
    if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in "\"'":
        flags.append("quoted")
        tok = tok[1:-1].strip()
    if not tok:
        return "empty", tok, flags
    if tok.lower() in MISSING_TOKENS:
        return "missing", tok, flags
    if RE_DATE_D_MON.match(tok) or RE_DATE_MON_D.match(tok) or RE_DATE_NUMERIC.match(tok):
        return "date", tok, flags
    if RE_SCIENTIFIC.match(tok) or RE_DECIMAL.match(tok):
        return "number", tok, flags
    if RE_PROBE.match(tok):
        return "probe", tok, flags
    if RE_INTEGER.match(tok):
        return "integer", tok, flags
    if RE_ENSEMBL_GENE.match(tok):
        return "ensembl_gene", tok, flags
    if RE_ENSEMBL_OTHER.match(tok):
        return "ensembl_other", tok, flags
    if RE_REFSEQ.match(tok):
        return "refseq", tok, flags
    if re.search(r"[\s/,;|]", tok):
        return "multi", tok, flags
    if RE_SYMBOL.match(tok) and re.search(r"[A-Za-z]", tok):
        return "symbol", tok, flags
    return "invalid", tok, flags


def case_style(token: str) -> str:
    """upper (HGNC style, C#orf# included), sentence (mouse style), lower, mixed or none."""
    if RE_ORF.match(token):
        return "upper"
    letters = [c for c in token if c.isalpha()]
    if not letters:
        return "none"
    if all(c.isupper() for c in letters):
        return "upper"
    if all(c.islower() for c in letters):
        return "lower"
    if letters[0].isupper() and all(c.islower() for c in letters[1:]):
        return "sentence"
    return "mixed"


def date_candidates(token: str) -> list:
    """Symbols that a spreadsheet turns into this date, limited to documented cases: SEPT2 to
    2-Sep and MARCH1 to 1-Mar [Ziemann2016], DEC1 [Abeysooriya2021]; current symbols from
    [HGNC2026], every approved gene listed when HGNC links the old symbol to several."""
    m = RE_DATE_D_MON.match(token)
    if m:
        day, mon = int(m.group(1)), m.group(2)[:3].lower()
    else:
        m = RE_DATE_MON_D.match(token)
        if not m:
            return []
        mon, day = m.group(1)[:3].lower(), int(m.group(2))
    found = []
    if mon == "sep" and 1 <= day <= 14:
        found.append(f"SEPT{day}")
    if mon == "mar" and 1 <= day <= 11:
        found.append(f"MARCH{day}")
    if mon == "dec" and day == 1:
        found.append("DEC1")
    return [{"original": o, "current": " or ".join(map_symbol(o)[1])} for o in found]


MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


def serial_note(token: str) -> str:
    """Decode a five-digit spreadsheet serial (1900 date system) and name documented candidates."""
    d = datetime.date(1899, 12, 30) + datetime.timedelta(days=int(token))
    cands = date_candidates(f"{d.day}-{MONTHS[d.month - 1]}")
    guess = " or ".join(f"{c['original']} (now {c['current']})" for c in cands)
    return f"{token} ({d.isoformat()}): " + (f"was {guess}?" if cands else "origin unknown")


def _split_multi(value: str) -> list:
    value = value.strip().strip('"')
    return [v.strip().strip('"') for v in re.split(r"[|,]", value) if v.strip().strip('"')]


def load_hgnc_table(text: str) -> dict:
    """Read a tab-separated HGNC table: the complete-set download (symbol, status, prev_symbol,
    alias_symbol; several values separated by '|') or a custom download (Approved symbol,
    Previous symbols, Alias symbols; separated by ', '). Rows whose status is not Approved are
    ignored."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise InputError("the HGNC table is empty")
    header = [h.strip().strip('"').strip().lower().replace(" ", "_") for h in lines[0].split("\t")]
    idx = {}
    for key, names in HGNC_COLUMNS.items():
        for n in names:
            if n in header:
                idx[key] = header.index(n)
                break
    if "symbol" not in idx:
        raise InputError("the HGNC table needs a 'symbol' or 'Approved symbol' column (tab-separated)")
    if "prev" not in idx and "alias" not in idx:
        raise InputError("the HGNC table needs a prev_symbol and/or alias_symbol column")
    approved, prev, alias = set(), defaultdict(set), defaultdict(set)
    for line in lines[1:]:
        fields = line.split("\t")

        def get(key):
            i = idx.get(key)
            return fields[i] if i is not None and i < len(fields) else ""

        symbol = get("symbol").strip().strip('"')
        status = get("status").strip().strip('"').lower()
        if not symbol or (status and not status.startswith("approved")):
            continue
        approved.add(symbol)
        for p in _split_multi(get("prev")):
            prev[p.upper()].add(symbol)
        for a in _split_multi(get("alias")):
            alias[a.upper()].add(symbol)
    return {"approved": approved,
            "approved_upper": {s.upper(): s for s in approved},
            "prev": {k: sorted(v) for k, v in prev.items()},
            "alias": {k: sorted(v) for k, v in alias.items()}}


def map_symbol(token: str, hgnc: dict | None = None):
    """Return (status, targets). Status: approved, case_variant, retired, ambiguous_previous,
    alias, ambiguous_alias, unknown (these need an HGNC table) or not_checked. The token is
    looked up upper-cased as written, then without hyphens. Previous-symbol and alias targets
    are pooled, previous first: one approved gene gives retired or alias, several give an
    ambiguous status (DEC1: previous symbol of DELEC1, alias of BHLHE40 [HGNC2026])."""
    up = token.upper()
    keys = [up] if "-" not in up else [up, up.replace("-", "")]
    if hgnc is not None:
        if token in hgnc["approved"]:
            return "approved", [token]
        if up in hgnc["approved_upper"]:
            return "case_variant", [hgnc["approved_upper"][up]]
        prev_table, alias_table = hgnc["prev"], hgnc["alias"]
    else:
        prev_table, alias_table = RETIRED_BUILTIN, ALIAS_BUILTIN
    for key in keys:
        prev = list(prev_table.get(key, ()))
        targets = prev + [t for t in alias_table.get(key, ()) if t not in prev]
        if len(targets) == 1:
            return ("retired" if prev else "alias"), targets
        if targets:
            return ("ambiguous_previous" if prev else "ambiguous_alias"), targets
    return ("unknown", []) if hgnc is not None else ("not_checked", [token])


def description_kind(desc: str, hgnc: dict | None = None) -> str:
    d = desc.strip()
    low = d.lower()
    if not d or low in MISSING_TOKENS or low == ".":
        return "empty"
    if RE_DOI.search(d):
        return "doi"
    if RE_PMID.match(d):
        return "pmid"
    if low.startswith(("http://", "https://")):
        return "url"
    if low == "published":
        return "published"
    if low.startswith("curated"):
        return "curated"
    if low.startswith("derived"):
        return "derived"
    if low == "user":
        return "user"
    if RE_GENE_LIKE.match(d) and d not in NOT_GENE_WORDS:
        if hgnc is None or map_symbol(d, hgnc)[0] != "unknown":
            return "gene_like"
    return "text"


# --------------------------------------------------------------------------- linting

def lint_gmt(text: str, hgnc: dict | None = None, species: str = "auto",
             min_size: int = 10, max_size: int = 500, max_items: int = 20) -> dict:
    """Lint GMT text and return the report as a dictionary."""
    if species not in ("auto", "human", "mouse"):
        raise ValueError("species must be auto, human or mouse")
    if min_size < 1 or max_size < min_size:
        raise ValueError("need 1 <= --min-size <= --max-size")
    issues = []

    def add(code, severity, message, line=None, set_name=None, items=None):
        items = [str(i) for i in (items or [])]
        issues.append({"code": code, "severity": severity, "line": line, "set": set_name,
                       "message": message, "items": items[:max_items], "n_items": len(items)})

    bom = text.startswith("﻿")
    if bom:
        text = text[1:]
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    crlf = [i for i, ln in enumerate(lines, 1) if ln.endswith("\r")]
    blank, comments, parsed = [], [], []
    for lineno, raw in enumerate(lines, 1):
        line = raw[:-1] if raw.endswith("\r") else raw
        if not line.strip():
            blank.append(lineno)
            continue
        if line.startswith("#"):
            comments.append(lineno)
            continue
        parts = line.split("\t")
        if len(parts) == 1:
            if re.search(r"\S[ ,;]+\S", line.strip()):
                add("WRONG_DELIMITER", "error",
                    "No tab on this line, but spaces, commas or semicolons between fields: GMT fields "
                    "are tab-separated, so a reader takes the whole line as one name. Re-export with tabs.",
                    lineno, None, [line.strip()[:80]])
            else:
                add("TOO_FEW_FIELDS", "error",
                    "One field only: a set needs a name, a description and at least one gene.", lineno)
            continue
        if len(parts) < 3:
            add("TOO_FEW_FIELDS", "error",
                f"{len(parts)} tab-separated fields: a set needs a name, a description and at least one "
                "gene. tmelens read_gmt skips such a line without warning [Wang2026].",
                lineno, parts[0].strip() or None)
            continue
        name = parts[0].strip()
        if not name:
            add("EMPTY_SET_NAME", "error", "The first field (the set name) is empty.", lineno)
            continue
        if name != parts[0]:
            add("SET_NAME_WHITESPACE", "warning",
                "Set name with leading or trailing spaces: matching by name will fail in some tools.",
                lineno, name)
        rec = {"name": name, "line": lineno, "description": parts[1].strip(), "n_fields": len(parts) - 2,
               "tokens": [], "empty_fields": 0, "whitespace": [], "quoted": []}
        for field in parts[2:]:
            kind, tok, flags = classify_token(field)
            if kind == "empty":
                rec["empty_fields"] += 1
                continue
            if "whitespace" in flags:
                rec["whitespace"].append(tok)
            if "quoted" in flags:
                rec["quoted"].append(tok)
            rec["tokens"].append((kind, tok))
        parsed.append(rec)

    id_counts = Counter()
    for rec in parsed:
        for kind, _ in rec["tokens"]:
            if kind in ID_KINDS:
                id_counts["entrez" if kind == "integer" else kind] += 1
    file_type = max(sorted(id_counts), key=lambda k: id_counts[k]) if id_counts else "none"

    # Mouse MHC symbols (H2-K1, H2-D1) are upper case by convention: they do not vote on species.
    styles = Counter(case_style(t) for rec in parsed for k, t in rec["tokens"]
                     if k == "symbol" and not RE_MOUSE_MHC.match(t))
    n_sym = sum(styles.values())
    if n_sym == 0:
        style_call = "none"
    elif styles["upper"] >= 0.9 * n_sym:
        style_call = "human-style"
    elif styles["sentence"] >= 0.9 * n_sym:
        style_call = "mouse-style"
    else:
        style_call = "mixed"
    if species != "auto":
        effective = species
    else:
        effective = {"human-style": "human", "mouse-style": "mouse"}.get(style_call, "unknown")
    if species == "auto" and style_call == "mixed":
        add("MIXED_CASE_STYLES", "warning",
            f"Symbol case is mixed across the file ({dict(sorted(styles.items()))}): human HGNC symbols "
            "are upper case [Bruford2020] and mouse symbols sentence case [MGI2026]. Find out which "
            "species each set comes from, map orthologues where needed, and rerun with --species.")

    report_sets, name_lines, content = [], defaultdict(list), defaultdict(list)
    for rec in parsed:
        name, ln = rec["name"], rec["line"]
        first_issue = len(issues)
        by_kind = defaultdict(list)
        for kind, tok in rec["tokens"]:
            by_kind[kind].append(tok)
        if rec["empty_fields"]:
            add("EMPTY_FIELDS", "info",
                f"{rec['empty_fields']} empty gene field(s) from doubled or trailing tabs; most readers skip "
                "them, some count them.", ln, name)
        if rec["whitespace"]:
            add("GENE_WHITESPACE", "warning",
                "Genes with leading or trailing spaces: tmelens strips them [Wang2026], but many readers do "
                "not and a padded symbol matches nothing.", ln, name, rec["whitespace"])
        if rec["quoted"]:
            add("QUOTED_GENES", "warning",
                "Genes wrapped in quotes (typical of CSV or write.table export): remove the quotes.",
                ln, name, rec["quoted"])
        if by_kind["date"]:
            items = []
            for tok in by_kind["date"]:
                cands = date_candidates(tok)
                guess = " or ".join(f"{c['original']} (now {c['current']})" for c in cands)
                items.append(f"{tok}: was {guess}?" if cands else f"{tok}: origin unknown")
            add("SPREADSHEET_DATE", "error",
                "Gene fields that are dates: spreadsheet software converted gene symbols [Ziemann2016]. "
                "Recover the symbols from the original source; the candidates listed are the documented "
                "cases only, so confirm before replacing.", ln, name, items)
        if by_kind["number"]:
            add("SPREADSHEET_NUMBER", "error",
                "Floating-point values in gene fields: spreadsheet conversion of identifiers "
                "[Ziemann2016]. Recover them from the original source.", ln, name, by_kind["number"])
        mixed = []
        if file_type != "entrez":
            serial = [t for t in by_kind["integer"] if len(t) == 5]
            mixed += [t for t in by_kind["integer"] if len(t) != 5]
            if serial:
                add("SPREADSHEET_SERIAL", "error",
                    "Five-digit numbers in a file of symbols: spreadsheets store dates as five-digit "
                    "serial numbers [Abeysooriya2021] (or Entrez IDs were mixed in). Each is shown as the "
                    "date it encodes (1900 date system) with documented candidates; recover the "
                    "identifiers from the original source.", ln, name, [serial_note(t) for t in serial])
        for kind in ("symbol", "ensembl_gene", "ensembl_other", "refseq", "probe"):
            if kind != file_type:
                mixed += by_kind[kind]
        if mixed:
            add("MIXED_IDENTIFIERS", "warning",
                f"Identifiers of another type in a file of {file_type} identifiers: they will not match "
                "the expression data. Map every gene to one identifier type.", ln, name, mixed)
        versioned = [t for t in by_kind["ensembl_gene"] + by_kind["ensembl_other"] if "." in t]
        if versioned:
            add("ENSEMBL_VERSION", "info",
                "Ensembl identifiers with version suffixes: strip the suffix on both sides before matching.",
                ln, name, versioned)
        if by_kind["multi"]:
            add("MULTI_GENE_FIELD", "warning",
                "Fields holding several identifiers (separated by '/', ',', ';', '|' or spaces): split "
                "them into separate fields or decide which gene is meant.", ln, name, by_kind["multi"])
        if by_kind["missing"]:
            add("MISSING_TOKEN", "warning",
                "Missing-value tokens (NA, NaN, #N/A) in gene fields: delete them; some readers keep "
                "them as genes.", ln, name, by_kind["missing"])
        if by_kind["invalid"]:
            add("INVALID_CHARACTERS", "warning",
                "Gene fields with characters no symbol or identifier uses: check the encoding and export.",
                ln, name, by_kind["invalid"])

        genes = [t for k, t in rec["tokens"] if k != "missing"]
        unique = list(dict.fromkeys(genes))
        dups = sorted(g for g, c in Counter(genes).items() if c > 1)
        if dups:
            add("DUPLICATE_GENE", "warning",
                "Genes listed more than once: mean and z-score scores weight them twice, rank-based "
                "tools may drop the copy. Keep one.", ln, name, dups)
        by_upper = defaultdict(set)
        for g in unique:
            by_upper[g.upper()].add(g)
        case_dups = sorted("/".join(sorted(v)) for v in by_upper.values() if len(v) > 1)
        if case_dups:
            add("CASE_DUPLICATE", "warning",
                "The same gene written in different case (for example CD8A and Cd8a): case-insensitive "
                "readers merge them, case-sensitive ones keep both. Settle the species and keep one.",
                ln, name, case_dups)

        symbols = sorted(set(by_kind["symbol"]))
        set_styles = Counter(case_style(t) for t in by_kind["symbol"])

        # Previous symbols and aliases (human nomenclature; not applied to mouse-style sets).
        renamed, mapped = {}, set()
        if effective != "mouse":
            found = defaultdict(list)
            for t in symbols:
                status, targets = map_symbol(t, hgnc)
                if status in ("retired", "alias"):
                    renamed[t] = targets[0]
                    mapped.add(t)
                    found[status].append(f"{t} -> {targets[0]}")
                elif status in ("ambiguous_previous", "ambiguous_alias"):
                    mapped.add(t)
                    found["ambiguous"].append(f"{t} -> {' or '.join(targets)}")
                elif status == "unknown" and case_style(t) == "upper":
                    found["unknown"].append(t)
            if found["retired"]:
                add("RETIRED_SYMBOL", "warning",
                    "Previous HGNC symbols: data annotated with current symbols will not match them. Use "
                    "the current symbol and keep the original in the provenance record (tmelens: "
                    "original_symbols) [HGNC2026] [Wang2026].", ln, name, found["retired"])
            if found["alias"]:
                add("ALIAS_SYMBOL", "warning",
                    "Alias symbols (often protein names), not approved HGNC symbols: use the approved "
                    "symbol [HGNC2026].", ln, name, found["alias"])
            if found["ambiguous"]:
                add("AMBIGUOUS_ALIAS", "warning",
                    "Previous symbols or aliases that HGNC links to several approved genes (for example "
                    "DEC1: previous symbol of DELEC1, alias of BHLHE40) [HGNC2026]: the symbol alone does "
                    "not say which gene was meant. Check the source; never take the first target.",
                    ln, name, found["ambiguous"])
            if found["unknown"]:
                add("UNKNOWN_SYMBOL", "info",
                    "Not an approved, previous or alias symbol in the HGNC table given: check for typos, "
                    "withdrawn symbols, non-human genes or non-gene identifiers.", ln, name, found["unknown"])
        if renamed:
            after = Counter(renamed.get(g, g) for g in unique)
            clash = sorted(g for g, c in after.items() if c > 1)
            if clash:
                add("DUPLICATE_AFTER_MAPPING", "warning",
                    "Entries that become the same gene once previous or alias symbols are renamed (for "
                    "example CTGF and CCN2): keep one.", ln, name, clash)

        lower = [t for t in symbols if case_style(t) == "lower" and t not in mapped]
        if lower:
            add("LOWERCASE_SYMBOLS", "warning",
                "All-lower-case symbols follow neither the human (upper case) nor the mouse (sentence "
                "case) convention [Bruford2020] [MGI2026]: a tool probably lower-cased them. Restore the "
                "case from the source before matching.", ln, name, lower)
        if effective == "human":
            bad = [t for t in symbols if (case_style(t) in ("sentence", "mixed") or RE_MOUSE_MHC.match(t))
                   and t not in mapped]
            if bad:
                add("SPECIES_CASE_MISMATCH", "warning",
                    "Symbols not in HGNC upper case, or mouse MHC symbols (H2-...), in a human file; sentence "
                    "case is the mouse convention [Bruford2020] [MGI2026]. If these are mouse genes, map "
                    "orthologues: upper-casing is not orthologue mapping.", ln, name, bad)
        elif effective == "mouse":
            bad = [t for t in symbols if case_style(t) == "upper" and sum(c.isalpha() for c in t) >= 2
                   and not RE_MOUSE_MHC.match(t)]
            if bad:
                add("SPECIES_CASE_MISMATCH", "warning",
                    "Upper-case symbols in a mouse file look like human symbols [MGI2026]: map orthologues "
                    "before scoring mouse data.", ln, name, bad)

        kind = description_kind(rec["description"], hgnc)
        if kind == "gene_like":
            add("DESCRIPTION_LOOKS_LIKE_GENE", "warning",
                f"The description field '{rec['description']}' looks like a gene symbol. If this line has "
                "no description column, its first gene is being read as the description and lost: insert "
                "a description (a DOI, 'curated' or a short note).", ln, name)
        elif kind == "empty":
            add("NO_PROVENANCE", "info",
                "Empty description: record where the set comes from (a DOI for a published list, "
                "'curated' for an assembled panel).", ln, name)
        elif kind == "published":
            add("PROVENANCE_WITHOUT_REFERENCE", "info",
                "Labelled 'published' without a DOI or PMID: add the reference so the list can be checked "
                "against the paper.", ln, name)

        n = len(unique)
        if n == 0:
            add("EMPTY_SET", "error", "No genes left after removing empty and missing fields.", ln, name)
        elif n == 1:
            add("SINGLE_GENE_SET", "warning",
                "A one-gene set is a single-gene readout, not a signature: report it as that gene.",
                ln, name, unique)
        elif n < min_size:
            add("SMALL_SET", "info",
                f"{n} genes, below --min-size {min_size}. Normal for a marker panel scored by mean or "
                "z-score, but one missing or mis-mapped gene moves the score, and Reimand et al. often "
                "recommend excluding sets of fewer than 10 to 15 genes from pathway enrichment analysis "
                "[Reimand2019].", ln, name)
        if n > max_size:
            add("LARGE_SET", "warning",
                f"{n} genes, above --max-size {max_size} [Reimand2019]. A large set averages many "
                "processes: in breast cancer more than 90% of random signatures of over 100 genes were "
                "significant outcome predictors, largely through proliferation [Venet2011]. Check what "
                "drives the score.", ln, name)

        name_lines[name].append(ln)
        usable = sorted({t.upper() for k, t in rec["tokens"] if k in ID_KINDS})
        if usable:
            content[tuple(usable)].append(name)
        report_sets.append({
            "name": name, "line": ln, "description": rec["description"], "provenance": kind,
            "n_fields": rec["n_fields"], "n_genes": n,
            "n_usable": len({t for k, t in rec["tokens"] if k in ID_KINDS}),
            "case_styles": dict(sorted(set_styles.items())),
            "issues": sorted({i["code"] for i in issues[first_issue:]}),
        })

    for nm, lns in name_lines.items():
        if len(lns) > 1:
            add("DUPLICATE_SET_NAME", "error",
                f"Set name used on {len(lns)} lines: dictionary-based readers such as tmelens read_gmt keep "
                "only the last definition [Wang2026]. Rename or merge deliberately.", lns[0], nm, lns)
    by_lower = defaultdict(set)
    for nm in name_lines:
        by_lower[nm.lower()].add(nm)
    for group in by_lower.values():
        if len(group) > 1:
            add("SET_NAME_CASE_DUPLICATE", "warning",
                "Set names that differ only in case: case-insensitive tools treat them as one set.",
                None, None, sorted(group))
    for names in content.values():
        distinct = list(dict.fromkeys(names))
        if len(distinct) > 1:
            add("IDENTICAL_SETS", "warning",
                "Different names, identical genes: report one, or say why both are kept.", None, None, distinct)
    if bom:
        add("BOM", "warning",
            "The file starts with a byte-order mark: readers that do not strip it put an invisible "
            "character into the first set name.", 1)
    if crlf:
        add("CRLF_LINE_ENDINGS", "info",
            "Windows line endings: readers that split on newline only leave a carriage return on the "
            "last gene of each line.", None, None, crlf)
    if blank:
        add("BLANK_LINES", "info", "Blank lines: harmless for most readers; remove them.", None, None, blank)
    if comments:
        add("COMMENT_LINE", "warning",
            "Lines starting with '#': in GMT every line defines a gene set [Reimand2019], and tmelens "
            "read_gmt loads such a line as a set when it has three or more tab-separated fields "
            "[Wang2026]. Delete these lines.", comments[0], None, comments)
    if not parsed:
        add("NO_SETS", "error", "No gene set could be read from this file.")

    issues.sort(key=lambda i: (SEVERITY_RANK[i["severity"]], i["line"] or 0, i["code"], i["set"] or ""))
    sev = Counter(i["severity"] for i in issues)
    sizes = [s["n_genes"] for s in report_sets]
    if effective == "mouse":
        symbol_check = "skipped_mouse"
    else:
        symbol_check = "hgnc_table" if hgnc is not None else "builtin_table"
    return {
        "settings": {"species": species, "min_size": min_size, "max_size": max_size,
                     "symbol_check": symbol_check},
        "format": {"n_lines": len(lines), "n_sets": len(report_sets), "bom": bom, "crlf": bool(crlf),
                   "blank_lines": len(blank), "comment_lines": len(comments)},
        "identifiers": {"file_type": file_type, "counts": dict(sorted(id_counts.items())),
                        "case_style": style_call, "case_counts": dict(sorted(styles.items())),
                        "species_assumed": effective},
        "provenance": dict(sorted(Counter(s["provenance"] for s in report_sets).items())),
        "sizes": ({"min": min(sizes), "median": statistics.median(sizes), "max": max(sizes),
                   "n_single_gene": sizes.count(1), "n_below_min": sum(1 for x in sizes if 1 < x < min_size),
                   "n_above_max": sum(1 for x in sizes if x > max_size)} if sizes else {}),
        "summary": {"errors": sev["error"], "warnings": sev["warning"], "info": sev["info"],
                    "ok": sev["error"] == 0,
                    "by_code": dict(sorted(Counter(i["code"] for i in issues).items()))},
        "issues": issues,
        "sets": report_sets,
    }


# --------------------------------------------------------------------------- command line

def read_text(path: str) -> str:
    p = Path(path)
    if not p.is_file():
        raise InputError(f"{path}: no such file")
    try:
        data = gzip.open(p, "rb").read() if p.suffix == ".gz" else p.read_bytes()
    except OSError as exc:
        raise InputError(f"{path}: cannot read ({exc})")
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InputError(f"{path}: not UTF-8 text (first bad byte at offset {exc.start}); "
                         "re-save the file as UTF-8 and lint again")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Lint a GMT gene-set file (format, names, provenance, spreadsheet damage, symbol "
                    "hygiene, species case, set size). Prints a JSON report.")
    p.add_argument("gmt", help="GMT file (tab-separated: name, description, genes); .gz accepted")
    p.add_argument("--hgnc", metavar="TSV",
                   help="HGNC table with symbol plus prev_symbol and/or alias_symbol columns (for example "
                        "hgnc_complete_set.txt); checks every symbol instead of the small built-in table")
    p.add_argument("--species", choices=("auto", "human", "mouse"), default="auto",
                   help="species of the expression data the sets will be matched to (default: inferred)")
    p.add_argument("--min-size", type=int, default=10, help="sets smaller than this are reported (default 10)")
    p.add_argument("--max-size", type=int, default=500, help="sets larger than this are reported (default 500)")
    p.add_argument("--max-items", type=int, default=20, help="examples listed per issue (default 20)")
    p.add_argument("--strict", action="store_true", help="exit with status 1 when the report has errors")
    p.add_argument("--out", help="write the JSON report here instead of standard output")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")  # legacy Windows code pages cannot print every character
    try:
        text = read_text(args.gmt)
        hgnc = load_hgnc_table(read_text(args.hgnc)) if args.hgnc else None
        report = lint_gmt(text, hgnc=hgnc, species=args.species, min_size=args.min_size,
                          max_size=args.max_size, max_items=max(1, args.max_items))
    except (InputError, ValueError) as exc:
        print(f"gmt_lint.py: error: {exc}", file=sys.stderr)
        return 2
    report = {"file": Path(args.gmt).name, **report}
    out = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        try:
            Path(args.out).write_text(out, encoding="utf-8")
        except OSError as exc:
            print(f"gmt_lint.py: error: cannot write {args.out} ({exc})", file=sys.stderr)
            return 2
    else:
        sys.stdout.write(out)
    return 1 if (args.strict and report["summary"]["errors"]) else 0


if __name__ == "__main__":
    sys.exit(main())
