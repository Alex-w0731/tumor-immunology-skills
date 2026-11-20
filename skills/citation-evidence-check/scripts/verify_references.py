#!/usr/bin/env python3
"""Check that cited references exist, match their citation details and carry no retraction notice.

Input is a JSON, CSV or TSV list of references with any of: doi, pmid, title, year, first_author
(aliases: author; or authors, from which the first name is taken), plus optional id, journal and
claim. Each reference is looked up at Crossref (api.crossref.org) and PubMed (NCBI E-utilities).
A reference without a working identifier, or whose identifier resolves to a work with a different
title, is also searched at Crossref by its title, so that a wrong identifier on a real paper can be
told apart from a citation whose title exists nowhere. With --abstracts the abstract is fetched
from Europe PMC and placed next to the claim.

What is computed (pure functions, unit-tested on saved responses):
  * title agreement, the Dice coefficient on title word sets [Dice1945]:
        Dice = 2 * |A & B| / (|A| + |B|)
    where A and B are the words of the cited and the registered title after folding case and
    accents, spelling out Greek letters, removing markup and a "RETRACTED:" prefix, and dropping
    stop words. >= 0.9 is a match, 0.6 to < 0.9 partial, < 0.6 a different work (this skill's
    heuristic thresholds);
  * year agreement with any publication date Crossref or PubMed holds (exact, off by one, or
    different; print and online years often differ by one);
  * first-author surname agreement, and journal agreement allowing NLM-style abbreviations;
  * notices: Crossref `updated-by` and `update-to` entries (source publisher or Retraction Watch)
    [Crossref2026b], PubMed publication types such as "Retracted Publication" or "Retraction Notice"
    and linked notices such as "Retraction in" or "Expression of concern in" [NLM2026], and
    "RETRACTED:" title prefixes. notice_status is "unknown" when the cited work was not identified
    (mismatch, not_found, lookup_failed, or a DOI registered outside Crossref);
  * evidence flags: preprint (with its published version), review, editorial or comment, trial,
    trial protocol, case report;
  * with --abstracts: numbers in the claim that do not occur in the abstract (both ends of a range
    such as 26-73% are read; numbers written as words are not compared; numbers glued to letters,
    as in PD-1, CD8 or TGF-beta1, are skipped), and the species or model-system words the abstract
    uses, including a hand-made list of murine tumour lines and strains (B16, MC38, CT26, 4T1,
    GL261, ID8, C57BL/6, BALB/c ...). These prompt a reading; they are not verdicts.

Exit status: 0 when the run completes, whatever the verdicts; 2 for input or option errors.

The script does not decide whether a paper supports a claim. "none_reported" means the services
report no notice, not that the paper is sound.

Network hosts: api.crossref.org, doi.org, eutils.ncbi.nlm.nih.gov, www.ebi.ac.uk. The User-Agent
is generic and carries no e-mail address. Requests are spaced by --delay seconds (default 0.4,
below NCBI's limit of three requests per second without an API key [NCBI2017]) and back off on
HTTP 429 [Crossref2026a]. --save DIR keeps every raw response; --offline DIR replays saved
responses without network, so a report can be reproduced and archived.

usage examples:
  python3 verify_references.py refs.csv --abstracts --save responses/ --out report.json
  python3 verify_references.py refs.json --offline responses/ --format tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

USER_AGENT = "verify-references/0.1 (reference checking script; python-urllib)"
CROSSREF_WORK = "https://api.crossref.org/works/"
CROSSREF_SEARCH = "https://api.crossref.org/works?"
DOI_HANDLE = "https://doi.org/api/handles/"
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"

TITLE_MATCH = 0.9      # heuristic
TITLE_PARTIAL = 0.6    # heuristic; tools/check_citations.py in this repository uses 0.6 as well
SEARCH_ROWS = 5

VERDICTS = ("verified", "check_title", "metadata_errors", "wrong_identifier", "identifier_conflict",
            "mismatch", "unconfirmed", "not_found", "lookup_failed")
NOTICE_STATUSES = ("retracted", "expression_of_concern", "is_notice", "corrected", "none_reported", "unknown")

STOPWORDS = frozenset("a an and the of in to for with on by from is are as at or via".split())
GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "ζ": "zeta",
         "η": "eta", "θ": "theta", "ι": "iota", "κ": "kappa", "λ": "lambda", "μ": "mu", "ν": "nu",
         "ξ": "xi", "π": "pi", "ρ": "rho", "σ": "sigma", "ς": "sigma", "τ": "tau", "φ": "phi",
         "χ": "chi", "ψ": "psi", "ω": "omega"}

COLUMN_ALIASES = {
    "id": ("id", "ref", "ref_id", "key", "label", "number"),
    "doi": ("doi",),
    "pmid": ("pmid", "pubmed", "pubmed_id"),
    "title": ("title", "article_title"),
    "year": ("year", "pub_year", "publication_year"),
    "first_author": ("first_author", "firstauthor", "author"),
    "authors": ("authors", "author_list"),
    "journal": ("journal", "journal_title", "container_title"),
    "claim": ("claim", "statement", "cited_for"),
}

CROSSREF_NOTICE_KIND = {
    "retraction": "retraction", "partial_retraction": "retraction", "withdrawal": "retraction",
    "removal": "retraction", "expression_of_concern": "concern", "correction": "correction",
    "corrigendum": "correction", "erratum": "correction", "addendum": "correction",
    "clarification": "correction",
}
# PubMed esummary `references[].reftype`, reduced to lower-case letters.
PUBMED_REF_NOTICE = {"retractionin": "retraction", "partialretractionin": "retraction",
                     "expressionofconcernin": "concern", "erratumin": "correction",
                     "correctedandrepublishedin": "correction"}
PUBMED_REF_IS_NOTICE = {"retractionof": "retraction", "partialretractionof": "retraction",
                        "expressionofconcernfor": "concern", "erratumfor": "correction"}
PUBMED_TYPE_IS_NOTICE = {"retraction of publication": "retraction", "retraction notice": "retraction",
                         "expression of concern": "concern", "published erratum": "correction"}
REVIEW_TYPES = {"review", "systematic review", "meta-analysis", "scoping review"}
COMMENT_TYPES = {"editorial", "comment", "letter", "news", "newspaper article", "interview", "biography"}
TRIAL_TYPES = {"randomized controlled trial", "controlled clinical trial", "pragmatic clinical trial",
               "equivalence trial", "adaptive clinical trial"}

RETRACTED_TITLE_RE = re.compile(r"^\s*(retracted|withdrawn)(\s+article)?\s*[:.\-]\s*", re.I)
NOTICE_TITLE_RE = re.compile(r"^\s*(?:(?:retraction(?:\s+(?:note|notice|statement))?|notice of retraction|"
                             r"expression of concern|correction|erratum|corrigendum)\s*(?::|\bto\b|\bregarding\b)"
                             r"|(?:retraction|correction|erratum|expression of concern)\s+for\b.{0,80}?\bet al\b)", re.I)
# A number may follow a hyphen only when the hyphen follows a digit (a range such as 26-73%), not a
# letter (PD-1, IL-2, COVID-19); numbers glued to letters, Greek included (CD8, TGF-β1), are skipped.
NUMBER_RE = re.compile(r"(?<![\w./])(?<![^\W\d_]-)(\d{1,3}(?:,\d{3})+|\d+(?:\.\d+)?)(\s?%|\s?percent)?(?!\w)")
MODEL_CUES = (
    ("mouse", re.compile(r"\b(mice|mouse|murine|syngeneic|C57BL/6|BALB/c|B16(?:[-.]?[A-Z0-9]+)?|MC-?38|CT-?26|4T1|EMT6|"
                         r"MB49|Renca|E0771|YUMM[0-9.]*|GL261|ID8|Pan02|Hepa1-6|Lewis lung)(?![A-Za-z0-9])", re.I)),
    ("rat", re.compile(r"\brats?\b", re.I)),
    ("zebrafish", re.compile(r"\bzebrafish\b", re.I)),
    ("non-human primate", re.compile(r"\b(macaques?|non-?human primates?|rhesus)\b", re.I)),
    ("human", re.compile(r"\b(patients?|humans?|participants?|donors?)\b", re.I)),
    ("in vitro", re.compile(r"\b(in vitro|cell lines?|cultured cells)\b", re.I)),
    ("xenograft or organoid", re.compile(r"\b(xenografts?|PDX|organoids?)\b", re.I)),
)


class InputError(Exception):
    """Problem with the input file or options (exit code 2)."""


class LookupFailed(Exception):
    """A service could not be reached or answered with an error."""


NOT_FOUND = object()   # a definite "no such record" (HTTP 404, empty search, unknown PMID)


# ----------------------------------------------------------------------------- text helpers

def clean_text(text: str) -> str:
    """Strip markup and entities, collapse whitespace."""
    text = re.sub(r"<[^>]+>", " ", text or "")
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def fold(text: str) -> str:
    text = clean_text(text).lower()
    text = "".join(GREEK.get(ch, ch) for ch in text)
    text = unicodedata.normalize("NFKD", text)
    return "".join(ch for ch in text if not unicodedata.combining(ch))


def title_words(title: str) -> set:
    title = RETRACTED_TITLE_RE.sub("", clean_text(title))
    return {w for w in re.findall(r"[a-z0-9]+", fold(title)) if w not in STOPWORDS and (len(w) > 1 or w.isdigit())}


def dice(title_a: str, title_b: str) -> float:
    """Dice coefficient on title word sets [Dice1945]; 0.0 when either title has no words."""
    a, b = title_words(title_a), title_words(title_b)
    if not a or not b:
        return 0.0
    return 2.0 * len(a & b) / (len(a) + len(b))


def _is_initials(token: str) -> bool:
    letters = token.replace(".", "").replace("-", "")
    return 0 < len(letters) <= 3 and letters.isalpha() and letters.isupper()


def surname_key(name: str) -> str:
    """Folded surname from 'Ayers M', 'M. Ayers', 'Ayers, Mark', 'van der Burg SH'."""
    s = (name or "").strip()
    if ";" in s:
        s = s.split(";")[0]
    if "," in s:
        s = s.split(",")[0]
    tokens = [t for t in s.replace(".", ". ").split() if t.lower().strip(".") not in ("et", "al")]
    while len(tokens) > 1 and _is_initials(tokens[-1]):
        tokens.pop()
    while len(tokens) > 1 and _is_initials(tokens[0]):
        tokens.pop(0)
    return " ".join(re.sub(r"[^a-z\- ]", "", fold(t)) for t in tokens).strip()


def surnames_agree(a: str, b: str) -> bool:
    if not a or not b:
        return False
    joined_a, joined_b = re.sub(r"[^a-z]", "", a), re.sub(r"[^a-z]", "", b)
    if joined_a == joined_b:
        return True
    last_a, last_b = re.split(r"[\s\-]+", a)[-1], re.split(r"[\s\-]+", b)[-1]
    return bool(last_a) and last_a == last_b


def first_from_list(authors: str) -> str:
    authors = (authors or "").strip()
    if not authors:
        return ""
    return authors.split(";")[0] if ";" in authors else authors.split(",")[0]


def compare_year(cited, years) -> str:
    if cited is None or not years:
        return "not_compared"
    diff = min(abs(int(cited) - int(y)) for y in years)
    return "exact" if diff == 0 else ("off_by_one" if diff == 1 else "disagree")


def compare_author(cited: str, names) -> str:
    key = surname_key(cited)
    names = [n for n in names if n]
    if not key or not names:
        return "not_compared"
    return "agree" if any(surnames_agree(key, surname_key(n)) for n in names) else "disagree"


JOURNAL_STOP = frozenset("the of and for in on".split())


def journal_tokens(name: str) -> list:
    return [w for w in re.findall(r"[a-z0-9]+", fold(name).replace("&", " and ")) if w not in JOURNAL_STOP]


def journals_agree(cited: str, registered: str) -> bool:
    """Same journal, allowing NLM-style abbreviations ('J Clin Invest' = 'Journal of Clinical Investigation')."""
    a, b = journal_tokens(cited), journal_tokens(registered)
    if not a or not b:
        return False
    if "".join(a) == "".join(b):
        return True
    return len(a) == len(b) and all(y.startswith(x) for x, y in zip(a, b))


def compare_journal(cited: str, names) -> str:
    names = [n for n in names if n]
    if not cited or not names:
        return "not_compared"
    return "agree" if any(journals_agree(cited, n) for n in names) else "disagree"


def numbers_in(text: str) -> list:
    """Numbers written in the text, skipping bare years; '1,234' is read as 1234."""
    found = []
    for m in NUMBER_RE.finditer(text or ""):
        value = float(m.group(1).replace(",", ""))
        if not m.group(2) and "," not in m.group(1) and value.is_integer() and 1900 <= value <= 2100:
            continue
        found.append((m.group(0).strip(), value))
    return found


def claim_numbers_missing(claim: str, abstract: str) -> list:
    """Numbers stated in the claim that do not occur anywhere in the abstract."""
    present = {v for _, v in numbers_in(abstract)}
    return [text for text, v in numbers_in(claim) if not any(abs(v - p) < 1e-9 for p in present)]


def model_cues(text: str) -> list:
    return [label for label, rx in MODEL_CUES if rx.search(text or "")]


# ----------------------------------------------------------------------------- input

def normalise_doi(raw: str):
    s = re.sub(r"^((https?://)?(dx\.|www\.)?doi\.org/|doi:\s*|doi\s+)", "", (raw or "").strip(), flags=re.I)
    s = s.strip().rstrip(".,;")
    if not s:
        return "", None
    if not re.match(r"^10\.\d{4,9}/\S+$", s):
        return "", f"DOI value {raw!r} is not a valid DOI"
    return s.lower(), None


def normalise_pmid(raw: str):
    s = re.sub(r"^pmid:?\s*", "", (raw or "").strip(), flags=re.I).rstrip(".,;")
    if not s:
        return "", None
    if not re.fullmatch(r"\d{1,9}", s):
        return "", f"PMID value {raw!r} is not a number"
    return s.lstrip("0") or s, None


def normalise_year(raw: str):
    s = (raw or "").strip()
    if not s:
        return None, None
    m = re.match(r"(1[89]\d\d|20\d\d|2100)(?!\d)", s)     # '2017', '2017a', '2017 Jun 26', '2017/06/26'
    if not m:
        return None, f"year value {raw!r} not understood; year not compared"
    return int(m.group(1)), None


def normalise_reference(row: dict, index: int) -> dict:
    lowered = {str(k).strip().lower().replace(" ", "_"): v for k, v in row.items() if k is not None}

    def pick(field):
        for alias in COLUMN_ALIASES[field]:
            value = lowered.get(alias)
            if isinstance(value, list):        # JSON author lists
                value = "; ".join(str(v) for v in value)
            if value not in (None, "") and str(value).strip():
                return str(value).strip()
        return ""

    raw_doi, raw_pmid, title = pick("doi"), pick("pmid"), clean_text(pick("title"))
    ref_id = pick("id") or str(index)
    if not (raw_doi or raw_pmid or title):
        raise InputError(f"reference {index} (id {ref_id}) needs at least one of doi, pmid or title")
    issues = []
    doi, issue = normalise_doi(raw_doi)
    issues += [issue] if issue else []
    pmid, issue = normalise_pmid(raw_pmid)
    issues += [issue] if issue else []
    year, issue = normalise_year(pick("year"))
    issues += [issue] if issue else []
    return {"id": ref_id, "doi": doi, "pmid": pmid, "title": title, "year": year,
            "first_author": pick("first_author") or first_from_list(pick("authors")),
            "journal": pick("journal"), "claim": pick("claim"), "input_issues": issues}


def read_references(path: str) -> list:
    if not os.path.isfile(path):
        raise InputError(f"input file not found: {path}")
    try:
        with open(path, encoding="utf-8-sig") as handle:
            text = handle.read()
    except UnicodeDecodeError as exc:
        raise InputError(f"{path} is not UTF-8 text ({exc.reason} at byte {exc.start}); save it as UTF-8") from None
    ext = os.path.splitext(path)[1].lower()
    if ext == ".json" or text.lstrip().startswith(("[", "{")):
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise InputError(f"{path} is not valid JSON: {exc}") from None
        if isinstance(data, dict):
            data = data.get("references")
        if not isinstance(data, list):
            raise InputError('JSON input must be a list of reference objects or {"references": [...]}')
        for i, item in enumerate(data, 1):
            if not isinstance(item, dict):
                raise InputError(f"JSON item {i} is not an object")
        rows = data
    else:
        first_line = text.splitlines()[0] if text.strip() else ""
        delimiter = "\t" if ext in (".tsv", ".tab") or ("\t" in first_line and "," not in first_line) else ","
        reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
        rows = list(reader)
        header = [h.strip().lower().replace(" ", "_") for h in (reader.fieldnames or [])]
        known = set(COLUMN_ALIASES["doi"] + COLUMN_ALIASES["pmid"] + COLUMN_ALIASES["title"])
        if rows and not known & set(header):
            raise InputError(f"no doi, pmid or title column in {path}; found columns: {', '.join(header)}")
    refs = [normalise_reference(row, i) for i, row in enumerate(rows, 1)]
    if not refs:
        raise InputError(f"no references found in {path}")
    return refs


# ----------------------------------------------------------------------------- record parsing

def _date(entry) -> str:
    parts = ((entry or {}).get("date-parts") or [[None]])[0]
    return "-".join(f"{int(p):02d}" if i else str(int(p)) for i, p in enumerate(parts) if p) if parts and parts[0] else ""


def crossref_years(msg: dict) -> list:
    years = set()
    for field in ("published-print", "published-online", "published", "issued"):
        parts = ((msg.get(field) or {}).get("date-parts") or [[None]])[0]
        if parts and parts[0]:
            years.add(int(parts[0]))
    return sorted(years)


def crossref_record(msg: dict) -> dict:
    doi = (msg.get("DOI") or "").lower()
    authors = msg.get("author") or []
    first = next((a for a in authors if a.get("sequence") == "first"), authors[0] if authors else {})
    notices, notice_for = [], []
    for entry in msg.get("updated-by") or []:
        kind = CROSSREF_NOTICE_KIND.get((entry.get("type") or "").lower())
        if kind:
            notices.append({"kind": kind, "type": entry.get("type"), "notice": (entry.get("DOI") or "").lower(),
                            "source": "crossref/" + (entry.get("source") or "publisher"), "date": _date(entry.get("updated"))})
    for entry in msg.get("update-to") or []:
        kind = CROSSREF_NOTICE_KIND.get((entry.get("type") or "").lower())
        target = (entry.get("DOI") or "").lower()
        if not kind:
            continue
        item = {"kind": kind, "type": entry.get("type"), "source": "crossref/" + (entry.get("source") or "publisher"),
                "date": _date(entry.get("updated"))}
        if target == doi:      # the article itself was replaced by its notice
            notices.append(dict(item, notice=doi))
        else:
            notice_for.append(dict(item, target=target))
    seen, unique = set(), []
    for n in notices:
        key = (n["kind"], n["notice"], n["source"])
        if key not in seen:
            seen.add(key)
            unique.append(n)
    relation = msg.get("relation") or {}
    return {
        "source": "crossref", "doi": doi, "pmid": "",
        "title": clean_text((msg.get("title") or [""])[0]),
        "journal": clean_text((msg.get("container-title") or [""])[0]),
        "journals": [clean_text(j) for j in (msg.get("container-title") or []) + (msg.get("short-container-title") or []) if j],
        "years": crossref_years(msg), "first_author": first.get("family") or first.get("name") or "",
        "type": msg.get("type") or "", "subtype": msg.get("subtype") or "", "publication_types": [],
        "notices": unique, "notice_for": notice_for,
        "published_versions": [r.get("id") for r in relation.get("is-preprint-of", []) if r.get("id")],
    }


def pubmed_record(rec: dict):
    if not rec or rec.get("error"):
        return None
    types = list(rec.get("pubtype") or [])
    notices, notice_for = [], []
    for ref in rec.get("references") or []:
        key = re.sub(r"[^a-z]", "", (ref.get("reftype") or "").lower())
        item = {"type": ref.get("reftype"), "source": "pubmed", "citation": ref.get("refsource") or ""}
        label = f"PMID {ref['pmid']}" if ref.get("pmid") else item["citation"]
        if key in PUBMED_REF_NOTICE:
            notices.append(dict(item, kind=PUBMED_REF_NOTICE[key], notice=label))
        elif key in PUBMED_REF_IS_NOTICE:
            notice_for.append(dict(item, kind=PUBMED_REF_IS_NOTICE[key], target=label))
    lowered = [t.lower() for t in types]
    if "retracted publication" in lowered and not any(n["kind"] == "retraction" for n in notices):
        notices.append({"kind": "retraction", "type": "Retracted Publication", "source": "pubmed", "notice": "", "citation": ""})
    for t in lowered:
        if t in PUBMED_TYPE_IS_NOTICE and not notice_for:
            notice_for.append({"kind": PUBMED_TYPE_IS_NOTICE[t], "type": t, "source": "pubmed", "target": "", "citation": ""})
    years = set()
    for field in ("pubdate", "epubdate"):
        m = re.match(r"\s*(\d{4})", rec.get(field) or "")
        if m:
            years.add(int(m.group(1)))
    authors = [a.get("name", "") for a in rec.get("authors") or [] if a.get("authtype", "Author") == "Author"]
    doi = next((a.get("value", "").lower() for a in rec.get("articleids") or [] if a.get("idtype") == "doi"), "")
    return {
        "source": "pubmed", "doi": doi, "pmid": str(rec.get("uid") or ""), "title": clean_text(rec.get("title") or ""),
        "journal": rec.get("fulljournalname") or rec.get("source") or "", "years": sorted(years),
        "journals": [j for j in (rec.get("fulljournalname"), rec.get("source")) if j],
        "first_author": authors[0] if authors else (rec.get("sortfirstauthor") or ""),
        "type": "", "subtype": "", "publication_types": types, "notices": notices, "notice_for": notice_for,
        "published_versions": [],
    }


def europepmc_record(data: dict):
    results = ((data or {}).get("resultList") or {}).get("result") or []
    if not results:
        return None
    r = results[0]
    abstract = re.sub(r"<h4>(.*?)</h4>", r" \1: ", r.get("abstractText") or "")
    return {"abstract": clean_text(abstract), "pmcid": r.get("pmcid") or "", "pmid": r.get("pmid") or "",
            "open_access": r.get("isOpenAccess") == "Y",
            "publication_types": list(((r.get("pubTypeList") or {}).get("pubType")) or [])}


def evidence_flags(cr, pm) -> list:
    flags = set()
    if cr:
        if cr["type"] == "posted-content" and cr["subtype"] == "preprint":
            flags.add("preprint")
        if cr["published_versions"]:
            flags.add("published_version_available")
        if cr["type"] and cr["type"] not in ("journal-article", "posted-content", "proceedings-article", "book-chapter", "book"):
            flags.add("not_a_journal_article:" + cr["type"])
    if pm:
        for t in (x.lower() for x in pm["publication_types"]):
            if t == "preprint":
                flags.add("preprint")
            elif t in REVIEW_TYPES:
                flags.add("review")
            elif t in COMMENT_TYPES:
                flags.add("editorial_or_comment")
            elif t == "clinical trial protocol":
                flags.add("trial_protocol")
            elif t.startswith("clinical trial") or t in TRIAL_TYPES:
                flags.add("clinical_trial")
            elif t == "case reports":
                flags.add("case_report")
    for rec in (cr, pm):
        if rec and RETRACTED_TITLE_RE.match(rec["title"]):
            flags.add("title_says_retracted")
        if rec and NOTICE_TITLE_RE.match(rec["title"]):
            flags.add("title_says_notice")
    return sorted(flags)


def notice_status(notices, notice_for, flags) -> str:
    kinds = {n["kind"] for n in notices}
    if "retraction" in kinds or "title_says_retracted" in flags:
        return "retracted"
    if "concern" in kinds:
        return "expression_of_concern"
    if notice_for or "title_says_notice" in flags:
        return "is_notice"
    if "correction" in kinds:
        return "corrected"
    return "none_reported"


def search_hit_acceptable(ref: dict, score: float, rec: dict) -> bool:
    """Whether a title-search hit may be treated as the cited work (this skill's heuristic).

    A hit with title agreement >= TITLE_MATCH is accepted (its year and author are then compared and
    reported). A hit between TITLE_PARTIAL and TITLE_MATCH is accepted only when the cited year is
    within one of the record's and the cited first author agrees; otherwise it is reported as the
    closest title and never substituted, because a loosely similar real paper is how an invented
    citation gets 'repaired' into a wrong one.
    """
    if score >= TITLE_MATCH:
        return True
    if score < TITLE_PARTIAL:
        return False
    return (compare_year(ref["year"], rec["years"]) in ("exact", "off_by_one")
            and compare_author(ref["first_author"], [rec["first_author"]]) == "agree")


def identifiers_conflict(ref: dict, cr, pm) -> bool:
    """DOI and PMID were both supplied, both resolved, and point at different works."""
    if not (ref["doi"] and ref["pmid"] and cr and pm):
        return False
    if pm["doi"]:
        return pm["doi"] != cr["doi"]
    return bool(cr["title"] and pm["title"]) and dice(cr["title"], pm["title"]) < TITLE_PARTIAL


def assess(ref: dict, cr=None, pm=None, *, matched_by="", doi_status="", pmid_status="",
           search_best=None, errors=(), epmc=None, other_work=None) -> dict:
    """Combine looked-up records into one verdict per reference (pure; see the verdict table in SKILL.md).

    other_work: the record a supplied identifier resolved to when a title search then found the cited
    work under another identifier; the reference is then judged against the work it describes.
    """
    issues = list(ref.get("input_issues") or [])
    if other_work:
        ident = f"DOI {ref['doi']}" if ref["doi"] else f"PMID {ref['pmid']}"
        issues.append(f"{ident} belongs to a different work: '{other_work['title']}'; the cited title was found under "
                      f"another identifier (see record)")
    if doi_status == "not_found":
        issues.append(f"DOI {ref['doi']} does not resolve (not at Crossref, not at doi.org)")
    elif doi_status == "elsewhere":
        issues.append(f"DOI {ref['doi']} resolves but is not registered with Crossref; its metadata was not compared there")
    if pmid_status == "not_found":
        issues.append(f"PMID {ref['pmid']} not found at PubMed")
    issues += ["lookup error: " + e for e in errors]
    records = [r for r in (cr, pm) if r]
    title_scores = {}
    if ref["title"]:
        title_scores = {r["source"]: round(dice(ref["title"], r["title"]), 3) for r in records if r["title"]}
    title_score = max(title_scores.values()) if title_scores else None
    years = sorted({y for r in records for y in r["years"]})
    year_check = compare_year(ref["year"], years)
    author_check = compare_author(ref["first_author"], [r["first_author"] for r in records])
    journal_check = compare_journal(ref["journal"], [j for r in records for j in r["journals"]])
    checks = (year_check, author_check, journal_check)
    conflict = identifiers_conflict(ref, cr, pm)
    wrong_id = bool(records) and ("not_found" in (doi_status, pmid_status) or bool(other_work))

    if not records:
        verdict = "unconfirmed" if doi_status == "elsewhere" else ("lookup_failed" if errors else "not_found")
    elif conflict:
        verdict = "identifier_conflict"
        issues.append(f"DOI {ref['doi']} and PMID {ref['pmid']} belong to different records")
    elif title_score is not None and title_score < TITLE_PARTIAL:
        verdict = "mismatch"
    elif title_score is None and sum(c == "disagree" for c in checks) >= 2:
        verdict = "mismatch"
    elif wrong_id:
        verdict = "wrong_identifier"
    elif "disagree" in checks:
        verdict = "metadata_errors"
    elif title_score is not None and title_score < TITLE_MATCH:
        verdict = "check_title"
    elif title_score is None and all(c == "not_compared" for c in checks):
        verdict = "unconfirmed"
        issues.append("only an identifier was supplied; compare the registered title with what the text cites")
    else:
        verdict = "verified"

    main = cr or pm
    if title_score is not None and title_score < TITLE_MATCH and main:
        issues.append(f"title agreement {title_score:.2f}; registered title: {main['title']}")
    if year_check in ("disagree", "off_by_one"):
        issues.append(f"year {ref['year']} cited; registered {', '.join(map(str, years))}")
    if author_check == "disagree":
        issues.append(f"first author '{ref['first_author']}' cited; registered '{main['first_author'] if main else ''}'")
    if journal_check == "disagree":
        registered = sorted({j for r in records for j in r["journals"]})
        issues.append(f"journal '{ref['journal']}' cited; registered {' / '.join(registered)}")
    if matched_by == "title_search" and records:
        issues.append("record found by title search, not through a supplied identifier; confirm it is the intended work")
    if search_best and (not records or verdict == "mismatch"):
        issues.append(f"no record with the cited title at Crossref; closest title found: '{search_best['title']}' "
                      f"(doi {search_best['doi']}, agreement {search_best['score']:.2f})")

    notices = [n for r in records for n in r["notices"]]
    notice_for = [n for r in records for n in r["notice_for"]]
    flags = evidence_flags(cr, pm)
    # A mismatch means the record is some other work: its notices say nothing about the cited one.
    status = notice_status(notices, notice_for, flags) if records and verdict != "mismatch" else "unknown"
    result = {
        "id": ref["id"], "verdict": verdict, "notice_status": status, "flags": flags,
        "title_score": title_score, "title_scores": title_scores, "year_check": year_check,
        "author_check": author_check, "journal_check": journal_check, "matched_by": matched_by or ("identifier" if records else ""),
        "record": None if not main else {
            "doi": (cr or {}).get("doi") or (pm or {}).get("doi") or "",
            "pmid": (pm or {}).get("pmid") or "", "title": main["title"],
            "journal": main["journal"] or (pm or {}).get("journal", ""), "years": years,
            "first_author": main["first_author"], "crossref_type": (cr or {}).get("type", ""),
            "publication_types": (pm or {}).get("publication_types", []),
        },
        "notices": notices, "notice_for": notice_for,
        "published_versions": (cr or {}).get("published_versions", []),
        "issues": issues, "input": {k: ref[k] for k in ("doi", "pmid", "title", "year", "first_author", "journal")},
    }
    if ref.get("claim"):
        result["claim"] = ref["claim"]
    if epmc is not None:
        result["abstract"] = epmc["abstract"]
        result["pmcid"] = epmc["pmcid"]
        result["open_access"] = epmc["open_access"]
        result["model_cues"] = model_cues(epmc["abstract"])
        if ref.get("claim"):
            result["claim_numbers_not_in_abstract"] = claim_numbers_missing(ref["claim"], epmc["abstract"])
    return result


# ----------------------------------------------------------------------------- network layer

def cache_name(kind: str, key: str) -> str:
    """File name under which a raw response is saved (--save) and replayed (--offline)."""
    slug = re.sub(r"[^a-z0-9._-]+", "_", key.lower()).strip("_")
    if len(slug) > 80:
        slug = slug[:60] + "-" + hashlib.sha1(key.encode("utf-8")).hexdigest()[:12]
    return f"{kind}--{slug}.json"


def http_get_json(url: str, timeout: float, retries: int = 3):
    host = urllib.parse.urlsplit(url).netloc
    last = ""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return NOT_FOUND
            last = f"HTTP {exc.code} from {host}"
            if exc.code not in (429, 500, 502, 503, 504):
                break
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            wait = float(retry_after) if retry_after and retry_after.isdigit() else 1.5 * 2 ** attempt
        except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
            last = f"{type(exc).__name__} from {host}: {exc}"
            wait = 1.5 * 2 ** attempt
        if attempt < retries - 1:
            time.sleep(min(wait, 30.0))
    raise LookupFailed(last or f"request to {host} failed")


class Fetcher:
    def __init__(self, offline_dir=None, save_dir=None, delay=0.4, timeout=30.0):
        self.offline_dir, self.save_dir, self.delay, self.timeout = offline_dir, save_dir, delay, timeout
        self._last = 0.0

    def get(self, kind: str, key: str, url: str):
        name = cache_name(kind, key)
        if self.offline_dir:
            path = os.path.join(self.offline_dir, name)
            if not os.path.isfile(path):
                raise LookupFailed(f"offline: no saved response for {kind} {key} (expected {name})")
            with open(path, encoding="utf-8") as handle:
                data = json.load(handle)
            return NOT_FOUND if isinstance(data, dict) and data.get("__status") == 404 else data
        wait = self.delay - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        try:
            data = http_get_json(url, self.timeout)
        finally:
            self._last = time.monotonic()
        if self.save_dir:
            os.makedirs(self.save_dir, exist_ok=True)
            with open(os.path.join(self.save_dir, name), "w", encoding="utf-8") as handle:
                json.dump({"__status": 404} if data is NOT_FOUND else data, handle, ensure_ascii=False)
        return data


def fetch_crossref(fetcher, doi):
    data = fetcher.get("crossref-work", doi, CROSSREF_WORK + urllib.parse.quote(doi, safe="/"))
    return data if data is NOT_FOUND else crossref_record(data.get("message") or {})


def fetch_doi_exists(fetcher, doi) -> bool:
    data = fetcher.get("doi-handle", doi, DOI_HANDLE + urllib.parse.quote(doi, safe="/"))
    return data is not NOT_FOUND and data.get("responseCode") == 1


def fetch_pubmed(fetcher, pmid):
    data = fetcher.get("pubmed-esummary", pmid, ESUMMARY + urllib.parse.urlencode(
        {"db": "pubmed", "id": pmid, "retmode": "json", "tool": "verify_references"}))
    result = (data.get("result") if data is not NOT_FOUND else None) or {}
    rec = pubmed_record(result.get(pmid))
    return rec if rec else NOT_FOUND


def fetch_pmid_for_doi(fetcher, doi):
    data = fetcher.get("pubmed-esearch", doi, ESEARCH + urllib.parse.urlencode(
        {"db": "pubmed", "term": f"{doi}[doi]", "retmode": "json", "tool": "verify_references"}))
    ids = ((data if data is not NOT_FOUND else {}).get("esearchresult") or {}).get("idlist") or []
    return ids[0] if len(ids) == 1 else NOT_FOUND


def search_crossref(fetcher, ref):
    query = " ".join(x for x in (ref["title"], ref["first_author"], str(ref["year"] or "")) if x)
    data = fetcher.get("crossref-search", query, CROSSREF_SEARCH + urllib.parse.urlencode(
        {"query.bibliographic": query, "rows": SEARCH_ROWS}))
    items = ((data if data is not NOT_FOUND else {}).get("message") or {}).get("items") or []
    best = None
    for item in items:
        rec = crossref_record(item)
        score = dice(ref["title"], rec["title"])
        if best is None or score > best[0]:
            best = (score, rec)
    return best


def fetch_europepmc(fetcher, doi, pmid):
    key, query = (f"doi:{doi}", f'DOI:"{doi}"') if doi else (f"pmid:{pmid}", f"EXT_ID:{pmid} AND SRC:MED")
    data = fetcher.get("europepmc", key, EUROPEPMC + urllib.parse.urlencode(
        {"query": query, "resultType": "core", "format": "json", "pageSize": 1}))
    return europepmc_record(data if data is not NOT_FOUND else {}) or {"abstract": "", "pmcid": "", "pmid": "",
                                                                       "open_access": False, "publication_types": []}


def lookup(ref: dict, fetcher: Fetcher, with_abstracts: bool = False) -> dict:
    errors, cr, pm = [], None, None
    doi_status = pmid_status = matched_by = ""
    search_best = None

    def guarded(fn, *args):
        try:
            return fn(*args)
        except LookupFailed as exc:
            errors.append(str(exc))
            return None

    if ref["doi"]:
        got = guarded(fetch_crossref, fetcher, ref["doi"])
        if got is NOT_FOUND:
            exists = guarded(fetch_doi_exists, fetcher, ref["doi"])
            doi_status = "elsewhere" if exists else ("not_found" if exists is not None else "")
        elif got:
            cr, doi_status = got, "found"
    if ref["pmid"]:
        got = guarded(fetch_pubmed, fetcher, ref["pmid"])
        if got is NOT_FOUND:
            pmid_status = "not_found"
        elif got:
            pm, pmid_status = got, "found"
    elif ref["doi"] and doi_status in ("found", "elsewhere"):
        pmid = guarded(fetch_pmid_for_doi, fetcher, ref["doi"])
        if pmid and pmid is not NOT_FOUND:
            got = guarded(fetch_pubmed, fetcher, pmid)
            pm = got if got and got is not NOT_FOUND else None
    if pm and not cr and pm["doi"] and pm["doi"] != ref["doi"]:
        got = guarded(fetch_crossref, fetcher, pm["doi"])
        cr = got if got and got is not NOT_FOUND else None
    other_work = None
    resolved = [r for r in (cr, pm) if r and r["title"]]
    if resolved and ref["title"] and max(dice(ref["title"], r["title"]) for r in resolved) < TITLE_PARTIAL:
        # The identifier resolves, but to a work with another title: does the cited title exist anywhere?
        best = guarded(search_crossref, fetcher, ref)
        if best and search_hit_acceptable(ref, *best) and best[1]["doi"] not in {r["doi"] for r in resolved}:
            other_work = {"title": resolved[0]["title"], "doi": resolved[0]["doi"]}
            cr, pm, matched_by = best[1], None, "title_search"
            pmid = guarded(fetch_pmid_for_doi, fetcher, cr["doi"]) if cr["doi"] else None
            if pmid and pmid is not NOT_FOUND:
                got = guarded(fetch_pubmed, fetcher, pmid)
                pm = got if got and got is not NOT_FOUND else None
        elif best:
            search_best = {"title": best[1]["title"], "doi": best[1]["doi"], "score": round(best[0], 3)}
    if not cr and not pm and ref["title"]:
        best = guarded(search_crossref, fetcher, ref)
        if best:
            score, rec = best
            if search_hit_acceptable(ref, score, rec):
                cr, matched_by = rec, "title_search"
                pmid = guarded(fetch_pmid_for_doi, fetcher, rec["doi"]) if rec["doi"] else None
                if pmid and pmid is not NOT_FOUND:
                    got = guarded(fetch_pubmed, fetcher, pmid)
                    pm = got if got and got is not NOT_FOUND else None
            else:
                search_best = {"title": rec["title"], "doi": rec["doi"], "score": round(score, 3)}
    epmc = None
    if with_abstracts and (cr or pm):
        doi = (cr or {}).get("doi") or (pm or {}).get("doi")
        epmc = guarded(fetch_europepmc, fetcher, doi, (pm or {}).get("pmid"))
    return assess(ref, cr, pm, matched_by=matched_by, doi_status=doi_status, pmid_status=pmid_status,
                  search_best=search_best, errors=errors, epmc=epmc, other_work=other_work)


# ----------------------------------------------------------------------------- output

def summarise(results: list) -> dict:
    return {"references": len(results),
            "verdicts": {v: sum(r["verdict"] == v for r in results) for v in VERDICTS if any(r["verdict"] == v for r in results)},
            "notice_status": {s: sum(r["notice_status"] == s for r in results) for s in NOTICE_STATUSES
                              if any(r["notice_status"] == s for r in results)}}


TSV_COLUMNS = ("id", "verdict", "notice_status", "title_score", "year_check", "author_check", "journal_check", "record_doi",
               "record_pmid", "record_years", "record_first_author", "record_title", "flags", "notices", "issues")


def to_tsv(results: list) -> str:
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_COLUMNS)
    for r in results:
        rec = r["record"] or {}
        notices = [f"{n['kind']} ({n['source']}{', ' + n['date'] if n.get('date') else ''}{', ' + n['notice'] if n.get('notice') else ''})"
                   for n in r["notices"]] + [f"this item is a {n['kind']} notice for {n.get('target') or '?'}" for n in r["notice_for"]]
        writer.writerow([r["id"], r["verdict"], r["notice_status"], "" if r["title_score"] is None else f"{r['title_score']:.3f}",
                         r["year_check"], r["author_check"], r["journal_check"], rec.get("doi", ""), rec.get("pmid", ""),
                         ",".join(map(str, rec.get("years", []))), rec.get("first_author", ""), rec.get("title", ""),
                         "; ".join(r["flags"]), "; ".join(notices), "; ".join(r["issues"])])
    return out.getvalue()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Check that references resolve at Crossref/PubMed, match their citation details, and carry no "
                    "retraction or correction notice. Does not judge whether a paper supports a claim.",
        epilog="Input columns (case-insensitive): doi, pmid, title, year, first_author (or author/authors), "
               "optional id, journal, claim. JSON: a list of objects with the same keys, or {\"references\": [...]}. "
               "Exit status: 0 when the run completes, whatever the verdicts; 2 for input or option errors.")
    parser.add_argument("input", help="JSON, CSV or TSV file of references")
    parser.add_argument("--format", choices=("json", "tsv"), default="json", help="report format (default json)")
    parser.add_argument("--out", help="write the report to this file instead of stdout")
    parser.add_argument("--abstracts", action="store_true", help="also fetch abstracts from Europe PMC for claim checking")
    parser.add_argument("--save", metavar="DIR", help="save every raw API response in DIR (provenance, later --offline runs)")
    parser.add_argument("--offline", metavar="DIR", help="read saved responses from DIR instead of the network")
    parser.add_argument("--delay", type=float, default=0.4, help="seconds between requests (default 0.4)")
    parser.add_argument("--timeout", type=float, default=30.0, help="seconds per request (default 30)")
    args = parser.parse_args(argv)
    try:
        if args.offline and args.save:
            raise InputError("--offline and --save cannot be combined")
        if args.offline and not os.path.isdir(args.offline):
            raise InputError(f"--offline folder not found: {args.offline}")
        if args.delay < 0 or args.timeout <= 0:
            raise InputError("--delay must be >= 0 and --timeout > 0")
        if args.out and not os.path.isdir(os.path.dirname(os.path.abspath(args.out))):
            raise InputError(f"--out folder does not exist: {os.path.dirname(os.path.abspath(args.out))}")
        if args.save and os.path.exists(args.save) and not os.path.isdir(args.save):
            raise InputError(f"--save must be a folder, but {args.save} is a file")
        refs = read_references(args.input)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    fetcher = Fetcher(offline_dir=args.offline, save_dir=args.save, delay=args.delay, timeout=args.timeout)
    results = [lookup(ref, fetcher, args.abstracts) for ref in refs]
    report = {"tool": "verify_references.py", "mode": "offline" if args.offline else "online",
              "thresholds": {"title_match": TITLE_MATCH, "title_partial": TITLE_PARTIAL},
              "summary": summarise(results), "results": results}
    text = to_tsv(results) if args.format == "tsv" else json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)
    s = report["summary"]
    print(f"{s['references']} references; verdicts: {s['verdicts']}; notices: {s['notice_status']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
