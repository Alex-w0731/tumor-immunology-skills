#!/usr/bin/env python3
"""Check every citation in skills/*/references/sources.md against Crossref and PubMed.

For each source entry (a list item that starts with a bold key such as **[Ayers2017]**):

  * the DOI must resolve at Crossref, and the record's year, first author and title must
    agree with what the entry says;
  * the PMID must exist at PubMed, and when PubMed lists a DOI it must equal the entry's DOI.

This catches the two citation failures that matter most in a methods repository: references
that do not exist, and identifiers that point at a different paper than the one described.
It cannot tell whether a paper really supports the claim it is cited for; that is what human
and agent review are for (see CONTRIBUTING.md).

Needs network access, so it is not part of the default test run. A separate workflow runs it
whenever a bibliography changes.

usage:
  python3 tools/check_citations.py                  # all skills
  python3 tools/check_citations.py --skill NAME
  python3 tools/check_citations.py --json report.json

On python.org builds for macOS without a CA bundle, run with SSL_CERT_FILE=/etc/ssl/cert.pem.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
USER_AGENT = "citation-check/1.0 (python-urllib)"

ENTRY_RE = re.compile(r"^\s*[-*]\s+\*\*\[([A-Z][A-Za-z'’\-]*((?:19|20)\d{2})[a-z]?)\]\*\*\s*(.*)$")
# A DOI may contain balanced parentheses (Elsevier SICI DOIs such as 10.1016/0895-4356(95)00048-8),
# so a parenthesised group is accepted while an unbalanced closing parenthesis ends the match.
DOI_RE = re.compile(r"\b(10\.\d{4,9}/(?:[^\s\"<>()\[\]]|\([^\s\"<>()]*\))+)")
PMID_RE = re.compile(r"\bPMID:?\s*(\d{6,9})\b")
STOPWORDS = {"a", "an", "and", "the", "of", "in", "to", "for", "with", "on", "by", "from", "is", "are", "as", "at", "or", "via"}


def fold(text: str) -> str:
    text = unicodedata.normalize("NFKD", text)
    text = "".join(c for c in text if not unicodedata.combining(c))
    return text.lower()


def words(text: str) -> set[str]:
    text = re.sub(r"<[^>]+>", " ", text)            # Crossref titles may carry markup
    text = re.sub(r"&[a-z]+;", " ", text)
    return {w for w in re.findall(r"[a-z0-9]+", fold(text)) if w not in STOPWORDS and len(w) > 1}


def title_overlap(record_title: str, entry_text: str) -> float:
    """Share of the record's title words that also occur in the entry text."""
    title_words = words(record_title)
    if not title_words:
        return 0.0
    return len(title_words & words(entry_text)) / len(title_words)


def parse_sources(path: Path) -> list[dict]:
    entries: list[dict] = []
    current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        m = ENTRY_RE.match(line)
        if m:
            current = {"key": m.group(1), "year": int(m.group(2)), "text": m.group(3), "file": str(path.relative_to(ROOT))}
            entries.append(current)
        elif current is not None and line.strip() and not line.lstrip().startswith("#"):
            current["text"] += " " + line.strip()
        else:
            current = None
    for e in entries:
        doi = DOI_RE.search(e["text"])
        e["doi"] = doi.group(1).rstrip(".,;") if doi else None
        pmid = PMID_RE.search(e["text"])
        e["pmid"] = pmid.group(1) if pmid else None
    return entries


def fetch_json(url: str, retries: int = 3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=40) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            last = exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"request failed after {retries} attempts: {url} ({last})")


def crossref_year(msg: dict) -> list[int]:
    years = []
    for field in ("published-print", "published-online", "published", "issued", "created"):
        parts = (msg.get(field) or {}).get("date-parts") or []
        if parts and parts[0] and parts[0][0]:
            years.append(int(parts[0][0]))
    return years


def check_doi(entry: dict) -> list[str]:
    problems: list[str] = []
    data = fetch_json("https://api.crossref.org/works/" + urllib.parse.quote(entry["doi"], safe="/"))
    if data is None:
        # DataCite and other agencies are not in Crossref: fall back to the DOI resolver.
        other = fetch_json("https://doi.org/api/handles/" + urllib.parse.quote(entry["doi"], safe="/"))
        if other is None or other.get("responseCode") != 1:
            problems.append(f"DOI {entry['doi']} does not resolve")
        else:
            entry["note"] = "DOI resolves but is not registered with Crossref; metadata not compared"
        return problems
    msg = data.get("message", {})
    title = (msg.get("title") or [""])[0]
    years = crossref_year(msg)
    if years and all(abs(entry["year"] - y) > 1 for y in years):
        problems.append(f"year {entry['year']} in the key vs {sorted(set(years))} at Crossref")
    authors = msg.get("author") or []
    family = fold(authors[0].get("family") or authors[0].get("name") or "") if authors else ""
    if family and family.replace(" ", "") not in fold(entry["text"]).replace(" ", "") and family.split()[-1] not in fold(entry["text"]):
        problems.append(f"first author at Crossref is '{authors[0].get('family') or authors[0].get('name')}', not found in the entry")
    overlap = title_overlap(title, entry["text"])
    if title and overlap < 0.6:
        problems.append(f"title mismatch (overlap {overlap:.2f}); Crossref title: {title[:110]}")
    entry["crossref_title"] = title
    return problems


def check_pmids(entries: list[dict]) -> dict[str, list[str]]:
    problems: dict[str, list[str]] = {}
    with_pmid = [e for e in entries if e["pmid"]]
    for start in range(0, len(with_pmid), 150):
        chunk = with_pmid[start:start + 150]
        ids = ",".join(e["pmid"] for e in chunk)
        data = fetch_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={ids}&retmode=json&tool=citation-check")
        result = (data or {}).get("result", {})
        for e in chunk:
            rec = result.get(e["pmid"])
            issues: list[str] = []
            if not rec or rec.get("error"):
                issues.append(f"PMID {e['pmid']} not found at PubMed")
            else:
                dois = [a.get("value", "").lower() for a in rec.get("articleids", []) if a.get("idtype") == "doi"]
                if e["doi"] and dois and e["doi"].lower() not in dois:
                    issues.append(f"PMID {e['pmid']} belongs to DOI {dois[0]}, but the entry cites {e['doi']}")
                overlap = title_overlap(rec.get("title", ""), e["text"])
                if rec.get("title") and overlap < 0.6:
                    issues.append(f"PMID {e['pmid']} title mismatch (overlap {overlap:.2f}); PubMed title: {rec.get('title', '')[:110]}")
            if issues:
                problems[e["file"] + "::" + e["key"]] = issues
        time.sleep(0.4)
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify DOIs and PMIDs cited by the skills.")
    parser.add_argument("--skill", help="check a single skill by folder name")
    parser.add_argument("--json", metavar="FILE", help="also write a JSON report")
    parser.add_argument("--delay", type=float, default=0.15, help="seconds between Crossref requests")
    args = parser.parse_args(argv)

    files = sorted((ROOT / "skills").glob("*/references/sources.md"))
    if args.skill:
        files = [f for f in files if f.parts[-3] == args.skill]
    if not files:
        print("no references/sources.md files found", file=sys.stderr)
        return 2
    entries: list[dict] = []
    for f in files:
        entries.extend(parse_sources(f))

    report: dict[str, list[str]] = {}
    cache: dict[str, list[str]] = {}
    no_identifier = []
    for e in entries:
        ident = e["file"] + "::" + e["key"]
        if not e["doi"] and not e["pmid"]:
            no_identifier.append(ident)
            continue
        if e["doi"]:
            try:
                cache_key = e["doi"].lower() + "|" + e["key"]
                if cache_key not in cache:
                    cache[cache_key] = check_doi(e)
                    time.sleep(args.delay)
                if cache[cache_key]:
                    report.setdefault(ident, []).extend(cache[cache_key])
            except RuntimeError as exc:
                report.setdefault(ident, []).append(f"could not check DOI: {exc}")
    try:
        for ident, issues in check_pmids(entries).items():
            report.setdefault(ident, []).extend(issues)
    except RuntimeError as exc:
        report.setdefault("pubmed", []).append(f"could not check PMIDs: {exc}")

    checked = len(entries) - len(no_identifier)
    print(f"{len(entries)} source entries in {len(files)} file(s); {checked} with a DOI or PMID")
    for ident in no_identifier:
        print(f"  note   {ident}: no DOI or PMID (URL-only source, not checked)")
    for ident, issues in sorted(report.items()):
        for issue in issues:
            print(f"  FAIL   {ident}: {issue}")
    print(f"\n{len(report)} entr{'y' if len(report) == 1 else 'ies'} with problems")
    if args.json:
        Path(args.json).write_text(json.dumps({"entries": len(entries), "checked": checked, "problems": report}, indent=2, ensure_ascii=False), encoding="utf-8")
    return 1 if report else 0


if __name__ == "__main__":
    sys.exit(main())
