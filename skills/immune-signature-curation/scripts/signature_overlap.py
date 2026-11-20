#!/usr/bin/env python3
"""Pairwise overlap between gene sets: Jaccard index and overlap coefficient, with flags for
identical, nested, near-duplicate and overlapping sets.

Formulas, for gene sets A and B after duplicate genes are removed (both coefficients are used to
measure gene-set redundancy in EnrichmentMap [Merico2010]):

  Jaccard index        J  = |A and B| / |A or B|          (shared genes / union)
  overlap coefficient  OC = |A and B| / min(|A|, |B|)     (shared genes / smaller set)

OC is 1 whenever the smaller set lies inside the larger one; J stays small when the sizes differ,
even for nested sets. Report both.

Flags (the first that applies):
  empty           a set has no genes (for example after --universe)
  identical       same genes (J = 1)
  nested          the smaller set lies entirely inside the larger (OC = 1)
  near_duplicate  J >= --jaccard (default 0.5; a heuristic of this script)
  overlapping     OC >= --overlap (default 0.5, inclusive; borrowed from EnrichmentMap, whose
                  use cases set the overlap coefficient cut-off to 0.5 [Merico2010] and whose
                  slider preset is written 'overlap > 0.5' [Reimand2019]; a setting for linking
                  related gene sets in a visualisation, not a test)
  shared          at least one shared gene, below both thresholds
  distinct        no shared gene (listed only with --all-pairs)

Without --against every pair within one file is compared; with --against each set of the first
file is compared with each set of the second (for example a local copy against the published
list). --universe restricts every set to the genes measured in the expression data before
comparing, because genes that are not measured cannot make two scores agree. Gene identifiers
are compared exactly; --ignore-case upper-cases them first (upper-casing is not orthologue
mapping). Lines starting with '#' and lines with fewer than three tab-separated fields are
skipped and counted. Duplicate set names within a file stop the run: lint the file first.

Usage examples:
  python3 signature_overlap.py panel.gmt --out overlap.json
  python3 signature_overlap.py my_sets.gmt --against published.gmt --universe measured_genes.txt --format tsv

Exit status: 0 report written; 2 unreadable input, duplicate set names, unwritable --out or bad
arguments.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

FLAG_ORDER = ("identical", "nested", "near_duplicate", "overlapping", "shared", "empty", "distinct")
TSV_COLUMNS = ("set_a", "set_b", "n_a", "n_b", "n_shared", "jaccard", "overlap_coefficient",
               "flag", "shared_genes")


class InputError(ValueError):
    """The input cannot be used."""


# --------------------------------------------------------------------------- pure functions

def parse_gmt(text: str, source: str = "gmt") -> tuple[dict, list]:
    """Return ({name: [genes]} in file order, skipped line numbers)."""
    sets: dict = {}
    skipped: list = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.rstrip("\r")
        if lineno == 1:
            line = line.lstrip("﻿")
        if not line.strip():
            continue
        if line.startswith("#"):
            skipped.append(lineno)
            continue
        parts = line.split("\t")
        name = parts[0].strip()
        if len(parts) < 3 or not name:
            skipped.append(lineno)
            continue
        if name in sets:
            raise InputError(f"{source}: set name {name!r} appears more than once (again on line {lineno}); "
                             "run gmt_lint.py and make the names unique")
        sets[name] = [g.strip() for g in parts[2:] if g.strip()]
    return sets, skipped


def parse_universe(text: str) -> set:
    """One identifier per line (the first tab-separated field); blank and '#' lines ignored."""
    genes = set()
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        first = line.split("\t")[0].strip()
        if first:
            genes.add(first)
    return genes


def prepare(sets: dict, universe: set | None = None, ignore_case: bool = False) -> dict:
    """Deduplicate, optionally upper-case and restrict to the universe.
    Returns {name: {"genes": set, "n_input": int, "dropped": sorted list}}."""
    uni = None if universe is None else ({g.upper() for g in universe} if ignore_case else set(universe))
    out = {}
    for name, genes in sets.items():
        g = {x.upper() for x in genes} if ignore_case else set(genes)
        kept = g if uni is None else g & uni
        out[name] = {"genes": kept, "n_input": len(g), "dropped": sorted(g - kept)}
    return out


def overlap_stats(a: set, b: set) -> dict:
    """Jaccard index and overlap coefficient [Merico2010]; None when undefined (an empty set)."""
    shared = a & b
    n_a, n_b, n_s = len(a), len(b), len(shared)
    union = n_a + n_b - n_s
    small = min(n_a, n_b)
    return {"n_a": n_a, "n_b": n_b, "n_shared": n_s,
            "jaccard": (n_s / union) if (union and small) else None,
            "overlap_coefficient": (n_s / small) if small else None,
            "shared_genes": sorted(shared)}


def classify(stats: dict, jaccard_min: float = 0.5, overlap_min: float = 0.5) -> str:
    """Flag a pair; identity and nesting are decided on integer counts, not floats."""
    n_a, n_b, n_s = stats["n_a"], stats["n_b"], stats["n_shared"]
    if n_a == 0 or n_b == 0:
        return "empty"
    if n_s == n_a == n_b:
        return "identical"
    if n_s == min(n_a, n_b):
        return "nested"
    if stats["jaccard"] >= jaccard_min:
        return "near_duplicate"
    if stats["overlap_coefficient"] >= overlap_min:
        return "overlapping"
    return "shared" if n_s else "distinct"


def _clusters(names: list, links: list) -> list:
    """Connected components (size >= 2) of identical / near-duplicate links, in input order."""
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in links:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    groups: dict = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return [g for g in groups.values() if len(g) > 1]


def compare(sets_a: dict, sets_b: dict | None = None, jaccard_min: float = 0.5,
            overlap_min: float = 0.5, all_pairs: bool = False) -> dict:
    """Compare prepared sets (output of prepare). Within one collection when sets_b is None."""
    if not (0 < jaccard_min <= 1 and 0 < overlap_min <= 1):
        raise ValueError("--jaccard and --overlap must be in (0, 1]")
    names_a = list(sets_a)
    if sets_b is None:
        pairs_idx = [(names_a[i], names_a[j]) for i in range(len(names_a)) for j in range(i + 1, len(names_a))]
        pool_b = sets_a
    else:
        pairs_idx = [(a, b) for a in names_a for b in sets_b]
        pool_b = sets_b
    pairs, best, links = [], {}, []
    for a, b in pairs_idx:
        st = overlap_stats(sets_a[a]["genes"], pool_b[b]["genes"])
        flag = classify(st, jaccard_min, overlap_min)
        row = {"set_a": a, "set_b": b, **st, "flag": flag}
        for me, other, side in ((a, b, "a"), (b, a, "b")):
            if sets_b is not None and side == "b":
                continue
            key = (st["jaccard"] or 0.0, st["overlap_coefficient"] or 0.0)
            cur = best.get(me)
            if st["n_shared"] and (cur is None or key > cur[0] or (key == cur[0] and other < cur[1])):
                best[me] = (key, other)
        if sets_b is None and flag in ("identical", "near_duplicate"):
            links.append((a, b))
        if flag != "distinct" or all_pairs:
            pairs.append(row)
    pairs.sort(key=lambda r: (FLAG_ORDER.index(r["flag"]), -(r["jaccard"] or 0.0),
                              -(r["overlap_coefficient"] or 0.0), r["set_a"], r["set_b"]))
    for r in pairs:
        for k in ("jaccard", "overlap_coefficient"):
            if r[k] is not None:
                r[k] = round(r[k], 4)
    set_rows = []
    for name in names_a:
        info = sets_a[name]
        b = best.get(name)
        set_rows.append({"name": name, "n_input": info["n_input"], "n_used": len(info["genes"]),
                         "dropped_by_universe": info["dropped"][:20],
                         "best_match": None if b is None else {"set": b[1], "jaccard": round(b[0][0], 4),
                                                               "overlap_coefficient": round(b[0][1], 4)}})
    by_flag: dict = {}
    for r in pairs:
        by_flag[r["flag"]] = by_flag.get(r["flag"], 0) + 1
    return {"sets": set_rows, "pairs": pairs,
            "clusters": _clusters(names_a, links) if sets_b is None else [],
            "summary": {"n_sets_a": len(names_a), "n_sets_b": None if sets_b is None else len(sets_b),
                        "n_pairs_compared": len(pairs_idx), "n_pairs_listed": len(pairs),
                        "by_flag": dict(sorted(by_flag.items()))}}


def to_tsv(report: dict) -> str:
    lines = ["\t".join(TSV_COLUMNS)]
    for r in report["pairs"]:
        vals = [r["set_a"], r["set_b"], r["n_a"], r["n_b"], r["n_shared"],
                "" if r["jaccard"] is None else f"{r['jaccard']:.4f}",
                "" if r["overlap_coefficient"] is None else f"{r['overlap_coefficient']:.4f}",
                r["flag"], ",".join(r["shared_genes"])]
        lines.append("\t".join(str(v) for v in vals))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- command line

def read_text(path: str) -> str:
    p = Path(path)
    if not p.is_file():
        raise InputError(f"{path}: no such file")
    try:
        data = gzip.open(p, "rb").read() if p.suffix == ".gz" else p.read_bytes()
        return data.decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InputError(f"{path}: cannot read as UTF-8 text ({exc})")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pairwise Jaccard index and overlap coefficient between gene sets in GMT files; "
                    "flags identical, nested, near-duplicate and overlapping sets.")
    p.add_argument("gmt", help="GMT file (tab-separated: name, description, genes); .gz accepted")
    p.add_argument("--against", metavar="GMT", help="compare each set of GMT with each set of this file")
    p.add_argument("--universe", metavar="FILE",
                   help="measured genes, one per line (first tab-separated field); sets are restricted to them")
    p.add_argument("--ignore-case", action="store_true", help="upper-case all identifiers before comparing")
    p.add_argument("--jaccard", type=float, default=0.5, help="near-duplicate threshold on J (default 0.5)")
    p.add_argument("--overlap", type=float, default=0.5, help="overlapping threshold on OC (default 0.5)")
    p.add_argument("--all-pairs", action="store_true", help="also list pairs without shared genes")
    p.add_argument("--format", choices=("json", "tsv"), default="json", help="output format (default json)")
    p.add_argument("--out", help="write the report here instead of standard output")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")  # legacy Windows code pages cannot print every character
    try:
        sets_a, skipped_a = parse_gmt(read_text(args.gmt), Path(args.gmt).name)
        sets_b, skipped_b = (parse_gmt(read_text(args.against), Path(args.against).name)
                             if args.against else (None, []))
        universe = parse_universe(read_text(args.universe)) if args.universe else None
        if universe is not None and not universe:
            raise InputError(f"{args.universe}: no gene identifiers found")
        prep_a = prepare(sets_a, universe, args.ignore_case)
        prep_b = prepare(sets_b, universe, args.ignore_case) if sets_b is not None else None
        report = compare(prep_a, prep_b, args.jaccard, args.overlap, args.all_pairs)
    except (InputError, ValueError) as exc:
        print(f"signature_overlap.py: error: {exc}", file=sys.stderr)
        return 2
    report = {"files": {"a": Path(args.gmt).name, "b": Path(args.against).name if args.against else None},
              "settings": {"jaccard_min": args.jaccard, "overlap_min": args.overlap,
                           "ignore_case": args.ignore_case, "all_pairs": args.all_pairs,
                           "universe_size": None if universe is None else len(universe)},
              "skipped_lines": {"a": skipped_a, "b": skipped_b}, **report}
    out = to_tsv(report) if args.format == "tsv" else json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        try:
            Path(args.out).write_text(out, encoding="utf-8")
        except OSError as exc:
            print(f"signature_overlap.py: error: cannot write {args.out} ({exc})", file=sys.stderr)
            return 2
    else:
        sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
