#!/usr/bin/env python3
r"""Diversity and clonality metrics for T-cell receptor clonotype count tables.

The input already lists clonotypes (output of a bulk TCR-seq or single-cell V(D)J pipeline):
either one row per clonotype per sample with a count column (--count-col), or one row per
cell or contig with a cell identifier (--cell-col; clone size = number of distinct cells).
Rows that share a clonotype key within a sample are summed. When a `productive` column
exists only productive rows are kept. A `locus` or `chain` column holding more than one
receptor chain is refused unless --chain selects one, because pooling chains double-counts
cells and mixes two different distributions.

Metrics (natural logarithms; c_i = count of clonotype i, N = sum of c_i, p_i = c_i / N,
over the S clonotypes observed in the sample):
  shannon           H = -sum p_i ln p_i                                   [Shannon1948]
  pielou_evenness   J = H / ln S, undefined when S = 1                    [Pielou1966]
  clonality         1 - J (Shannon entropy normalised by ln S, inverted)  [Tumeh2014]
  simpson           lambda = sum p_i^2                                    [Simpson1949]
  inverse_simpson   1 / lambda, the Hill number of order 2                [Hill1973]
  exp_shannon       exp(H), the Hill number of order 1                    [Hill1973]
  gini_simpson      1 - lambda
  gini              Gini coefficient of the S observed clone sizes x sorted ascending,
                    sum_k (2k - S - 1) x_k / (S sum x)                    (convention)
  topK_fraction     share of N held by the K largest clonotypes
  chao1             S + f1 (f1 - 1) / (2 (f2 + 1)), f1 and f2 = clonotypes seen once and
                    twice; bias-corrected form of vegan estimateR; integer counts only [Oksanen2026]
  expanded_*        clonotypes with c_i >= 2 and their share of N (meaningful for cell counts)

Downsampling (--downsample D, or "min" for the smallest sample depth) draws D units without
replacement from each sample, repeats it --iterations times from --seed and reports the
mean, SD and 2.5 / 97.5 percentiles of every metric, plus the exact expected richness
  E[S_D] = sum_i [1 - C(N - c_i, D) / C(N, D)]                            [Hurlbert1971]
Samples below D are excluded, never upsampled. Draws use only random.random(), seeded per
sample from "<seed>:<sample name>", so results do not depend on sample order.
Without a sample column each file is one sample, named after the file without '.gz' and its
table extension; rows that end up in one sample from several files are pooled and noted.
Source keys resolve in references/sources.md of the tcr-repertoire-metrics skill.

Examples:
  python3 repertoire_metrics.py --input clones.tsv --sample-col sample_id \
      --clone-col junction_aa,v_call --count-col duplicate_count --chain TRB
  python3 repertoire_metrics.py --input contigs.csv --sample-col sample --cell-col barcode \
      --clone-col cdr3 --chain TRB --downsample min --iterations 50 --seed 7 --format tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import random
import statistics
import sys

TRUE_VALUES = {"t", "true", "1", "yes", "y"}
CHAIN_COLUMNS = ("locus", "chain")
DEFAULT_TOP = (1, 10, 100)


class InputError(ValueError):
    """Invalid input or arguments; the command line prints it and exits with code 2."""


# ----------------------------------------------------------------------------- reading

def read_table(path, sep=None):
    """Return (header, rows) of a TSV or CSV file (optionally gzipped); cells are stripped strings."""
    if not os.path.isfile(path):
        raise InputError(f"input file not found: {path}")
    gz = path.lower().endswith(".gz")
    base = path[:-3] if gz else path
    delimiter = sep or ("," if base.lower().endswith(".csv") else "\t")
    with (gzip.open if gz else open)(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, None)
        if not header or not any(h.strip() for h in header):
            raise InputError(f"{path}: the file is empty")
        header = [h.strip() for h in header]
        rows = [
            {h: (rec[i].strip() if i < len(rec) else "") for i, h in enumerate(header)}
            for rec in reader if any(cell.strip() for cell in rec)
        ]
    if len(header) == 1 and delimiter == "\t" and "," in header[0]:
        raise InputError(f"{path}: only one column found; if the file is comma-separated pass --sep ,")
    return header, rows


def _require(header, column, path, flag):
    if column not in header:
        raise InputError(f"{path}: column '{column}' ({flag}) not found; columns are: {', '.join(header)}")


TABLE_EXTENSIONS = (".tsv", ".csv", ".txt", ".tab")


def sample_name_from_path(path):
    """Sample name for a file without a sample column: drop '.gz' and one table extension only,
    so 'P1.pre.tsv' and 'P1.post.tsv' stay two samples ('P1.pre', 'P1.post')."""
    name = os.path.basename(path)
    if name.lower().endswith(".gz"):
        name = name[:-3]
    root, ext = os.path.splitext(name)
    return root if root and ext.lower() in TABLE_EXTENSIONS else name


def load_repertoires(paths, clone_cols, count_col=None, cell_col=None, sample_col=None,
                     chain_col=None, chain=None, productive_col=None, keep_nonproductive=False,
                     sep=None):
    """Read clonotype tables into ({sample: {clonotype_key: count}}, {sample: qc}, notes)."""
    if (count_col is None) == (cell_col is None):
        raise InputError("give exactly one of --count-col (a count per row) or --cell-col (one row per cell)")
    if not clone_cols:
        raise InputError("--clone-col needs at least one column name")
    samples, qc, cells, notes, sources = {}, {}, {}, [], {}
    for path in paths:
        header, rows = read_table(path, sep)
        if not rows:
            raise InputError(f"{path}: the file has a header but no data rows")
        for col in clone_cols:
            _require(header, col, path, "--clone-col")
        for col, flag in ((count_col, "--count-col"), (cell_col, "--cell-col"), (sample_col, "--sample-col"),
                          (chain_col, "--chain-col"), (productive_col, "--productive-col")):
            if col:
                _require(header, col, path, flag)
        s_col = sample_col or ("sample" if "sample" in header else None)
        default_name = sample_name_from_path(path)
        c_col = chain_col or next((c for c in CHAIN_COLUMNS if c in header), None)
        p_col = None if keep_nonproductive else (productive_col or ("productive" if "productive" in header else None))
        if c_col is None:
            if chain:
                raise InputError(f"{path}: --chain {chain} was given but there is no locus/chain column; name it with --chain-col")
            notes.append(f"{path}: no locus/chain column, so mixing of receptor chains could not be checked")
        else:
            present = {}
            for row in rows:
                present[row[c_col].upper()] = present.get(row[c_col].upper(), 0) + 1
            listing = ", ".join(f"{v or '(empty)'}: {n} rows" for v, n in sorted(present.items()))
            if chain is None and len(present) > 1:
                raise InputError(
                    f"{path}: column '{c_col}' holds several receptor chains ({listing}). Pass --chain "
                    "(for example --chain TRB) and run once per chain: pooling chains double-counts cells "
                    "and mixes two different distributions")
            if chain is not None and chain.upper() not in present:
                raise InputError(f"{path}: no rows with {c_col} = {chain}; values present: {listing}")
        if p_col is None and not keep_nonproductive:
            notes.append(f"{path}: no productive column; assuming non-productive rearrangements were removed upstream")
        for number, row in enumerate(rows, start=2):
            name = row[s_col] if s_col else default_name
            if not name:
                raise InputError(f"{path}: empty sample name in column '{s_col}' (line {number})")
            q = qc.setdefault(name, {"rows_read": 0, "rows_used": 0, "dropped_other_chain": 0,
                                     "dropped_nonproductive": 0, "dropped_missing_clonotype": 0,
                                     "dropped_zero_count": 0})
            q["rows_read"] += 1
            files = sources.setdefault(name, [])
            if path not in files:
                files.append(path)
            if c_col:
                q["chain"] = chain.upper() if chain else row[c_col].upper()
                if chain and row[c_col].upper() != chain.upper():
                    q["dropped_other_chain"] += 1
                    continue
            if p_col and row[p_col].lower() not in TRUE_VALUES:
                q["dropped_nonproductive"] += 1
                continue
            parts = [row[c] for c in clone_cols]
            if not all(parts):
                q["dropped_missing_clonotype"] += 1
                continue
            key = "|".join(parts)
            table = samples.setdefault(name, {})
            if cell_col:
                cell = row[cell_col]
                if not cell:
                    raise InputError(f"{path}: empty cell identifier in column '{cell_col}' (line {number})")
                seen = cells.setdefault(name, {}).setdefault(cell, set())
                q["rows_used"] += 1
                if key not in seen:          # a repeated contig of the same cell counts once
                    seen.add(key)
                    table[key] = table.get(key, 0) + 1
                continue
            raw = row[count_col]
            try:
                value = float(raw)
            except ValueError:
                raise InputError(f"{path}: non-numeric count '{raw}' in column '{count_col}' (line {number})") from None
            if math.isnan(value) or math.isinf(value) or value < 0:
                raise InputError(f"{path}: invalid count '{raw}' in column '{count_col}' (line {number}); counts must be finite and >= 0")
            if value == 0:
                q["dropped_zero_count"] += 1
                continue
            table[key] = table.get(key, 0) + value
            q["rows_used"] += 1
    for name, q in qc.items():
        table = samples.get(name)
        if not table:
            raise InputError(f"sample '{name}': no rows left after filtering ({json.dumps(q)})")
        integer = all(float(v).is_integer() for v in table.values())
        if integer:
            for key in table:
                table[key] = int(table[key])
        q["integer_counts"] = integer
        q["clonotypes"] = len(table)
        if cell_col:
            per_cell = cells.get(name, {})
            q["cells"] = len(per_cell)
            q["cells_in_several_clonotypes"] = sum(1 for keys in per_cell.values() if len(keys) > 1)
        else:
            q["rows_collapsed"] = q["rows_used"] - len(table)
    for name, files in sources.items():
        if len(files) > 1:
            notes.append(f"sample '{name}' pools rows from {len(files)} files ({', '.join(files)}); "
                         "check that these files really are one sample")
    return samples, qc, notes


# ----------------------------------------------------------------------------- metrics

def gini(values):
    """Gini coefficient of positive values (0 = all equal); population form, no small-sample factor."""
    xs = sorted(values)
    s, total = len(xs), float(sum(xs))
    if s == 0 or total <= 0:
        raise InputError("Gini needs at least one positive value")
    return sum((2 * k - s - 1) * x for k, x in enumerate(xs, start=1)) / (s * total)


def diversity_metrics(counts, top=DEFAULT_TOP):
    """All single-sample metrics for one list of clonotype counts (zeros are ignored)."""
    values = [c for c in counts if c > 0]
    if not values:
        raise InputError("no clonotype has a positive count")
    n = float(sum(values))
    s = len(values)
    integer = all(float(c).is_integer() for c in values)
    h = max(0.0, math.log(n) - sum(c * math.log(c) for c in values) / n)   # = -sum p ln p
    lam = sum((c / n) ** 2 for c in values)
    evenness = min(1.0, max(0.0, h / math.log(s))) if s > 1 else None
    out = {
        "depth": int(n) if integer else n,
        "richness": s,
        "shannon": h,
        "pielou_evenness": evenness,
        "clonality": None if evenness is None else 1.0 - evenness,
        "simpson": lam,
        "inverse_simpson": 1.0 / lam,
        "gini_simpson": 1.0 - lam,
        "exp_shannon": math.exp(h),
        "gini": gini(values),
    }
    ranked = sorted(values, reverse=True)
    for k in top:
        out[f"top{k}_fraction"] = sum(ranked[:k]) / n
    if integer:
        f1 = sum(1 for c in values if c == 1)
        f2 = sum(1 for c in values if c == 2)
        out.update(singletons=f1, doubletons=f2, chao1=s + f1 * (f1 - 1) / (2.0 * (f2 + 1)),
                   expanded_clonotypes=sum(1 for c in values if c >= 2),
                   expanded_fraction=sum(c for c in values if c >= 2) / n)
    else:
        out.update(singletons=None, doubletons=None, chao1=None, expanded_clonotypes=None, expanded_fraction=None)
    return out


def log_comb(a, b):
    return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)


def expected_richness(counts, depth):
    """Exact expected number of clonotypes in a random subsample of `depth` units [Hurlbert1971]."""
    counts = [int(c) for c in counts if c > 0]
    n = sum(counts)
    if not 0 < depth <= n:
        raise InputError(f"rarefaction depth {depth} must be between 1 and the sample depth {n}")
    log_all = log_comb(n, depth)
    return sum(1.0 if n - c < depth else 1.0 - math.exp(log_comb(n - c, depth) - log_all) for c in counts)


def floyd_sample(total, m, rng):
    """m distinct integers from range(total) (Floyd's algorithm, using only rng.random())."""
    chosen = set()
    for j in range(total - m, total):
        t = int(rng.random() * (j + 1))
        chosen.add(j if t in chosen else t)
    return chosen


def subsample(counts, depth, rng):
    """Draw `depth` units without replacement; return counts in the same order (zeros kept)."""
    counts = [int(c) for c in counts]
    total = sum(counts)
    if not 0 < depth <= total:
        raise InputError(f"cannot draw {depth} units from a sample of depth {total}")
    if depth == total:
        return list(counts)
    out = [0] * len(counts)
    j, upper = 0, counts[0]
    for x in sorted(floyd_sample(total, depth, rng)):
        while x >= upper:
            j += 1
            upper += counts[j]
        out[j] += 1
    return out


def quantile(sorted_values, q):
    """Linear-interpolation quantile (the default 'type 7' definition)."""
    pos = (len(sorted_values) - 1) * q
    lo = math.floor(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (pos - lo) * (sorted_values[hi] - sorted_values[lo])


def summarize(values):
    """Mean, SD, 2.5 and 97.5 percentiles of repeated values; None entries are counted, not used."""
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return None
    return {"mean": sum(vals) / len(vals), "sd": statistics.stdev(vals) if len(vals) > 1 else 0.0,
            "p2_5": quantile(vals, 0.025), "p97_5": quantile(vals, 0.975),
            "n_undefined": len(values) - len(vals)}


def resolve_depth(downsample, samples):
    """Target depth from --downsample ('min' or a positive integer); None when not downsampling."""
    if downsample is None:
        return None
    bad = [name for name, table in samples.items() if not all(isinstance(v, int) for v in table.values())]
    if bad:
        raise InputError("downsampling needs integer counts (reads, UMIs or cells); non-integer values in: "
                         + ", ".join(bad) + ". Frequencies cannot be downsampled; obtain the counts")
    depths = {name: sum(table.values()) for name, table in samples.items()}
    if str(downsample).lower() == "min":
        return min(depths.values())
    try:
        target = int(str(downsample))
    except ValueError:
        raise InputError(f"--downsample must be a positive integer or 'min', not '{downsample}'") from None
    if target < 1:
        raise InputError("--downsample must be at least 1")
    return target


def downsampled_metrics(counts, depth, iterations, rng, top=DEFAULT_TOP):
    runs = [diversity_metrics(subsample(counts, depth, rng), top) for _ in range(iterations)]
    keys = [k for k in runs[0] if k != "depth"]
    return {"depth": depth, "iterations": iterations,
            "expected_richness": expected_richness(counts, depth),
            "metrics": {k: summarize([r[k] for r in runs]) for k in keys}}


def depth_note(depths):
    if len(depths) > 1 and min(depths.values()) > 0:
        ratio = max(depths.values()) / min(depths.values())
        if ratio > 2:
            return (f"sample depths differ {ratio:.1f}-fold; richness, Chao1, Gini, clonality and top-K "
                    "fractions depend on depth, so compare samples on downsampled values")
    return None


def shallow_note(depths):
    """Heuristic flag for --downsample min: with three or more samples, the smallest depth is
    below a fifth of the median depth, so one shallow sample sets the depth for all others."""
    if len(depths) < 3:
        return None
    smallest = min(depths, key=depths.get)
    median = statistics.median(depths.values())
    if depths[smallest] < median / 5:
        return (f"the smallest depth ({depths[smallest]}, sample '{smallest}') is below a fifth of the median "
                f"depth ({median:g}); also report the comparison at a fixed higher depth that excludes the "
                "shallow sample(s) (heuristic)")
    return None


def analyse(samples, qc, notes, top=DEFAULT_TOP, downsample=None, iterations=20, seed=1):
    if not samples:
        raise InputError("no samples to analyse")
    if iterations < 1:
        raise InputError("--iterations must be at least 1")
    target = resolve_depth(downsample, samples)
    notes = list(notes)
    if target is not None and str(downsample).lower() == "min":
        note = shallow_note({name: sum(t.values()) for name, t in samples.items()})
        if note:
            notes.append(note)
    results = []
    for name, table in samples.items():
        counts = list(table.values())
        entry = {"sample": name, "qc": qc[name], "observed": diversity_metrics(counts, top)}
        if target is not None:
            depth = sum(counts)
            if depth < target:
                entry["downsampled"] = {"depth": target, "skipped": f"depth {depth} is below the target; excluded, not upsampled"}
                notes.append(f"sample '{name}' (depth {depth}) is below the downsampling depth {target} and was excluded")
            else:
                rng = random.Random(f"{seed}:{name}")
                entry["downsampled"] = downsampled_metrics(counts, target, iterations, rng, top)
        results.append(entry)
    if target is None:
        note = depth_note({name: sum(t.values()) for name, t in samples.items()})
        if note:
            notes.append(note)
    return {"tool": "repertoire_metrics.py",
            "settings": {"top": list(top), "downsample_depth": target,
                         "iterations": iterations if target is not None else None,
                         "seed": seed if target is not None else None, "log_base": "e"},
            "samples": results, "notes": notes}


# ----------------------------------------------------------------------------- output

def _round(value):
    if isinstance(value, float):
        return float(f"{value:.10g}")
    if isinstance(value, dict):
        return {k: _round(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_round(v) for v in value]
    return value


def _cell(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def to_tsv(result):
    metric_keys = list(result["samples"][0]["observed"])
    ds_keys = [k for k in metric_keys if k != "depth"]
    downsampled = any("metrics" in e.get("downsampled", {}) for e in result["samples"])
    header = ["sample", "chain"] + metric_keys
    if downsampled:
        header += ["ds_depth", "ds_iterations", "ds_expected_richness"]
        header += [f"ds_{k}_{stat}" for k in ds_keys for stat in ("mean", "sd")]
    lines = ["\t".join(header)]
    for e in result["samples"]:
        row = [e["sample"], e["qc"].get("chain")] + [e["observed"][k] for k in metric_keys]
        if downsampled:
            ds = e.get("downsampled", {})
            m = ds.get("metrics")
            row += [ds.get("depth"), ds.get("iterations"), ds.get("expected_richness")]
            for k in ds_keys:
                summary = m.get(k) if m else None
                row += [summary["mean"] if summary else None, summary["sd"] if summary else None]
        lines.append("\t".join(_cell(v) for v in row))
    return "\n".join(lines) + "\n"


def add_input_arguments(parser):
    """Arguments shared with repertoire_overlap.py."""
    parser.add_argument("--input", nargs="+", required=True, help="TSV/CSV file(s), optionally .gz; one file per sample is fine")
    parser.add_argument("--clone-col", required=True,
                        help="column(s) defining a clonotype, comma-separated (e.g. junction_aa or junction_aa,v_call)")
    unit = parser.add_mutually_exclusive_group(required=True)
    unit.add_argument("--count-col", help="column with the count per row (reads, UMIs, cells or templates)")
    unit.add_argument("--cell-col", help="column with a cell barcode; clone size = number of distinct cells")
    parser.add_argument("--sample-col", help="sample column (default: 'sample' if present, else the file name)")
    parser.add_argument("--chain-col", help="receptor chain column (default: 'locus' or 'chain' if present)")
    parser.add_argument("--chain", help="keep only this chain, e.g. TRB; required when the chain column has several values")
    parser.add_argument("--productive-col", help="productive flag column (default: 'productive' if present)")
    parser.add_argument("--keep-nonproductive", action="store_true", help="do not filter on a productive column")
    parser.add_argument("--sep", help="field separator (default: ',' for .csv, tab otherwise)")
    parser.add_argument("--downsample", help="common depth: a positive integer, or 'min' for the smallest sample depth")
    parser.add_argument("--iterations", type=int, default=20, help="downsampling repeats (default 20)")
    parser.add_argument("--seed", type=int, default=1, help="random seed for downsampling (default 1)")
    parser.add_argument("--format", choices=("json", "tsv"), default="json", help="output format (default json)")
    parser.add_argument("--out", help="write here instead of standard output")


def load_from_args(args):
    sep = "\t" if (args.sep or "").lower() in ("tab", "\\t") else args.sep
    return load_repertoires(args.input, [c.strip() for c in args.clone_col.split(",") if c.strip()],
                            count_col=args.count_col, cell_col=args.cell_col, sample_col=args.sample_col,
                            chain_col=args.chain_col, chain=args.chain, productive_col=args.productive_col,
                            keep_nonproductive=args.keep_nonproductive, sep=sep)


def write_output(text, out):
    if out:
        with open(out, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)


def parse_top(text):
    try:
        values = [int(v) for v in text.split(",") if v.strip()]
    except ValueError:
        raise InputError(f"--top must be comma-separated positive integers, not '{text}'") from None
    if not values or any(v < 1 for v in values):
        raise InputError(f"--top must be comma-separated positive integers, not '{text}'")
    return tuple(values)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Diversity and clonality metrics (Shannon, Pielou, clonality, Simpson, Gini, top-K, Chao1) "
                    "for TCR clonotype count tables, with optional seeded downsampling to a common depth.",
        epilog="Exit code 2 means the input or the arguments need fixing; the message says which.")
    add_input_arguments(parser)
    parser.add_argument("--top", default="1,10,100", help="K values for topK_fraction (default 1,10,100)")
    args = parser.parse_args(argv)
    try:
        samples, qc, notes = load_from_args(args)
        result = analyse(samples, qc, notes, top=parse_top(args.top), downsample=args.downsample,
                         iterations=args.iterations, seed=args.seed)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    text = to_tsv(result) if args.format == "tsv" else json.dumps(_round(result), indent=2) + "\n"
    if args.format == "tsv" and result["notes"]:
        for note in result["notes"]:
            print(f"note: {note}", file=sys.stderr)
    write_output(text, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
