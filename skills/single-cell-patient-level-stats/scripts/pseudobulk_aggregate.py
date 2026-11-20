#!/usr/bin/env python3
"""Sum raw counts per sample x cell type (pseudobulk) and count the cells behind each sum.

Pseudobulk aggregation as in [Lun2017] and [Crowell2020]: the raw integer counts (UMIs or
reads) of all cells that share a sample and a cell type are summed gene by gene, giving one
count profile per sample x cell type for a negative-binomial or limma-voom model in which
the sample (patient) is the replicate. Pseudobulks built from fewer than --min-cells cells
are left out of the count matrix and kept, flagged, in the metadata table; a sample without
cells of a type gets no pseudobulk for it (none is invented with zero counts). Sample-level
labels (--group-col, --patient-col, --covariates) must be constant within a sample;
otherwise the sample column is not the unit of replication or a label is wrong, and the
script stops. Non-integer or negative counts stop the script: they mean normalised, scaled
or integrated values, which a count model must not receive.

The JSON summary reports, per cell type, how many independent units (patients when
--patient-col is given, otherwise samples) pass the filter in each group, flags a cell type
as testable only when every group reaches --min-units, and lists patients that contribute
more than one sample to the same group.

Input formats
  wide (default)  one row per cell: sample, cell type, optional cell id, labels, then one
                  column per gene with raw counts. Columns that are neither labels nor genes
                  must be named in --exclude-cols. Integer QC columns would otherwise be summed
                  as genes, so the script stops when a gene column has a common per-cell
                  metadata name (nCount_RNA, nFeature_RNA, n_genes, total_counts, n_counts,
                  nUMI, seurat_clusters, leiden, ...); the JSON lists the first gene columns
                  so that other metadata columns can be spotted.
  long            one row per cell x gene (zero rows may be omitted): sample, cell type,
                  cell id, gene and count columns. Cells without any row are not counted.

Outputs
  --out        gene x pseudobulk count matrix (TSV); columns named <sample>::<cell type>,
               sorted by cell type, then sample; only pseudobulks passing --min-cells.
  --meta-out   one row per sample x cell type: pseudobulk_id, sample, cell_type, n_cells,
               total_counts, kept, then the label columns. Its n_cells column is the input
               that composition_test.py expects with --count-col n_cells.
  stdout       JSON summary.

Usage:
  python3 pseudobulk_aggregate.py cells.tsv --sample-col sample_id --celltype-col cell_type \\
      --cell-col cell_id --group-col condition --patient-col patient_id --covariates batch \\
      --min-cells 10 --out pb_counts.tsv --meta-out pb_meta.tsv > pb_summary.json
  python3 pseudobulk_aggregate.py umi_long.tsv --format long --cell-col barcode \\
      --gene-col gene --count-col umi --group-col condition --out pb_counts.tsv

Exit codes: 0 success, 2 input problem (the message names the line, column or value).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys

ID_SEP = "::"
# Per-cell metadata columns from Seurat and scanpy exports that hold integers and would
# silently pass as genes (compared in lower case). None of these is an HGNC or MGI symbol.
QC_NAMES = {"ncount_rna", "nfeature_rna", "n_genes", "n_genes_by_counts", "total_counts", "n_counts",
            "numi", "ngene", "seurat_clusters", "leiden", "louvain", "orig.ident", "phase"}
QC_PREFIXES = ("ncount_", "nfeature_", "total_counts_", "pct_counts_", "percent.", "percent_", "n_genes_")


class InputError(ValueError):
    """An input problem that the user must fix (exit code 2)."""


def metadata_like(columns):
    """Column names that look like per-cell QC or clustering metadata rather than genes."""
    return [c for c in columns if c.lower() in QC_NAMES or c.lower().startswith(QC_PREFIXES)]


def parse_count(raw, line, column):
    raw = raw.strip()
    if raw == "0":
        return 0
    try:
        value = float(raw)
    except ValueError:
        value = float("nan")
    if not math.isfinite(value) or value < 0 or value != int(value):
        raise InputError(
            f"line {line}, column '{column}': {raw!r} is not a non-negative integer count. Pseudobulk "
            "needs raw counts (UMIs or reads), not normalised, log-transformed, scaled or integrated "
            f"values. If '{column}' is a metadata column, list it in --exclude-cols or --covariates.")
    return int(value)


def open_table(path):
    try:
        fh = open(path, newline="", encoding="utf-8-sig")
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc}") from exc
    first = fh.readline()
    if not first.strip():
        fh.close()
        raise InputError(f"{path} is empty")
    delimiter = "\t" if "\t" in first else ","
    header = [h.strip() for h in next(csv.reader([first], delimiter=delimiter))]
    dupes = sorted({h for h in header if header.count(h) > 1})
    if dupes:
        fh.close()
        raise InputError(f"duplicate column name(s) in the header: {', '.join(dupes)}")
    return fh, header, csv.reader(fh, delimiter=delimiter)


def _require(header, columns):
    missing = [c for c in columns if c and c not in header]
    if missing:
        raise InputError(f"column(s) not found: {', '.join(missing)}; available: {', '.join(header[:30])}"
                         + (" ..." if len(header) > 30 else ""))


def _records(reader, header):
    for line, rec in enumerate(reader, start=2):
        if not rec or all(not f.strip() for f in rec):
            continue
        if len(rec) != len(header):
            raise InputError(f"line {line}: {len(rec)} fields but the header has {len(header)}")
        yield line, rec


def _check_labels(labels, sample, values, line, required=()):
    for col in required:
        if not values.get(col):
            raise InputError(f"line {line}: empty value in '{col}'")
    previous = labels.setdefault(sample, values)
    if previous != values:
        col = next(c for c in values if values[c] != previous[c])
        raise InputError(
            f"line {line}: sample '{sample}' has more than one value in '{col}' ({previous[col]!r} and "
            f"{values[col]!r}). Labels must be constant within a sample; if a patient's samples are being "
            "merged, they must share these labels.")


def aggregate_wide(header, reader, sample_col, celltype_col, cell_col=None, label_cols=(),
                   exclude_cols=(), required_labels=()):
    """Return (genes, sums, n_cells, labels, n_input_cells) for a one-row-per-cell table."""
    _require(header, [sample_col, celltype_col, cell_col, *label_cols, *exclude_cols])
    idx = {name: i for i, name in enumerate(header)}
    meta = {c for c in (sample_col, celltype_col, cell_col, *label_cols, *exclude_cols) if c}
    genes = [h for h in header if h not in meta]
    suspicious = metadata_like(genes)
    if suspicious:
        raise InputError(
            f"column(s) {', '.join(suspicious)} look like per-cell metadata, not genes; summing them as "
            "genes would distort library sizes. List them in --exclude-cols (or --covariates for "
            "sample-level labels).")
    gene_idx = [idx[g] for g in genes]
    sums, n_cells, labels, seen = {}, {}, {}, set()
    n_input = 0
    for line, rec in _records(reader, header):
        sample, celltype = rec[idx[sample_col]].strip(), rec[idx[celltype_col]].strip()
        if not sample or not celltype:
            raise InputError(f"line {line}: empty value in '{sample_col}' or '{celltype_col}'")
        if cell_col:
            cid = rec[idx[cell_col]].strip()
            if cid in seen:
                raise InputError(f"line {line}: cell id '{cid}' appears more than once")
            seen.add(cid)
        _check_labels(labels, sample, {c: rec[idx[c]].strip() for c in label_cols}, line, required_labels)
        key = (sample, celltype)
        vec = sums.get(key)
        if vec is None:
            vec = sums[key] = [0] * len(genes)
            n_cells[key] = 0
        n_cells[key] += 1
        n_input += 1
        for j, gi in enumerate(gene_idx):
            vec[j] += parse_count(rec[gi], line, genes[j])
    return genes, sums, n_cells, labels, n_input


def aggregate_long(header, reader, sample_col, celltype_col, cell_col, gene_col, count_col,
                   label_cols=(), required_labels=()):
    """Return (genes, sums, n_cells, labels, n_input_cells) for a cell x gene long table."""
    if not cell_col:
        raise InputError("--format long needs --cell-col so that cells can be counted")
    _require(header, [sample_col, celltype_col, cell_col, gene_col, count_col, *label_cols])
    idx = {name: i for i, name in enumerate(header)}
    gene_pos, cell_of, pairs, labels, partial = {}, {}, set(), {}, {}
    for line, rec in _records(reader, header):
        sample, celltype = rec[idx[sample_col]].strip(), rec[idx[celltype_col]].strip()
        cid, gene = rec[idx[cell_col]].strip(), rec[idx[gene_col]].strip()
        if not sample or not celltype or not cid or not gene:
            raise InputError(f"line {line}: empty sample, cell type, cell id or gene value")
        key = (sample, celltype)
        previous = cell_of.setdefault(cid, key)
        if previous != key:
            raise InputError(f"line {line}: cell '{cid}' is in sample/cell type {previous} earlier and {key} here")
        if (cid, gene) in pairs:
            raise InputError(f"line {line}: cell '{cid}' has more than one row for gene '{gene}'")
        pairs.add((cid, gene))
        _check_labels(labels, sample, {c: rec[idx[c]].strip() for c in label_cols}, line, required_labels)
        gene_pos.setdefault(gene, len(gene_pos))
        per_key = partial.setdefault(key, {})
        per_key[gene] = per_key.get(gene, 0) + parse_count(rec[idx[count_col]], line, count_col)
    genes = sorted(gene_pos, key=gene_pos.get)
    n_cells = {}
    for key in cell_of.values():
        n_cells[key] = n_cells.get(key, 0) + 1
    sums = {key: [partial.get(key, {}).get(g, 0) for g in genes] for key in n_cells}
    return genes, sums, n_cells, labels, len(cell_of)


def summarise(genes, sums, n_cells, labels, n_input, *, min_cells, min_units, group_col=None,
              patient_col=None, label_cols=(), sample_col="sample"):
    """Build the metadata rows, the kept pseudobulk ids and the JSON summary."""
    keys = sorted(n_cells, key=lambda k: (k[1], k[0]))
    meta_rows, kept_ids = [], []
    for sample, celltype in keys:
        n = n_cells[(sample, celltype)]
        kept = n >= min_cells
        pid = f"{sample}{ID_SEP}{celltype}"
        row = {"pseudobulk_id": pid, "sample": sample, "cell_type": celltype, "n_cells": n,
               "total_counts": sum(sums[(sample, celltype)]), "kept": "yes" if kept else "no"}
        row.update({c: labels[sample].get(c, "") for c in label_cols})
        meta_rows.append(row)
        if kept:
            kept_ids.append(pid)
    groups = sorted({labels[s][group_col] for s in labels}) if group_col else []
    unit_name = patient_col or sample_col
    by_type = []
    for celltype in sorted({k[1] for k in keys}):
        kept_samples = sorted(s for (s, c) in keys if c == celltype and n_cells[(s, c)] >= min_cells)
        entry = {"cell_type": celltype, "samples_kept": len(kept_samples),
                 "samples_dropped": sum(1 for (s, c) in keys if c == celltype) - len(kept_samples)}
        if group_col:
            units = {g: {} for g in groups}
            for s in kept_samples:
                unit = labels[s][patient_col] if patient_col else s
                units[labels[s][group_col]].setdefault(unit, []).append(s)
            entry["units_per_group"] = {g: len(units[g]) for g in groups}
            repeated = [{"unit": u, "group": g, "samples": ss}
                        for g in groups for u, ss in sorted(units[g].items()) if len(ss) > 1]
            reasons = [f"group '{g}' has {len(units[g])} {unit_name} value(s) after filtering, fewer than {min_units}"
                       for g in groups if len(units[g]) < min_units]
            if len(groups) < 2:
                reasons.append("fewer than two groups in the data")
            if repeated:
                entry["repeated_units"] = repeated
            entry["testable"] = not reasons
            entry["reasons"] = reasons
        by_type.append(entry)
    notes = []
    if any(e.get("repeated_units") for e in by_type):
        notes.append(f"Some {unit_name} values contribute more than one pseudobulk to the same group and cell "
                     "type. They are not independent replicates: keep one per patient by a prespecified rule, "
                     "merge technical replicates of one specimen, or use a mixed model with the patient as a "
                     "random effect.")
    if not genes:
        notes.append("No gene columns found: only cell numbers were computed.")
    summary = {
        "input_cells": n_input,
        "samples": len(labels),
        "cell_types": len({k[1] for k in keys}),
        "genes": len(genes),
        "first_genes": genes[:5],
        "min_cells": min_cells,
        "min_units_per_group": min_units,
        "unit_of_replication": unit_name,
        "pseudobulks": {"total": len(keys), "kept": len(kept_ids), "dropped": len(keys) - len(kept_ids)},
        "dropped": [{"sample": r["sample"], "cell_type": r["cell_type"], "n_cells": r["n_cells"]}
                    for r in meta_rows if r["kept"] == "no"],
        "by_cell_type": by_type,
        "notes": notes,
    }
    return meta_rows, kept_ids, summary


def write_counts(path, genes, sums, kept_ids):
    lookup = {f"{s}{ID_SEP}{c}": vec for (s, c), vec in sums.items()}
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(["gene"] + kept_ids) + "\n")
        for j, gene in enumerate(genes):
            fh.write("\t".join([gene] + [str(lookup[pid][j]) for pid in kept_ids]) + "\n")


def write_meta(path, meta_rows, label_cols):
    cols = ["pseudobulk_id", "sample", "cell_type", "n_cells", "total_counts", "kept", *label_cols]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in meta_rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")


def build_parser():
    p = argparse.ArgumentParser(
        description="Sum raw counts per sample x cell type (pseudobulk), count cells, apply a minimum-cells "
                    "filter and report how many independent units per group remain for each cell type.",
        epilog="Example: pseudobulk_aggregate.py cells.tsv --group-col condition --patient-col patient_id "
               "--out pb_counts.tsv --meta-out pb_meta.tsv",
    )
    p.add_argument("input", help="TSV or CSV; delimiter taken from the header line")
    p.add_argument("--format", choices=["wide", "long"], default="wide",
                   help="wide: one row per cell with gene columns (default); long: one row per cell x gene")
    p.add_argument("--sample-col", default="sample_id", help="sample column (default: sample_id)")
    p.add_argument("--celltype-col", default="cell_type", help="cell-type column (default: cell_type)")
    p.add_argument("--cell-col", help="cell id column (required for --format long; checked for duplicates)")
    p.add_argument("--gene-col", default="gene", help="long format: gene column (default: gene)")
    p.add_argument("--count-col", default="count", help="long format: count column (default: count)")
    p.add_argument("--group-col", help="sample-level group label (condition, response, subtype)")
    p.add_argument("--patient-col", help="patient label; units are counted per patient when given")
    p.add_argument("--covariates", help="comma-separated sample-level columns to carry into the metadata")
    p.add_argument("--exclude-cols", help="wide format: comma-separated per-cell columns that are not genes")
    p.add_argument("--min-cells", type=int, default=10, help="minimum cells per pseudobulk (default: 10)")
    p.add_argument("--min-units", type=int, default=4,
                   help="units per group needed to flag a cell type as testable (default: 4)")
    p.add_argument("--out", help="write the gene x pseudobulk count matrix (TSV) here")
    p.add_argument("--meta-out", help="write the per-pseudobulk metadata (TSV) here")
    return p


def _split(text):
    return [t.strip() for t in text.split(",") if t.strip()] if text else []


def main(argv=None):
    args = build_parser().parse_args(argv)
    label_cols = list(dict.fromkeys(c for c in [args.group_col, args.patient_col, *_split(args.covariates)] if c))
    required = [c for c in (args.group_col, args.patient_col) if c]
    try:
        if args.min_cells < 1 or args.min_units < 1:
            raise InputError("--min-cells and --min-units must be at least 1")
        fh, header, reader = open_table(args.input)
        with fh:
            if args.format == "wide":
                genes, sums, n_cells, labels, n_input = aggregate_wide(
                    header, reader, args.sample_col, args.celltype_col, args.cell_col, label_cols,
                    _split(args.exclude_cols), required)
            else:
                genes, sums, n_cells, labels, n_input = aggregate_long(
                    header, reader, args.sample_col, args.celltype_col, args.cell_col, args.gene_col,
                    args.count_col, label_cols, required)
        if not n_cells:
            raise InputError(f"{args.input} has a header but no data rows")
        meta_rows, kept_ids, summary = summarise(
            genes, sums, n_cells, labels, n_input, min_cells=args.min_cells, min_units=args.min_units,
            group_col=args.group_col, patient_col=args.patient_col, label_cols=label_cols,
            sample_col=args.sample_col)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.out and genes:
        write_counts(args.out, genes, sums, kept_ids)
    elif genes:
        summary["notes"].append("Count matrix not written (no --out given).")
    if args.meta_out:
        write_meta(args.meta_out, meta_rows, label_cols)
    sys.stdout.write(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
