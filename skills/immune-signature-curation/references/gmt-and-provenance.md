# GMT files and provenance records

## Contents

- The GMT format
- How tmelens reads and writes GMT
- The provenance record
- Choosing the label
- Version stamps

## The GMT format

A GMT file is plain text in which each line defines one gene set: its name, a second field (a
description or identifier) and its genes, tab-separated; the gene identifiers must match those of
the data the sets are applied to [Reimand2019]. Lines differ in length and there is no header.

Conventions of this skill that keep every reader happy:

- UTF-8 without a byte-order mark, LF line endings, no quotes, no empty fields.
- Set names unique, also when case is ignored; no spaces at either end.
- No lines starting with `#` and no header line: every line is read as a set.
- One identifier type and one species per file.
- The second field holds provenance: the DOI for a `published` set, `curated`, or
  `derived:<cohort>`. `gmt_lint.py` records the kind it finds (`doi`, `pmid`, `url`, `published`,
  `curated`, `derived`, `user`, `empty`, `gene_like`, `text`).

## How tmelens reads and writes GMT

Version 0.1.0, read in November 2026 [Wang2026]; re-check later versions.

- `read_gmt` splits each line on tabs, skips lines with fewer than three fields, strips spaces
  from genes and drops empty fields, and stores sets in a dictionary by name: a repeated name keeps
  only its last definition, and a `#` line with three or more fields becomes a set.
- `write_gmt` writes name, description and genes; the description is the library DOI when there
  is one, otherwise the source label (`published` or `curated`), and `user` for sets that are not
  in the library.
- `symbol_aliases()` returns the original-to-current symbol map built from `original_symbols`.
- Export the built-in library for linting or overlap checks:

```bash
python3 -c "from tmelens.signatures import load_signatures, write_gmt; write_gmt(load_signatures(), 'tmelens.gmt')"
```

Other readers handle duplicates, `#` lines and spaces in their own ways. Before trusting a reader,
feed it a test file with a duplicated name and a `#` line and see what it loads (heuristic).

## The provenance record

Keep the record as JSON beside the GMT, with the fields of the tmelens library
(`src/tmelens/data/signatures.json`: `name`, `category`, `genes`, `source`, `reference`, `doi`,
`note`, `original_symbols`) [Wang2026] and four fields this skill adds (`location`, `species`,
`identifier_type`, `retrieved`). The example reproduces the tmelens entry for Pan_F_TBRS, whose
original paper used the retired symbols CTGF and FAM101B; the added fields are left as
placeholders.

```json
{
  "version": "2026.11",
  "hgnc_table": "hgnc_complete_set.txt, downloaded <date>",
  "signatures": [
    {
      "name": "Pan_F_TBRS",
      "category": "stroma",
      "genes": ["ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1", "RFLNB",
                "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", "TGFBI",
                "TNS1", "TPM1"],
      "original_symbols": {"CCN2": "CTGF", "RFLNB": "FAM101B"},
      "source": "published",
      "reference": "Mariathasan S et al. Nature 2018;554:544-548.",
      "doi": "10.1038/nature25501",
      "location": "<table or figure of the paper>",
      "species": "human",
      "identifier_type": "HGNC symbol",
      "retrieved": "<date>",
      "note": "Gene list only; scored here as <method>."
    }
  ]
}
```

`original_symbols` maps current to original symbol, the direction tmelens uses; reverse it when
matching data annotated with old symbols.

## Choosing the label

Ask in this order:

1. Was the list learned from data that this project analyses (differential expression,
   clustering, a model)? Then `derived:<cohort>`, even if the genes look familiar.
2. Is it, apart from recorded symbol renames, exactly the list printed in a paper or released in
   a database? Then `published`. tmelens uses `published` for sets that reproduce the gene list
   of the cited paper, including per-cell-type lists taken from a larger marker panel
   ("adapted from" in `note`) [Wang2026].
3. Otherwise `curated`: canonical marker or pathway panels, published lists with genes added or
   removed, and lists whose origin nobody can show. tmelens uses `curated` for panels assembled
   for tmelens that do not reproduce a single publication [Wang2026].

tmelens asks users to verify any list against its original publication before reporting it in a
paper [Wang2026]; do that for every `published` set and record the table or figure.

MSigDB sets: label them `published`, give the collection paper (for the hallmarks
[Liberzon2015]) and the set's own source as listed by MSigDB, and record the release. Mouse
collections: record whether a set is mouse-native or ortholog-converted; the mouse Hallmark
collection (MH) is ortholog-converted [Castanza2023].

## Version stamps

Write into the record and the methods: the HGNC table date, the MSigDB release, the tmelens
version, and the date each list was retrieved, all "as of <month year>". Recheck symbols when any
of these change.
