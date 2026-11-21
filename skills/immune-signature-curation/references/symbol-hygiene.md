# Symbol hygiene

## Human symbols and the HGNC table

- The guidelines require upper-case Latin letters and Arabic numerals in human symbols and keep
  `orf` lower case in C#orf# symbols such as C3orf18; some symbols carry a hyphen (readthrough
  transcripts such as INS-IGF2), and legacy exceptions remain approved [Bruford2020].
- Download the complete HGNC set (as of November 2026 at
  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt; if it
  has moved, start from https://www.genenames.org). `gmt_lint.py --hgnc` reads its `symbol`,
  `status`, `prev_symbol` and `alias_symbol` columns, and also the column names of the HGNC custom
  download (`Approved symbol`, `Previous symbols`, `Alias symbols`) [HGNC2026].
- With the table, each symbol gets one status: `approved`; `case_variant` (approved once upper
  cased; reported by the case checks); `retired` (a previous symbol of one gene); `alias` (an
  alias of one gene); `ambiguous_previous` or `ambiguous_alias` (several genes, counting previous
  symbols and aliases together, so DEC1, the previous symbol of DELEC1 and an alias of BHLHE40,
  is ambiguous [HGNC2026]); `unknown`.
- Act on them in this order. Approved symbols stay, even if they are also an old name of another
  gene. Retired symbols and single-target aliases are replaced by the approved symbol, with the
  original kept in `original_symbols` [Wang2026]. Ambiguous ones are settled from the source
  paper, never by taking the first target. Unknown ones are typos, withdrawn symbols, non-human
  genes or not gene symbols at all: fix or drop, and say which.
- Apply the same mapping to the expression matrix. If the matrix carries older symbols than the
  GMT (or the reverse), genes fail to match without any error. Two old symbols that map to one
  current symbol create a duplicate in the matrix too.

## The built-in table of gmt_lint.py

Used when no HGNC table is given. Every entry was queried at rest.genenames.org in November 2026
[HGNC2026]. A symbol that is not in the table is "not checked", not "approved".

- Previous symbols: CTGF to CCN2, FAM101B to RFLNB, KIAA0125 to FAM30A (the three renames tmelens
  applies [Wang2026]); SEPT1 to SEPT14 to SEPTIN1 to SEPTIN14, except SEPT13 to SEPTIN7P2; MARCH1
  to MARCH11 to MARCHF1 to MARCHF11; MARC1 and MARC2 to MTARC1 and MTARC2; DEC1 to DELEC1; CD8 to
  CD8A; CD45 to PTPRC; CD39 to ENTPD1; IDO to IDO1; CD11B to ITGAM; CD11C to ITGAX; FCGR3 to
  FCGR3A or FCGR3B (ambiguous).
- Aliases: DEC1 to BHLHE40 and SEPT2 to SEPTIN6, which make DEC1 and SEPT2 ambiguous; TIM-3 and
  TIM3 to HAVCR2; PD-L1 and PDL1 to CD274; PD-1 to PDCD1, but PD1 to PDCD1,
  SNCA or SPATA2; PD-L2 to PDCD1LG2; B7-H3 and B7H3 to CD276; CTLA-4 to CTLA4; VISTA to VSIR;
  CD16 to FCGR3A or FCGR3B; CD16a to FCGR3A; CD32 to FCGR2A or FCGR2B; CD56 and NCAM to NCAM1;
  CD103 to ITGAE; CD25 to IL2RA or ISG20; CD127 to IL7R; CD31 to PECAM1; CD73 to NT5E; CD206 to
  MRC1; GITR to TNFRSF18; OX40 to TNFRSF4; CD137 and 4-1BB to TNFRSF9; NKp46 to NCR1; NKG2D to
  KLRK1; TIM-1 to HAVCR1.
- Lookup is case-insensitive and tries the symbol as written, then without hyphens.

## Spreadsheet damage

What happens:

- Excel with default settings converts gene symbols to dates and floating-point numbers: SEPT2
  became `2-Sep` or `2006/09/02`, MARCH1 `1-Mar`, and the RIKEN accession 2310009E13 became
  `2.31E+13` [Ziemann2016].
- Dates are also stored as five-digit serial numbers, which a later export writes out as numbers;
  counting these raised the share of affected articles to 30.9% [Abeysooriya2021].
- In tests published in 2021, Excel and Google Sheets converted gene names whether a text file
  was opened, data were pasted or names were typed; LibreOffice and Gnumeric did not; formatting
  the destination cells as plain text before pasting or typing avoided the conversion
  [Abeysooriya2021]. Software changes, so lint whatever any program produced.
- HGNC renamed the human symbols that turned into dates (SEPT1 to SEPTIN1, MARCH1 to MARCHF1)
  [Bruford2020]. Lists copied from older papers, older annotation files and non-HGNC identifiers
  still carry the risk.

Recovery, in order:

1. Find the original: the authors' text or GMT file, the table in the paper, MSigDB, or the
   platform annotation. Prefer it over any repair.
2. Without the original, take candidates from the linter only: they cover the documented patterns
   (SEPT1-14 to `#-Sep`, MARCH1-11 to `#-Mar`, DEC1 to `1-Dec`) and list every gene HGNC links
   to the old symbol (`1-Dec` gives DELEC1 or BHLHE40). Five-digit serials are decoded in the
   1900 date system, so `2-Sep` typed in 2019 is 43710 and in 2020 is 44076. Read in the 1904
   date system, the same number is a date four years and one day later (arithmetic), so `2-Sep`
   becomes `3-Sep`: a serial from a workbook of unknown date system has two candidates.
3. A candidate is a lookup, not a recovered symbol, even when it fits the set's biology. Keep the
   field out of any file that will be scored, record the field, its candidates and the source to
   check, and report it as an open lookup. Write the current symbol only once the source confirms
   it, and mark that gene as repaired in the record.
4. If no candidate fits, drop the field and report the loss; never guess a symbol.

Prevention: read spreadsheets programmatically with every column as text (for example
`readxl::read_excel(path, col_types = "text")` in R), keep gene lists in text formats, and lint
after every conversion. Reading as text cannot repair cells that were already saved as dates.

## Species

- Mouse gene symbols start with an upper-case letter followed by lower-case letters or numbers,
  with stated exceptions, and protein symbols are written all upper case in both species
  [MGI2026]. A list of protein names can therefore look human even when the data are mouse.
- Mouse MHC symbols (H2-K1, H2-D1, H2-Aa, H2-Ab1, H2-Eb1) keep an upper-case letter after `H2-`
  and fit neither case pattern; `gmt_lint.py` keeps them out of the species vote, accepts them in
  mouse files and flags them in human files.
- To use human sets on mouse data, or the reverse, map orthologues with a named resource and
  version, and report genes without an orthologue and one-to-many mappings. Upper-casing is not
  orthologue mapping.
- Prefer mouse-native gene sets where they exist: MSigDB added them because orthologue mapping
  is uncertain and orthologous genes can differ in function [Castanza2023].
- Run `gmt_lint.py --species mouse` on mouse sets: the HGNC checks are skipped
  (`symbol_check: skipped_mouse`) and upper-case symbols are flagged as probably human.
