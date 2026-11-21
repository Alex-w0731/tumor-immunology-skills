---
name: immune-signature-curation
description: >-
  Builds, imports and audits gene signatures (gene sets) for immune scoring and enrichment
  analysis: provenance labels (published, curated, derived), HGNC symbol hygiene (previous
  symbols, aliases, gene names turned into dates such as 1-Mar or 2-Sep by spreadsheets), human
  versus mouse nomenclature, signature size and coverage, overlap and redundancy (Jaccard index,
  overlap coefficient), proliferation and tumor-purity contamination, random-signature
  benchmarks, non-circular validation of a new signature, and the GMT format. Bundles a GMT
  linter and an overlap calculator. Use when a researcher imports signatures from papers or
  MSigDB, writes or cleans a GMT file, builds or validates a new immune signature, or asks
  whether signatures are redundant or a survival association is specific, or writes 基因集,
  免疫特征基因集, 基因签名, GMT 文件, 基因名变成日期, 随机基因集 or 新签名验证.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.1"
  domain: tumor-immunology
  summary: Builds, imports and audits immune gene signatures and GMT files for provenance, symbol hygiene, redundancy and confounding.
  summary-zh: 构建、导入并审核免疫基因签名与GMT文件，核查来源、基因符号、冗余与混杂。
---

# Immune signature curation

A researcher wants to use, build or audit gene signatures for immune scoring or enrichment.
Deliver a gene-set collection in which every set carries a provenance label and a reference, uses
current symbols of the right species with the originals recorded, has a known size and coverage,
and has its overlaps with the other sets documented; then say what each score can and cannot
show. For a new signature, deliver a validation plan in which no patient informs both the gene
selection and the evaluation.

## Scope

- Covers: provenance labels; importing from papers, supplements and MSigDB; symbol and species
  hygiene; what size and coverage do to scores; overlap and redundancy; proliferation, purity and
  random-signature checks; validating a new signature; GMT files.
- Not covered: scoring a matrix and per-matrix coverage (use `tme-bulk-profiling`, which bundles
  `signature_coverage.py`); deconvolution (`immune-deconvolution-review`); fitting survival models
  (`immune-survival-analysis`); response classifiers, AUC and cut-points
  (`response-biomarker-evaluation`); whether a paper supports what it is cited for
  (`citation-evidence-check`).
- Research use only: a signature score is never a statement about one patient.

## Workflow

### 1. Record provenance before touching genes

Keep one record per set with the fields of the tmelens library (`name`, `category`, `genes`,
`source`, `reference`, `doi`, `note`, `original_symbols`) [Wang2026], plus `species`,
`identifier_type`, `location` (table or figure, or MSigDB collection and release) and `retrieved`
(date). The template is in [references/gmt-and-provenance.md](references/gmt-and-provenance.md).
Give every set one label:

| Label | Use when | Consequence |
|---|---|---|
| `published` | The genes are the list printed in the cited paper; only symbol renames differ, and they are recorded in `original_symbols` [Wang2026] | Cite paper and table. The label covers the gene list, not a validated score |
| `curated` | A panel assembled for this project, a published list with genes added or removed, or a list nobody can trace [Wang2026] | Say "curated panel" in methods and figures; it has no validation of its own |
| `derived` | Learned from data: differential expression, clusters, model weights | Name the cohort and never evaluate the set on it (step 7). This label is this skill's addition to the tmelens scheme |

A published gene list is not a validated measurement. The T cell-inflamed GEP score is a weighted
sum of housekeeping-normalised values of its 18 genes and became a clinical-grade assay
[Ayers2017]; the unweighted mean tmelens computes from the same genes [Wang2026] is a different
score. Write "genes from [paper], scored by [method]" and do not carry the paper's cut-offs or
performance over to that score.

### 2. Import without spreadsheet damage

- Take each list from the most machine-readable original: the authors' GMT or text file, the
  MSigDB file (record the release), or the paper's table. Never open or save a GMT in a
  spreadsheet program.
- Excel with default settings turns gene symbols into dates and identifiers into floating-point
  numbers [Ziemann2016], and so did Google Sheets in tests; counting five-digit date serials, 30.9%
  of 11,117 articles with supplementary Excel gene lists (2014-2020) had such errors
  [Abeysooriya2021]. When the only source is an Excel supplement, read it programmatically with
  every column as text. A file saved after conversion holds dates, not genes: reading it as text
  does not undo the damage.
- Run `scripts/gmt_lint.py` on every incoming file. Resolve every `error`, decide every
  `warning`, and record each decision for the report.
- A date (`1-Mar`, `2006/09/02`), a float (`2.31E+13`) or a five-digit serial (`43710`) in a gene
  field means the symbol was lost. The linter decodes serials and names the documented candidates
  (SEPT1-14, MARCH1-11, DEC1; `1-Dec` has two, DELEC1 and BHLHE40). A candidate is a lookup, not a
  repair: write the current symbol (`SEPTIN2`, `MARCHF1`) [Bruford2020] once the source confirms
  it; until then keep the field out of the scoring file and list it, with its candidates, as a
  lookup to settle before scoring, because changing the collection later changes ssGSEA scores.

### 3. Bring symbols to current nomenclature for the right species

The procedure, the HGNC download and the alias table are in
[references/symbol-hygiene.md](references/symbol-hygiene.md).

- Human symbols are upper case except the `orf` of C#orf# symbols [Bruford2020]; mouse symbols
  start with an upper-case letter followed by lower case [MGI2026]. Protein names in a list
  (PD-L1, TIM-3, CD45) mean it was typed by hand: check every gene of that list.
- Map previous symbols and aliases with a current HGNC table (`--hgnc`), use the approved symbol
  and keep the original in `original_symbols`, as tmelens does for CTGF to CCN2, FAM101B to RFLNB
  and KIAA0125 to FAM30A [Wang2026].
- Ambiguous symbols need the source: PD1 is an HGNC alias of PDCD1, SNCA and SPATA2; CD16 of
  FCGR3A and FCGR3B; CD25 of IL2RA and ISG20; DEC1 is the previous symbol of DELEC1 and an alias
  of BHLHE40 [HGNC2026]. Never take the first target.
- Annotate the matrix and the GMT with the same HGNC version; renaming one side drops genes
  silently. After renaming, run `signature_overlap.py --universe` with the matrix's genes (or give
  the user that command): a renamed gene in `dropped_by_universe` whose original symbol is in the
  matrix means the matrix carries the older annotation; map both sides with one table.
- Upper-casing a mouse symbol is not orthologue mapping. Map orthologues with a named resource
  and report unmapped genes, or use mouse-native sets, which MSigDB added because orthologue
  mapping is uncertain and orthologues can differ in function [Castanza2023].
- Lint again after renaming: a set that held both CTGF and CCN2 now has a duplicate.

### 4. Know what size and coverage do to a score

- Report for each set the number of genes and the number measured in the data (`--universe` lists
  both, step 3); coverage rules for scoring are in `tme-bulk-profiling`.
- For the mean of k standardised genes with average pairwise correlation r, the variance of the
  score is [1 + (k - 1)r]/k (algebra; equal variances assumed). Small sets carry gene-specific
  noise; large sets track whatever program their genes share, which in bulk tumours is often
  proliferation, purity or total infiltration rather than the named biology.
- One gene: report it as that gene. Two to nine genes: normal for marker panels scored by mean
  or z-score, but one missing or mis-mapped gene moves the score, so name the genes measured.
- Enrichment tools drop sets outside their size filters (Reimand et al. often recommend excluding
  sets of fewer than 10-15 or more than 200-500 genes for human data; GSEA's default upper limit
  is 500 [Reimand2019]). Check which sets were dropped before reporting a set as not enriched.
- In breast cancer more than 90% of random signatures of over 100 genes predicted outcome
  [Venet2011]: every outcome claim needs the step 6 benchmark, and for large sets it decides.
- ssGSEA [Barbie2009], as described in [Hanzelmann2013], normalises by the range of scores
  across all gene sets and samples: adding or removing a set rescales every score. Check whether
  your implementation does this, fix the collection before scoring, and compare scores only
  within one run.

### 5. Measure overlap and redundancy

- Run `scripts/signature_overlap.py` on every set, with `--universe` set to the measured genes.
- `identical`, `nested`, `near_duplicate`: keep one, or keep both and say why; their scores agree
  partly by construction. In the tmelens library CYT (GZMA, PRF1) lies inside
  Danaher_Cytotoxic_cells (J = 0.2, OC = 1), and Immune_Checkpoints and T_cell_Exhaustion share
  five genes (J = 0.42, OC = 0.625) [Wang2026].
- `overlapping`: before interpreting a correlation or a shared association, recompute with the
  shared genes removed from one set (sensitivity analysis; heuristic).
- Low gene overlap does not make two scores independent: sets with only partial overlap can be
  partial readouts of the same process [Liberzon2015], and immune scores co-vary through total
  infiltration (reasoning). Report the score correlation in the data next to the gene overlap.
- MSigDB holds more than 10,000 sets, with redundancy across them [Liberzon2015]: run the overlap
  check on any collection drawn from it, prefer a collection built to reduce redundancy, such as
  the 50 hallmarks, for broad screens [Liberzon2015], and correct for multiple testing across
  every set tested.

### 6. Rule out proliferation, purity and "any gene set works"

Procedures and pseudo-code are in
[references/confounding-and-validation.md](references/confounding-and-validation.md).

- Correlate each score with a proliferation score and with a purity estimate that does not use
  the signature's genes; purity confounds correlations between tumour transcriptomes [Aran2015].
  Purity adjustment of scores is covered in `tme-bulk-profiling`.
- Benchmark any outcome or response association against random gene sets of identical size,
  drawn from the measured genes and run through the same scoring, model, endpoint and covariates;
  report where the signature falls. In breast cancer 28 of 47 published outcome signatures (60%)
  did no better than random sets of identical size [Venet2011].
- Report results with and without adjustment for proliferation and purity. Adjusting for a
  proliferation metagene almost abolished the outcome association of published and random
  signatures in breast cancer, and removing cell-cycle genes from a signature did not remove the
  confounding [Venet2011]: adjust for a proliferation score; do not just delete genes.
- Test the association with patients as the sampling unit (one score per patient in the outcome
  model): p-values from models that take genes as the sampling unit, such as gene permutation,
  assume independent genes and can be wildly anti-conservative [Goeman2007]. The random-set
  benchmark samples genes to answer another question (is this set better than any set of its
  size?), so report its percentile beside the patient-level estimate, not instead of it.

### 7. Validate a new signature without circularity

- No patient may inform both the gene list and its evaluation. With cross-validation, repeat the
  gene selection inside every fold. Error estimated after selecting genes on all samples is
  biased: published few-gene rules reported negligible error that was no longer zero once the
  selection bias was corrected. Use 10-fold cross-validation or the .632+ bootstrap external to
  selection [Ambroise2002], or better an independent cohort with genes, weights and scoring frozen.
- Batch correction that protects the outcome groups can create group differences when groups are
  unbalanced across batches [Nygaard2016]: never pass the outcome to a correction applied before
  evaluation.
- Validate on the endpoint that is claimed. A survival association in TCGA does not validate a
  response signature, and it needs the step 6 checks itself.
- Show incremental value in the validation data: beat random sets of identical size
  [Venet2011] and established gene lists for the same biology, in the same model.
- Check coherence in independent data: the genes should co-vary. Danaher et al. used
  co-expression in 9986 TCGA samples to discard false-positive markers [Danaher2017]; hallmark
  sets had to show coherent expression and were checked in data not used to refine them
  [Liberzon2015].
- For a cell-type signature, compare the score with an orthogonal measurement on the same
  samples, such as flow cytometry or immunohistochemistry [Danaher2017], or single-cell data.
- Until all of this holds, the set stays `derived`, with its cohort named.

### 8. Export

Write UTF-8 text with LF line endings, one set per line, tab-separated: name, description, genes;
unique names; no header or comment lines. Put the DOI, `curated` or `derived:<cohort>` in the
description (tmelens writes the DOI, else the source label [Wang2026]). Ship the provenance
records beside the GMT, run `gmt_lint.py --strict` on the exported file, and version the
collection with the HGNC table date and MSigDB release.

## Defaults and decision rules

| Decision | Default | Source |
|---|---|---|
| Set copied from a paper, only renames | `published`, DOI in the GMT description | [Wang2026] |
| Team panel, modified or untraceable list | `curated` | [Wang2026] |
| Set learned from the data | `derived:<cohort>`; not evaluated on that cohort | heuristic |
| Previous symbol or alias | Approved HGNC symbol; original kept in `original_symbols` | [HGNC2026] [Wang2026] |
| Symbol HGNC links to several genes (PD1, CD16, CD25, DEC1) | Settle from the source paper; never the first target | [HGNC2026] |
| Date, float or five-digit serial in a gene field | Error; out of the scoring file until the source confirms a candidate | [Ziemann2016] [Abeysooriya2021]; hold-out: heuristic |
| Case | Human upper case (C#orf# excepted); mouse upper-case initial, then lower case | [Bruford2020] [MGI2026] |
| Human-style or mouse-style call in `gmt_lint.py` | At least 90% of symbols in that case | heuristic |
| Size flags in `gmt_lint.py` | Below 10 (info), above 500 (warning): permissive ends of the recommended 10-15 and 200-500 bounds | [Reimand2019] |
| Near-duplicate pair | Jaccard index at least 0.5 | heuristic |
| Overlapping pair | Overlap coefficient at least 0.5, borrowed from EnrichmentMap (cut-off 0.5; preset "overlap > 0.5"): a display setting, not a test | [Merico2010] [Reimand2019] |
| Random-set benchmark | Identical size, from measured genes, same model; at least 1,000 draws | [Venet2011]; draws: convention |
| Proliferation and purity | Report adjusted and unadjusted results | [Venet2011] [Aran2015] |
| Gene selection with cross-validation | Inside every fold; 10-fold or .632+ bootstrap | [Ambroise2002] |

## Pitfalls

- Saving a GMT or a gene table from Excel or Google Sheets, which convert symbols on opening,
  pasting and typing [Abeysooriya2021]. Use text tools; lint after every conversion.
- Calling a published list "validated" and quoting the paper's performance for a different score.
  Name the scoring method separately (step 1).
- Upper-casing mouse genes to use human sets: mouse Cd8b1 becomes CD8B1, a previous symbol of
  human CD8B [HGNC2026], and genes without a same-named orthologue are lost. Map orthologues
  (step 3).
- Duplicated set names: dictionary-based readers such as tmelens `read_gmt` keep the last
  definition without warning [Wang2026]. The linter reports `DUPLICATE_SET_NAME` as an error.
- Comment or header lines: tmelens loads a `#` line with three or more tab-separated fields as a
  gene set [Wang2026]. Delete them.
- Reading the joint associations of overlapping or co-varying signatures as independent findings.
  Quantify overlap (step 5) and test the signatures in one model.
- Taking an outcome association as evidence for the named biology: signatures of postprandial
  laughter and of mouse social defeat predicted breast cancer outcome [Venet2011].
- Selecting genes on all samples and cross-validating only the classifier [Ambroise2002].
- Comparing ssGSEA scores from runs with different gene-set collections [Hanzelmann2013].
- Renaming symbols in the GMT but not in the matrix, or the reverse: the renamed genes drop out
  silently. Map both sides with one HGNC table and recheck coverage (step 3).

## Scripts

Both scripts use the Python standard library only; paths are relative to this skill's folder.
Run them; do not read them.

```bash
python3 scripts/gmt_lint.py sets.gmt --out lint.json
python3 scripts/gmt_lint.py sets.gmt --hgnc hgnc_complete_set.txt --species human --strict
```

Input: a GMT file (`.gz` accepted) and optionally an HGNC table with columns `symbol`, `status`,
`prev_symbol`, `alias_symbol`. Without `--hgnc` only a small built-in table of previous symbols
and immune aliases (checked November 2026) is applied, so a missing warning does not show that
every symbol is current. Output (JSON): `summary` (errors, warnings, info, `ok`, counts by code),
`issues` (code, severity, line, set, message, items), `sets` (size, provenance kind, codes),
`identifiers` (identifier type, case style, species assumed), `sizes`. Exit 0; 1 with `--strict`
when there are errors; 2 for unreadable input.

| Codes | Action |
|---|---|
| `SPREADSHEET_DATE`, `SPREADSHEET_NUMBER`, `SPREADSHEET_SERIAL` (errors) | Recover from the original source (step 2) |
| `DUPLICATE_SET_NAME`, `TOO_FEW_FIELDS`, `WRONG_DELIMITER`, `EMPTY_SET`, `NO_SETS` (errors) | Fix the file before any reader sees it |
| `RETIRED_SYMBOL`, `ALIAS_SYMBOL`, `AMBIGUOUS_ALIAS`, `DUPLICATE_AFTER_MAPPING`, `UNKNOWN_SYMBOL` | Step 3 |
| `SPECIES_CASE_MISMATCH`, `MIXED_CASE_STYLES`, `LOWERCASE_SYMBOLS`, `CASE_DUPLICATE` | Settle the species of each set; rerun with `--species` |
| `DESCRIPTION_LOOKS_LIKE_GENE` | The description column is missing and the first gene is lost |
| `NO_PROVENANCE`, `PROVENANCE_WITHOUT_REFERENCE` | Step 1 |
| `SINGLE_GENE_SET`, `SMALL_SET`, `LARGE_SET` | Step 4 |
| `IDENTICAL_SETS`, `DUPLICATE_GENE`, `MIXED_IDENTIFIERS`, other format codes | Remove copies, map to one identifier type, re-export clean text |

```bash
python3 scripts/signature_overlap.py panel.gmt --universe measured_genes.txt --out overlap.json
python3 scripts/signature_overlap.py mine.gmt --against published.gmt --format tsv
```

Input: one GMT, or two with `--against` (each set of the first against each set of the second,
for example a local copy against the published list); `--universe` is one gene per line.
Output: `pairs` (`n_a`, `n_b`, `n_shared`, `jaccard`, `overlap_coefficient`, `flag`,
`shared_genes`; pairs without shared genes only with `--all-pairs`), `sets` (genes before and
after the universe, best match), `clusters` (sets linked by identical or near-duplicate pairs),
`summary`. Flags: `identical`, `nested`, `near_duplicate` (`--jaccard`, default 0.5),
`overlapping` (`--overlap`, default 0.5), `shared`, `empty`. Identifiers are compared exactly
unless `--ignore-case` is given; duplicate set names stop the run with exit 2.

## Report

The answer to the researcher contains:

1. Assumptions: species, identifier type, HGNC table date, MSigDB release, scoring method.
2. A provenance table: set, label, reference and location, genes, genes measured, renamed
   symbols (original to current), repaired or dropped genes with the reason.
3. Checks run: the lint summary and the decision on every warning, overlap flags next to score
   correlations, coverage.
4. Contamination checks: correlation with proliferation and purity, the random-set percentile,
   estimates with 95% confidence intervals with and without adjustment, and the multiple-testing
   correction across all sets tested.
5. For a new signature: how selection and evaluation were separated, the validation cohort and
   endpoint, the comparators.
6. Limits: what each score measures and relative to what, what remains unvalidated, and the
   orthogonal checks still to do (flow cytometry, immunohistochemistry or single-cell data).

```text
Collection: <n> sets (<n> published, <n> curated, <n> derived); <species>; HGNC table of <month year>.
Repairs: <n> spreadsheet-damaged fields recovered from <source>; renamed <old -> new>; dropped <gene: why>.
Coverage: <set> <measured>/<total> genes measured; missing <genes>.
Overlap: <A> vs <B>: <n> shared, J = <..>, OC = <..>, score rho = <..>; they agree partly by construction.
Confounding: rho with proliferation <..>, with purity (<method>) <..>; <statistic> at the <..>th
  percentile of <m> random sets of <k> genes (seed <..>); HR per SD <..> (95% CI <..>) unadjusted,
  <..> (<..>) adjusted for <covariates>; <FDR method> across <n> sets.
Status: <set> is <label>; if derived, from <cohort> and untested in independent patients. Research use only.
```

## References

- [references/gmt-and-provenance.md](references/gmt-and-provenance.md): open when reading or
  writing GMT files, building the provenance record, or labelling MSigDB sets.
- [references/symbol-hygiene.md](references/symbol-hygiene.md): open when symbols need mapping,
  a list came through a spreadsheet, or the species of sets and data differ.
- [references/confounding-and-validation.md](references/confounding-and-validation.md): open
  before testing a signature against outcome or response, or when designing a validation.
- [references/sources.md](references/sources.md): the bibliography; open it to see what each
  source key supports and where that was checked.
- `scripts/gmt_lint.py` and `scripts/signature_overlap.py`: run as shown under Scripts.
