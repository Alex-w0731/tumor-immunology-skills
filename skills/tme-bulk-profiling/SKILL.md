---
name: tme-bulk-profiling
description: >-
  Profiles the immune context of tumors from a bulk expression matrix (RNA-seq or microarray).
  Checks orientation, identifier type, linear versus log scale, duplicated symbols and missing
  values; chooses between z-score, mean and rank-based ssGSEA signature scoring; checks signature
  gene coverage; handles tumor-purity confounding of immune correlations; and states what
  inflamed, excluded and desert calls from bulk data can and cannot support. Uses the tmelens
  Python package when it is installed and computes the same scores from their definitions
  otherwise. Use when a researcher asks
  for immune signature scores, hot versus cold tumors, immune infiltration of a cohort, CYT or TIS
  scores, or the immune context of a gene, or writes 免疫浸润, 免疫评分, 冷热肿瘤, 免疫表型,
  ssGSEA 打分 or 肿瘤纯度校正.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.0"
  domain: tumor-immunology
  summary: Checks a bulk tumor expression matrix, scores immune signatures and reports immune phenotypes with their limits.
  summary-zh: 检查肿瘤 bulk 表达矩阵，计算免疫特征评分，并在说明局限的前提下给出免疫表型。
---

# Bulk tumor immune profiling

A researcher has a bulk expression matrix and wants immune signature scores, inflamed / excluded /
desert calls, or the immune correlates of a gene. Deliver scores only after the matrix has passed
the checks below, name the frame of reference of every score (relative to these samples, or per
sample within one pipeline), and state coverage, purity handling and limits. Use the tmelens
package when `python3 -c "import tmelens"` succeeds; otherwise compute the same quantities from the
definitions in [references/scoring-methods.md](references/scoring-methods.md) and say so.

## Scope

- Covers: matrix checks, choice of scoring method, signature coverage, tumor-purity confounding,
  phenotype calls from bulk data, and the write-up.
- Not covered: cell-type deconvolution and reference matrices (use `immune-deconvolution-review`);
  Kaplan-Meier, log-rank or Cox models of scores (use `immune-survival-analysis`); single-cell or
  spatial data; testing whether a score separates immunotherapy responders (use
  `response-biomarker-evaluation`).
- Research use only: never turn a score into a statement about one patient's prognosis or treatment.

## Workflow

### 1. Establish what the matrix is

Before computing, find out from the user, the GEO / TCGA record or the file: platform, unit (raw
counts, TPM / FPKM, log2(x + 1), RMA, M-values), species, sample types (primary, metastasis and its
site, normal adjacent) and processing batches. The unit decides the transform, the species decides
whether human signatures apply, and mixed sample types or batches shift every cohort-relative score.
A lymph-node metastasis contains nodal lymphocytes, so its immune scores do not mean the same as a
primary tumor's (reasoning). Score jointly only samples you would compare.

### 2. Check the matrix

Run `scripts/check_expression_matrix.py`. Resolve every `error` before scoring; decide every
`warning` and record the decision for the report.

| Issue code | Action |
|---|---|
| `ORIENTATION_SAMPLES_BY_GENES` | Transpose to genes x samples; tmelens expects genes in rows [Wang2026] |
| `IDENTIFIERS_ENSEMBL`, `IDENTIFIERS_NOT_SYMBOLS` | Map to HGNC symbols first (strip Ensembl version suffixes; for arrays choose and report a probe-collapsing rule). tmelens 0.1.0 does not map identifiers [Wang2026] |
| `SCALE_LOG` | Do not log again: `--log no` on the tmelens command line, `log=False` in Python [Wang2026] |
| `SCALE_LINEAR`, `SUM_TO_ONE_MILLION` | Apply log2(x + 1) before z-score or mean scoring |
| `COUNT_LIKE` | Normalize for depth (CPM, or TPM with gene lengths), then log2(x + 1); tmelens logs but never normalizes library size [Wang2026]. Integer data whose samples each sum to one million are rounded TPM / CPM and get `SUM_TO_ONE_MILLION` instead |
| `SCALE_AUTO_LOG_WOULD_NOT_TRIGGER` | Linear data the auto rule misses: request the transform explicitly |
| `SCALE_CENTERED_PER_GENE`, `SCALE_CENTERED_PER_SAMPLE` | Obtain the uncentered matrix: centered values cannot be logged, and per-gene centering makes within-sample ranks reflect the cohort instead of expression |
| `DUPLICATE_GENE_IDS*` | Choose a collapsing rule and report it; tmelens upper-cases symbols and averages duplicates by default [Wang2026] |
| `MISSING_VALUES` | tmelens drops all-missing genes and fills other gaps with the gene median [Wang2026]; drop genes with many gaps instead when imputation would drive their scores |
| `SYMBOL_CASE_NOT_HUMAN` | Map mouse genes to human orthologs; upper-casing is not an ortholog mapping |
| `EXCEL_DATE_SYMBOLS` | Recover the identifiers from the original file |
| `FEW_SAMPLES` (fewer than 10) | Prefer mean scores and report continuous values; z-scores and quantile-based phenotype labels are unstable at this size (heuristic) |
| `FEW_GENES` (fewer than 5000) | Do not use ssGSEA, which ranks a set only against the genes present; use mean or z-scores (heuristic) |

The scale call is a documented heuristic. Confirm the unit with the data source whenever the call
changes what you do: a column labelled "TPM" is often log2(TPM + 1) already. Once the unit is
known, state the transform instead of relying on `--log auto`: `tmelens score` maps `--log auto`,
`yes` and `no` to `"auto"`, `True` and `False` before calling `prepare_expression` [Wang2026].
Statements about tmelens in this skill describe version 0.1.0 (as of November 2026); check
`tmelens --version` and re-read the changelog if it differs.

### 3. Choose signatures and check their coverage

Export what you will score (`tmelens signatures --export sets.gmt`, or the user's GMT) and run
`scripts/signature_coverage.py` on it with the matrix that will be scored.

- tmelens scores a set as NaN when fewer than 50% of its genes are measured [Wang2026]; the script
  applies the same rule and the same three symbol renamings, so it flags the same sets as long as
  the GMT genes are upper case. The 50% is a software default, not a validated threshold, so also
  read which genes are missing in any set you will interpret.
- Danaher_Treg (FOXP3) and Danaher_Th1 (TBX21) are single genes, and Danaher_CD8_T_cells (CD8A,
  CD8B) and CYT (GZMA, PRF1) still pass the 50% rule with one gene missing: report such scores as
  single-gene readouts.
- Composite axes lose failed members silently: `combine_scores` skips all-NaN signatures, so the
  stromal axis rests on CAF_core alone when Pan_F_TBRS fails coverage [Wang2026]. Report which
  members entered each axis.
- Near-zero coverage for every set means an identifier mismatch, not biology.
- tmelens upper-cases matrix symbols but not GMT genes: upper-case a custom GMT first. It renames
  only three retired symbols (CTGF to CCN2, FAM101B to RFLNB, KIAA0125 to FAM30A) [Wang2026].
- Constant genes count as covered but carry no between-sample information; they only shrink the
  differences between samples in z-score and mean scores (tmelens sets a zero-variance gene's z to
  0) [Wang2026]. The script reports `informative_coverage`.
- tmelens "published" sets reproduce gene lists only, not the published scoring algorithm or its
  cut-offs; "curated" sets are tmelens' own panels [Wang2026].

### 4. Choose the scoring method

First decide the comparison: samples within one jointly processed cohort, a new sample against a
reference, or one cohort against another. Then pick:

| Method | Frame of reference | Valid for | Not valid for |
|---|---|---|---|
| `zscore` (tmelens default) | Each gene standardized across the samples in the run; every signature's cohort mean is 0 | Ranking samples within one jointly processed cohort | A single sample (tmelens returns 0 for every signature); cohorts scored separately (each is centered at 0); scores pooled across runs |
| `mean` | Per-sample level on the log scale; Danaher's cell-score definition [Danaher2017] | Samples from the same pipeline and units, including new ones | Different platforms, units or normalizations |
| `ssgsea` | Within-sample ranks [Barbie2009] [Hanzelmann2013]; tmelens then divides by the range of the run's score matrix [Wang2026] | Per-sample scoring of genome-wide matrices; robust to monotone transforms | Cut-down or panel matrices (ranks among the genes present); raw counts or CPM (their ranks carry transcript length); absolute values compared across runs or platforms |

Rules that hold for every method:

- Compare samples within a signature, never signatures within a sample: a mean score is log
  abundance plus an unknown set-specific constant [Danaher2017], z-scores are centered per
  signature, and an ssGSEA baseline depends on how highly the set's genes are expressed.
- When two cohorts differ in platform or pipeline, cohort and platform are confounded and no
  rescaling, rank-based method or batch correction separates them. Compare within-cohort
  associations side by side, or use an orthogonal assay applied to both.
- A single new sample: use the mean score placed in the reference distribution, or z-scores with
  frozen reference means and SDs, only if the pipeline is identical (heuristics; see
  [references/scoring-methods.md](references/scoring-methods.md)).
- GSVA's own method is cohort-relative and absent from tmelens; do not call tmelens scores GSVA
  scores. GSVA's z-score (Lee's combined z-score) divides by the square root of the gene count,
  tmelens' by the count [Hanzelmann2013] [Lee2008] [Wang2026], so never mix values from the two tools.

### 5. Account for purity before reading any immune correlation

Malignant-cell transcripts scale with purity and immune transcripts with the non-malignant fraction,
so a gene made mainly by tumor cells anticorrelates with immune scores through purity alone
(reasoning; purity confounds correlation and clustering of tumor transcriptomes [Aran2015]).

1. Get a per-sample purity estimate that does not come from expression: ABSOLUTE from copy-number
   and mutation data [Carter2012], or pathology. ESTIMATE purity is built from ssGSEA scores of
   stromal and immune genes [Yoshihara2013], so adjusting an immune score for it partly adjusts the
   score for itself; use it only as a labelled sensitivity analysis.
2. Report raw and purity-adjusted partial Spearman correlations side by side; tmelens
   `gene_immune_context(expr, gene, scores, covariates=clin[["purity"]])` computes the adjusted
   version with BH FDR, and `covariates=None` the raw one [Wang2026].
3. Control the false discovery rate across every signature tested [Benjamini1995].
4. Word results as associations. An association that survives adjustment is still correlational;
   name the spatial, perturbation or orthogonal data that would test a mechanism.

Recipe and interpretation: [references/purity-and-phenotypes.md](references/purity-and-phenotypes.md).

### 6. Phenotype calls, only when asked

The tmelens rule [Wang2026]: infiltration axis = z-standardized mean of Danaher_T_cells,
Danaher_CD8_T_cells, Danaher_Cytotoxic_cells and CYT; stromal axis = the same for Pan_F_TBRS and
CAF_core; inflamed = infiltration at or above its 2/3 quantile; excluded = not inflamed with stroma
at or above its median over all samples; desert = the rest. What follows from it:

- About one third of the samples are inflamed by construction. The proportions are thresholds, not
  findings: never report them as prevalence or compare them between cohorts. For contrast, CD8
  immunohistochemistry in a metastatic urothelial cohort found 26% inflamed, 47% excluded and 27%
  desert [Mariathasan2018].
- The phenotypes are spatial: in excluded tumors CD8 T cells sit in the peritumoral stroma rather
  than the parenchyma [Mariathasan2018], within the immune-profile framework of [Chen2017]. Bulk RNA
  cannot locate cells, so a bulk "excluded" call means low T-cell signal with an active stromal /
  TGF-beta program. Test spatial statements with CD8 IHC or multiplex immunofluorescence.
- Labels are cohort-relative: they change when samples are added or removed and do not transfer
  between cohorts.
- Use the continuous axes, not the labels, for association tests, and never tune the quantiles to
  separate outcomes.

## Defaults and decision rules

| Rule | Value | Source |
|---|---|---|
| Auto-log rule (tmelens `log="auto"`) | Linear if q99 > 100, or max - min > 50 with q25 > 0; then log2(x + 1) | [Wang2026]; convention, documented by tmelens as the NCBI GEO2R rule |
| Further linear calls in the checker | max > 50; or all values non-negative integers with max > 1 | heuristic |
| Centered matrix | Negative values and at least 90% of features with mean within 0.05 SD of zero | heuristic |
| Minimum signature coverage | At least 50% of genes and at least 1 gene, else NaN | [Wang2026]; convention |
| Duplicate symbols | Mean after upper-casing | [Wang2026]; convention |
| Missing values | Drop all-missing genes; median-impute the rest | [Wang2026]; convention |
| ssGSEA | alpha = 0.25; scores divided by the score-matrix range | [Wang2026]; alpha = 1/4 as in [Barbie2009]; range normalization as described in [Hanzelmann2013] |
| Inflamed | Infiltration axis at or above its 2/3 quantile | [Wang2026]; convention |
| Excluded | Not inflamed, stromal axis at or above its median over all samples | [Wang2026]; convention |
| Too few samples for cohort-relative scores or quantiles | Fewer than 10 | heuristic |
| Too few genes for rank-based scores | Fewer than 5000 | heuristic |
| Multiple testing | Benjamini-Hochberg FDR across signatures | [Benjamini1995] |

## Pitfalls

- **Logging twice.** Values labelled TPM are often log2(TPM + 1) already; a second log compresses
  the scale and distorts every score. Trust the checker plus the data source, not the label, then
  state the transform explicitly (`--log no` or `log=False`; step 2).
- **Scoring raw counts.** Library size then drives z-score and mean scores, and count (or CPM)
  ranks carry transcript length into ssGSEA. Normalize first; use TPM for ssGSEA.
- **ssGSEA on a centered or cut-down matrix.** Per-gene centering turns within-sample ranks into
  cohort-relative ranks; a 50-gene export ranks genes only against each other. Score the uncentered
  genome-wide matrix, or use mean or z-scores.
- **Comparing z-scores across cohorts, or adding samples later.** Each run is centered at zero and
  re-scoring moves every value. Score jointly processed samples together.
- **Reading a published gene list as the published score.** A tmelens TIS_Ayers18 value is an
  unweighted summary of the gene list; the published profile was developed into a clinical-grade
  assay [Ayers2017], and its cut-offs do not transfer.
- **Identifier mismatch read as biology.** Ensembl or Entrez rows give all-NaN scores, or scores from
  the few genes that happen to match. Map first and read coverage.
- **Mouse data scored with human signatures.** Upper-cased mouse symbols match only where names
  coincide; genes without a same-named ortholog (the human HLA genes, for instance) drop out
  silently. Map orthologs and report coverage.
- **Purity ignored, or adjusted with ESTIMATE without saying it is circular.** Report the purity
  source and both raw and adjusted estimates.
- **Threshold-defined labels reported as prevalence or as spatial fact.** See step 6.
- **Optimal cut-points or raw p-values across many signatures.** Use continuous scores or
  pre-specified splits and FDR; survival modelling belongs to `immune-survival-analysis`.

## Scripts

Both scripts use only the Python standard library (3.10 or later), are deterministic, and are run,
not read. Exit status 0 means a report was produced; 2 means unreadable input or a bad argument.

```bash
python3 scripts/check_expression_matrix.py expression.tsv --out check.json
python3 scripts/check_expression_matrix.py cohort.csv.gz --orientation samples-by-genes
python3 scripts/signature_coverage.py sets.gmt expression.tsv --out coverage.json
python3 scripts/signature_coverage.py sets.gmt expression.tsv --min-coverage 0.8 --format tsv
```

`scripts/check_expression_matrix.py` reads a TSV or CSV matrix (optionally `.gz`; identifiers in
the first row and first column; the R `write.table` header without a corner cell and leading `#`
lines are accepted). It prints JSON with `orientation` (call, confidence, basis), `identifiers`
(type, symbol case, Ensembl version suffixes, date-like symbols), `scale` (call, quantiles, reasons,
`tmelens_auto_log_would_transform`), `sample_totals`, `values` (missing, negative, non-numeric),
`constant`, `duplicates` and `issues`. Read `issues` first: each has a `code`, a `severity`
(`error`, `warning`, `info`) and the action it calls for. Pass `--orientation` when the
orientation confidence is `low` or the call is `undetermined`.

`scripts/signature_coverage.py` reads a GMT file (name, description, genes; tab-separated) and the
matrix (genes in the first column, or `--orientation samples-by-genes`). Per set it reports
`n_genes`, `n_found`, `coverage`, `n_constant`, `informative_coverage`, `flagged` and the `missing`
genes; `renamed_symbols` lists the tmelens renamings it applied to the matrix; `warnings` cover
identifier mismatch, lower-case GMT genes, duplicated set names and sets that pass only because of
constant genes. With the defaults (`--min-coverage 0.5`, `--min-genes 1`) and an upper-case GMT,
`flagged` marks exactly the sets tmelens 0.1.0 would score as NaN.

## Report

The answer to the researcher contains:

```text
Data      platform; unit as received -> transform applied; species; n samples and types;
          identifier mapping; tool and version (tmelens x.y.z or own code)
Checks    checker issues and the decision for each (scale, duplicates, missing values)
Signatures  name, published gene list or curated panel, coverage, missing genes, sets
          dropped, members that entered each composite axis
Method    zscore / mean / ssgsea and its frame: "relative to these N samples" or
          "per sample within this pipeline"
Results   scores or associations with effect sizes and FDR; raw and purity-adjusted side by
          side; phenotype counts described as threshold-defined
Limits    no spatial information; cohort-relative labels; purity source and circularity;
          platform or batch confounding; research use only
Next      orthogonal checks: CD8 IHC or multiplex IF, expression-independent purity,
          replication in an independent cohort
```

Write in the user's language; keep gene symbols, signature names and statistics untranslated.

## References

- [references/sources.md](references/sources.md): the bibliography with what each source supports;
  open it to cite a claim or check a number.
- [references/scoring-methods.md](references/scoring-methods.md): formulas as tmelens implements
  them, method properties, matching CYT, Danaher and Ayers definitions, GSVA translation, scoring
  without tmelens, placing a new sample. Open it when choosing or implementing a method.
- [references/purity-and-phenotypes.md](references/purity-and-phenotypes.md): purity estimates and
  their circularity, the partial-correlation recipe, what phenotype calls support, wording
  templates. Open it for any gene-immune correlation or phenotype question.
