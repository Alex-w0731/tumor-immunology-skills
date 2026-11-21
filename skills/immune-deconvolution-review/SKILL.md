---
name: immune-deconvolution-review
description: >-
  Guides choosing, running and auditing methods that estimate immune cell composition from
  bulk tumor transcriptomes: marker-based scores (MCP-counter, xCell) and reference-based
  deconvolution (CIBERSORT, CIBERSORTx, EPIC, quanTIseq, NNLS). Covers what each output can be
  compared across (samples, cell types or both), relative versus absolute fractions, reference
  and data mismatch (platform, tissue, normalization, linear scale), collinear cell types such
  as M1/M2 macrophages, tumor purity, validation on mixtures of known composition or by flow
  cytometry and IHC, and licence limits on matrices such as LM22. Bundles a pure-Python NNLS
  deconvolver with diagnostics and a mixture benchmark. Use when asked to estimate immune cell
  fractions, run or interpret CIBERSORT or xCell output, compare immune infiltration between
  groups, or validate a deconvolution; also for 免疫浸润, 去卷积, 免疫细胞比例, 细胞组成估计.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.1"
  domain: tumor-immunology
  summary: Chooses, runs and sanity-checks immune deconvolution of bulk tumor transcriptomes, with an NNLS script and a mixture benchmark.
  summary-zh: 选择、运行并核查肿瘤 bulk 转录组的免疫去卷积，附 NNLS 脚本与混合样本基准评估。
---

# Immune deconvolution review

A researcher wants immune composition from bulk tumor RNA-seq or microarray data, or brings
output from CIBERSORT, xCell, MCP-counter, EPIC or quanTIseq and wants it interpreted. Deliver
a method matched to the comparison the question needs, inputs checked for scale, platform and
gene coverage, a run or audit with collinearity, fit and purity diagnostics, the validation
evidence that exists, and a report that says what each number can be compared with.

## Scope

- Covers: bulk deconvolution and cell-type enrichment for immune and stromal composition;
  method choice, output meaning, reference and data mismatch, collinearity, tumor purity,
  validation, and licences of reference matrices.
- Not here: gene-set scores and inflamed / excluded / desert phenotypes (use
  `tme-bulk-profiling`); survival models that use the estimates (use `immune-survival-analysis`).
  Single-cell annotation and spatial spot deconvolution are outside this skill.
- Research use only; nothing here is a basis for decisions about an individual patient.

## Workflow

1. **Name the comparison the question needs; it fixes which methods and claims are allowed.**
   - Same cell type between samples or groups ("more CD8 T cells in responders?"): use an
     output that is comparable between samples: MCP-counter or xCell scores, or EPIC,
     quanTIseq or CIBERSORT absolute mode [Sturm2019]. CIBERSORT relative mode alone cannot
     answer it, because its fractions are shares of total immune content [Sturm2019]: a higher
     CD8 share can mean fewer other immune cells, not more CD8 T cells.
   - Different cell types within a sample ("macrophages dominate"): use deconvolution output
     (CIBERSORT relative or absolute, EPIC, quanTIseq). Marker-based scores do not allow it
     [Sturm2019]; MCP-counter was designed to compare one population across samples
     [Becht2016]. xCell is contested; see Defaults.
   - Share of all cells in the tumor: use a method that estimates fractions of all cells and
     models uncharacterized (cancer) cells, EPIC or quanTIseq [Racle2017] [Finotello2019].
     CIBERSORT absolute-mode scores are comparable but are not fractions [Sturm2019].
   If existing output does not support the comparison the researcher made, say so first.

2. **Check the expression input.**
   - Scale: regression methods need linear, non-log values. CIBERSORT works best in non-log
     linear space [Newman2015]; in a benchmark, log or VST data gave two- to four-fold higher
     median RMSE than linear data [AvilaCobos2020]. Ask how the matrix was made. Linearize
     log2(x + 1) with 2^x - 1; `nnls_deconvolve.py` stops on log-looking input until the scale
     is declared.
   - Normalization: match what the reference was built on. EPIC and quanTIseq expect TPM
     [Racle2017] [Finotello2019]; quantile normalization was always sub-optimal for bulk
     deconvolution methods [AvilaCobos2020]. Z-scores, log ratios and log-scale batch-corrected
     values (ComBat output, for example) are not linear mixtures and can be negative; start
     from the uncorrected linear matrix.
   - Identifiers and species: HGNC symbols in both matrices, duplicates collapsed, then marker
     coverage per cell type. A cell type whose markers are absent from the platform or were
     filtered out is estimated from almost nothing. The script reports coverage and flags
     genes that match only when letter case is ignored: human CD8A against mouse Cd8a is not
     an ortholog mapping.
   - Platform: MCP-counter units depend on the platform, so only samples profiled on the same
     platform can be compared [Becht2016]. Treat cohorts pooled from several platforms the
     same way for any method until a batch check shows otherwise (heuristic).

3. **Check the reference against the data.**
   - Tissue and state: reference profiles of earlier methods came mainly from circulating
     immune cells and were generally microarray-based [Racle2017]. CIBERSORT with LM22 (547
     genes, 22 hematopoietic phenotypes) was not explicitly tested on RNA-seq data in the
     original paper [Newman2015].
   - Platform: a single-cell or other-platform reference needs explicit correction; CIBERSORTx
     was built to minimize platform-specific variation for this purpose [Newman2019].
   - Missing cell types: omitting a cell type that is present in the mixture made results
     substantially worse regardless of the other pipeline choices [AvilaCobos2020]. In tumors
     this always includes malignant cells, and often fibroblasts and endothelium. When the
     share of all cells matters, prefer EPIC or quanTIseq; otherwise inspect residuals for the
     missing signal.
   - Tumor-expressed signature genes inflate immune estimates. quanTIseq removed signature
     genes highly expressed in cancer cell lines for this reason [Finotello2019]; with any
     other matrix, look at the genes with the largest positive residuals.
   - Licence: record where the matrix came from and its terms. As of November 2026 the
     CIBERSORTx service is free for non-commercial use only, use by any commercial entity is
     prohibited, and the terms grant a non-transferable right to use the service while
     reserving all other rights [Stanford2026]; re-check the site. Nothing there licenses
     redistribution, so do not commit LM22 or other matrices downloaded from it to a public
     repository (tmelens does not ship LM22). Cite [Newman2015] and point readers to the
     source, or build a matrix from reference profiles you may share (for example with
     tmelens `make_signature_matrix`).

4. **Run.** For a published method use the authors' implementation and record version, mode,
   reference, input scale and normalization. For a custom or shareable matrix, run
   `nnls_deconvolve.py` (see Scripts); it fits the same model as tmelens `deconvolve`.

5. **Sanity-check the output before interpreting it.**
   - Collinearity: a lower condition number means more stable estimates; Abbas et al. and
     CIBERSORT both chose signature genes partly by minimizing it [Abbas2009] [Newman2015],
     and LM22's value of 11.4 [Newman2015] is a reference point, not a cut-off. The number
     summarizes the whole matrix, so one inseparable pair can hide behind an ordinary value;
     read the pair flags too (heuristic). For each flagged pair, typically closely related
     subsets such as M1 and M2 macrophages or naive and memory CD4 T cells, report the sum, or
     merge the columns and refit. In simulated samples every benchmarked method showed
     spillover between CD8 and CD4 T cells, between NK and CD8 T cells, and from DCs to B cells;
     the first two were confirmed in validation data [Sturm2019].
   - Fit: compare per-sample `pearson_r` and `rel_residual` across samples; a sample that fits
     much worse than the rest has content the reference does not describe. Genes with large
     positive mean residuals are expressed by something absent from the reference.
   - Total explained: when signature and mixture share units (both TPM, for example),
     `coef_sum` falls as unexplained content such as tumor cells rises, and fractions hide
     this (heuristic). Compare it between groups before reading differences in fractions.
   - Purity: purity confounds correlation and clustering of tumor transcriptomes [Aran2015].
     Correlate the estimates and `coef_sum` with an independent purity estimate (ESTIMATE
     [Yoshihara2013], or copy-number purity when available). If the groups differ in purity,
     adjust or stratify before attributing immune differences to the grouping. quanTIseq's
     "other" fraction is not purity, since it can include stromal cells [Finotello2019].
   - Background: deconvolution-based methods mostly predicted small amounts of immune cells
     that were absent [Sturm2019]. Read estimates near that background level as "not
     detected"; measure the level on mixtures that lack the cell type (step 6).
   - Composition: relative fractions sum to 1, so a rise in one forces falls elsewhere and
     negative correlations among them are partly arithmetic. Use log ratios [Aitchison1982]
     or absolute-scale estimates for association analyses. NNLS output often contains exact
     zeros, which log ratios cannot take: merge sparse types or state the pseudocount used
     (heuristic).

6. **Validate, or label the estimates unvalidated.**
   - Mixtures of known composition from scRNA-seq of the same tumor type: build the reference
     and the mixtures from different cells. Keeping each patient's cells in one split gave
     slightly higher error than splitting cells [AvilaCobos2020], so split by patient when
     there are enough patients. Include malignant cells at a realistic purity, and mixtures
     that lack each cell type to measure background [Sturm2019].
   - Orthogonal measurements on a subset (flow cytometry or IHC), as used to validate
     MCP-counter, EPIC and quanTIseq [Becht2016] [Racle2017] [Finotello2019]. mRNA fractions
     are not cell fractions; methods that renormalize by mRNA content match cell counts better
     [Racle2017] [Monaco2019].
   - Score with `mixture_benchmark.py` and read per-cell-type r, RMSE, bias and slope, and the
     per-sample r across cell types. Put both tables on the same denominator first.
   - Coverage: a benchmark is evidence only for what it spans. Read from the truth the range of
     the non-reference share (tumor and other cells), each type's lowest true share, and which
     types are ever absent (only those have a measured background); note the platform of the
     reference and the mixtures. Name the target samples that fall outside, such as lower
     immune content, rarer cells, or bulk data scored with a single-cell reference without the
     correction in step 3: there the benchmark gives no evidence (heuristic).
   - With no validation data, say so and restrict claims to the comparison type the method
     supports (step 1). Design details: `references/validation-design.md`.

7. **Report** with the template below; state the comparison type next to every result.

## Defaults and decision rules

What common outputs mean (comparisons from [Sturm2019], Table 1, unless another key is given):

| Method | Family | Output | Between samples | Between cell types | Denominator |
|---|---|---|---|---|---|
| CIBERSORT, relative | regression (nu-SVR) | fractions summing to 1 | no: a change in share is not a change in abundance | yes | immune cells in the reference |
| CIBERSORT, absolute mode | regression | arbitrary units | yes | yes | none, not a fraction |
| EPIC | constrained least squares | fractions with an uncharacterized share | yes | yes | all cells |
| quanTIseq | constrained least squares | fractions with an "other" share | yes | yes | all cells |
| MCP-counter | marker means | arbitrary units, platform-specific [Becht2016] | yes | no | none |
| xCell | marker enrichment | transformed enrichment scores | yes | contested: yes in the authors' simulations, though they say the scores cannot be reliably used as proportions [Aran2017]; no in [Sturm2019] | none |
| TIMER | least squares | arbitrary units | yes, within a cancer type | no | none |
| `nnls_deconvolve.py` | NNLS | fractions plus raw coefficients | coefficients only, if units match (heuristic) | fractions, yes | cell types in the signature |

| Setting | Default | Label |
|---|---|---|
| Input scale for regression | linear, non-log | [Newman2015] [AvilaCobos2020] |
| Log detection | linear if 99th percentile > 100, or max - min > 50 with 25th percentile > 0 | convention (NCBI GEO2R rule, as in tmelens) |
| Minimum shared genes | 2 x number of cell types | convention (tmelens) |
| Row scaling | each gene of signature and mixture divided by the same factor, its mean across signature columns; this re-weights genes, unlike the per-matrix row scaling that [AvilaCobos2020] advise against | convention (tmelens) |
| Collinear pair flag | Pearson r >= 0.9 between row-scaled signature columns | heuristic |
| Condition number | compare alternatives; no published cut-off; LM22 = 11.4 | [Newman2015] |
| Marker coverage warning | fewer than half of a cell type's marker genes present | heuristic |
| Correlation in a benchmark | at least 3 pairs, always beside RMSE and bias | heuristic |
| Denominator mismatch flag | median row sums: one within 0.02 of 1, the other below 0.9 | heuristic |
| High r but wrong scale flag | r >= 0.8 with abs(bias) > 0.5 x mean truth | heuristic |
| Tests across cell types | Benjamini-Hochberg FDR over every cell type tested | convention |
| CIBERSORT P value | a global per-sample test of the whole deconvolution, so it cannot vouch for any single cell type | [Newman2015] |

## Pitfalls

- **Score read as a fraction.** MCP-counter values are arbitrary units [Becht2016]; xCell's
  are enrichment scores that its authors say cannot be reliably used as proportions
  [Aran2017]; marker-based scores cannot rank cell types within a sample [Sturm2019]. Fix:
  compare one cell type across samples, or switch to a deconvolution method.
- **Relative fractions compared as abundance.** A group difference in CIBERSORT relative
  fractions can come from a change in other immune cells or in total infiltration. Fix: pair
  with an absolute-scale method or an immune-content estimate, and word results as shares.
- **Log data into a linear model.** Mixing is additive only on the linear scale; fitting log
  values degrades the estimates [AvilaCobos2020]. Fix: linearize and refit.
- **Unmodelled tumor cells.** Signature genes expressed by cancer cells are credited to
  immune cells, typically the type whose markers they share (the reason quanTIseq filters
  them [Finotello2019]). Fix: check positive residuals, remove or flag those genes, prefer
  methods with an uncharacterized fraction.
- **Reporting collinear subsets separately.** An M1/M2 ratio from near-identical columns is
  dominated by noise and by tumor expression of either type's markers. Fix: report the sum,
  or validate the split on mixtures before using it.
- **Ignoring purity.** Groups that differ in purity can differ in immune estimates for that
  reason alone [Aran2015]. Fix: purity estimate, adjustment or stratification.
- **Validating with correlation only, or on circular mixtures.** A high r coexists with a
  large bias or a wrong denominator, and mixtures built from the reference's own cells
  overstate accuracy. Fix: RMSE, bias and slope on held-out cells or patients.
- **Over-reading small values.** Estimates near the background level are noise [Sturm2019].
- **Filtering on the CIBERSORT P value as quality control.** It tests the whole fit, not each
  cell type [Newman2015]. If samples are removed, report how many and whether they differ in
  purity or group, since the filter can select. (heuristic)
- **Redistributing licensed matrices.** Fix: link to the distributor and document the version.
- **No multiple-testing control.** Testing 10 to 22 cell types per comparison is a family of
  tests. Fix: FDR across all cell types tested, reported beside the raw P values (convention).

## Scripts

Run the scripts; there is no need to read them. Both use the Python standard library only and
exit with code 2 and a message naming the bad file, gene, sample or column.

```bash
# NNLS deconvolution with diagnostics (JSON by default; --format tsv gives one row per sample)
python3 scripts/nnls_deconvolve.py --signature signature.tsv --mixture expression.tsv \
    --mixture-scale log2 --out deconvolution.json --residuals-out residuals.tsv
# Benchmark against known fractions; merge collinear types; put truth on the reference denominator
python3 scripts/mixture_benchmark.py --estimated fractions.tsv --truth truth.tsv \
    --merge "Macrophage=Macrophage_M1+Macrophage_M2" --renormalize truth
```

`scripts/nnls_deconvolve.py` inputs: two TSV files (CSV if the name ends in `.csv`), genes in
the first column; the signature has cell types as columns and linear values, the mixture has
samples as columns. Scale flags take `auto` (default: refuse log-looking data), `linear` or
`log2` (log2(x + 1), linearized). Output fields to read: `fractions` (shares of the signature
cell types), `coefficients` and `coef_sum` (raw scale), `pearson_r`, `rmse`, `rel_residual`
(fit on the row-scaled system), `condition_number_fit` and `condition_number_unscaled` (null
means numerically singular), `collinear_pairs`, `marker_coverage`, `missing_genes`,
`n_case_mismatched_genes`, `top_residual_genes` (positive `mean_residual` means more signal
than the fit explains), and `warnings`. `--format tsv` output can be passed straight to the
benchmark.

`scripts/mixture_benchmark.py` inputs: two samples x cell types tables with the sample ID in
the first column; `NA` is missing. Only shared samples and cell types are compared; others
are listed. Read `per_cell_type` (n, pearson_r, spearman_rho, rmse, bias = mean(est - true),
slope and intercept of est on true, mean_est_when_true_zero), `per_sample` (r across cell
types, sum_true, sum_est), `overall`, and `warnings`, which flag a denominator mismatch, cell
types with high r but large bias, and rows renormalized over fewer cell types because of
missing values. Pooled r over all pairs mixes differences between cell types with accuracy
within each, so lead with the per-cell-type values.

## Report

Include, in this order:

```text
Question and comparison type: <between samples | within sample | share of all cells>
Method: <tool, version, mode, reference matrix and its source and licence>
Input checks: <platform, normalization, scale used; shared genes n/N; cell types with poor marker coverage>
Diagnostics: <condition number; collinear pairs and how they were merged; fit r range and
  outlier samples; genes with large positive residuals; coef_sum by group; background level>
Purity: <estimate used; correlation with immune estimates; adjustment or stratification>
Validation: <none | mixtures (design, held-out cells or patients) | orthogonal (measure, n)>;
  per cell type r, RMSE, bias; what it covers (purity range, lowest shares, platform)
Results: <estimates with spread; group differences with effect size and FDR across cell types>
Limits: <relative vs absolute; unvalidated or merged cell types; values near background>
To verify: <orthogonal measurement (IHC, flow cytometry) of the cell types that drive conclusions>
```

## References

- `references/sources.md`: full citations for every key above; open to cite or check a claim.
- `references/method-notes.md`: inputs, reference origin, modes and caveats per method
  (CIBERSORT/LM22, CIBERSORTx, EPIC, quanTIseq, MCP-counter, xCell, TIMER, tmelens); open when
  the researcher names a tool or you must check an input requirement.
- `references/validation-design.md`: how to build mixture benchmarks and orthogonal
  validation and how to read the metrics; open when designing or reviewing a validation.
