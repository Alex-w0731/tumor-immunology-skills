# Purity adjustment and immune phenotypes from bulk data

## Contents

- Why purity confounds immune correlations
- Purity estimates and their circularity
- Recipe: purity-adjusted association
- Inflamed, excluded and desert: definitions and the tmelens rule
- What bulk phenotype calls can and cannot support
- Wording templates

## Why purity confounds immune correlations

A bulk sample mixes malignant and non-malignant cells. Transcripts made only by malignant cells
scale with purity; immune and stromal transcripts scale with the non-malignant fraction. So a
gene expressed mainly by tumor cells correlates negatively with every immune and stromal score,
and a gene expressed by immune or stromal cells correlates positively, with no biological link to
T cells (reasoning from mixture arithmetic). Purity has been shown to confound correlation and
clustering analyses of tumor transcriptomes [Aran2015]. Expect raw gene-immune correlations to
shrink, vanish or change sign after adjustment.

## Purity estimates and their circularity

| Estimate | Based on | Use when adjusting immune scores |
|---|---|---|
| ABSOLUTE [Carter2012] | somatic copy-number and mutation data | preferred: independent of expression |
| Pathology review | tumor-cell fraction on a slide | independent but coarse; the section may not match the RNA aliquot (heuristic) |
| CPE [Aran2015] | median of four purity estimates after normalizing them to equal means and SDs | not expression-independent: ESTIMATE is one of its four inputs [Aran2015]; treat it like ESTIMATE when adjusting immune scores |
| ESTIMATE [Yoshihara2013] | ssGSEA scores of 141 stromal and 141 immune genes, combined into the ESTIMATE score | circular: adjusting an immune score for it removes part of the immune signal under study (reasoning from how it is built); use only as a sensitivity analysis and say so |

If only ESTIMATE purity exists, report unadjusted results as the primary analysis, ESTIMATE-adjusted
results as a sensitivity analysis, and state that the adjustment is conservative because it
partly adjusts the immune score for itself.

## Recipe: purity-adjusted association

This is what tmelens `gene_immune_context` computes [Wang2026]:

1. Rank-transform the gene, each signature score and each covariate.
2. Residualize the ranked gene and each ranked score on the ranked covariates by least squares.
3. Correlate the residuals: partial Spearman r, with t = r * sqrt(df / (1 - r^2)) and
   df = n - 2 - k for k covariates.
4. Adjust p-values by Benjamini-Hochberg across all signatures tested [Benjamini1995].

```python
import tmelens as tme
ctx = tme.gene_immune_context(expr, "GENEX", scores, covariates=clin[["purity"]])
```

`expr` is the prepared genes x samples matrix, `scores` the samples x signatures table and `clin`
indexed by the same sample IDs. Report the unadjusted correlation (no covariates) next to the
adjusted one. Interpret: an association that disappears after adjustment was largely purity; one
that persists is an association beyond purity, still not a mechanism. Do not add covariates that
are themselves immune readouts (another immune score, ESTIMATE immune score) unless the question
is explicitly conditional on them.

## Inflamed, excluded and desert: definitions and the tmelens rule

The phenotypes are tissue-level, spatial categories. Most solid tumors show one of three: immune
inflamed, immune excluded or immune desert [Mariathasan2018]; in excluded tumors CD8 T cells are
held in the fibroblast- and collagen-rich peritumoral stroma instead of the parenchyma, and this
pattern was associated with a fibroblast TGF-beta signature [Mariathasan2018]. Chen and Mellman
discuss such immune profiles within the cancer-immune set point framework [Chen2017].

tmelens `classify_immune_phenotype` [Wang2026]:

- infiltration axis: `combine_scores` of Danaher_T_cells, Danaher_CD8_T_cells,
  Danaher_Cytotoxic_cells and CYT;
- stromal axis: `combine_scores` of Pan_F_TBRS and CAF_core;
- inflamed: infiltration at or above its 2/3 quantile; excluded: not inflamed, stroma at or above
  its median over all samples; desert: the rest. Both quantiles are arguments.

Consequences: about one third of samples are inflamed by construction (with n samples and no ties
on the infiltration axis, exactly n - ceil(2(n - 1)/3): 4 of 12, 14 of 40, because tmelens uses
linear-interpolation quantiles); the excluded and desert shares depend on how the two axes
correlate; a signature that fails coverage drops out of its axis silently; and the thresholds
move whenever samples are added or removed.

## What bulk phenotype calls can and cannot support

| Supported, with care | Not supported |
|---|---|
| Ranking tumors within one jointly processed cohort on the continuous axes | Where T cells sit in the tissue (bulk RNA has no spatial information) |
| Association of the continuous axes with a gene or outcome, purity-adjusted, with FDR | The prevalence of inflamed, excluded or desert tumors (fixed by the thresholds) |
| Choosing tumors for spatial validation | Comparing label proportions between cohorts or studies |
| A hypothesis that a stromal / TGF-beta program accompanies low T-cell signal | A claim that TGF-beta or fibroblasts cause exclusion |
| | Any statement about an individual patient's prognosis or treatment |

For scale: by CD8 immunohistochemistry a metastatic urothelial cohort was 47% excluded, 27% desert
and 26% inflamed, and the mean CD8 T-effector signature rose from desert to excluded to inflamed
[Mariathasan2018]. A bulk T-cell signal tracks the histological phenotype on average, which is why
it can rank tumors but cannot place an individual tumor in a spatial category.

## Wording templates

- "Using cohort-relative thresholds (top third of a T-cell infiltration score; stromal score at or
  above the cohort median), k of N tumors were labelled inflamed, a proportion fixed by the
  threshold. Labels are hypotheses about spatial phenotype and were not confirmed histologically."
- "GENEX expression was inversely associated with CYT (Spearman rho = -0.41, FDR = 0.002). After
  adjustment for ABSOLUTE purity the association was attenuated (partial rho = -0.08, FDR = 0.61),
  consistent with confounding by tumor purity." (numbers illustrative)
- Chinese reports: "该分析基于批量转录组的相对评分，只能提示关联，不能说明免疫细胞的空间分布或因果关系。"
