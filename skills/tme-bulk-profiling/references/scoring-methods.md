# Signature scoring methods in detail

## Contents

- What each tmelens method computes
- Properties that decide validity
- Matching published score definitions
- Translating to and from the GSVA package
- Scoring without tmelens
- Placing a new sample against an existing cohort
- Matrix anchor genes used by the checker

## What each tmelens method computes

tmelens 0.1.0 `score_signatures` [Wang2026]. Let x(g,s) be log2-scale expression of gene g in
sample s and G the measured genes of a set (a set with fewer than 50% of its genes, or fewer than
`min_genes`, is returned as NaN).

- **zscore** (default): z(g,s) = (x(g,s) - mean of x(g,.)) / SD of x(g,.), SD with n - 1;
  genes whose SD is zero or undefined get z = 0; score = mean of z over G.
- **mean**: mean of x(g,s) over G, the log of the geometric mean of the linear values.
- **ssgsea**: genes are ranked within each sample (average ranks for ties, truncated to integers;
  r = 1 for the lowest). Walking down from the highest-ranked gene,
  P_hit(i) = sum of r^alpha over set genes ranked at or above i / sum of r^alpha over the set and
  P_miss(i) = genes outside the set among the first i / (N - |G|); the score is the sum over i of
  P_hit - P_miss, with alpha = 0.25. This is the enrichment score of [Barbie2009], who set alpha
  to 1/4. With `normalize=True` (the default, also used by the command line) every score is
  divided by max - min of the whole samples x signatures score matrix.

`combine_scores` (used by the phenotype axes) z-standardises each requested signature across samples
and averages them, silently skipping signatures that are entirely NaN [Wang2026].

## Properties that decide validity

| Property | zscore | mean | ssgsea (tmelens) |
|---|---|---|---|
| Depends on the other samples in the run | yes: gene means and SDs come from them | no | only through the range normalisation, one scalar per run |
| Cohort mean of each signature | 0 by construction | free | free |
| Changed by forgetting or repeating the log | yes | yes | no: ranks survive monotone transforms |
| Depends on which other genes are in the matrix | no | no | yes: ranks are among all N genes and P_miss uses N - \|G\| |
| Units (counts versus TPM) | depth must be normalised | depth must be normalised | count and CPM ranks also carry transcript length; use TPM |
| One sample on its own | meaningless: every z is 0 | defined | defined; the normalisation then uses that sample's own range |
| Compare two signatures within a sample | no | no [Danaher2017] | no |

The last row holds for every method: a mean score equals log abundance plus an unknown
set-specific constant [Danaher2017], z-scores are centred per signature, and an ssGSEA baseline
depends on how highly the set's genes are expressed. Compare samples within one signature.

Rank-based scores remove monotone scale differences but not gene-specific platform effects, so
array and RNA-seq scores are not on a common scale even with ssGSEA (reasoning). When cohort and
platform coincide, no rescaling can separate biology from platform.

## Matching published score definitions

- **CYT**: Rooney et al. used the geometric mean of GZMA and PRF1 in TPM with a 0.01 offset
  [Rooney2015]. tmelens `mean` on log2(TPM + 1) gives log2 of the geometric mean of (TPM + 1): close
  but not identical at low expression. State which one you computed.
- **Danaher cell types**: the published score is the simple average log2 expression of the markers
  [Danaher2017], which is tmelens `mean` on log2 data. A z-score of the same genes is a different,
  cohort-relative quantity; one log2 unit means a doubling only for the mean score.
- **Ayers profiles**: tmelens `IFNG_Ayers6` and `TIS_Ayers18` reproduce the published gene lists
  [Wang2026]. The published T cell-inflamed profile was developed into a clinical-grade assay
  [Ayers2017]; an unweighted tmelens score of the same genes on another platform is not that assay,
  and the paper's cut-offs and validation do not transfer to it.
- **Curated sets** (`CAF_core`, `Immune_Checkpoints`, `M2_Macrophages` and others) are tmelens panels,
  not published signatures; cite them as tmelens curated panels [Wang2026].

## Translating to and from the GSVA package

- GSVA's own method estimates where each gene's expression sits within the sample population
  [Hanzelmann2013], so it is cohort-relative like the z-score. tmelens has no GSVA method: do not
  call tmelens scores "GSVA scores".
- GSVA's ssGSEA (attributed to [Barbie2009]) ranks expression within each sample and normalises by
  the range across all gene sets and samples [Hanzelmann2013], the behaviour of tmelens
  `normalize=True`.
- GSVA's z-score is the combined z-score of Lee et al., the sum of member z-scores divided by the
  square root of the number of genes [Lee2008] [Hanzelmann2013]. tmelens divides by the number of
  genes. Sample order within a signature is the same; magnitudes differ by the square root of |G|,
  so never mix values from the two tools.
- Package defaults and interfaces change between releases (as of November 2026): read the installed
  version's documentation before reproducing published numbers.

## Scoring without tmelens

When tmelens cannot be imported, compute the same quantities with pandas or plain Python from the
definitions above, on the checked, log2-scale, genes x samples matrix:

```python
import pandas as pd
x = pd.read_csv("expression_log2.tsv", sep="\t", index_col=0)
x.index = x.index.str.strip().str.upper()
wanted = list(dict.fromkeys(g.strip().upper() for g in gene_list))
genes = [g for g in wanted if g in x.index]                      # report coverage = len(genes) / len(wanted)
mean_score = x.loc[genes].mean(axis=0)                           # per sample
sd = x.loc[genes].std(axis=1, ddof=1)
z = x.loc[genes].sub(x.loc[genes].mean(axis=1), axis=0).div(sd.where(sd > 0), axis=0).fillna(0.0)
z_score = z.mean(axis=0)                                         # cohort-relative
```

Implement ssGSEA only if you can test it against a reference implementation; otherwise use mean
scores and say why.

## Placing a new sample against an existing cohort

tmelens has no option for this. If, and only if, the new sample went through the identical
pipeline (same platform, annotation, normalisation and transform), two heuristics of this skill
are defensible for research description:

1. Compute the mean score and report the sample's percentile within the cohort's mean scores.
2. Compute z-scores with the cohort's frozen gene means and SDs instead of re-estimating them.

Neither turns the result into a classification of that patient.

## Matrix anchor genes used by the checker

`check_expression_matrix.py` recognises the gene axis when at least half of the labels match an
identifier pattern or at least three labels are among 48 widely measured genes (ACTB, GAPDH, PTPRC,
CD8A and others), all approved HGNC symbols when checked [HGNC2026]. This is a heuristic for
orientation only; pass `--orientation` when the call is "low" or "undetermined".
