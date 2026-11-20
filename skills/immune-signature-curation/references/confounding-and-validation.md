# Confounding checks and non-circular validation

## Contents

- Random gene-set benchmark
- Proliferation
- Tumor purity
- Circularity checklist
- Coherence and incremental value

## Random gene-set benchmark

Why: signatures unrelated to cancer were associated with breast cancer outcome, 60% of 47
published outcome signatures did no better than random signatures of identical size, and more
than 90% of random signatures of over 100 genes were significant outcome predictors [Venet2011].
A significant association alone therefore says nothing about the genes in the set.

```text
measured  = genes present and expressed in the matrix after QC
k         = number of signature genes in `measured`
observed  = statistic(score(signature))       # same scoring, model, covariates and endpoint
for i in 1..m                                 # m >= 1000 (convention); fixed, reported seed
    R       = k genes drawn from `measured` without replacement
              (optionally within mean-expression bins matched to the signature: heuristic)
    null[i] = statistic(score(R))
b          = count of null[i] at least as extreme as observed (direction fixed beforehand)
p          = (b + 1) / (m + 1)                # convention; b / m understates p by about 1/m [Phipson2010]
percentile = 100 * count(null[i] < observed) / m
```

- Report the percentile, p, m, the seed, the median and 2.5-97.5% range of the random statistics,
  and the same benchmark for any published comparator signature.
- Random genes drawn from the whole matrix include many weakly expressed genes, whose scores are
  noisier than those of a curated set; matching on expression bins avoids an unfair comparison
  (heuristic). Say which draw you used.
- Failing the benchmark shows that the association is not specific to the genes. Passing it does
  not prove the named biology; it removes one alternative explanation.
- The benchmark samples genes, so its p treats the signature as exchangeable with random sets,
  which a coherent set is not. It complements, and does not replace, the patient-level test of
  the association, whose p-value must not come from a model that takes genes as the sampling
  unit [Goeman2007].

## Proliferation

- Score a proliferation set: the curated tmelens panel Proliferation_core (MKI67, TOP2A, BIRC5,
  CCNB1, CDK1, MCM2, PCNA, TYMS, AURKA) [Wang2026], or a meta-PCNA-style metagene built as Venet
  et al. did, from the genes most correlated with PCNA [Venet2011].
- Run `signature_overlap.py` between the signature and the proliferation set first: shared genes
  make any adjustment remove signal by construction.
- Do not treat deleting cell-cycle genes from the signature as the fix: purging cell cycle genes
  failed to rule out the confounding effect of proliferation, because more than half of the breast
  cancer transcriptome correlated with meta-PCNA [Venet2011]. Adjust instead.
- Report the correlation, then the association with and without proliferation as a covariate.
  Adjusting the data for meta-PCNA almost abolished the outcome association of published and
  random signatures in breast cancer, while adjustment with permuted meta-PCNA values left it
  unchanged [Venet2011]; a permuted covariate is a useful control in your own analysis.

## Tumor purity

- Purity confounds correlations between tumour transcriptomes [Aran2015]. Immune and stromal
  scores measure non-tumour content partly by construction (reasoning).
- Use a purity estimate that does not depend on the signature's genes, for example one from DNA
  copy number or pathology. An expression-based estimate built from immune genes is circular for
  an immune score.
- Adjusting an immune score for purity asks about immune composition at a given amount of
  non-tumour tissue, which may remove the effect of interest (reasoning). Report adjusted and
  unadjusted results and say which question each answers. How to adjust scores is covered in the
  sibling skill `tme-bulk-profiling`.

## Circularity checklist

| Design | Verdict | Basis |
|---|---|---|
| Genes selected on all samples, then the classifier cross-validated | Biased error estimate | [Ambroise2002] |
| Gene selection repeated inside each fold (10-fold or .632+ bootstrap) | Acceptable internal estimate | [Ambroise2002] |
| Derived in cohort A, genes, weights and scoring frozen, tested in independent cohort B on the same endpoint | Validation | reasoning |
| Cohort B shares patients with A (a dataset reprocessed under another accession, one trial reported twice) | Not independent; compare sample and patient identifiers | reasoning |
| Batch correction protecting the outcome groups, applied before evaluation | Can create group differences when groups are unbalanced across batches | [Nygaard2016] |
| Cut-off or weights tuned in cohort B | Not validation; see `response-biomarker-evaluation` | reasoning |
| TCGA survival used to validate a response signature | Different endpoint | reasoning |
| No random or established comparator | Incomplete | [Venet2011] |

## Coherence and incremental value

- Coherence in the validation data: compute the median pairwise Spearman correlation among the
  set's genes and the share of variance on their first principal component, and compare both with
  random sets of the same size (heuristic). Danaher et al. evaluated candidate immune markers by
  their co-expression in TCGA [Danaher2017], and each hallmark set had to show coherent expression
  [Liberzon2015]. An incoherent set averages several programs, and its score follows whichever
  dominates a given sample.
- Incremental value: fit the established signatures alone and together with the new one, in the
  same model and patients, and report the change with its uncertainty. Model metrics are covered in
  `immune-survival-analysis` and `response-biomarker-evaluation`.
