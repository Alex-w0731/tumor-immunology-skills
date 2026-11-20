# Validation, leakage and reporting

Open this file for the leakage audit, when no independent cohort exists and cross-validation
is needed, when deciding whether cohorts may be pooled, and when writing up.

## Leakage audit

For every step, ask: did it see the response labels of the patients whose performance is
being reported? If yes, the reported performance is apparent, not validated.

| Step | Leaks when done on the evaluation patients | Remedy |
|---|---|---|
| Selecting genes or features by responder vs non-responder difference | yes | independent cohort, or repeat the selection inside each training fold [Ambroise2002] |
| Fitting weights or a classifier | yes | same |
| Tuning hyperparameters, choosing the number of genes | yes | nested cross-validation [Varma2006] |
| Choosing the marker's direction | yes | prespecify; report the AUC unflipped |
| Choosing a threshold (Youden, minimum p, maximum chi-square) | yes | prespecify, or learn it in the training folds only [Leeflang2008] |
| Picking the best of many candidate signatures | yes | report all candidates with FDR [Benjamini1995]; validate the winner elsewhere |
| Batch correction told the response labels as a covariate to preserve | yes, by construction | correct without the outcome, or within training folds |
| z-scoring or ranking without labels | no label leakage | fine within a cohort; scores become cohort-relative, so thresholds do not transfer |

## Nested cross-validation when no independent cohort exists

1. Split patients into outer folds stratified by response (and by cohort, if several);
   10-fold is preferred to leave-one-out [Ambroise2002].
2. Inside each outer training set, run every label-using step: feature selection, fitting,
   tuning by an inner cross-validation, threshold choice. Freeze the result [Varma2006].
3. Score the outer test fold with the frozen model.
4. Estimate the AUC from the out-of-fold scores. Either average the per-fold AUCs or compute
   one AUC on pooled out-of-fold scores when the fold models share a scale; say which
   (convention). Repeat the whole procedure over several random splits and report the spread.
5. The bootstrap interval from `scripts/roc_auc.py` on out-of-fold scores ignores the
   variability of model fitting and is too narrow for this use (heuristic). Present the
   result as internal validation; an independent cohort from another centre or platform
   remains the stronger test.

## External validation

- Freeze genes, weights, preprocessing recipe, direction and threshold before opening the
  validation data, and record them with a date.
- Apply them once and report the result whatever it is.
- A threshold on a cohort-relative score (within-cohort z-score, ssGSEA divided by the score
  range, a top-tertile call) transfers only as a rule, and a quantile rule (for example the
  tmelens "inflamed" call at the 2/3 quantile) fixes the share of patients called positive by
  design, so its sensitivity and PPV depend on the cohort's response rate as well as on the
  marker.

## Pooling cohorts

| Situation | What to do |
|---|---|
| One cohort | standard analysis |
| Several cohorts, same assay and scale | per-cohort AUCs, within-cohort AUC as the summary, pooled AUC descriptive only [Janes2008] |
| Different assays or cohort-relative scales | per-cohort and within-cohort AUC only; never pool raw scores or thresholds |
| Per-cohort AUCs disagree | report the spread; look for differences in tumour type, line of therapy, label criteria, assay and timing before any single summary |
| Two markers compared across cohorts | paired comparison per cohort and the within-cohort difference (`roc_auc.py --cohort --score2`, this skill's stratified extension of [DeLong1988]); do not rely on a pooled comparison |

`scripts/roc_auc.py --cohort` warns when the pooled and within-cohort AUCs differ by more than
0.05 (heuristic flag). A formal random-effects meta-analysis of per-cohort AUCs is beyond this
skill; involve a statistician if one is needed.

## Reporting checklists

As of November 2026; check each for a newer version before use.

| Study | Checklist | Items that address this skill's pitfalls |
|---|---|---|
| Marker or threshold assessed as a classifier of response | STARD 2015 [Bossuyt2015] | 12a cut-offs, distinguishing prespecified from exploratory; 15 indeterminate results; 16 missing data; 18 intended sample size; 19 participant flow diagram; 23 cross-tabulation; 24 accuracy estimates with their precision |
| Tumour marker related to outcome (prognostic) | REMARK [McShane2005] | use the full checklist |
| Multivariable model developed or validated | TRIPOD [Collins2015], updated as TRIPOD+AI [Collins2024] | use the full checklist, including how the model was validated |

## Wording apparent versus validated results

- Apparent: "In the discovery cohort, where the signature and its cut-off were derived, the
  apparent AUC was 0.xx (95% CI ...); this estimate is optimistic and requires validation."
- Internal validation: "Nested 10-fold cross-validation, repeating gene selection and tuning
  within each training fold, gave an AUC of 0.xx (range over repeats ...)."
- External validation: "With all parameters frozen before analysis, the AUC in the
  independent cohort was 0.xx (95% CI ...)."
- Single-arm data: "These data show an association with response under treatment; whether the
  marker predicts benefit from the treatment requires a comparator arm [Ballman2015]."
