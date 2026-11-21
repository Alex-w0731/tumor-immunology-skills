---
name: response-biomarker-evaluation
description: >-
  Evaluates whether a biomarker (gene signature score, TMB, PD-L1, any single measurement)
  separates responders from non-responders in immunotherapy cohorts. Defines response labels
  from RECIST 1.1 or iRECIST best overall response and tests how stable disease and
  non-evaluable patients change the result; reports ROC AUC with a bootstrap confidence
  interval and permutation p-value, compares two markers on the same patients (DeLong), shows
  small-cohort uncertainty, evaluates a prespecified threshold with exact intervals, audits
  cut-point and feature-selection leakage, separates per-cohort from pooled analysis, and maps
  the report to STARD, REMARK and TRIPOD. Use when asked whether a score predicts response to
  immune checkpoint blockade, for AUC, ROC curve, sensitivity, specificity, PPV, cut-off,
  responder versus non-responder comparisons, or terms such as 疗效预测, 生物标志物, 应答者,
  ROC曲线, 截断值.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.1"
  domain: tumor-immunology
  summary: Tests whether a biomarker separates immunotherapy responders from non-responders, with honest uncertainty and leakage checks.
  summary-zh: 评估生物标志物能否区分免疫治疗应答者与非应答者，含不确定性与信息泄漏检查。
---

# Response biomarker evaluation

A researcher has one value per patient (a signature score, TMB, PD-L1, a cell fraction) and
response outcomes from one or more immune checkpoint blockade cohorts, and asks whether the
marker "predicts response". Deliver an estimate that survives review: an explicit response
label with sensitivity analyses, an AUC with its interval and the counts behind it, a paired
test when markers compete, threshold metrics only for a cut-off fixed in advance, an audit of
everything learned from the same patients, and a report that separates apparent from
validated performance. Research use only: nothing here supports a decision for an individual
patient.

## Scope

- Covers binary response endpoints built from RECIST 1.1 or iRECIST best overall response, or
  a study's clinical-benefit label, for one marker or a fixed score; one or several cohorts.
- Survival endpoints (PFS, OS, Cox models, time-dependent ROC): use `immune-survival-analysis`.
  Computing signature scores from expression: `tme-bulk-profiling`. Checking deconvolution
  estimates used as markers: `immune-deconvolution-review`.
- Building a new multivariable classifier is out of scope beyond the leakage rules below.
- Definitions match the owner's tmelens `biomarker_auc` (label 1 = responder, AUC = U/(n1 n0)
  with ties counted 1/2, stratified percentile bootstrap, AUC unflipped, at least 3 patients
  per class), so results are comparable. Two deliberate differences: 2000 resamples instead of
  1000, and a permutation p-value instead of scipy's Mann-Whitney p. Both test the same null;
  they differ by Monte Carlo error and, where scipy uses its normal approximation, by that.

## Workflow

1. **Pin down the claim before computing.** Write down the marker and its prespecified
   direction (higher or lower in responders), the endpoint, the cohort and treatment, the
   sample timing, and whether the marker's genes, weights, direction or threshold were derived
   using any of these patients. If anything was derived here, every number you produce is
   apparent (in-sample) performance; carry that label into the report.

2. **Define the response label, then plan sensitivity analyses.** Primary analysis:
   responder = CR or PR as best overall response; non-responder = SD or PD; NE (not
   evaluable, no post-baseline scan) excluded but counted. This is how the tmelens
   IMvigor210 example codes `binaryResponse`. Check how the labels were produced: confirmed
   or unconfirmed responses (RECIST 1.1 requires confirmation when response is the primary
   endpoint of a non-randomized trial [Eisenhauer2009]), and RECIST 1.1 or iRECIST (an iUPD
   needs confirmation, iCPD, at the next assessment 4-8 weeks later, and does not override a
   later iSD, iPR or iCR [Seymour2017]). Then repeat the analysis with
   (a) NE counted as non-responders and (b) a clinical-benefit label when durations exist,
   for example durable clinical benefit = PR or SD lasting more than 6 months, with patients
   followed for less than 6 months set aside rather than called non-responders [Rizvi2015].
   Do not headline CR/PR versus PD with SD dropped, even when a collaborator asks for the
   "cleaner" number: diagnostic studies that sample diseased patients and a separate control
   group overestimate accuracy (relative diagnostic odds ratio 3.0) [Lijmer1999], and dropping
   SD builds a similar contrast of extremes, so expect the same bias by analogy. Show it, if at
   all, as a labelled sensitivity analysis. If the conclusion changes between definitions
   (heuristic: the AUC moves by more than the half-width of its own interval, or the interval
   starts or stops including 0.5), say so in the first sentence of the result. Mapping tables
   and dataset quirks: `references/response-labels.md`.

3. **Check the data.** One row per patient: the patient is the unit, and several samples
   per patient shrink the interval falsely, so keep the prespecified pre-treatment sample.
   Use pre-treatment samples for the primary analysis; on-treatment biopsies can reflect the
   response itself (analyse them separately, if at all). Count responders, non-responders
   and exclusions by reason, and save a per-patient inclusion table (patient, label, sample
   timing, included or the reason for exclusion): it lets the researcher audit every count
   and draw the participant flow diagram that STARD 2015 expects (item 19) [Bossuyt2015].
   Check that missing marker values are not concentrated in one class
   (`rows.missing_score_by_class` in the output); note heavy ties (IHC categories, integer
   counts).

4. **Estimate discrimination** with `scripts/roc_auc.py`. The AUC is the probability that a
   random responder scores on the response side of a random non-responder, ties counting
   one half; 0.5 is no discrimination [Hanley1982]. Report the AUC, its 95% interval, n1,
   n0 and the permutation p-value. An AUC below 0.5 in the prespecified direction is
   reported as it is, never flipped to 1 - AUC: choosing the direction after seeing the
   data is itself a form of leakage. An interval that includes 0.5 means discrimination
   was not shown, not that there is no association. If many markers were screened, report
   how many and control the FDR across them [Benjamini1995]; the best of many is optimistic
   on the cohort that picked it.

5. **Size the uncertainty.** Asymptotic intervals for a single AUC under-cover in small
   samples (for high AUCs they needed more than 200 patients), and no small-sample
   alternative is uniformly best [Obuchowski1998]. When the smaller class has fewer than 10
   patients (heuristic; the script warns, also per cohort), quote `se_delong` beside the
   bootstrap interval and call the result exploratory. When an AUC is 0 or 1 (complete
   separation), both collapse to a point: report the counts instead. Before analysing, or
   when someone plans a cohort, run `--plan` to show the precision available; the table
   under Defaults gives typical values.

6. **Compare two markers on the same patients** with `--score2`: paired DeLong test and a
   paired bootstrap interval for the AUC difference [DeLong1988]. Never compare overlapping
   intervals or two separate p-values; the two AUCs are correlated (the output reports the
   correlation), which the paired test uses. Asymptotic intervals for a difference of two
   AUCs kept their coverage even with 20 patients in simulation, unlike single-AUC intervals
   [Obuchowski1998]. If one marker is a model that contains the other and both were fitted
   on these data, do not use DeLong; test the added component's association and then
   estimate the AUC change with an interval [Demler2012]. With several cohorts, see step 9.

7. **Evaluate a threshold only if it was fixed in advance** (a published cut-off, or one
   learned in a separate training cohort) with `scripts/threshold_metrics.py`, which gives
   exact Clopper-Pearson intervals [Clopper1934]. A cut-off optimized on the same patients
   (Youden index, minimum p, maximum chi-square) inflates sensitivity and specificity: true
   values of 84% are estimated near 90% in studies of 40 patients [Leeflang2008]; see also
   [Altman1994]. If only an in-sample cut-off exists, report it as apparent, make the
   threshold-free AUC the primary measure, and ask for validation. PPV and NPV depend on the
   response rate; recompute them with `--prevalence` for the population the claim targets.

8. **Audit leakage.** List every step that used response labels: gene selection, weights,
   direction, threshold, hyperparameters, choice among candidate signatures. Each must be
   learned without the evaluation patients, either in an independent cohort or inside every
   training fold of a cross-validation (nested when anything is tuned) [Ambroise2002]
   [Varma2006]. Unsupervised steps such as within-cohort z-scores do not leak labels but
   make scores cohort-relative. Checklist and nested design: `references/validation-and-reporting.md`.

9. **With several cohorts, analyse per cohort first.** Run with `--cohort` to get per-cohort
   AUCs, the within-cohort AUC (each responder compared only with non-responders of its own
   cohort; this skill's summary, weighting cohorts by their number of pairs) and the pooled
   AUC. Accuracy estimates need covariate adjustment just as associations do [Janes2008]:
   cohorts differ in platform, scale, tumour type and response rate, so a pooled AUC can be
   higher or lower than every cohort's. Treat the pooled AUC as descriptive; when cohorts
   disagree, report the spread instead of one number. Compare markers the same way: with
   `--cohort --score2`, base the conclusion on the per-cohort paired comparisons and the
   within-cohort difference (`comparison.within_cohort`: per-cohort differences weighted by
   pairs, variance from the per-cohort DeLong variances; this skill's extension of
   [DeLong1988], not a published test), never on a pooled comparison that disagrees with them.
   Cohort-relative scores (z-scored within cohort, or ssGSEA divided by the score range as in
   tmelens) never carry a threshold from one cohort to another.

10. **Prognostic or predictive.** In a single-arm treated cohort an association with response
    may be prognostic. A claim that the marker predicts benefit from immunotherapy needs
    treated and untreated patients, preferably randomized, and a treatment-by-marker
    interaction test [Ballman2015]. State this whenever the data are single-arm.

11. **Report** with the template below and the matching checklist (as of November 2026; check
    for newer versions): STARD 2015 when a marker or threshold is assessed as a classifier of
    response [Bossuyt2015]; REMARK for tumour-marker prognostic studies [McShane2005]; TRIPOD
    for multivariable prediction models [Collins2015], updated as TRIPOD+AI [Collins2024].

## Defaults and decision rules

| Question | Default | Source |
|---|---|---|
| Responder | CR or PR as best overall response, confirmed when the protocol requires it | [Eisenhauer2009]; coding as in tmelens |
| Non-responder | SD or PD | convention |
| NE, no post-baseline scan | excluded from the primary analysis and counted; sensitivity analysis as non-responders | convention |
| iRECIST data | use iBOR; iUPD needs confirmation (iCPD) at the next assessment, 4-8 weeks later; a later iSD, iPR or iCR overrides it | [Seymour2017] |
| Clinical-benefit label | the source study's definition, e.g. PR or SD lasting more than 6 months; shorter follow-up set aside | [Rizvi2015] |
| SD dropped (CR/PR vs PD) | sensitivity analysis only, labelled as a contrast of extremes | [Lijmer1999] by analogy |
| AUC | U/(n1 n0), ties 1/2, direction prespecified, never flipped | [Hanley1982] |
| Interval for one AUC | stratified percentile bootstrap; DeLong SE alongside; 2000 resamples (tmelens uses 1000) | [Efron1986]; [DeLong1988]; resample count: convention |
| p-value | label permutation, p = (b + 1)/(m + 1), never 0; 5000 permutations | [Phipson2010]; permutation count: convention |
| Patients per class | at least 3 to compute; fewer than 10 in the smaller class = exploratory | convention (tmelens); heuristic |
| Two markers, same patients | paired DeLong test and paired bootstrap interval of the difference | [DeLong1988] |
| Two markers, several cohorts | per-cohort paired comparisons plus the within-cohort difference; pooled comparison descriptive only | [DeLong1988] per cohort; stratified combination is this skill's extension |
| Nested models fitted on these data | no DeLong; test the added term, then an interval for the AUC change | [Demler2012] |
| Many markers screened | Benjamini-Hochberg FDR; report the number tested | [Benjamini1995] |
| Threshold | prespecified only; exact intervals; no in-sample optimum | [Clopper1934]; [Leeflang2008]; [Altman1994] |
| Gene selection, tuning | outside the evaluation patients or inside every training fold (nested) | [Ambroise2002]; [Varma2006] |
| Several cohorts | per-cohort and within-cohort AUC; pooled AUC descriptive only | [Janes2008]; within-cohort summary is this skill's convention |
| "Predictive" claim | comparator arm and treatment-by-marker interaction | [Ballman2015] |

Approximate 95% intervals an observed AUC can carry (Hanley-McNeil standard error
[Hanley1982], AUC +/- 1.96 SE, from `--plan`; planning only):

| Responders / non-responders | AUC 0.60 | AUC 0.70 | AUC 0.80 |
|---|---|---|---|
| 10 / 30 | 0.39-0.81 | 0.50-0.90 | 0.62-0.98 |
| 15 / 45 | 0.43-0.77 | 0.54-0.86 | 0.66-0.94 |
| 20 / 40 | 0.44-0.76 | 0.55-0.85 | 0.67-0.93 |
| 30 / 90 | 0.48-0.72 | 0.58-0.82 | 0.70-0.90 |
| 50 / 150 | 0.51-0.69 | 0.61-0.79 | 0.72-0.88 |

## Pitfalls

- **Silent label choices.** NE folded into non-responders, or SD dropped, without saying so.
  The AUC changes and cannot be compared with other studies. State the primary definition,
  count every category, report the sensitivity analyses.
- **On-treatment samples mixed with baseline.** A marker induced by response inflates the
  AUC. Keep baseline samples for the primary analysis.
- **Direction chosen after the fact.** Reporting 1 - AUC, or whichever direction looks
  better, is selection on the outcome. Prespecify; report as observed.
- **A point AUC from 25 patients called "good discrimination".** Give n1, n0 and the interval;
  the precision table shows how wide it is.
- **Markers compared by overlapping intervals or separate p-values.** Both miss the
  correlation between AUCs on the same patients. Use the paired comparison.
- **Youden cut-off reported as sensitivity and specificity.** Optimistic, most of all in
  small cohorts [Leeflang2008]. Use a prespecified cut-off or validate elsewhere.
- **Signature genes picked by response in the cohort used to evaluate them**, or
  cross-validation that does not redo the selection inside each fold. The cross-validated
  error can look like zero and is no longer zero once selection is repeated inside the folds
  [Ambroise2002]; tuning by cross-validation on data with no signal gave error estimates well
  below chance [Varma2006]. Use an independent cohort or nested cross-validation.
- **Pooling cohorts with different scales or response rates**, or moving a threshold across
  cohort-relative scores. Analyse per cohort; use the within-cohort AUC as the summary.
- **PPV and NPV quoted from a cohort whose response rate was fixed by sampling** (for example
  equal numbers of responders and non-responders). Recompute at a realistic response rate.
- **"Predictive" claimed from a single-arm cohort** [Ballman2015].

## Scripts

Run the scripts; there is no need to read them. Both use the standard library only, are
deterministic, take TSV (or `.csv`) with one row per patient, print JSON (or write `--out`),
and exit with code 2 and a message naming the bad column or value.

`scripts/roc_auc.py`: AUC, interval, p-value, per-cohort, within-cohort and paired analyses.

```bash
# primary analysis: categories mapped explicitly; unmapped labels (e.g. NE) excluded and counted
python3 scripts/roc_auc.py cohort.tsv --score ifng_score --label bor --positive CR,PR --negative SD,PD --seed 1
# sensitivity analysis: NE as non-responders
python3 scripts/roc_auc.py cohort.tsv --score ifng_score --label bor --positive CR,PR --negative SD,PD,NE --seed 1
# several cohorts, second marker on the same patients, lower values expected in responders
python3 scripts/roc_auc.py all.tsv --score tide --direction lower --label response --cohort study --score2 cd274 --seed 1
# planning: precision of an AUC of 0.70 with 15 responders and 45 non-responders
python3 scripts/roc_auc.py --plan --auc 0.70 --n-pos 15 --n-neg 45
```

The label column holds 1/0, or categories mapped with `--positive` and `--negative` (case is
ignored). Output: `rows` (read, used, `excluded_by_label`, `excluded_missing_score`,
`missing_score_by_class`); `overall` (`n_pos`, `n_neg`, `auc`, `ci_low`, `ci_high`, `p_perm`,
`se_delong`, `observed_direction`); with `--cohort`, `by_cohort` (not estimated below 3 per
class) and `within_cohort`; with `--score2`, `comparison` with `overall` (`auc1`, `auc2`,
`diff`, `diff_ci_low`, `diff_ci_high`, `p_delong`, `corr_auc`) and, with `--cohort`,
`by_cohort` and `within_cohort` (`cohorts_used`: cohorts with at least 3 per class);
`warnings` (small classes, complete separation, missing scores, pooling gaps, excluded
labels), which you should read and pass on. With `--cohort`, `overall` is the pooled AUC.
Keep the same `--seed` to reproduce a run; defaults are 2000 resamples and 5000 permutations.

`scripts/threshold_metrics.py`: one prespecified threshold, no search.

```bash
python3 scripts/threshold_metrics.py cohort.tsv --score tmb --label bor --positive CR,PR --negative SD,PD --threshold 10 --prevalence 0.2
```

Test-positive means score >= threshold (`--direction lower`: score <= threshold). Output:
`tp`, `fp`, `fn`, `tn`, and `sensitivity`, `specificity`, `ppv`, `npv`, `response_rate`,
`test_positive_rate`, each as `x`, `n`, `estimate`, `ci_low`, `ci_high`; `by_cohort` with
`--cohort`; `at_prevalence` with `--prevalence` (point estimates only).

## Report

Give the researcher, in this order:

```text
Question and claim type: <marker>, prespecified direction, endpoint, cohort(s); discovery
  (apparent) or validation.
Labels: responder = ..., non-responder = ...; excluded: NE = n, other = n; criteria
  (RECIST 1.1 / iRECIST; confirmed or not); sample timing.
Data checks: one sample per patient, baseline only, missing values by class, ties; the
  per-patient inclusion table (file name).
Result: AUC = 0.xx (95% CI 0.xx-0.xx; n = a responders, b non-responders; permutation p = ...).
  Sensitivity analyses, each with its own responder and non-responder counts: NE as
  non-responders AUC = ...; clinical benefit AUC = ...
Comparison (if any): difference = 0.xx (95% CI ...), DeLong p = ...; with cohorts: per
  cohort ..., within-cohort ..., pooled (descriptive) ...
Threshold (only if prespecified): sensitivity x/n = 0.xx (exact 95% CI ...), specificity ...,
  PPV and NPV at the cohort's response rate and at <rate>.
Cohorts: per-cohort AUCs, within-cohort AUC, pooled AUC (descriptive).
Leakage audit: what was learned from these patients; consequence.
Script warnings: each one, and what it means for this result.
Limits: single-arm (prognostic vs predictive), sample size, label definition, research use only.
Next step: independent validation cohort, prespecified analysis plan, checklist (STARD / REMARK / TRIPOD+AI).
```

Mark synthetic or simulated data as such, and never present an apparent estimate as
validated performance.

## References

- `references/response-labels.md`: open when mapping best overall response, iRECIST or
  clinical-benefit fields to labels, or when a dataset's response column is unclear.
- `references/validation-and-reporting.md`: open for the leakage audit, the nested
  cross-validation design, pooling decisions and the checklist items to report.
- `references/sources.md`: full citations for every source key, with what each supports.
