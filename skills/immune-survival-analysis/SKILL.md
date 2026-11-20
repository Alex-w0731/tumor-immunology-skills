---
name: immune-survival-analysis
description: >-
  Tests whether an immune feature (signature score, estimated cell fraction or single gene) is
  associated with overall survival, progression-free interval or another time-to-event outcome
  in a patient cohort: endpoint choice including the TCGA recommendations, a prespecified Cox
  model with the hazard ratio per 1 SD, events-per-variable and proportional-hazards checks,
  adjustment for stage and tumour purity, Benjamini-Hochberg FDR across features, immortal-time
  bias, and tested scripts for Kaplan-Meier, log-rank, FDR and the optimal cut-point problem.
  Use when the user asks for survival analysis, Kaplan-Meier curves, Cox regression, hazard
  ratios, prognostic value of immune infiltration, high versus low groups, a best cut-off or
  surv_cutpoint, or responder versus non-responder survival; also for 生存分析, KM曲线, Cox回归,
  预后, 免疫浸润与生存, 最佳截断值, 风险比.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.0"
  domain: tumor-immunology
  summary: Tests immune features against survival with continuous Cox models, sound endpoints, FDR and no cut-point fishing.
  summary-zh: 用连续Cox模型、合理终点与FDR检验免疫特征与生存的关联，不做截断值挖掘。
---

# Immune feature survival analysis

A researcher wants to know whether an immune feature measured at baseline (a signature score,
an estimated cell fraction, a gene) is associated with a time-to-event outcome in a cohort
such as TCGA, a GEO series or a trial. Deliver a prespecified analysis: a justified endpoint,
the Cox hazard ratio per standard deviation with its confidence interval, unadjusted and
adjusted for the main confounders, checks of events per variable and proportional hazards,
FDR control when several features are tested, a Kaplan-Meier figure from a prespecified split,
and a report that says what was not checked. The result is an association in one cohort, for
research use only, never a prognosis for an individual patient.

## Scope

- Covers: baseline immune features against OS, PFI, PFS, DSS or DFI in retrospective or trial
  cohorts; one feature or a screen of many; reviewing someone else's survival claim (an
  "optimal" cut-off, a responder analysis, an unadjusted hazard ratio).
- Not covered: computing the scores or fractions (use `tme-bulk-profiling` for signature scores
  and `immune-deconvolution-review` for cell-type fractions), developing and validating a
  multivariable prognostic model or risk score, or testing whether a marker classifies
  immunotherapy response (use `response-biomarker-evaluation`).

## Workflow

1. **Write the plan before looking at outcomes.** Fix the feature and how it is computed, the
   endpoint with its event definition and time origin, the population, the covariates, the
   model form (continuous, per 1 SD), every feature and endpoint that will be tested (this is
   the FDR family) and the Kaplan-Meier display split. REMARK asks for prespecified hypotheses,
   precisely defined endpoints, the list of candidate variables and how marker values were
   handled (items 1, 7, 8, 11) [Altman2012]. Each choice made after seeing p-values is a hidden
   extra test; if the user has already explored, say so and treat the result as exploratory.

2. **Build one row per patient and check it.**
   - One sample per patient: choose a rule (primary tumour, first vial or aliquot) and report how
     many samples and patients it removes. Duplicate samples are pseudo-replicates; pass `--id`
     to km_logrank.py and it stops when an ID occurs twice. In TCGA the barcode begins with the
     tissue source site and the patient [Liu2018], so `TCGA-XX-XXXX` (12 characters) is the
     patient and the next two digits are the sample type (01 primary solid tumour, 06
     metastatic, 11 solid tissue normal in the GDC sample-type code table; confirm there, for
     example for blood cancers). If a primary-only rule removes many patients, decide
     explicitly whether metastatic samples can answer the question.
   - Event coded 1 = event, 0 = censored (recode labels such as `1:DECEASED`); one time unit;
     decide what to do with time 0 rows (km_logrank.py lists them).
   - Count patients, events and censored, and the median follow-up by the reverse Kaplan-Meier
     method (censoring treated as the event) [Schemper1996]; km_logrank.py prints all of these.
   - Complete-case analysis drops patients with any missing covariate: report n at each step
     (REMARK item 12 asks for patients and events overall and per subgroup) [Altman2012].

3. **Choose the endpoint.**
   - TCGA: read [references/tcga-endpoints.md](references/tcga-endpoints.md). Table 3 of the
     TCGA Pan-Cancer Clinical Data Resource recommends PFI without reservation for 29 of the 33
     cancer types (its Discussion says 27); OS carries a caution or is not recommended where
     deaths are few (for example PRAD, THCA, PCPG), and DSS is mostly approximated [Liu2018].
   - Other cohorts: use the cohort's own definitions and state them (for immunotherapy cohorts,
     time from treatment start; the progression criteria). OS is the least ambiguous endpoint
     but needs enough deaths.
   - Decide before analysis. If several endpoints are reported, report all of them and include
     all of them in the FDR family.

4. **Check that the events support the model.** Count model parameters, not variables: a
   continuous feature is 1, stage I to IV as categories is 3, each spline term is 1 (the usual
   reading of the rule; convention). Aim for at least 10 events per parameter: in simulations at
   2 to 25 events per variable, Cox coefficients became more biased, confidence intervals lost
   coverage and Wald tests lost validity as the ratio fell, and the authors judged 10 the most
   prudent boundary [Peduzzi1995]. Fewer can be acceptable for a sensitivity analysis that shows
   confounding is controlled [Vittinghoff2007]. If the adjusted model falls short, collapse
   categories (stage I-II versus III-IV), keep the prespecified core covariates and present the
   fuller model as a sensitivity analysis, rather than selecting covariates by p-value.

5. **Fit the continuous Cox model; this is the primary analysis.** Standardise the feature in
   the analysed cohort so the hazard ratio is per 1 SD, and report the SD in original units:
   an HR per SD depends on the cohort's spread and is not comparable across cohorts without it.
   Dichotomising loses power and leaves residual confounding [Royston2006]. Report HR, 95% CI,
   p, n and events, unadjusted and adjusted. Check the functional form with a spline term or a
   martingale-residual plot (convention); if the relation is clearly non-linear, show its shape
   instead of a single HR.
   ```r
   library(survival)
   d$feature_z <- as.numeric(scale(d$feature))
   fit <- coxph(Surv(time, event) ~ feature_z + age + stage_group + purity, data = d)
   summary(fit); cox.zph(fit); plot(cox.zph(fit))   # pan-cancer: add + strata(cancer_type)
   ```
   In Python use `lifelines.CoxPHFitter().fit(df, "time", "event")` on dummy-coded covariates and
   `check_assumptions(df)`. With the tmelens package, `hazard_ratio_table(scores, time, event,
   covariates)` fits one Cox model per feature (HR per 1 SD) and adds a BH `fdr` column over the
   features of that call only, skipping features with fewer than 10 values; when several
   endpoints or models are run, adjust all their p-values together with fdr_adjust.py (step 9).
   If no Cox software is available, the permutation mode of cutpoint_inflation.py reports the
   continuous score test (the Cox score test at HR = 1: p and direction, unadjusted, no HR); say
   that the HR and the adjusted model are missing.

6. **Adjust for confounders, not for consequences.**
   - Include stage and age. REMARK asks for estimates from a model with the standard prognostic
     variables, regardless of their significance (item 17) [Altman2012].
   - Immune cells are a main part of the non-cancer fraction of a bulk sample, so immune scores
     and tumour purity are entangled, and purity confounds transcriptomic correlations
     [Aran2015]: adjust for purity or show the HR with and without it. Prefer a
     purity estimate from DNA (copy number) to one computed from immune and stromal genes in the
     same expression data, which partly adjusts the immune signal away (heuristic).
   - Pan-cancer: stratify the baseline hazard by cancer type and look at per-type estimates;
     pooling lets cancer type masquerade as an immune effect (convention).
   - TCGA: the TCGA-CDR authors compared outcomes from the largest tissue source sites with
     those from the other sites and concluded that the site may need to be taken into account,
     as a proxy for differences in age, stage, treatment and completeness of annotation
     [Liu2018]. When one site dominates, add a sensitivity analysis stratified by site
     (heuristic).
   - Never adjust for anything measured after baseline (response, later treatment).
   - Positive control: in the TCGA-CDR, the HR for stage III-IV versus I-II was significantly
     above 1 for most of the 14 cancer types tested, but not for MESO, PAAD and UVM [Liu2018]. An
     inverted stage HR in your data suggests reversed event coding or a broken join; an
     imprecise one in a small cohort does not (heuristic).

7. **Check proportional hazards** with the scaled Schoenfeld residual test and plot
   [Grambsch1994] (`cox.zph` in R). If the feature violates it, report HRs in prespecified time
   windows or a time-varying coefficient, and say that a single HR is an average over follow-up;
   if a covariate violates it, stratify on that covariate. Read the plot as well as the p-value:
   large cohorts flag trivial departures and small ones miss real ones (heuristic).

8. **Show Kaplan-Meier curves from a prespecified split only.** Run km_logrank.py with
   `--split median` or `--split tertile`; show numbers at risk and do not interpret the tail
   beyond `last_time_with_10_at_risk` (heuristic). The log-rank p answers a cruder question than
   the continuous HR and is secondary. Do not test a cut that was found by searching for the
   largest log-rank statistic or smallest p-value (survminer `surv_cutpoint`, X-tile, a manual
   scan) with an ordinary log-rank test [Altman1994]. If it has been done, run
   cutpoint_inflation.py in permutation mode on the data for the search-adjusted p, run its
   simulation mode with the cohort's n and censored fraction to quote the false-positive rate
   of the search, and treat the cut as a hypothesis for an independent cohort. The maximally
   selected rank statistic (R package maxstat) is an analytic alternative to the permutation
   test [Lausen1992].

9. **Control multiplicity across features.** The family is every feature x endpoint x model
   actually run, not the rows that will be shown; pass `--n-tests` to fdr_adjust.py when the
   file is a subset. Use BH for a screen [Benjamini1995] and Holm when each named feature is a
   claim on its own [Holm1979]. Fix the q threshold in the plan (0.05, or 0.10 for an openly
   exploratory screen; convention). Correlated signatures (CD8 T cell, cytolytic activity,
   IFN-gamma) that all pass are one finding, not several independent confirmations (heuristic).

10. **Watch for groups defined after baseline.** Comparing responders with non-responders, or
    patients with and without immune-related adverse events, from the start of treatment
    favours the group whose definition requires surviving long enough [Anderson1983]. A feature
    measured on treatment (an on-treatment biopsy) is such a classifying event too
    [GiobbieHurder2013]. Use a landmark analysis (fix the landmark before analysis, for example
    a scheduled response assessment, and report how many responses come after it; drop
    patients with an event or censoring before it; classify by status at the landmark; measure
    time from the landmark), a Cox model with a time-dependent covariate, or inverse
    probability weighting [GiobbieHurder2013]. A pre-treatment immune feature analysed from
    treatment start is not affected.

11. **Replicate, then report.** Test the same prespecified model in an independent cohort when
    one exists (same direction, similar HR) and fill in the Report template. Without
    replication, call the result exploratory.

## Defaults and decision rules

| Decision | Default | Basis |
|---|---|---|
| Primary model | Cox, feature continuous and standardised: HR per 1 SD, SD reported in units | [Royston2006]; convention |
| Display split | median (ties go low) or tertiles, fixed in the plan | convention; [Altman1994] |
| Data-derived cut-point | never for inference; if done: search-adjusted p and independent validation | [Altman1994] [Lausen1992] |
| Cut search range in the permutation test | at least 10% of patients in each group | convention |
| Events per model parameter | 10 or more, counting parameters rather than variables | [Peduzzi1995]; counting is convention |
| Below 10 events per parameter | acceptable for a sensitivity analysis of confounding | [Vittinghoff2007] |
| TCGA endpoint | from Table 3; PFI is recommended without reservation for 29 of 33 types | [Liu2018] |
| Proportional hazards | scaled Schoenfeld residual test and plot | [Grambsch1994] |
| Many features | BH over all tests run; q threshold fixed in advance | [Benjamini1995]; threshold is convention |
| Claims about named features | Holm | [Holm1979] |
| Median follow-up | reverse Kaplan-Meier | [Schemper1996] |
| Curve tail | interpret only while 10 or more remain at risk | heuristic |
| Post-baseline groups | landmark or time-dependent Cox, never the naive comparison | [Anderson1983] [GiobbieHurder2013] |
| Reporting | REMARK checklist, in particular items 1, 7, 8, 11, 12 and 15 to 18 | [McShane2005] [Altman2012] |

## Pitfalls

- **Optimal cut-point fishing.** With 150 patients, 40% censoring and cuts searched between the
  10th and 90th percentiles, the best cut of a pure-noise marker gave p < 0.05 in 42% of 1,000
  simulated cohorts, against 5.5% for a prespecified median split and 5.6% for the continuous
  test (`cutpoint_inflation.py --seed 1`; synthetic, this skill's simulation). Fix: continuous
  HR first; search-adjusted p if a searched cut must be discussed.
- **Selective reporting.** Showing the one significant endpoint or feature hides the others.
  Fix: report every test in the plan and adjust over all of them.
- **FDR over the shown rows only.** Adjusting the top 10 of 64 tests is far too lenient. Fix:
  `--n-tests 64`.
- **Responder versus non-responder survival from baseline.** Immortal time makes response look
  protective even when it is not. Fix: landmark or time-dependent analysis (step 10).
- **OS in TCGA types with few deaths.** PRAD has 10 OS events in the whole TCGA-CDR [Liu2018];
  any HR is unstable. Fix: use PFI and recount events in the analysed subset.
- **Reversed event coding.** 0 = dead inverts every HR. Fix: check the event fraction against
  the source and use the stage positive control.
- **Univariable HR presented as independent prognostic value.** Fix: report the adjusted model
  with stage, age and purity next to it.
- **Comparing HR per SD across cohorts.** Different spreads give different HRs for the same
  biology. Fix: report the SD in units, or the HR per fixed unit, next to it.
- **Pooling cancer types.** Immune scores and survival both differ by cancer type. Fix:
  `strata(cancer_type)` and per-type estimates.

## Scripts

Run the scripts; do not read them into context. All take TSV (or CSV if the name ends in
.csv), need only Python 3.10+, exit with code 2 and a message naming the bad column or line, and
print `--help`.

- **scripts/km_logrank.py**: Kaplan-Meier [Kaplan1958] with Greenwood log-scale confidence
  intervals, median survival with its interval, survival at fixed times, the reverse
  Kaplan-Meier median follow-up and the log-rank test for two or more groups [Mantel1966]
  [Peto1972] (same chi-square as R `survdiff`).
  ```bash
  python3 scripts/km_logrank.py --input cohort.tsv --time pfi_months --event pfi_event --group stage_group --id patient_id
  python3 scripts/km_logrank.py --input cohort.tsv --time os_months --event os_event \
      --score cd8_score --split median --times 12,24,60 --curve-out km_curve.tsv
  ```
  Columns: one row per patient, a time column (>= 0), an event column (1/0), and either a group
  column or a continuous `--score` with `--split median|tertile` or a `--cut` fixed in advance;
  `--id` names the patient column to check for duplicates. It offers no cut search, by design.
  Output JSON: `n`, `events`, `follow_up`, per-group `n`, `events`, `median_survival`,
  `median_ci`, `last_time_with_10_at_risk`, `survival_at`; `logrank` (`chi2`, `df`, `p`,
  observed and expected per group, or a `reason` when undefined); `warnings` (groups with fewer
  than 10 events, time 0 rows, dropped rows). `--curve-out` writes the step function for
  plotting. When S(t) is exactly 0.5 over an interval, the median is the start of that interval
  (R prints its midpoint).
- **scripts/fdr_adjust.py**: appends `p_bh` and `p_holm` to a table of p-values (identical to R
  `p.adjust` "BH" and "holm"); missing p-values stay missing.
  ```bash
  python3 scripts/fdr_adjust.py --input cox_per_feature.tsv --column p --sort > cox_per_feature_adj.tsv
  python3 scripts/fdr_adjust.py --input top_hits.tsv --column p --n-tests 64 --out top_hits_adj.tsv
  ```
  The summary line on stderr gives the number of p-values, N and how many pass `--alpha`.
- **scripts/cutpoint_inflation.py**: two modes, both seeded.
  ```bash
  python3 scripts/cutpoint_inflation.py --n 150 --censoring 0.4 --n-sim 1000 --seed 1
  python3 scripts/cutpoint_inflation.py --input cohort.tsv --time os_months --event os_event \
      --score cd8_score --n-perm 1000 --seed 1
  ```
  Simulation mode generates null cohorts (marker independent of survival) and reports
  `false_positive_rate` for `best_cut_point`, `median_split` and `continuous_score_test`, each
  with its Monte Carlo SE. Set `--n` and `--censoring` to match the user's cohort. Permutation
  mode reports `best_cut` (cut, group sizes, `p_naive`), `search_adjusted_p` (permutation p for
  the best cut after the search), the `median_split` test and the continuous score test with its
  direction. A small adjusted p supports some association; it does not validate the cut-point.

The Cox model itself is fitted in R `survival` or Python `lifelines` (step 5); no Cox script is
bundled because the scripts use the standard library only.

## Report

Give the researcher, in this order:

```text
Question and plan   feature (how computed), endpoint (event, time origin), cohort and inclusion
                    rule, model and covariates, number of tests in the FDR family; what was
                    prespecified and what was not.
Data                patients (excluded and why), events, censored, median follow-up (reverse
                    KM) and maximum; events per model parameter.
Primary result      HR per 1 SD (SD = x units) with 95% CI and p, unadjusted and adjusted for
                    stage, age, purity; BH q over N tests.
Checks              proportional hazards (test and plot), functional form, purity sensitivity,
                    stage positive control, cut-point search-adjusted p if a cut was searched.
Display             Kaplan-Meier by prespecified median split: log-rank p, survival at fixed
                    times with numbers at risk.
Limits              retrospective, residual confounding, one cohort, endpoint caveats,
                    association not causation, not for individual prognosis; what was not
                    checked.
Next steps          same model in an independent cohort; measure the same cell population
                    with an orthogonal method (for example IHC or multiplex imaging).
```

## References

- [references/sources.md](references/sources.md): open to check what a source key supports and
  how closely it was read (abstract, full text or metadata) before quoting it.
- [references/tcga-endpoints.md](references/tcga-endpoints.md): open whenever the cohort is
  TCGA, to pick the endpoint for the cancer type and see its event count.
