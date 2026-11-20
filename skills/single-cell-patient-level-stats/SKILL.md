---
name: single-cell-patient-level-stats
description: >-
  Plans, runs and audits group comparisons in single-cell RNA-seq where cells are nested in
  patients or samples, keeping the patient as the unit of replication. Covers pseudoreplication
  in cell-level tests, pseudobulk differential expression per cell type, cell-type proportion
  tests with the compositional caveat, paired pre/post or tumor/normal designs, the p-value
  floor of small designs, and honest reporting with few patients. Bundles a pseudobulk
  aggregator and an exact small-sample composition test. Use when comparing responders and
  non-responders, treatment arms, timepoints or tumor subtypes at single-cell level, when
  FindMarkers or rank_genes_groups returns thousands of condition DE genes, when testing
  cell-type abundance, or with only a few patients. Terms include pseudobulk,
  pseudoreplication, differential abundance, muscat, propeller, scCODA, 伪重复, 伪bulk, 细胞比例,
  患者水平差异分析.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.0"
  domain: tumor-immunology
  summary: "Compares groups in single-cell data at the patient level: pseudobulk DE, composition tests, paired designs and few-patient limits."
  summary-zh: 单细胞组间比较的患者水平统计：伪bulk差异表达、细胞比例检验、配对设计与小样本边界
---

# Patient-level statistics for single-cell group comparisons

A single-cell study that compares conditions (immunotherapy response, treatment arm, timepoint,
tumor subtype) usually has tens of thousands of cells but only a handful of patients, and any
statement about the condition rests on the patients. Deliver an analysis in which the patient (or
an independent sample from a different patient) is the unit of replication, effect sizes reported
against the real number of patients, and wording that claims no more than that number supports.

## Scope

- Covers: condition comparisons of expression within a cell type (pseudobulk), of cell-type
  proportions, paired designs, minimum numbers of patients, reporting.
- Not covered: cluster markers and annotation (cell-level tests are acceptable there as
  descriptive rankings, not as inference about conditions), integration, trajectories,
  cell-cell communication, spatial statistics.
- Use `immune-survival-analysis` for outcome models, `immune-deconvolution-review` for
  composition estimated from bulk RNA-seq, `tme-bulk-profiling` for bulk immune scores,
  `spatial-neighborhood-statistics` for spatial statistics.

## Workflow

### 1. Fix the unit of replication before any test

Build a sample sheet, one row per sample: sample, patient, group, pairing variable, tissue site,
capture or batch, processing protocol. Count independent patients per group, then choose:

| Situation | Design |
|---|---|
| Each patient gives one sample to one group | Unpaired; patient = sample |
| Same patients in both groups (pre/on treatment, tumor/adjacent normal) | Paired, step 5 |
| Some patients give several samples to one group (regions, primary plus metastasis) | Keep one per patient by a rule fixed before seeing results (for example primary tumors only), sum only technical replicates of one specimen, or fit a mixed model with a patient random effect [Hoffman2021] |
| Group fully confounded with capture batch or processing date | Group and batch cannot be separated statistically; say so and stop |

Why: cells of one patient share genotype, history, sampling and processing, so they are not
independent [Zimmerman2021]; testing a treatment effect with non-independent subsamples as
replicates is pseudoreplication [Hurlbert1984]. Adding patients raises power more efficiently
than adding cells per patient [Killip2004].

### 2. Audit existing cell-level results

If conditions were already compared across cells (Seurat `FindMarkers` or scanpy
`rank_genes_groups` with the condition as grouping, Wilcoxon or t-tests over cells, chi-square or
Fisher tests on cell counts pooled across patients), those p-values use the number of cells as the
sample size and ignore patient-to-patient variation. The most widely used single-cell methods can
report hundreds of DE genes when there is no biological difference, and they favour highly
expressed genes [Squair2021]. Type 1 error is highly inflated, including when differences between
individuals are first removed by batch correction (ComBat with individuals as batches) and the
cells are then tested as independent [Zimmerman2021]. Do not carry these p-values into a results
section; redo the analysis at patient level (steps 3 and 4) and explain why the number of findings
shrinks.

### 3. Expression within a cell type: pseudobulk

1. Start from raw integer counts (the RNA assay counts in Seurat, a raw `counts` layer in
   AnnData). Log-normalised, scaled, integrated or SCTransform-corrected values are not raw
   counts, and a count model must not receive them (convention); `scripts/pseudobulk_aggregate.py`
   stops on non-integers, but cannot detect corrected values that are still integers.
2. Sum counts over the cells of each sample and cell type [Lun2017] [Crowell2020]. Summation
   restores type 1 error control when samples (plates) carry their own effects [Lun2017].
3. Drop pseudobulks built from fewer than 10 cells. The muscat case study used a 10-cell
   threshold, and 10 is the default `min_cells` of muscat `pbDS()` (package source, as of November
   2026; re-check) [Crowell2020]. Report how many were dropped and rerun at 30 as a sensitivity
   check (convention). Ten is a floor, not a target: in the muscat simulations, aggregation methods
   gained markedly from 20 to 100 cells per cell type and sample, and only moderately at 200 or
   400 [Crowell2020]. Do not create zero-count pseudobulks for samples that lack the cell type:
   absence of cells is not absence of expression.
4. Per cell type, count patients per group after the filter. Below 4 per group, label the cell
   type exploratory (convention; step 6 says what fewer can support).
5. Fit a count model with library-size normalisation: edgeR quasi-likelihood, DESeq2 or
   limma-voom. Normalisation is required because cell numbers per pseudobulk differ; summed counts
   without it perform poorly when cell numbers are unbalanced [Murphy2022]. Code:
   [references/model-recipes.md](references/model-recipes.md).
6. Design formula. Unpaired: `~ group` plus prespecified covariates that vary between patients;
   never a patient fixed effect, which is aliased with the group. Paired: `~ patient + group`
   [McCarthy2012]. Several samples per patient with a between-patient contrast: a mixed model with
   a patient random effect (dream) [Hoffman2021].
7. Declare the multiple-testing family before looking: BH across all gene and cell-type tests
   (convention), or within each cell type if stated [Benjamini1995].

When asked about the debate: GLMMs with a patient random effect were reported to be more powerful
than pseudobulk [Zimmerman2021]; a reanalysis with a balanced performance measure found pseudobulk
best [Murphy2022]; the original authors replied that type 1 error has to be controlled before
power is compared [Zimmerman2022]. In other benchmarks pseudobulk methods outperformed
single-cell methods [Squair2021] and were among the best performers, alongside cell-level mixed
models, while being far faster [Crowell2020]. Use pseudobulk by default. A cell-level GLMM with
a patient random effect is a defensible alternative; a cell-level test without one is not.

### 4. Cell-type proportions

1. Fix the denominator: all cells, immune cells (PTPRC [CD45] positive) or a parent lineage
   (`--within`). Samples sorted or enriched for a population are comparable with others only
   inside that population.
2. Use one value per sample. A sample without cells of a type has proportion 0, not missing (the
   opposite of the pseudobulk rule: here absence is an abundance measurement).
3. Test with `scripts/composition_test.py`: exact Wilcoxon rank-sum (unpaired) [Mann1947]
   [Wilcoxon1945] or signed-rank (paired) on per-sample values, with rank-biserial correlation
   [Kerby2014], Hodges-Lehmann shift [Hodges1963] and BH across cell types [Benjamini1995].
4. Do not pool cells across patients into a contingency table. Chi-square or Fisher tests on pooled
   counts are the composition version of the pseudoreplication in step 2.
5. Compositional caveat: proportions sum to one, so a real increase in one type produces apparent
   decreases in the others, and testing each type separately can report those as changes;
   Wilcoxon and Poisson approaches do not address this [Buttner2021]. Word decreases as relative.
   Check with `--transform alr --reference TYPE`, choosing an abundant type with no expected
   biological change; a shift that persists is a shift relative to that reference. This rank test
   on a log-ratio is this skill's heuristic, not a published compositional method. For a formal
   compositional model use scCODA, whose results are also relative to a reference type
   [Buttner2021]. With covariates or more than two groups, use propeller (logit or arcsine square
   root transform, empirical Bayes linear models); in its simulations the logit transform
   controlled type 1 error better with 3 or 5 samples per group, arcsine with 10 or 20
   [Phipson2022].
6. Proportions describe the cells that survived dissociation, storage and capture, and protocols
   shift them [Denisenko2020]. Compare only samples processed the same way, call the result
   "proportion of captured cells", and use an orthogonal measure (cell densities from IHC or
   multiplex imaging) for absolute abundance.

### 5. Paired designs

- Compare within patients: signed-rank on per-patient differences for proportions;
  `~ patient + group` for pseudobulk [McCarthy2012].
- Leave incomplete pairs out of the paired test and list them; do not mix paired and unpaired
  samples in one rank test. With many incomplete pairs, a mixed model can use every sample
  [Hoffman2021].
- Report direction counts ("higher on treatment in 5 of 5 patients") beside the effect size.
- Record whether repeat biopsies came from the same lesion or site; a site change is confounded
  with time.

### 6. Match claims to the number of patients

An exact rank test has a floor: without ties, the smallest two-sided p-value is
2 / C(n1 + n2, n1) for two groups and 2 / 2^n for n pairs (arithmetic on the exact null
distributions [Mann1947] [Wilcoxon1945]). Ties, such as several zero proportions, raise it, and BH
across cell types raises the bar further. `composition_test.py` reports it as `p_floor`.

| Design | 2 vs 2 | 3 vs 3 | 3 vs 4 | 4 vs 4 | 5 vs 5 | 4 pairs | 5 pairs | 6 pairs |
|---|---|---|---|---|---|---|---|---|
| Smallest exact p | 0.333 | 0.100 | 0.057 | 0.029 | 0.008 | 0.125 | 0.0625 | 0.031 |

| Patients per group (or pairs) | Defensible | Not defensible |
|---|---|---|
| 1 | Describing these individuals | Any statement or p-value about the condition |
| 2 to 3 | Per-patient values, direction counts, effect sizes, hypotheses | "Significant"; "no difference" (the test could not reject) |
| 4 to 5 | Unpaired: exact p below 0.05 only with complete (4 vs 4) or near-complete (5 vs 5) separation. Paired: 4 or 5 pairs cannot reach a two-sided exact p below 0.05 at all | Claims about modest effects, or of absence |
| 6 or more | Conventional inference with effect sizes and intervals; 6 pairs still reach p below 0.05 only if all six move the same way | Ignoring the floor after multiplicity correction |

The row labels are this skill's heuristic; the p-value statements are arithmetic on the exact null
distributions. For pseudobulk DE the evidence is from bulk RNA-seq: with 3 replicates, 9 of 11
tools found only 20 to 40% of the DE genes detected with 42 replicates, and the authors recommend
at least 6 replicates, rising to at least 12 when DE genes of all fold changes matter
[Schurch2016]. Those were yeast cultures; patients differ far more, so read these numbers as
optimistic (this skill's reasoning, not a tested result). With fewer than 6 patients per group,
add a leave-one-patient-out check (`--exclude-samples`; heuristic) and report whether the
direction holds.

## Defaults and decision rules

| Item | Default | Basis |
|---|---|---|
| Unit of replication | Patient; samples only when from different patients | [Hurlbert1984] [Zimmerman2021] |
| Pseudobulk | Sum of raw counts per sample and cell type | [Lun2017] [Crowell2020] |
| Minimum cells per pseudobulk | 10; sensitivity run at 30 | [Crowell2020] (muscat default, as of Nov 2026); 30 is convention |
| Patients per group before calling a cell type testable | 4 (an execution guard, not a power guarantee) | convention |
| Normalisation | TMM or DESeq2 size factors | need: [Murphy2022]; method: package defaults (convention) |
| Paired DE formula | `~ patient + group` | [McCarthy2012] |
| Several samples per patient, between-patient contrast | Mixed model, patient random effect (dream) | [Hoffman2021] |
| Composition test with few patients | Exact rank-sum or signed-rank on per-sample values | tests: [Mann1947] [Wilcoxon1945]; choice: heuristic |
| Exact or normal approximation | Exact up to 50 samples or non-zero pairs | heuristic (run time) |
| Effect sizes | Rank-biserial r; Hodges-Lehmann shift in proportion units | [Kerby2014] [Hodges1963] |
| Multiplicity | BH within a declared family; BY for arbitrary dependence | [Benjamini1995] [Benjamini2001] |
| Log-ratio check | log2((n + 0.5) / (n_ref + 0.5)), rank test | heuristic; offset as in propeller's logit [Phipson2022] |
| Composition with covariates or several groups | propeller; logit with 3 to 5 samples per group | [Phipson2022] |
| Formal compositional model | scCODA, read against a reference type | [Buttner2021] |
| Replicates for DE | At least 6 per group; 12 when DE genes of all fold changes matter | [Schurch2016] (bulk, yeast) |

The 30-cell sensitivity level, the 4-patient guard and BH across all gene and cell-type tests are
conventions taken from the owner's breast-cancer-multiomics analysis plan, not published thresholds.

## Pitfalls

- Testing conditions across cells. Hundreds of false DE genes, biased to highly expressed ones
  [Squair2021]. Use pseudobulk.
- Batch-correcting the patient away (ComBat, or by the same logic integrated values, with patients
  as batches) and then testing cells. Type 1 error rose markedly with ComBat [Zimmerman2021]: the
  between-patient variance is removed while every cell still counts as a replicate. A patient fixed effect in a
  cell-level model is no fix either: for a between-patient contrast it is aliased with the group.
  Use pseudobulk or a GLMM with a patient random effect.
- Fisher or chi-square tests on pooled cell counts. The p-value tracks cell numbers. Compare
  per-sample proportions.
- Counting regions, biopsies or captures of one patient as replicates. The real n is smaller than
  reported. Keep one per patient or model the patient.
- Zero-filling missing pseudobulks, or treating absent cell types as missing in composition.
  Both bias results; the rules are opposite for the two analyses.
- Reading a falling proportion as depletion. It may be compositional [Buttner2021]. Word it as
  relative and run the log-ratio check.
- Reading p above 0.05 at 3 vs 3 as "no difference". The floor is 0.1. Report the floor, effect
  sizes and direction counts.
- Changing the cell threshold, denominator, test or cell-type merging after seeing results. The
  p-values lose their meaning. Fix choices first and report every sensitivity run.
- Comparing proportions between samples dissociated, stored or sorted differently. Differences may
  be protocol artefacts [Denisenko2020]. Restrict or stratify by protocol.

## Scripts

Run the scripts; do not read them into context. Both exit with code 2 and a message naming the
column, line or value when the input is wrong. Python 3.9 or later, standard library only.

`scripts/pseudobulk_aggregate.py` sums raw counts per sample and cell type, applies
`--min-cells`, and counts units per group.

```bash
python3 scripts/pseudobulk_aggregate.py cells.tsv --sample-col sample_id --celltype-col cell_type \
  --cell-col cell_id --group-col condition --patient-col patient_id --covariates batch \
  --exclude-cols percent_mt,n_genes --min-cells 10 --out pb_counts.tsv --meta-out pb_meta.tsv > pb_summary.json
```

- Input, wide format: one row per cell; sample, cell type, optional cell id, labels, then one
  integer column per gene. Name every other per-cell column in `--exclude-cols`: an integer QC
  column left in (nCount_RNA, n_genes) would be summed as a gene and distort library sizes. The
  script stops on common QC column names and on non-integers, but cannot recognise every metadata
  column, so check `genes` and `first_genes` in the JSON. Long format (`--format long`): one row
  per cell and gene with `--cell-col`, `--gene-col`, `--count-col`.
- Output: `pb_counts.tsv` (genes by `<sample>::<cell type>`, kept pseudobulks only) and
  `pb_meta.tsv` (every sample and cell type with `n_cells`, `total_counts`, `kept` and labels).
- Read the JSON: `dropped`, then per cell type `units_per_group`, `testable`, `reasons`, and
  `repeated_units` (patients with more than one sample in a group; resolve them as in step 1).
- For a large AnnData or SingleCellExperiment, sum in place (recipes) and run the script on the
  exported cell metadata alone: with no gene columns it still reports cell numbers and units.

`scripts/composition_test.py` compares per-sample proportions between two groups.

```bash
python3 scripts/composition_test.py pb_meta.tsv --count-col n_cells --sample-col sample \
  --group-col condition --case treated --control control --patient-col patient_id --tsv-out comp.tsv > comp.json
python3 scripts/composition_test.py counts.tsv --count-col n_cells --group-col timepoint \
  --case on --control pre --pair-col patient_id --transform alr --reference B
```

- Input: one row per cell, or one row per sample and cell type with `--count-col`. Group, pair
  and patient labels must be constant within a sample. `--patient-col` stops the run if a patient
  is in both groups (use `--pair-col`) or has two samples in one group.
- Output per cell type: medians per group, `hl_shift` (case minus control, proportion units or
  log2 ratio), `rank_biserial` (from -1 to 1, positive means higher in case), `n_up`/`n_down`
  when paired, `p_value`, `p_floor`, `q_bh`, `method`. `notes` carries the floor, compositional
  and dropped-pair warnings; repeat them in the answer.
- Options: `--within` (denominator), `--min-sample-cells`, `--exclude-samples` (leave-one-out),
  `--method exact|normal`, `--proportions-out` (per-sample values for plotting).

## Report

Tell the researcher, in this order:

1. Design: patients per group, paired or not, unit of replication, samples or pairs excluded and why.
2. Data checks: raw counts confirmed, cell threshold and pseudobulks dropped, patients per group per
   cell type after filtering, proportion denominator, processing consistency.
3. Results with uncertainty: for composition, per cell type the medians, Hodges-Lehmann shift,
   rank-biserial r, direction count, p, p_floor and q; for DE, genes at the declared FDR per cell
   type with the patient n, log2 fold changes and the model formula.
4. Limits: patient numbers and the p floor, compositional reading, captured cells versus tissue,
   exploratory status where fewer than 4 patients per group remain. Research use only; nothing
   here concerns an individual patient's care.
5. What would confirm it: more patients or an independent cohort, an orthogonal abundance measure,
   targeted validation of DE genes.

Template:

```text
Design: <n_case> vs <n_control> patients (<paired|unpaired>); unit = patient. Excluded: <list, reason>.
Methods: counts summed per patient and cell type (min <10> cells); <edgeR QL|DESeq2|limma-voom>,
~ <formula>; BH across <family>. Proportions of <denominator> per patient, exact two-sided
Wilcoxon <rank-sum|signed-rank>, BH across <k> cell types.
Testable cell types (>= 4 patients per group): <list>; exploratory: <list>.
Result: <cell type> <median_case> vs <median_control> (shift <x> percentage points, r = <r>,
<k>/<n> patients in the same direction; p = <p>, smallest attainable p = <p_floor>, q = <q>).
Sensitivity: 30-cell threshold <same|changed: ...>; leaving out one patient at a time, direction
held in <k>/<n> runs; log-ratio to <reference> <same direction|not>.
Limits: <n> patients per group; relative proportions; exploratory, needs validation in <cohort/assay>.
```

For few patients, a sentence such as "exploratory: 3 vs 3 patients, smallest attainable exact
p = 0.1" (探索性结果：每组 3 名患者，精确检验最小可达 p = 0.1) is the honest summary.

## References

- [references/sources.md](references/sources.md): the bibliography; open it to check what a key
  supports before quoting a number to the researcher.
- [references/model-recipes.md](references/model-recipes.md): edgeR, DESeq2, limma-voom, dream,
  scanpy and muscat aggregation, and propeller code for the designs above; open it when writing
  the model-fitting code.
