# Eval results

Each skill ships three eval cases in `skills/<name>/evals/`. For this first run, one case per
skill was executed twice by fresh agents: once with the skill's `SKILL.md` to follow and once
without it, each in its own folder whose path revealed neither the skill nor the configuration.
A third agent graded both runs against the case's fixed statements without being told which run
had the skill, reading the replies and the files each run produced.

- Date: 2026-11-20
- Model for executors and graders: Claude Opus 5.5
- Iteration 1; one run per configuration, so a difference of one statement is within noise.

## Results

| Skill | Case | With skill | Without skill | Grader preferred overall |
|---|---|---|---|---|
| [`citation-evidence-check`](../skills/citation-evidence-check/SKILL.md) | `grant-significance-agent-draft` | 5/5 | 5/5 | with skill |
| [`immune-deconvolution-review`](../skills/immune-deconvolution-review/SKILL.md) | `benchmark-pseudobulk-zh` | 6/6 | 6/6 | with skill |
| [`immune-signature-curation`](../skills/immune-signature-curation/SKILL.md) | `excel-damaged-gmt-audit` | 6/6 | 5/6 | with skill |
| [`immune-survival-analysis`](../skills/immune-survival-analysis/SKILL.md) | `optimal-cutpoint-claim` | 5/5 | 4/5 | with skill |
| [`response-biomarker-evaluation`](../skills/response-biomarker-evaluation/SKILL.md) | `label-definition-sensitivity` | 5/5 | 5/5 | with skill |
| [`single-cell-patient-level-stats`](../skills/single-cell-patient-level-stats/SKILL.md) | `paired-composition-on-treatment` | 6/6 | 5/6 | with skill |
| [`spatial-neighborhood-statistics`](../skills/spatial-neighborhood-statistics/SKILL.md) | `cd8-tumor-compartments` | 5/5 | 5/5 | with skill |
| [`tcr-repertoire-metrics`](../skills/tcr-repertoire-metrics/SKILL.md) | `single-cell-pre-post` | 5/5 | 2/5 | with skill |
| [`tme-bulk-profiling`](../skills/tme-bulk-profiling/SKILL.md) | `melanoma-tpm-phenotypes` | 5/5 | 2/5 | with skill |
| **All** | 9 cases | **48/48 (100%)** | **39/48 (81%)** | with skill in 9 of 9 |

The last column is the grader's overall judgement of which reply would serve the researcher
better, taken from its written comparison; it was made without knowing which run had the skill.

## Where the two runs differed

**`immune-signature-curation` / `excel-damaged-gmt-audit`**: passed only with the skill.
> The result reports that the set name Checkpoints is defined twice with different genes (8 genes and 3 genes), explains that a reader would silently keep only one definition (dictionary-based readers keep the last, 3-gene one), and resolves it explicitly by renaming, merging or dropping one with a s…

- With the skill: Reply change #8: 'Checkpoints (7) | Second definition removed | Same name as line 3. Readers that store sets by name keep only the last definition (this 3-gene one) without saying so; other readers may merge or reject it. Its genes (PDCD1, CD274, CTLA4) are a…
- Without the skill: Reported and resolved explicitly: both sets renamed to Checkpoints_8gene and Checkpoints_3gene ('Two different sets had the same name'; log: 'renamed by size so none is overwritten. Decide which one you mean'). The explanation is only partial: 'Depending on t…

**`immune-survival-analysis` / `optimal-cutpoint-claim`**: passed only with the skill.
> In place of the searched cut-off, the answer presents a Kaplan-Meier comparison at a split fixed without reference to outcome, such as the median (60 versus 60 patients), and reports that it is not significant (log-rank p about 0.17).

- With the skill: km_median_split.png shows the median split at −0.098: CD8 low n = 60 (39 deaths) vs high n = 60 (34 deaths), subtitled "Median split chosen without reference to outcome; log-rank p = 0.17". Reply: "The figure for the paper should be the median split. Your 0.6…
- Without the skill: Only partly done. A reports a median-split log-rank p as non-significant: the reply says "A median split gives p = 0.17", and the results text gives "HR 0.72, 95% CI 0.45–1.15; p = 0.17". I confirmed this split is 60 vs 60 with p = 0.1725. But the median spli…

**`single-cell-patient-level-stats` / `paired-composition-on-treatment`**: passed only with the skill.
> The results sentence the response proposes does not call the CD8_T change (or any other cell-type change) statistically significant; it reports the CD8_T increase with its direction count or effect size and labels it exploratory or preliminary.

- With the skill: The sentence reads: "In an exploratory paired analysis of five patients, CD8+ T cells made up a larger share ... in all five (median 9.0% vs 16.4%; Hodges-Lehmann shift +8.1 percentage points; exact two-sided Wilcoxon signed-rank p = 0.0625, the smallest valu…
- Without the skill: The proposed sentence never calls the finding exploratory or preliminary. It opens with "increased in all five paired patients (median 9.0% to 16.4%; paired t-test on logit-transformed proportions, P = 0.016; FDR = 0.094...)" and goes on to quote P = 0.035 an…

**`tcr-repertoire-metrics` / `single-cell-pre-post`**: passed only with the skill.
> Diversity is compared at an equal number of cells (post subsampled to the 30 pre cells, with repeats), and the answer reports that at equal cell number post has fewer TRB clonotypes than pre (about 17 versus 26) even though its raw count is higher (37 versus 26).

- With the skill: The reply says "两个样本都无放回下采样到30个细胞…重复100次，seed = 1". Its table row is "克隆型数 | 26 | 17.0 [13–21] | 37", and the text says "全部数据下，治疗后的克隆型数（37）比治疗前（26）还多…在同样深度下，治疗后的克隆型其实更少（17对26）". The JSON gives post richness at 30 cells as mean 16.98 (2.5–97.5%: 13–21), Shanno…
- Without the skill: Only partly done. The subsampling is there: "从治疗后样本中随机抽 30 个细胞，重复 10,000 次" (the R script samples without replacement, N_RAREFY 10000, seed 20261120), and it reports Shannon 2.43 and clonality 0.142 at 30 cells. But the reply never gives the clonotype count a…

**`tcr-repertoire-metrics` / `single-cell-pre-post`**: passed only with the skill.
> The answer reports the share of cells in expanded clonotypes (at least 2 cells) at both time points: about 23% before treatment and about 55% to 63% after.

- With the skill: The reply table row is "≥2个细胞的克隆型所占细胞比例 | 23%（3个克隆型） | 54% [40–68%] | 63%（7个克隆型）". The JSON expanded_fraction values are 0.2333 pre, 0.5447 post at 30 cells and 0.625 post at full depth. My recount gives 7/30 = 23.3%, 50/80 = 62.5% and about 54.9% at 30 cells.
- Without the skill: The reply never states what share of post-treatment cells are in expanded clonotypes. For pre it only says "其中 77% 是只出现一次的克隆（Good's coverage 0.23）", which implies 23% by subtraction. For post it lists per-clone counts (25/80, 10/80, five clones of 3) that the…

**`tcr-repertoire-metrics` / `single-cell-pre-post`**: passed only with the skill.
> The answer states that with a single patient the pre/post change is descriptive and cannot by itself support a statistical test or a conclusion about response to treatment.

- With the skill: The reply says "不过只有一位患者，治疗前也只有30个细胞，所以这只是描述性结果，没有做统计检验" and, in the statistics section, "只有一位患者，每个时间点一个样本，而分析单位是患者，所以只做描述，不做检验。不能把细胞或下采样重复当作独立样本来算p值". It also says "仅供科研分析，不能用于个体患者的临床判断". Treatment response is covered by this clinical-judgment disclaimer, no…
- Without the skill: The reply runs cell-level tests and presents them as evidence: "置换检验…Shannon 差值为 −0.77，p = 0.004…克隆性差值为 +0.125，p = 0.012" and "CT001…Fisher 检验 q = 0.002，是唯一单独达到显著的克隆". Its caveat ("只有一位患者，每个时间点一份组织。结论只适用于这两份样本。p 值只反映抽到哪些细胞带来的随机性…也不能推广到其他患者") limits how far th…

**`tme-bulk-profiling` / `melanoma-tpm-phenotypes`**: passed only with the skill.
> The response explains that tmelens-style calls label the top third of the infiltration axis as inflamed, so the inflamed percentage (4 of 12) is set by the threshold and cannot be reported as the cohort's prevalence of inflamed tumours.

- With the skill: Section 3 says "It comes out at 4 of 12 (33%), but tmelens calls the top third of the cohort inflamed by definition. Any 12 tumours with distinct scores give exactly 4", and "the number says nothing about how often your melanomas are inflamed". The rule is gi…
- Without the skill: Only partly done. B never states that the tmelens rule labels the top third inflamed; it says "I had tmelens's exported gene lists but not its scoring code or cut-offs" and uses its own mean-split rule. Its headline is "4 to 6 of the 12 tumours (33–50%)" and…

**`tme-bulk-profiling` / `melanoma-tpm-phenotypes`**: passed only with the skill.
> The response reports that only 7 of the 19 Pan_F_TBRS genes are in the matrix, below the 50% coverage tmelens requires, so that signature is dropped and the stromal axis rests on CAF_core alone.

- With the skill: The reply says "Pan_F_TBRS (the TGF-beta signature): 7/19 genes (37%). That is below tmelens' 50% rule, so tmelens drops it without warning. The stromal score is therefore CAF_core alone" and lists the 12 missing genes. phenotype_run.json has stromal_axis_mem…
- Without the skill: B reports "Only 7 of 19 Pan-F-TBRS genes are there" and lists the missing genes. It never mentions the 50% tmelens coverage rule and does not drop the signature. tme_phenotype_scoring.R sets stroma_sets <- c("Pan_F_TBRS", "CAF_core") with no coverage filter.…

**`tme-bulk-profiling` / `melanoma-tpm-phenotypes`**: passed only with the skill.
> The response states that the phenotype labels (and any z-score based scores) are relative to these 12 samples and would change if samples were added or removed.

- With the skill: The reply says "Both scores are relative to these 12 samples only. *MEL11 becomes inflamed if any one of the four inflamed tumours is left out." and "Each gene is standardised across these 12 samples, so every score is relative to this cohort". Under Limits i…
- Without the skill: Only partly done. B says "the cut-offs are relative to this cohort, they largely decide the answer" and that scores are "z-scored across the 12 tumours". It never states that the labels or z-scores would change if samples were added or removed. There is no le…

## How to read this

- The baseline is strong. A current frontier model working without the skill already passes
  most statements, so the skills add value mainly at specific decision points: refusing an
  outcome-selected cut-off as the main figure, not calling a 5-pair comparison significant,
  not running cell-level tests on one patient, applying the documented tmelens rules. The
  differences above are those cases.
- Several statements passed in both configurations. The graders flagged the ones they judged
  non-discriminating or ambiguous; those notes are kept below and drive the next revision of
  the eval cases. A statement that every reasonable answer passes measures nothing.
- The statements were written together with the skills, so they reward what the skills teach.
  The grader did not know which run had the skill, but the cases were not designed by a party
  independent of the skills.
- Per-statement evidence for both runs is in [`evals/iteration-1.json`](evals/iteration-1.json).

## Grader notes on the statements

- `citation-evidence-check`: Statements 1, 2 and 4 would pass for any answer that actually looked up the references: Crossref and PubMed both flag the retraction, and the abstracts state the 18% (14 of 76) figure and the mouse-only design.
- `citation-evidence-check`: Statement 2 names a 'cumulative' rate. B gives the same 18% (14 of 76) all-dose figure without that word, so a strict grader and a lenient one could disagree.
- `citation-evidence-check`: Statement 3's clause 'unless it states that it looked that paper up and checked it against the sentence' does not say how explicit that statement must be. B shows it only through a PMID in a 'Checked against' column and the phrase 'rewrote the sentence to say what that paper found'.
- `citation-evidence-check`: Statement 5 penalizes only calling something supported on metadata alone. It does not catch claims that go beyond a stated abstract-only basis (B's '40% does not appear anywhere in the paper'), and it does not cover how replacement statements were checked, so it barely separates the two attempts.
- `immune-deconvolution-review`: Statement 2: it is unclear whether the parenthetical 'median about 0.47' is required. Neither reply states the median: A's benchmark_as_given.json warning gives 0.469, and B's files do not give it. The clause 'instead of reading bias or RMSE from the raw comparison' is also hard to grade, because both answers report raw bias and RMSE to show the denominator mismatch before switching to the rescal…
- `immune-deconvolution-review`: Statement 1: it is unclear whether per-type metrics that appear only in output files or figures count as 'the answer reports'. B's reply text names per-type error values only for B cells and CD4; the other three types are in metrics_by_celltype.tsv and fig1.
- `immune-deconvolution-review`: Statements 5 and 6 cover standard caveats (circular reference, transfer to TCGA, an orthogonal IHC or flow check) that any competent answer would raise, so they did not separate the attempts here.
- `immune-deconvolution-review`: No statement separated these two attempts: all six pass for both.
- `immune-signature-curation`: Statement 2 depends on whether its parenthetical ('dictionary-based readers keep the last, 3-gene one') is required. Under a lenient reading, B's 'kept only one of the two sets' would pass. It is the only statement that separates the two attempts.
- `immune-signature-curation`: Statement 3 would pass for almost any answer that ran a symbol check. Both attempts made the same two renames and logged the originals, so it does not discriminate.
- `immune-signature-curation`: Statement 6 is unclear about whether keeping the one-gene set in the cleaned GMT counts as 'presenting it as a signature score' when the text says to show FOXP3 expression instead (B's case).
- `immune-signature-curation`: Statement 5 accepts either orthologue mapping or setting the set aside, and does not require a named orthology resource or version. So an unsourced mapping (B) passes as easily as a documented set-aside (A).
- `immune-signature-curation`: Statement 1 accepts both restoring with a note and removing, by design. So it does not separate the more cautious choice (A removes the fields pending the source) from writing an inferred symbol into the scoring file with a warning only in the log and reply (B).
- `immune-survival-analysis`: Statement 5: 'In place of the searched cut-off' can be read two ways. It is unclear whether a median-split log-rank p reported next to a retained searched-cut KM figure (Run A) should count; I graded it strictly, as requiring the outcome-independent split to replace the searched one as the presented KM analysis.
- `immune-survival-analysis`: Statements 1, 3 and 4 would pass for nearly any competent answer that declines the claim and fits the obvious Cox models. Statements 3 and 4 overlap, since both come down to fitting the per-SD Cox model, so none of the three separated these two runs.
- `immune-survival-analysis`: Statement 2 offers two alternatives. Run A's 44% figure comes from permuting this dataset rather than simulating new cohorts, but both runs pass anyway on the permutation p of about 0.13.
- `response-biomarker-evaluation`: Statements 1 and 3 cover basic choices that any careful answer makes (CR/PR vs SD/PD on baseline samples; on-treatment biopsies set aside). Both attempts pass, so these two don't separate the runs.
- `response-biomarker-evaluation`: Statement 2: the parenthetical '(18 responders and 46 non-responders)' doesn't say whether the reply must give that split. Attempt B's reply gives only n=64; the 18/46 split appears only in its auc_results.tsv.
- `response-biomarker-evaluation`: Statement 4: 'labelled as a sensitivity analysis' is ambiguous when an attempt labels the SD-dropped AUC 'secondary', 'exploratory' or 'not for the headline' instead, as Attempt B does.
- `response-biomarker-evaluation`: Statement 5: it isn't clear whether the words 'predictive' and 'prognostic' are required. Attempt B explains the distinction in plain words without them.
- `response-biomarker-evaluation`: No statement checks the good-result point that the data file is synthetic. Attempt A flagged this and Attempt B did not, which is the clearest real difference between the runs, and it goes ungraded.
- `single-cell-patient-level-stats`: 1 (paired comparison, no pooled test): any reasonably careful answer passes it, and both attempts do.
- `single-cell-patient-level-stats`: 2 (P6 left out, 5 pairs): any careful answer passes it, and both attempts do.
- `single-cell-patient-level-stats`: 3 (0.0625 floor): it passes even when the floor is cited only as a reason to switch to a paired t-test that can reach p < 0.05, as Attempt A does. It doesn't check whether the attempt actually respected the limit.
- `single-cell-patient-level-stats`: 4 (CD8 5/5 plus effect size): 'effect size on the proportion scale' is ambiguous. It is unclear whether group medians (9.0% → 16.4%) or a fold change count, or whether a paired median or Hodges-Lehmann shift in pp is required. Attempt A passes only under the looser reading.
- `spatial-neighborhood-statistics`: Statement 5 does not separate the runs. B puts a pooled ROI-stratified p < 0.001 in its results paragraph, but the wording ('rather than one p-value obtained by pooling all cells') passes any response that reports per-patient effects alongside a pooled p.
- `spatial-neighborhood-statistics`: Statement 2 is ambiguous when a response states the sign-test floor but also reports a pooled cell-level (ROI-stratified) p in the paper paragraph. It does not say whether that breaks 'treats each of the four patients as one replicate'.
- `spatial-neighborhood-statistics`: Statement 4's last clause ('says what further evidence a mechanism would need') is ambiguous to grade in two cases: when a response names candidate mechanisms a snapshot cannot identify (B: ECM barrier, chemokine gradients, TGF-β) rather than evidence types, and when that text sits only in an output file rather than the reply.
- `spatial-neighborhood-statistics`: Statement 3 passes for almost any answer that converts to µm and tries a second radius. It does not test the more discriminating point from the rubric, that the radius must exceed the ~17 µm median nearest-neighbour distance; at 17 µm about half the cells are isolated.
- `spatial-neighborhood-statistics`: Statement 1 would pass for any answer that runs a region-stratified null on this table, because the effect disappears so completely. It does not separate competent runs.
- `tcr-repertoire-metrics`: Statements 1 and 2 pass for both attempts, so they do not tell this pair apart: both avoid pooling chains and both count productive-only cells, giving 30/26 and 80/37.
- `tcr-repertoire-metrics`: Statements 3 and 4 depend on whether a value that appears only in an output CSV counts as 'the answer reports'. Attempt A's diversity_metrics.csv has the clonotype count at 30 cells (17.06) and frac_cells_in_expanded (0.233/0.625/0.549), but its written reply does not state them. I graded on the written reply.
- `tcr-repertoire-metrics`: Statement 5 bundles two claims: the change is descriptive with no statistical test, and it cannot support a treatment-response conclusion. It is unclear whether a general 'not for clinical judgment' disclaimer satisfies the second part, and whether an attempt that reports p-values with a single-patient caveat deserves partial credit.
- `tme-bulk-profiling`: Statement 1 (log scale already applied) and statement 4 (MEL03 excluded call is not spatial; use CD8 IHC/mIF) passed for both runs, so they do not separate these two attempts. Both would still catch a naive answer that logs the data again or accepts the MEL03 sentence.
- `tme-bulk-profiling`: Statement 5 is borderline for Run B. B's "cut-offs are relative to this cohort" implies that the labels depend on the cohort, but it never states the added/removed-samples consequence that the statement requires. A lenient grader might pass it.

## Reproducing

Every case folder follows the `claude plugin eval` layout (`prompt.md`, `case.yaml`, `graders/`), so
the same cases can be run with Claude Code's own eval runner:

```bash
claude plugin eval . --allow-tools Write,Edit,Bash --scaffold
```
