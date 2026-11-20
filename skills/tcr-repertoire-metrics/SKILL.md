---
name: tcr-repertoire-metrics
description: >-
  Computes and interprets T-cell receptor repertoire diversity, clonality and overlap from
  clonotype count tables that already exist (bulk TCR-seq exports such as AIRR TSV, or
  single-cell V(D)J contig tables): Shannon entropy, Pielou evenness, clonality as 1 - Pielou,
  Simpson and inverse Simpson, Gini, top-clone fractions, Chao1, exact rarefaction and seeded
  downsampling to a common depth, and Jaccard and Morisita-Horn overlap. Use when someone asks
  to compare TCR diversity or clonal expansion between samples, time points or tissues, to
  normalise for sequencing depth or cell number, or to quantify tumour-blood clonotype sharing.
  Triggers include TCR repertoire, clonality, clonotype diversity, Shannon index, Chao1,
  rarefaction, downsampling, Morisita-Horn, Jaccard, clonal expansion, TRB, scTCR-seq,
  TCR克隆性, 克隆型多样性, 克隆扩增, 降采样, 重叠度.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.0"
  domain: tumor-immunology
  summary: Computes TCR repertoire diversity, clonality and overlap from clonotype counts, compared at equal sequencing depth.
  summary-zh: 基于克隆型计数表计算TCR库多样性、克隆性与样本间重叠，并在相同测序深度下比较。
---

# TCR repertoire diversity, clonality and overlap

A researcher has clonotype tables from a bulk TCR-seq or single-cell V(D)J pipeline and wants
to know which sample is more diverse or more clonal, whether expansion changed between time
points, or how much two compartments share. Deliver metrics with explicit definitions,
computed on comparable data (same chain, clonotype key, count unit and depth), with the
spread from repeated subsampling and an honest statement of what the numbers cannot show.

## Scope

- Covers: statistics on existing clonotype tables, one receptor chain at a time: diversity,
  clonality, richness estimation, depth normalisation, and two-sample overlap.
- Does not cover: read alignment, error correction or clonotype calling (do these upstream in
  the sequencing pipeline), TCR specificity or antigen prediction, sequence-similarity
  clustering, or BCR lineage analysis.
- To relate a metric to survival use `immune-survival-analysis`; to test it as a response
  biomarker use `response-biomarker-evaluation`; for patient-level tests on single-cell data
  use `single-cell-patient-level-stats`; for immune context from bulk RNA use
  `tme-bulk-profiling`.

## Workflow

1. **Establish what one row and one count mean.** Bulk tables have one row per clonotype and
   a count of reads, UMI-collapsed molecules or template estimates; use molecules
   (`umi_count`) rather than reads when the pipeline gives both (convention). Single-cell
   contig tables have one row per contig, so clone size is the number of distinct cells
   (barcodes), never the number of rows or UMIs. Write the unit down: every metric is in that
   unit, and values in different units are not comparable. If the table holds only
   frequencies, Shannon, Pielou, clonality, Simpson, Gini and top-K can still be computed, but
   Chao1, singletons and downsampling cannot; ask for counts and state the limitation.

2. **Fix the clonotype key and the chain.** Default key (convention): CDR3 amino acid plus V
   gene. A nucleotide key gives at least as many clonotypes as an amino-acid key, and adding
   the V gene splits clonotypes further, so use one key for every sample and for both metrics
   and overlap. Analyse one locus per run: TRB for most bulk data; for single-cell data either
   TRB or the paired alpha-beta clone ID the pipeline assigned. Pooling TRA and TRB counts each
   cell twice and merges two distributions; the script refuses a multi-chain table unless
   `--chain` picks one.

3. **Filter and read the QC block.** Keep productive rearrangements, the population on which
   clonality was defined [Tumeh2014] (the script does this whenever a `productive` column
   exists). Check the per-sample QC: rows dropped as non-productive, other chain, missing CDR3
   or zero count; rows collapsed into one key; and, for cells, `cells_in_several_clonotypes`.
   Such a cell adds one unit to each of its clonotypes, so the depth then exceeds the cell
   number. If more than a few percent of cells are affected (heuristic), remove those barcodes
   (possible doublets) or switch to the pipeline's clone ID, and report the choice. For low-input bulk
   samples compare the clonotype count with the number of input cells: in RNA aliquots
   equivalent to 500 to 50,000 sorted mouse CD4 T cells, [Chaara2018] found about 20-fold
   (500 cells) to 2-fold (50,000 cells) more clonotypes than cells, and in the small samples
   the excess remained after removing singletons. More clonotypes than cells means richness
   and Chao1 of that sample are unreliable.

4. **Compare depths before comparing anything.** Every metric here changes with depth.
   Richness and Chao1 change most. Pielou evenness, and so clonality, shifted with depth
   relative to the number of cells in [Chaara2018]. Gini and top-K fractions are computed over
   observed clonotypes only, so they shift too.
   - Describing one sample: report observed values together with the depth.
   - Comparing samples: subsample every sample without replacement to one common depth, by
     default the smallest sample depth, as VDJtools does [Shugay2015]. Rarefaction to a common
     depth is what makes richness comparable [Gotelli2001]. Choose the depth before looking at
     the results, never upsample, and name every excluded sample.
   - With three or more samples, if the smallest depth is below a fifth of the median depth
     (heuristic; `repertoire_metrics.py` adds a note), report the comparison twice: at the
     smallest depth with every sample, and at a fixed higher depth that excludes the shallow
     samples.
   - Check that a ranking does not depend on the depth chosen: repeat at a second, lower
     depth or draw rarefaction curves [Gotelli2001]. If the order changes, report that
     instead of a winner.
   - Repeat the subsampling (default 20; 50 to 100 for final figures), record the seed, and
     report mean with SD or 2.5 to 97.5 percentiles. For richness, the exact expected value at
     that depth [Hurlbert1971] is printed alongside; the simulated mean should agree with it
     within a few SD / √(repeats).
   - The spread is subsampling noise only. The shallowest sample is used whole and has none,
     and the spread is not a confidence interval for the repertoire. Never mix downsampled and
     full-depth values in one comparison.

5. **Pick the metric that answers the question.** Expansion or dominance: clonality, top-1
   and top-10 fractions, inverse Simpson (an effective number of clonotypes that weights
   abundant ones; lower means more dominated). Breadth: richness at equal depth and
   exp(Shannon). Unseen richness: Chao1, only as a lower bound (see Pitfalls). Metrics can
   disagree because richness is driven by rare clonotypes and Simpson or clonality by abundant
   ones: a sample can have both a long singleton tail and strong dominance. Report the
   disagreement and explain it instead of picking the convenient number. Before comparing
   with "clonality" values from another tool or paper, check how they were defined.

6. **Overlap, when asked.** Use the same chain, key, filters and a common depth (subsample to
   the smaller depth). Report Morisita-Horn (abundance-weighted and overly sensitive to the
   most abundant clonotypes [Rempala2013]), Jaccard (presence and absence), and the shared
   fraction in both directions, because the share of tumour counts found in blood and the
   share of blood counts found in tumour answer different questions.
   - Depth moves Jaccard in either direction. Subsampling removes rare clonotypes first, so
     Jaccard falls when the shared clonotypes are rare and rises when the private ones are
     mostly singletons (worked example in `references/formulas.md`). Only depth-matched values
     are comparable. A shared fraction falls as the other sample gets shallower.
   - Low indices next to a high shared fraction in one direction usually mean that one
     sample's largest clonotypes are private. Look at which clonotypes dominate each sample
     before calling the repertoires distinct.
   - A clonotype not detected in the shallower or more diverse sample may be present below
     detection: say "not detected at this depth", not "absent".

7. **Statistics across samples.** Use one value per sample per metric; the patient, not the
   cell, read, clonotype or subsampling repeat, is the unit of replication. Pair pre and post
   or tumour and blood samples from the same patient. Several metrics summarise the same
   distribution and move together, so name one primary metric before testing and adjust the
   rest for multiplicity. With one sample per group, describe and do not test. This is
   research analysis, not a decision rule for an individual patient.

8. **Report** with the template below.

## Defaults and decision rules

| Decision | Default | Source |
|---|---|---|
| Count unit | molecules (UMIs) over reads when both exist; distinct cells for single-cell data | convention |
| Clonotype key | CDR3 amino acid + V gene, identical for every sample and analysis | convention |
| Rows used | productive rearrangements only | [Tumeh2014] |
| Clonality | 1 − *H* / ln *S* over productive clonotypes (= 1 − Pielou evenness); log base cancels | [Tumeh2014], [Pielou1966], [Shannon1948] |
| Simpson naming | `simpson` = Σ*p*², `gini_simpson` = 1 − Σ*p*², `inverse_simpson` = 1/Σ*p*²; say which one you report | [Simpson1949], [Hill1973] |
| Effective numbers | richness, exp(*H*), inverse Simpson = Hill numbers of order 0, 1, 2 | [Hill1973] |
| Gini coefficient | over observed clone sizes, population form (formula in references) | convention |
| Chao1 | bias-corrected, *S* + *f*₁(*f*₁ − 1)/(2(*f*₂ + 1)) as in vegan `estimateR`; integer counts only; a lower bound | [Oksanen2026], [Chiu2014] |
| Chao1 for T cells | expect underestimation of total diversity and dependence on sample size | [Laydon2015] |
| Depth normalisation | random subsampling without replacement to a common depth | [Chaara2018], [Gotelli2001] |
| Common depth | smallest sample depth (`--downsample min`) | [Shugay2015] |
| Very shallow sample | smallest depth below a fifth of the median (three or more samples): report at the smallest depth and at a fixed higher depth without it | heuristic |
| Robust ranking | the same order at a second depth or along rarefaction curves | [Gotelli2001] |
| Exact richness at depth *D* | Hurlbert rarefaction expectation | [Hurlbert1971] |
| Repeats | 20 by default, 50 to 100 for final numbers (with 20, the 2.5 and 97.5 percentiles are nearly the extremes) | heuristic |
| Depth warning | scripts flag samples or pairs whose depths differ more than 2-fold | heuristic |
| Expanded clonotype | ≥ 2 cells (single-cell data only) | convention |
| Overlap indices | Morisita-Horn for abundance-weighted similarity, Jaccard for presence/absence, shared fractions in both directions | [Horn1966], [Jaccard1912], [Rempala2013] |
| AIRR column names | `junction_aa`, `v_call`, `locus`, `productive`, `duplicate_count` / `umi_count` (as of November 2026) | [AIRR2026] |

## Pitfalls

- **Richness or Chao1 compared at unequal depth.** The deeper sample looks more diverse by
  construction, and the ranking can reverse once depths match. Fix: subsample to a common
  depth (Workflow step 4) before ranking.
- **Clonality treated as depth-free.** In [Chaara2018], Pielou evenness of unfiltered data
  fell steeply as the ratio of clonotypes to cells rose, so clonality moves with depth too.
  Fix: compare clonality on downsampled values and give the depth.
- **"Clonality" from another tool compared with this one.** Some studies, for example those
  using the immunoSEQ Analyzer, report "productive Simpson clonality", the square root of
  Simpson's index Σ*p*² over productive rearrangements [Coghill2025]. Both run from 0 to 1 but
  they are different functions: counts 4, 2, 1, 1 give √0.34375 = 0.59 against 1 − Pielou =
  0.125. Fix: recompute from counts (the square root of the `simpson` output gives the
  former), or compare only values with the same definition.
- **Chains mixed.** One table with TRA and TRB rows counts each cell twice and blends two
  distributions. Fix: one locus per run, or the pipeline's paired clone ID per cell.
- **Contigs, UMIs or reads counted as cells.** In single-cell data, a cell with two contigs
  or many UMIs would get extra weight. Fix: `--cell-col` with the barcode column.
- **Different clonotype keys across samples or analyses.** An amino-acid key merges
  convergent nucleotide variants; a key without the V gene merges unrelated receptors with the
  same CDR3. Either choice is defensible, but switching between samples is not. Fix: one key,
  stated in the report.
- **"Simpson index" or "Gini" left ambiguous.** Simpson's λ, 1 − λ and 1/λ move in different
  directions, and the Gini coefficient is not the Gini-Simpson index. Fix: name the exact
  quantity and give its formula.
- **Chao1 reported as repertoire size.** It is a lower bound built from singletons and
  doubletons [Chiu2014], it is unstable when doubletons are few relative to singletons, every
  false singleton (a sequencing or PCR variant of a larger clonotype) pushes it up, and
  over-sequenced small samples carry such artefactual clonotypes [Chaara2018]; it still
  underestimates the diversity of a person's repertoire [Laydon2015]. Fix: call it a
  depth-dependent lower bound, give *f*₁ and *f*₂, and compare it only between samples
  processed the same way at the same depth.
- **Frequencies treated as counts.** Chao1, singletons and subsampling need integer counts;
  multiplying frequencies by a guessed total invents data. Fix: obtain the counts.
- **One random subsample, or shallow samples silently upsampled.** A single draw hides its
  own noise, and upsampling invents clonotype counts. Fix: repeat with a seed, report the
  spread, and exclude and name samples below the target depth.
- **Pseudo-replication, including subsampling repeats.** Treating cells or clonotypes as
  independent observations inflates significance. A test across subsampling repeats is worse:
  the repeats re-draw one library, so adding repeats makes any difference "significant".
  Fix: one value per sample, patients as the unit, paired tests for paired designs.
- **Overlap read as presence or absence.** Jaccard and shared fractions depend on depth, and a
  tumour clonotype missing from blood may sit below the blood's detection limit. Fix: compare
  overlap at a common depth and phrase results as "detected at this depth".
- **Clonotypes outnumbering input cells.** In low-input libraries this signals artefacts, and
  removing singletons does not cure it [Chaara2018]. Fix: report it and treat richness and
  Chao1 as unreliable.

## Scripts

Run these scripts; do not read or rewrite them. Paths are relative to this skill's folder.
Both read TSV or CSV (by extension, optionally `.gz`) and exit with code 2 and a message
naming the bad column, value or file when the input needs fixing. Downsampling is fast (well
under a second per repeat at 100,000 units), so 50 to 100 repeats cost little.

`scripts/repertoire_metrics.py`: run it for per-sample diversity, clonality and richness,
with or without downsampling.

```bash
# bulk AIRR table, all samples in one file, TRB only, downsampled to the smallest depth
python3 scripts/repertoire_metrics.py --input clones.tsv --sample-col sample \
  --clone-col junction_aa,v_call --count-col duplicate_count --chain TRB \
  --downsample min --iterations 50 --seed 1 > metrics.json

# single-cell contig table: clone size = distinct barcodes
python3 scripts/repertoire_metrics.py --input contigs.csv --sample-col sample \
  --cell-col barcode --clone-col cdr3 --chain TRB --downsample min --format tsv
```

For the pipeline's paired clone ID, pass that column as `--clone-col` together with
`--cell-col barcode --chain TRB`: each cell is then counted once through its TRB contig, and
cells without a TRB contig are dropped.

`scripts/repertoire_overlap.py`: run it for pairwise Jaccard, Morisita-Horn and shared
fractions. It loads `repertoire_metrics.py` from the same folder and takes the same input
options, plus `--pairs tumour:blood` to restrict the pairs and fix their order (`_a` is the
first name). With `--downsample min` the common depth is the smallest among the samples in
the requested pairs.

```bash
python3 scripts/repertoire_overlap.py --input tumour_blood.tsv --sample-col sample \
  --clone-col junction_aa,v_call --count-col duplicate_count --chain TRB \
  --pairs tumour:blood --downsample min --iterations 50 --seed 1
```

Input options shared by both scripts:

| Option | Meaning |
|---|---|
| `--clone-col` | required; column(s) forming the clonotype key, comma-separated |
| `--count-col` or `--cell-col` | exactly one: a count per row, or a cell barcode (clone size = distinct cells) |
| `--sample-col` | defaults to `sample` if present, else each file is one sample named after the file without `.gz` and its table extension (`P1.pre.tsv` becomes `P1.pre`) |
| `--chain-col`, `--chain` | chain column defaults to `locus` or `chain`; several values present means `--chain` is required |
| `--productive-col`, `--keep-nonproductive` | productive column defaults to `productive`; values T, TRUE, 1, yes are kept |
| `--downsample`, `--iterations`, `--seed` | common depth (integer or `min`), repeats (20), seed (1) |
| `--top` | metrics only; K values for top-K fractions (default 1,10,100) |
| `--format`, `--out` | `json` (default, full detail) or `tsv` (one row per sample or pair) |

Reading the output: JSON has, per sample, `qc` (rows read, used and dropped by reason, chain,
`integer_counts`, `rows_collapsed` or `cells`), `observed` (full-depth metrics) and, when
downsampling, `downsampled` with `depth`, `iterations`, `expected_richness` and, for each
metric, `mean`, `sd`, `p2_5`, `p97_5` and `n_undefined`. `notes` lists unchecked assumptions
(no chain or productive column), samples pooled from several files, excluded samples, depth
differences above 2-fold and, with `--downsample min`, a smallest depth below a fifth of the
median; carry each note into the report. Overlap JSON gives `observed` and
`downsampled` per pair, where `_a` and `_b` follow the `sample_a` / `sample_b` order. Null
means undefined (for example clonality with a single clonotype, or Chao1 from frequencies).

## Report

The answer to the researcher contains:

1. **Data and definitions**: chain, clonotype key, count unit, filters with rows dropped, and
   each sample's depth (and cell number for single-cell data).
2. **Depth handling**: common depth and how it was chosen, repeats, seed, excluded samples,
   and whether the ranking held at a second depth.
3. **Results table** at the common depth (mean and spread over repeats), with full-depth
   values only as context; say that the spread is subsampling noise only.
4. **Interpretation**: which metric answers the question, and why metrics disagree if they do.
5. **Statistics**: the unit of replication and the test, or why no test was run.
6. **Limits and checks**: depth and coverage, possible artefactual clonotypes, Chao1 as a lower
   bound, single chain, sample size; what would confirm the finding (replicate libraries,
   UMI-based counts, independent samples, more patients). Mark synthetic or example data as
   such.

```text
Data: TRB, productive only; clonotype = CDR3 aa + V gene; unit = UMIs.
      S1 12,400 UMIs; S2 3,050 UMIs (4.1-fold difference).
Depth: both subsampled to 3,050 UMIs, 50 repeats, seed 1; no sample excluded;
       same ranking at 1,500 UMIs.
| sample | richness at 3,050, mean (2.5-97.5%) | clonality | inverse Simpson | top-10 fraction | Chao1 (lower bound; f1, f2) |
Spread: subsampling noise only; S2 is used whole and has none.
Interpretation: ...
Statistics: one sample per patient, so descriptive only; no test.
Limits: one library per sample; Chao1 underestimates total diversity; research use only.
```

## References

- `references/formulas.md`: open it to write a methods section, to check an exact formula or
  edge case (single clonotype, frequencies, ties, Jaccard and depth), or to explain a number
  to a reviewer.
- `references/sources.md`: the bibliography behind every source key in this skill, with what
  each source supports; open it before quoting a source.
