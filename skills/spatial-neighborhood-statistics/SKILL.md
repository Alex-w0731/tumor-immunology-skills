---
name: spatial-neighborhood-statistics
description: >-
  Tests whether cell types are spatially associated (co-localized) or segregated in tissue from a
  table of cell or spot coordinates with labels, from multiplexed imaging (IMC, MIBI, CODEX, CyCIF,
  multiplex IF) or spatial transcriptomics. Covers neighbour-graph choice (k nearest neighbours or a
  radius in micrometres), the label-permutation neighbourhood enrichment test and its null, edge
  effects, confounding by tissue compartments and density, the patient as the unit of replication,
  combining slides, multiple testing across cell-type pairs and reporting. Bundles a seeded
  standard-library permutation script. Use when running or reviewing squidpy nhood_enrichment or
  histoCAT neighbourhood analysis, when a result says two cell types interact or avoid each other,
  or when ROIs are pooled. Terms include neighborhood enrichment, co-localization, cell proximity,
  空间邻域富集, 细胞共定位, 空间互作, 置换检验.
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.1"
  domain: tumor-immunology
  summary: Tests cell-type co-localization and segregation in tissue by label permutation, per slide, with the patient as the replicate.
  summary-zh: 基于标签置换检验细胞类型的空间共定位与分离：逐切片分析，以患者为重复单位。
---

# Spatial neighbourhood statistics for cell types

A researcher has cell (or spot) coordinates with cell-type labels from one or more tissue sections
and asks whether two types sit closer together, or further apart, than expected. Deliver a
per-sample neighbourhood test with a stated graph in micrometres and a stated null, effect sizes
that can be compared across samples, a patient-level summary that claims no more replication than
there are patients, and wording that keeps spatial association apart from cell-cell interaction.

## Scope

- Covers: neighbour graphs, the label-permutation enrichment test, restricted nulls for
  compartments, edge handling, combining slides and patients, multiple testing, reporting.
- Not covered: segmentation and phenotyping, spot deconvolution, spatially variable genes,
  ligand-receptor inference, cellular-neighbourhood clustering.
- Use `single-cell-patient-level-stats` for small-sample group tests on per-patient values,
  `immune-survival-analysis` or `response-biomarker-evaluation` when a spatial score is related to
  outcome or response.

## Workflow

### 1. Audit the table and the design before building any graph

Build a sample sheet: one row per sample (slide, ROI, core, field of view) with patient, group,
compartment annotation available or not, cell count per type. Then check:

- **Units.** Get micrometres per pixel from the image metadata (field size divided by pixel count),
  not from the stated resolution: [Keren2018] imaged 800 um fields as 2048 pixels "at 500nm
  resolution" (500 nm was the ion-beam probe size), so a pixel was 800/2048 = 0.39 um and 100
  pixels were 39 um. A radius typed in the wrong unit silently changes the biology being tested.
  If the coordinates are already in micrometres, pass `--um-per-unit 1` so the output says so.
- **One coordinate frame per sample.** Exported tables often stack slides with overlapping
  coordinates. A graph over the merged table links cells from different tissues. Exact duplicate
  coordinates are the usual symptom; the script reports them.
- **Replication.** Count patients, not cells or ROIs. Several ROIs or cores from one tumour are
  repeated measurements of one patient [Wilson2021]; treating them, or the cells, as independent
  is pseudoreplication [Hurlbert1984].

Decide what the data can support, before seeing results:

| Independent patients | What can be claimed |
|---|---|
| 1 | Description of this tissue only; per-sample p-values test random labelling inside it |
| 2 to 5 | Per-patient effects and direction counts; no two-sided sign test can reach p < 0.05 (floor 2^(1-n) = 0.0625 at n = 5; `derivation`) |
| 6 or more | Sign test across patients possible (6 of 6 in one direction gives p = 0.031; `derivation`) |
| Two groups of patients | Compare per-patient effects with a patient-level test or a mixed model (step 7) |

### 2. Fix the question and the pairs first

Write down which pairs are hypotheses (for example CD8 T cells next to tumour cells) and whether
the question is directional ("do CD8 T cells have more tumour neighbours than chance?") or symmetric.
All pairs of t types means t(t + 1)/2 unordered or t^2 ordered tests per sample (78 unordered pairs
for 12 types), and every one belongs to the multiple-testing family. Pre-specifying a few pairs is
the honest way to keep power with few patients; label everything else exploratory.

### 3. Build one neighbour graph per sample, in micrometres

Never let a graph cross samples. Pick the graph to match the question:

| Graph | Use when | Watch for |
|---|---|---|
| k nearest neighbours, union (symmetric) | Cell density varies and every cell should have neighbours; the default | Neighbour distance varies with density; in sparse areas edges can bridge gaps and tissue fragments, so cap edge length |
| Fixed radius (um) | The question has a physical scale (contact, short-range signalling) | Neighbour counts rise with density; sparse cells may be isolated (the script counts them) |
| Directed kNN | You want each cell's own k nearest, an asymmetric question | Pairs become ordered: A to B differs from B to A |

- State the neighbourhood in micrometres even for kNN (the script reports median edge length and
  the median nearest-neighbour distance).
- Keep a radius well above the median nearest-neighbour distance: below it, at least half the
  cells have no neighbour at all (`derivation`), and the tests run on a handful of edges.
- On a regular spot grid set the radius between the first and second ring of spots, not exactly at
  the spot spacing: rounding then decides at random whether neighbours count (`heuristic`).
- Inspect the diagnostics in `samples[].graph`: mean degree, isolated cells, longest edge. A
  longest kNN edge many times the median nearest-neighbour distance (the script's note fires
  above 5 times; `heuristic`) means bridged gaps; add `--max-dist` or `--max-dist-factor`.
- Plan a sensitivity analysis with a second neighbourhood definition (another k or radius);
  [Keren2018] reported that their enrichment results were robust to the distance chosen.

### 4. Choose the null: what the permutation keeps fixed

The test permutes labels over the fixed cell positions: a random-labelling null, which is a
different hypothesis from independence of two cell populations [Goreaud2003]. It asks whether
labels are exchangeable across the whole sample. Tumour nests, stroma, lymphoid aggregates and
necrosis make them non-exchangeable, so the unstratified test flags compartment structure itself:
two types that prefer different compartments come out "segregated" even if neither cell type
reacts to the other.

- Unstratified: describes tissue architecture. Report it as such.
- Stratified by compartment (`--strata-col`): labels move only within tumour or stroma, so the
  test asks whether the pair is closer than expected given compartment composition. This is
  usually the question behind "do CD8 T cells contact tumour cells?". [Keren2018] restricted
  their randomization in the same way (within the tumour compartment, or within lineage
  categories). Analysing each compartment as its own sample is the simpler alternative;
  [Wilson2022] restricted their TMA analyses to the tumour compartment.
- Block permutation (`--block-size`, um): labels move only within squares, which keeps large-scale
  density and composition gradients. Blocks must be several neighbourhoods wide; check
  `cells_in_single_type_strata`, because a block with one type has nothing to permute.
- Define compartments independently of the labels tested (pathologist annotation, a tissue
  classifier on other channels). A "tumour region" drawn around tumour cells removes the
  tumour-immune segregation by construction.

When compartments exist, report the unstratified and the stratified result side by side, per
sample: observed count, expected count under each null, and where the tested types sit (for
example the share of CD8 T cells in stroma). (E_unstratified - E_stratified) / (E_unstratified -
observed) is the share of the unstratified deficit or excess that compartment composition accounts
for (`heuristic`; near or above 100% means composition explains all of it; read it only when the
unstratified effect is clear).

Density acts the same way. With a radius graph, dense regions contribute most of the edges, so
their composition dominates the counts; kNN gives every cell about the same weight. If density
differs between compartments, run both graphs and say whether the conclusion changes.

### 5. Handle edges

Cells near the edge of the image or tissue lose neighbours that lie outside it [Wilson2022]. The
permutation null uses the same positions, so truncation affects observed and expected counts alike
and the test stays valid as a test of random labelling. What edges still cause: border artefacts
(crush, folds, segmentation failures) and a different composition at the tissue margin. Run a
sensitivity analysis with focal cells restricted to the interior:

- `--edge-buffer D --window bbox|hull`: only cells at least D from the bounding box or convex hull
  are focal; their neighbours anywhere still count. Use D of at least the radius (for kNN, the
  edge cap): every focal cell's whole neighbourhood then lies inside the window (the border or
  minus-sampling idea; `derivation`).
- For irregular tissue, holes or several fragments, the box and hull are wrong. Compute each
  cell's distance to the tissue mask boundary elsewhere and pass a 1/0 column with `--focal-col`,
  or analyse fragments separately.

### 6. Run the per-sample tests and read them correctly

Run `scripts/neighborhood_enrichment.py` (Scripts section). Read `notes` first. Then per sample:

- **Effect size.** Use `log2_ratio` = log2((observed + 0.5)/(expected + 0.5)). It does not grow
  with cell number, so it can be compared across samples analysed with the same graph definition;
  the pseudocount pulls it toward 0 when expected counts are small. `z` = (observed -
  expected)/SD as in [Palla2022] grows with the number of cells and edges: two copies of one
  tissue give the same log2 ratio and about 1.4 times the z (`derivation`, checked in the tests).
  A larger |z| in a bigger slide is not a stronger association; never compare or average z
  across samples.
- **p-values.** One-sided `p_enrich` and `p_deplete` as in histoCAT [Schapiro2017], computed as
  (b + 1)/(m + 1) so they are never zero [Phipson2010]; `p_two_sided` doubles the smaller; `q` is
  BH across the pairs of that sample [Benjamini1995].
- **Permutations.** With m permutations the smallest two-sided p is 2/(m + 1). For one pair to
  survive BH among P pairs when the others are null, 2/(m + 1) must be at most alpha/P, so
  m >= 2P/alpha roughly (`derivation`): 999 is enough for up to 25 pre-specified pairs; all 78
  pairs of 12 types need 2 x 78 / 0.05 = 3,120, so use 4,999 or more.
- **Rare types.** Avoidance of rare cell types is rarely detectable [Schapiro2017], because the
  expected number of contacts is already near zero. Types below `--min-cells` (default 10,
  `heuristic`) are not tested. Report "not assessable", never "no avoidance".

### 7. Combine evidence across samples at the patient level

The cohort question is whether the association holds across patients, and each patient is one
replicate. The script averages `log2_ratio` over a patient's samples (equal weight) to get one
effect per patient, then reports direction counts, the median and an exact two-sided sign test
[Dixon1946] with its floor 2^(1-n), and BH across pairs. The sign test treats each patient's
direction as a fair coin under random labelling. Counts are discrete and skewed, so that holds
only approximately, least well for small expected counts (`median_expected`; `derivation`):
read a borderline sign-test result for such pairs as approximate.

- For two groups (responders versus non-responders), compare per-patient effects with a
  rank-based test, or fit a linear mixed model with a random intercept per patient on per-sample
  scores, weighted for cell numbers, as spicyR does [Canete2022]. [Wilson2021] cautions that
  subject means can understate variability, inflating type I error, unless a subject's repeated
  values agree; report the spread of ROIs within patients. A test on one mean per patient still
  has the patient as its unit; the mixed model gains most when patients contribute unequal
  numbers of ROIs or cells (this skill's reasoning). An exact two-sided rank-sum test cannot go
  below 2/C(n1 + n2, n1): 0.1 for 3 versus 3 patients, 2/252 = 0.0079 for 5 versus 5
  (`derivation`).
- Do not pool all cells of all patients into one test, and that includes a permutation test that
  sums contacts over samples while permuting within each: its p still describes these cells, not
  the patients, so keep it out of a results paragraph. Do not combine per-slide p-values with
  Fisher or Stouffer (they test whether random labelling fails anywhere, which structured tissue
  almost guarantees), and do not count "significant in 30% of images" as a test; histoCAT used
  that as a display rule [Schapiro2017].
- Expect heterogeneity: [Keren2018] found immune-tumour mixing and compartmentalization in
  different patients of one cohort. Show every patient's effect, not only the summary.
- One small ROI may miss phenotypes and structures; how much area is needed depends on how
  segregated the tissue is [Bost2023]. Say so when a patient contributes one small field.

### 8. Interpret association, not interaction

A significant enrichment means the two labels are adjacent more often than random labelling of
these positions predicts. Shared compartments, shared recruitment cues, density gradients,
segmentation errors at cell borders (one cell's signal assigned to its neighbour) and sampling
all produce it. It does not show contact, signalling, suppression or causation. Name the
orthogonal evidence that would test a mechanism: membrane contact at higher resolution,
ligand-receptor expression as a hypothesis, perturbation in a model system, an independent cohort.

## Defaults and decision rules

| Choice | Default or rule | Source |
|---|---|---|
| Neighbour definition in histoCAT | Each segmented cell expanded by a user-set number of pixels (1 to 6 in most analyses), so measured from the cell outline, not a centroid radius | [Schapiro2017] |
| Neighbour definition in a MIBI TNBC cohort | 100 pixels = 39 um | [Keren2018] |
| kNN size | k = 6 | `convention`; Squidpy builds the hexagonal Visium grid with n_neigh = 6 [Palla2022] |
| kNN edge cap | 1.8 x the sample's median nearest-neighbour distance | `heuristic` |
| Spot grid radius | Between ring 1 and ring 2 distances | `heuristic` |
| Null | Labels permuted within sample, positions and graph fixed | [Schapiro2017], [Palla2022] |
| Permutations | 999 default; at least 2P/alpha for P pairs | `derivation` |
| Permutation p | (b + 1)/(m + 1), one-sided for enrichment and depletion | [Phipson2010], [Schapiro2017] |
| Significance in histoCAT | p < 0.01 per one-tailed test | [Schapiro2017]; this skill adjusts across pairs instead |
| Multiple testing | BH across pairs per sample and across pairs in the summary; BY when dependence worries you | [Benjamini1995], [Benjamini2001] |
| Effect size | log2((O + 0.5)/(E + 0.5)); 0.5 pseudocount | `convention` |
| z-score | Within-sample only | [Palla2022] for the definition |
| Minimum cells per type per sample | 10 | `heuristic` |
| Edge buffer | At least the radius, or the kNN edge cap | `derivation` |
| Across patients | One effect per patient; exact sign test; floor 2^(1-n) | [Dixon1946], [Hurlbert1984] |
| Several ROIs per patient, group comparison | Mixed model with random patient intercept | [Canete2022], [Wilson2021] |

## Pitfalls

- **Slides stacked in one coordinate frame.** Edges join unrelated tissues and every pair looks
  mixed. Fix: `--sample-col`, one graph per sample.
- **Cells or ROIs as replicates.** A pooled p of 1e-30 from three patients is a statement about
  cells. Fix: one effect per patient; report the floor of the test.
- **Comparing z-scores between slides or groups.** z tracks cell number. Fix: log2 ratio, or a
  model that weights by cell numbers [Canete2022].
- **Calling compartment segregation "exclusion" or "avoidance".** Unstratified depletion of CD8 T
  cells next to tumour cells often just says CD8 T cells sit in stroma. Fix: stratify by an
  independent compartment annotation and report both results.
- **Circular strata.** Compartments drawn from the tested labels erase the effect. Fix: annotate
  compartments from histology or other channels, and say how.
- **Pixels taken for micrometres.** Fix: convert with `--um-per-unit` from the metadata. A
  histoCAT distance of a few pixels is an expansion of the cell outline; used as a centroid
  radius it finds almost no neighbours.
- **kNN bridging gaps.** A small fragment with fewer than k cells borrows neighbours from across
  the gap. Fix: cap edges, use a radius, or drop or separate the fragment.
- **"No avoidance" for rare types.** Absence of a significant depletion with a few cells is absence
  of power. Fix: report "not assessable".
- **Too few permutations for the family.** With 999 permutations and 78 pairs no single pair can
  pass a Bonferroni-level threshold. Fix: more permutations or fewer, pre-specified pairs.
- **Spot data read as cells.** A dominant label per spot discards the mixture, and co-occurrence
  inside one spot is invisible to a neighbour graph. Fix: say so; analyse deconvolved abundances
  separately if that is the question.
- **Association written as interaction.** Fix: the step 8 wording.

## Scripts

Run the script; do not re-implement it. Python 3.10+, standard library only, deterministic for a
given `--seed`. Run time grows with cells x neighbours x permutations; tens of thousands of cells per
sample with 999 permutations is practical. Try `--n-perm 99` first on large data.

```bash
python3 scripts/neighborhood_enrichment.py cells.tsv \
  --x-col x_px --y-col y_px --type-col phenotype --sample-col roi --patient-col patient \
  --um-per-unit 0.5 --graph radius --radius 25 --n-perm 999 --seed 1 \
  --tsv-out pairs.tsv --summary-tsv summary.tsv --out result.json
# capped kNN, labels permuted within annotated compartments, interior cells focal (ordered pairs)
python3 scripts/neighborhood_enrichment.py cells.tsv \
  --x-col x_px --y-col y_px --type-col phenotype --sample-col roi --patient-col patient \
  --um-per-unit 0.5 --graph knn --k 6 --max-dist-factor 1.8 --strata-col region \
  --edge-buffer 25 --pairs CD8T:Tumor,Treg:CD8T --n-perm 999 --seed 1 --out strat.json
```

- **Input:** TSV or CSV (comma for `.csv`), one row per cell or spot. Columns named by options:
  `--x-col`, `--y-col` (numbers), `--type-col`; optional `--sample-col`, `--patient-col` (each
  sample must map to one patient), `--strata-col`, `--focal-col` (1/0). With `--um-per-unit`,
  coordinates are converted on input and every distance option (`--radius`, `--max-dist`,
  `--edge-buffer`, `--block-size`) is in micrometres. The graph is kNN unless `--graph radius` is
  given; `--radius` without it, or `--k` with it, stops the run rather than being ignored.
- **Output JSON:** `distance_unit`, `parameters`, `samples` (cell counts, graph diagnostics, focal
  cells, strata, untestable types, p floor), `pairs` (per sample: `n_a`, `n_b`, `observed`,
  `expected`, `sd`, `z`, `log2_ratio`, `p_enrich`, `p_deplete`, `p_two_sided`, `q`, `call`),
  `across_samples` (per pair: `n_samples`, `n_units`, `median_expected`, `median_log2_ratio`,
  `n_pos`, `n_neg`, `sign_test_p`, `sign_test_p_floor`, `q`, `unit_effects` per patient), `notes`.
  `--tsv-out` and `--summary-tsv` write the two tables.
- **Reading it:** `call` is enriched or depleted when that sample's q < `--alpha`. Pairs are
  unordered for a symmetric graph with all cells focal; with directed kNN, an edge buffer or a focal
  column they are ordered (A to B: neighbours of type B around focal A cells), in every sample.
- **Errors:** exit code 2 with the offending column, line or option named.
- **Not included:** Delaunay graphs, 3D coordinates, K or L functions, group comparisons and mixed
  models; use the per-patient `unit_effects` or per-sample `log2_ratio` as their input.

## Report

The answer to the researcher contains:

```text
Answer: one or two sentences: which pairs are associated or segregated, in how many patients, under which null
Question: pairs (pre-specified / exploratory), directional or symmetric; P pairs tested
Data: N patients, S samples (per patient), cells per type per sample (smallest), um per pixel and its source
Graph: kNN k = .. (union, edges capped at .. um) or radius .. um; median nearest-neighbour distance, mean degree, isolated cells, longest edge
Null: random labelling within each sample [unstratified | within <compartment>, defined by .. | blocks of .. um]
Edges: focal cells restricted by .. / sensitivity result
Per sample: observed, expected, log2(O/E), one-sided p, q (z only within a sample)
Across patients: per-patient effects, n positive / negative, sign test p (floor ..), q across pairs; a results paragraph cites these, never a p pooled over samples
Sensitivity: second graph; stratified versus unstratified (observed, expected under each null); edge buffer
Not assessable: types below --min-cells, pairs with too few patients for a cohort-level test
Interpretation: spatial association or segregation in these tissues; not interaction or causation
Limits and next checks: patients, ROI size, label quality; orthogonal evidence that would test a mechanism
```

## References

- `references/sources.md`: open to check what a cited source says or to cite it in a write-up.
- `references/null-models-and-alternatives.md`: read when choosing between random labelling and
  other nulls, when density gradients dominate, when a distance-dependent summary (cross-type K or
  L) or a mixing score is requested, or when setting up a mixed model across ROIs and patients.
