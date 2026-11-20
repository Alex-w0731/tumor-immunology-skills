# Null models, alternative statistics and across-patient models

Read this when the default (labels permuted within each sample) does not obviously match the
question, when density gradients dominate a tissue, when someone asks for a distance-dependent
summary or a mixing score, or when ROIs and patients must go into one model.

## Contents

1. Random labelling versus independence
2. When labels depend on location
3. Distance-dependent summaries
4. Composition-dependent scores
5. Models across ROIs and patients

## 1. Random labelling versus independence

- **Random labelling** (this script, histoCAT [Schapiro2017], Squidpy [Palla2022]): the cell
  positions are fixed and the labels are exchangeable among them. It fits the question "given
  where the cells are, are these labels arranged at random?", for example cell states or types
  assigned to the cells of one tissue.
- **Independence** of two populations: each pattern keeps its own structure and one is moved
  relative to the other. It fits two populations that arise from separate processes.
- The two hypotheses have different null distributions, and choosing the inappropriate one can
  lead to misinterpretation, in some cases to opposite conclusions about interaction
  [Goreaud2003]. State which one was tested.
- Asymmetric nulls exist: [Keren2018] randomized the locations of Y-positive cells (their number
  fixed) while X-positive cells stayed in place, so the X-Y and Y-X results need not agree. The
  bundled script permutes all labels together; with ordered pairs, A to B and B to A differ only
  through the focal set or a directed graph.

## 2. When labels depend on location

- A label-permutation null assumes exchangeability across the whole sample. When a type's density
  varies across the tissue, an apparent association can reflect that varying density rather than
  anything local (this skill's reasoning). [Baddeley2000] developed methods, including an
  analogue of the K-function, for studying interaction between points when the pattern is
  spatially inhomogeneous.
- Options, from simplest: (a) permute within compartments annotated independently of the tested
  labels (`--strata-col`); (b) permute within blocks several neighbourhoods wide (`--block-size`),
  which preserves gradients larger than a block; (c) model the covariate (for example distance to
  the tumour border) in a regression on per-cell neighbour counts, which the script does not do.
- Diagnostic: put the unstratified and the stratified `log2_ratio` of each patient side by side.
  If the effect disappears within compartments, the finding is about compartment composition,
  which may still be biologically interesting, but it is not a local association.
- Block size trades off: small blocks remove almost all structure and leave little to permute
  (watch `cells_in_single_type_strata`); large blocks remove little. Report the size in
  micrometres and, if results depend on it, show two sizes (`heuristic`).

## 3. Distance-dependent summaries

- A neighbour graph tests one scale. Intertype (cross-type) K12 or L12 functions summarise the
  association of two types over a range of distances, and need a stated null (independence or
  random labelling) [Goreaud2003]; Squidpy implements Ripley's L, F and G [Palla2022].
- They need edge correction (isotropic and translation corrections are common), and with holes or
  uneven cell density the theoretical value under complete spatial randomness is biased, more so
  as holes increase; use a permutation-based reference computed over the observed cells instead
  [Wilson2022].
- The replication rules do not change: one curve or one summary value per sample, then patients as
  the replicates.

## 4. Composition-dependent scores

- [Keren2018] defined a mixing score as the number of immune-tumour interactions divided by the
  number of immune-immune interactions (cells in direct contact) and used it to separate cold,
  mixed and compartmentalized tumours.
- Such a ratio depends on how many immune and tumour cells a sample has, so it describes patients
  rather than testing association. To ask whether mixing exceeds chance, compare it with its
  expectation under label permutation of the same sample, as the enrichment test does
  (this skill's reasoning).

## 5. Models across ROIs and patients

- **Per-patient summary** (the script): mean `log2_ratio` over a patient's samples, then a sign test
  across patients [Dixon1946]. Distribution-free and transparent; every patient weighs the same;
  within-patient precision and spread are ignored. [Wilson2021] caution that subject means can
  understate variability unless a subject's repeated values agree, so report the ROI-level values
  too. The test still has the patient as its unit (this skill's reasoning).
- **Mixed model** on one score per image: score_ij = a_i + b x_i (+ covariates) with a random
  intercept a_i per patient and condition x_i as a fixed effect, as in spicyR, which also weights
  images for very different cell numbers and computed mixed-model p-values with Satterthwaite's
  approximation [Canete2022]. Check convergence of the fit [Wilson2021].
- With few patients per group the variance between patients is poorly estimated whatever the
  model; prefer the per-patient summary with the exact floor stated, and describe rather than test
  (`heuristic`).
- Use `log2_ratio` (or a comparable, cell-number-aware score) as the response, not z.
