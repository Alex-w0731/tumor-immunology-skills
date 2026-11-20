# Formulas and conventions

Exactly what `repertoire_metrics.py` and `repertoire_overlap.py` compute, so results can be
described in a methods section and reproduced in other software.

## Contents

- Notation
- Diversity and clonality
- Richness
- Downsampling
- Overlap
- Methods paragraph template

## Notation

For one sample after filtering: *c*ᵢ is the count of clonotype *i* (reads, UMIs, cells or
templates, whatever unit the table holds), *N* = Σ *c*ᵢ the depth, *S* the number of observed
clonotypes, *p*ᵢ = *c*ᵢ / *N*, and *f*₁, *f*₂ the numbers of clonotypes seen exactly once and twice.
All sums run over observed clonotypes; unobserved clonotypes never enter any formula.

## Diversity and clonality

* **Shannon entropy** *H* = −Σ *p*ᵢ ln *p*ᵢ, computed as ln *N* − Σ *c*ᵢ ln *c*ᵢ / *N* [Shannon1948].
  Natural log. Entropy in bits = *H* / ln 2; the base cancels in evenness and clonality, but
  not in *H* itself, so state it when quoting *H*.
* **Pielou evenness** *J* = *H* / ln *S* [Pielou1966], clamped to [0, 1] against rounding.
* **Clonality** = 1 − *J* [Tumeh2014]. Undefined when *S* = 1 (0/0): the scripts return null
  and count such repeats as `n_undefined`; do not substitute 1 without saying so.
* **Simpson** λ = Σ *p*ᵢ² [Simpson1949], the probability that two units drawn with replacement
  belong to the same clonotype. **Gini-Simpson** = 1 − λ; **inverse Simpson** = 1 / λ.
  Hurlbert's finite-sample version of 1 − λ multiplies by *N* / (*N* − 1) [Hurlbert1971]; the
  difference is negligible at repertoire depths and is not reported.
* **Hill numbers** [Hill1973]: order 0 = *S*, order 1 = exp(*H*) (`exp_shannon`), order 2 = 1 / λ.
  They share one unit (effective number of clonotypes), which makes them easier to compare
  than *H* and λ, and they weight rare clonotypes less as the order rises.
* **Gini coefficient** of observed clone sizes *x*₍₁₎ ≤ … ≤ *x*₍S₎:
  *G* = Σₖ (2*k* − *S* − 1) *x*₍ₖ₎ / (*S* Σ *x*), which equals the mean absolute difference over
  twice the mean (population form, no *S* / (*S* − 1) factor). Convention (no verified
  repertoire source): implementations differ in that factor and in whether unobserved
  clonotypes count, so give the formula when reporting it.
* **top-K fraction** = (sum of the *K* largest *c*ᵢ) / *N*; ties at the boundary do not change it.
* **Expanded clonotypes**: clonotypes with *c*ᵢ ≥ 2 and their share of *N*. Convention for
  single-cell data, where the unit is a cell; with reads or UMIs it mainly reflects depth.

## Richness

* **Observed richness** *S*: the most depth-sensitive metric; compare only at equal depth.
* **Chao1, bias-corrected**: *S* + *f*₁(*f*₁ − 1) / (2(*f*₂ + 1)), as in vegan `estimateR`
  (`S.chao1`) [Oksanen2026]. Chao's estimator is a lower bound built from singletons and
  doubletons [Chiu2014]; the bias correction reduces but does not remove the downward bias
  [Oksanen2026]; for T cells it frequently underestimates total diversity and depends on
  sample size [Laydon2015]. It needs integer counts and becomes unstable when *f*₂ is small
  relative to *f*₁. Every false singleton raises it, and over-sequencing small samples adds
  artefactual clonotypes [Chaara2018].
* **Expected richness at depth *D*** (exact rarefaction) [Hurlbert1971]:
  E[*S*_D] = Σᵢ [1 − C(*N* − *c*ᵢ, *D*) / C(*N*, *D*)], computed with log-gamma; a term is 1
  when *N* − *c*ᵢ < *D*. It equals the mean richness of infinitely many random subsamples,
  so the simulated mean must agree with it within a few SD / √(iterations).

## Downsampling

* Each repeat draws *D* units without replacement from the *N* units of a sample
  (multivariate hypergeometric), the scheme used for TCR rarefaction in [Chaara2018]:
  *D* distinct unit indices by Floyd's algorithm, mapped to clonotypes by cumulative counts.
* Randomness comes only from `random.random()` of a generator seeded with the string
  "<seed>:<sample name>", so a sample's draws do not change when samples are added or
  reordered. Python documents that `random()` with the same seed gives the same sequence
  across versions; the scripts rely on that guarantee, which was tested on Python 3.10 only,
  so record the Python version with the seed.
* `--downsample min` uses the smallest sample depth, as VDJtools does [Shugay2015]. Samples
  below an explicit depth are excluded and listed in `notes`; nothing is upsampled.
* Summaries over repeats: mean, SD (*n* − 1 denominator), and the 2.5 and 97.5 percentiles by
  linear interpolation. These describe subsampling variability only, not biological or
  library-preparation variability.

## Overlap

For samples A and B with proportions *p*ᵢ and *q*ᵢ:

* **Jaccard** = |A ∩ B| / |A ∪ B| over observed clonotypes [Jaccard1912]; presence and
  absence only. Depth moves it in either direction, because subsampling loses rare clonotypes
  first (worked example below), so compare Jaccard values only at a common depth.
* **Morisita-Horn** = 2 Σ *p*ᵢ *q*ᵢ / (Σ *p*ᵢ² + Σ *q*ᵢ²) [Horn1966], the form analysed for TCR
  populations by [Rempala2013]. It is 1 for identical relative abundances and 0 with nothing
  shared, and it is driven mainly by the most abundant clonotypes [Rempala2013].
* **Shared fraction of A** = Σ over shared clonotypes of *c*ᵢ(A) / *N*_A, and likewise for B.
  These two values are not symmetric. Read them as the share of a sample in clonotypes
  detected in the other sample at that depth, not as proof that the rest is absent. The
  shared fraction of A falls as B gets shallower.

**Worked example: Jaccard and depth** (synthetic; checked in the unit tests). Each sample
has 110 units and is subsampled to 55, so each singleton survives with probability 0.5.

| Pair | Composition | Jaccard, full depth | Jaccard at 55 units (mean of repeats) |
|---|---|---|---|
| A, B | one shared clonotype of 100 units; 10 private singletons in each | 1/21 = 0.048 | about 0.09 (rises: private singletons drop out) |
| C, D | one private clonotype of 100 units in each; 10 shared singletons | 10/12 = 0.83 | about 0.26 (falls: a shared singleton is kept in both with probability 0.25) |

The shared fraction of A stays at 100/110 = 0.91 at both depths because the shared clonotype
is always drawn; that of C halves (0.09 to about 0.045).

## Methods paragraph template

> TCR repertoire metrics were computed for each sample from productive [TRB] rearrangements,
> defining a clonotype as [CDR3 amino-acid sequence plus V gene] and clone size as the number
> of [UMIs / reads / cells]. Shannon entropy (natural logarithm), Pielou evenness, clonality
> (1 − evenness), Simpson's index and its inverse, the Gini coefficient of clone sizes and the
> fraction of [units] in the [10] largest clonotypes were calculated. To remove the effect of
> sequencing depth, each sample was subsampled without replacement to [D] [units] (the
> smallest sample depth), [R] times with random seed [s], and the mean and 2.5–97.5
> percentiles over repeats are reported; the expected richness at depth D was computed
> exactly (Hurlbert 1971). Bias-corrected Chao1 is reported as a lower bound of richness.
> Repertoire overlap was quantified with the Morisita-Horn and Jaccard indices at a common
> depth. Calculations used repertoire_metrics.py and repertoire_overlap.py
> (tumor-immunology-skills, tcr-repertoire-metrics 0.1.0).
