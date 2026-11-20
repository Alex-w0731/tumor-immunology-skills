# Method notes

What the original papers and the independent benchmark say about each method, checked against
the cited text in November 2026. Software versions, web services and licence terms change; check
them at the source before relying on them.

## Contents

CIBERSORT and LM22; CIBERSORTx; EPIC; quanTIseq; MCP-counter; xCell; TIMER; tmelens and
nnls_deconvolve.py; what the benchmark covered.

## CIBERSORT and LM22

- nu-support vector regression on a signature matrix [Sturm2019].
- LM22 has 547 genes that distinguish 22 human hematopoietic phenotypes: seven T cell types,
  naive and memory B cells, plasma cells, NK cells and myeloid subsets. Marker numbers were
  chosen to give the lowest 2-norm condition number (11.4, computed with R's `kappa`)
  [Newman2015]. Other implementations give similar but not identical values.
- Works best on non-log linear values; the original comparison ran every method in non-log
  linear space [Newman2015].
- Reports a global empirical P value for each sample's deconvolution and no P values for
  individual cell types [Newman2015].
- Not explicitly tested on RNA-seq data in the original paper; the authors expected its
  linearity assumptions to hold there [Newman2015].
- Relative mode gives fractions of total immune content, for comparisons within a sample.
  Absolute mode gives a score in arbitrary units that reflects absolute proportion, comparable
  between samples and cell types, but it is not a fraction [Sturm2019].

## CIBERSORTx

- Extends CIBERSORT to infer cell-type-specific expression profiles from bulk tissue, and
  minimizes platform-specific variation so that single-cell references can be used on bulk
  data [Newman2019].
- Terms as of November 2026: free for non-commercial use only; use by any commercial entity,
  including for research, is prohibited; the user gets a revocable, non-transferable right to
  use the service, and all other rights are reserved [Stanford2026]. The terms do not mention
  redistributing downloaded files, so they give no permission for it: obtain LM22 and the
  software from the distributor and point others there.

## EPIC

- Constrained least squares; six immune cell types plus fibroblasts and endothelial cells in
  the benchmark's classification [Sturm2019].
- RNA-seq reference profiles from circulating immune cells or from tumor-infiltrating
  non-malignant cells [Racle2017].
- Works on TPM, which its full cell-proportion prediction relies on [Racle2017]; an
  independent benchmark found TPM gave EPIC its best results [AvilaCobos2020].
- Renormalizes predicted mRNA proportions by cell-type mRNA content and estimates an
  uncharacterized cell fraction such as cancer cells [Racle2017].
- Validated with flow cytometry, IHC and scRNA-seq of melanoma and colorectal tumors
  [Racle2017].

## quanTIseq

- Ten immune cell types; constrained least squares [Sturm2019]; signature built from RNA-seq of
  blood-derived immune cells [Finotello2019].
- TPM input (its pipeline can compute TPM from reads). Genes above 700 TPM were excluded
  because highly expressed genes can bias deconvolution, and signature genes with median log2
  expression above 7 in all non-hematopoietic cancer cell lines of CCLE were removed
  [Finotello2019].
- Scales fractions by per-cell-type mRNA content factors [Finotello2019].
- Estimates an "other" fraction of cells absent from the signature. It is not tumor purity,
  because it can include stromal cells [Finotello2019].
- Validated with simulated data, flow cytometry and IHC [Finotello2019].

## MCP-counter

- Eight immune and two stromal populations; the score is the mean expression of marker genes
  [Becht2016] [Sturm2019].
- Designed to compare the abundance of one population across samples; scores are arbitrary
  units, conceptually close to IHC densities [Becht2016]. Marker-based scores cannot be
  compared between populations [Sturm2019].
- The units depend on the expression platform: compare only samples from the same platform
  [Becht2016].
- Validated with in vitro mRNA mixtures and IHC [Becht2016].

## xCell

- 64 immune and stromal cell types, with signatures from 1822 pure cell transcriptomes; raw
  scores are averaged ssGSEA scores [Aran2017].
- Raw scores of RNA-seq samples are not linear in abundance and do not allow comparison across
  cell types; xCell transforms them using in silico mixtures and applies spillover
  compensation between closely related types [Aran2017].
- The authors report that transformed scores resembled known fractions in simulations, enabling
  comparison across cell types, and also state that the scores are enrichment scores that
  cannot be reliably used as proportions [Aran2017]. The independent benchmark classes xCell
  output as arbitrary units comparable between samples only [Sturm2019]. Treat comparisons
  across cell types as unvalidated unless your own mixtures support them.
- Results depend on the other samples analysed in the same run, so the same sample can score
  differently in another batch; xCell detects no signal in a few similar samples [Sturm2019].

## TIMER

- Linear least squares for six immune cell types; arbitrary units comparable between samples of
  one cancer type, not across cancer types or cell types [Sturm2019].

## tmelens and nnls_deconvolve.py

- tmelens `deconvolve` linearizes log2(x + 1) input, fits NNLS on the genes shared with the
  signature, rescales to sum 1, divides each gene by its signature mean by default, and reports
  per-sample Pearson r and RMSE. `make_signature_matrix` keeps the top marker genes per cell
  type by fold change (at least 2) from the user's own reference profiles. tmelens ships no
  licence-restricted matrix.
- `scripts/nnls_deconvolve.py` fits the same model with the standard library and adds the
  condition number, collinear pairs, marker coverage, residual genes, a flag for identifiers
  that match only when case is ignored, and a check of the input scale. NNLS is solved by the
  Lawson-Hanson active-set algorithm [Lawson1995], with each passive-set sub-problem solved
  through the normal equations rather than the original Householder QR.
- Row scaling divides each gene of the signature and of the mixture by the same factor, so it
  re-weights genes in the least-squares fit and leaves an exact fit unchanged. [AvilaCobos2020]
  advise against row scaling applied to the mixture and the reference separately, which
  changes the mixing relation; that is a different operation. `--no-row-scale` gives the
  unweighted fit for comparison.

## What the benchmark covered

[Sturm2019] compared seven methods on nine cell types, using bulk samples simulated from about
11,000 tumor microenvironment cells and independent gold-standard estimates. Accuracy was high
for well-defined cell-type signatures; spillover explained many errors (for example quanTIseq
predicted a considerable amount of macrophages and monocytes in a pure cancer-associated
fibroblast sample, traced to five genes expressed by both), and addressing it improved
"fuzzy" signatures. [AvilaCobos2020] compared preprocessing and methods on pseudo-bulk mixtures
from five scRNA-seq datasets: methods with scRNA-seq references performed comparably to the best
bulk methods.
