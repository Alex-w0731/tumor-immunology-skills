# Validation design

A published benchmark shows how a method behaved on someone else's data. The platform, tumor
type, purity range and reference decide how it behaves on yours, so validate the pipeline you
will run (method, reference matrix and preprocessing together), and re-validate after changing
any of them.

## In silico mixtures from single-cell data

1. Use scRNA-seq from the same tumor type, ideally from a protocol close to the bulk data.
2. Split cells into a reference part and a mixture part. [AvilaCobos2020] used 50:50 splits;
   keeping each individual's cells in a single split gave slightly higher RMSE and lower
   correlation with the same conclusions. Split by patient when there are enough patients,
   since that is the more realistic test (heuristic).
3. Build each mixture on the linear scale by summing or averaging the expression of sampled
   cells in known proportions. Include malignant cells across the purity range of the cohort,
   even when the reference lacks them, because real tumors contain them.
4. Add mixtures in which each cell type is absent (background prediction) and spike-in series of
   increasing abundance (minimal detection fraction), the two quantities defined in [Sturm2019].
5. Store the truth on both denominators: share of all cells, and share of the reference cell
   types. `mixture_benchmark.py --renormalize truth` derives the second from the first.
6. Summed pseudo-bulk weights each cell type by its RNA content per cell, so RNA shares and
   cell shares differ. Methods that correct for mRNA content target cell fractions
   [Racle2017] [Finotello2019] [Monaco2019]; compare them with a truth in cells. (heuristic)

## Physical mixtures and orthogonal measurements

- In vitro mRNA mixtures test a signature without tissue effects; MCP-counter used them
  together with IHC [Becht2016].
- Flow cytometry and IHC on a subset of the cohort's own samples are the closest check of real
  performance [Becht2016] [Racle2017] [Finotello2019]. Match what is measured: flow cytometry
  gives fractions of acquired cells and IHC gives densities in a section. quanTIseq converted
  fractions to cells per mm2 using total cell densities from images [Finotello2019]. Rank
  correlation is the fair comparison between arbitrary-unit scores and densities. (heuristic)
- With few samples, report n with every correlation and avoid claims about single cell types.

## Reading mixture_benchmark.py output

| Metric | Question it answers | If it is poor |
|---|---|---|
| per-cell-type r, rho | Are samples ranked correctly for this cell type? | Do not compare this cell type between samples |
| slope, intercept, bias | Is the scale right? | Report between-sample comparisons only, no absolute values |
| RMSE | Typical error in fraction units | Compare with the group difference you hope to detect |
| mean_est_when_true_zero | Background level | Read estimates at or below it as not detected |
| per-sample r | Is the composition within a sample right? | Do not rank cell types within samples |
| overall pooled r | Little on its own | It mixes differences in mean abundance between cell types |

Decision rules (heuristics):

- A cell type with low per-cell-type r in mixtures that resemble your data: leave it out of
  conclusions, or merge it with its spillover partner (`--merge`) and benchmark again.
- High r but slope far from 1: the cell type can be compared between samples but its values
  are not fractions to report.
- Good per-cell-type metrics with poor per-sample r: the method ranks samples but not cell
  types; restrict claims accordingly.
- A warning about row sums means truth and estimate use different denominators: renormalize
  before reading bias or RMSE.

## Mistakes to avoid

- Building mixtures from the same cells as the reference, which overstates accuracy.
- Reporting correlation alone, which hides bias and denominator errors.
- Leaving malignant cells out of mixtures meant to represent tumors.
- Treating pseudo-bulk accuracy as proof for a different platform or tissue; state where the
  benchmark data came from next to every number derived from it.
