# Model recipes for patient-level single-cell comparisons

Code for the designs in SKILL.md, written against the outputs of `pseudobulk_aggregate.py`
(`pb_counts.tsv`, `pb_meta.tsv`). Package interfaces change: these were written in November 2026,
so check the installed versions and their documentation before running. Replace `condition`,
`patient_id`, level names and the cell type with the study's own.

## Contents

- Aggregating large objects in place
- Loading one cell type in R
- edgeR quasi-likelihood (default)
- DESeq2 and limma-voom
- Several samples per patient: dream
- Composition with covariates: propeller

## Aggregating large objects in place

Sum raw counts, never normalised values, then apply the same rules as the script (minimum cells,
units per group). Running `pseudobulk_aggregate.py` on the exported cell metadata alone gives the
cell numbers and unit counts to check against.

```python
# scanpy >= 1.10: sum the raw-count layer per sample and cell type
import scanpy as sc
pb = sc.get.aggregate(adata, by=["sample_id", "cell_type"], func="sum", layer="counts")
counts = pb.layers["sum"]          # pseudobulks x genes; pb.obs holds sample_id and cell_type
```

```r
# muscat (Bioconductor): SingleCellExperiment with assay "counts"
library(muscat)
pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
```

## Loading one cell type in R

```r
counts <- as.matrix(read.delim("pb_counts.tsv", row.names = 1, check.names = FALSE))
meta <- read.delim("pb_meta.tsv", check.names = FALSE)
ct <- meta[meta$kept == "yes" & meta$cell_type == "CD8_T", ]
rownames(ct) <- ct$pseudobulk_id
ct$condition <- relevel(factor(ct$condition), ref = "control")
ct$patient_id <- factor(ct$patient_id)
y_counts <- counts[, ct$pseudobulk_id]
```

`check.names = FALSE` keeps the `sample::cell_type` column names intact.

## edgeR quasi-likelihood (default)

```r
library(edgeR)
y <- DGEList(y_counts, samples = ct)
design <- model.matrix(~ condition, data = ct)               # unpaired
# design <- model.matrix(~ patient_id + condition, data = ct) # paired: patients in both groups
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                                       # TMM; normLibSizes() in edgeR >= 4
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
res <- glmQLFTest(fit, coef = ncol(design))                   # last column = condition effect
tab <- topTags(res, n = Inf)$table
```

Collect the p-values of every cell type tested and apply BH across them all when the declared
family is all gene and cell-type tests: `p.adjust(p_all, method = "BH")`.

## DESeq2 and limma-voom

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(y_counts, colData = ct, design = ~ condition)  # paired: ~ patient_id + condition
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "control"))
```

```r
library(limma)
v <- voom(y, design)                     # y after filterByExpr and calcNormFactors
fit <- eBayes(lmFit(v, design))
topTable(fit, coef = ncol(design), number = Inf)
```

## Several samples per patient: dream

Use when patients contribute more than one sample and the contrast is between patients, or when
many pairs are incomplete.

```r
library(variancePartition)
form <- ~ condition + (1 | patient_id)
vobj <- voomWithDreamWeights(y, form, ct)
fit <- eBayes(dream(vobj, form, ct))
topTable(fit, coef = "conditiontreated", number = Inf)
```

## Composition with covariates: propeller

`composition_test.py` handles two groups, paired or not, without covariates. For covariates or
more groups:

```r
library(speckle)
propeller(clusters = cells$cell_type, sample = cells$sample_id, group = cells$condition,
          transform = "logit")
```

For a paired or covariate-adjusted design, build the transformed proportions with
`getTransformedProps()` and test them with `propeller.ttest()` and a design matrix that includes
the patient, following the speckle vignette. The results remain proportions and carry the same
compositional caveat as the rank tests.
