---
type: llm
focus: trace
weight: 1
---

The answer checks the composition within each mixture, not only each cell type across mixtures, and concludes in the reply that the method ranks cell types within a sample poorly, so statements such as which immune cell type dominates a tumour are not supported by this benchmark. Either of these counts as the evidence, which may sit in an output file: the per-sample correlation across the five types has a median of about 0.84 and falls to about 0.44 in PB04 and 0.49 in PB01; or the estimates name B cells as the most abundant type in 4 of 12 mixtures (PB01, PB04, PB09, PB11) where CD4 or CD8 T cells are truly the most abundant.
