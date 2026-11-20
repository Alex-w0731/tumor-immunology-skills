---
type: llm
focus: trace
weight: 1
---

The answer defines the multiple-testing family as every cell type, endpoint and model actually run (all 22 cell types, multiplied by the number of endpoints when more than one is tested), not only the results that will be reported, and applies Benjamini-Hochberg FDR across that whole family instead of selecting cell types with raw p < 0.05.
