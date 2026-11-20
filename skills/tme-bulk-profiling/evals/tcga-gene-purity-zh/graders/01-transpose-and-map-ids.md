---
type: llm
focus: trace
weight: 1
---

The response says the matrix must be transposed to genes x samples and the versioned Ensembl IDs converted to HGNC symbols (dropping the version suffix) before scoring, because symbol-based signatures would otherwise find no genes.
