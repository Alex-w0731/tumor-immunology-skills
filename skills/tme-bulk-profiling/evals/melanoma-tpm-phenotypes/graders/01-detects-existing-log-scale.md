---
type: llm
focus: trace
weight: 1
---

The response identifies that expression_tpm.tsv is already log-scaled (values top out near 13 and none exceed 50) even though it is called TPM, and it does not apply a second log2 transform before scoring.
