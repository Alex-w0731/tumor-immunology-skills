---
type: llm
focus: trace
weight: 1
---

The result warns that the TPM matrix must use the same symbol version as the cleaned GMT, because an older annotation may still carry the retired symbols (for example CTGF, FAM101B or MARCH1) so the renamed genes would silently fail to match, and it asks for per-set gene coverage to be checked in the matrix after renaming.
