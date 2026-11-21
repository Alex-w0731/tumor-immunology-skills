---
type: llm
focus: trace
weight: 1
---

The reply warns that ssGSEA scores depend on the whole gene-set collection (ssGSEA as implemented in GSVA normalises by the range of scores over all sets and samples), so changing the collection after scoring, for example by restoring a held-out gene or set, changes the scores, and it tells the user to settle the open items before scoring or to compare scores only within one run. The explanation in brackets is not required; the warning and the advice are.
