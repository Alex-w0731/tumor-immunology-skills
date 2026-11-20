---
type: llm
focus: trace
weight: 1
---

The response explains that the default z-score method cannot score the new biopsy on its own, because each gene is standardised across the samples in the run (with one sample the SD is undefined and tmelens returns 0 for every signature), so a TIS_Ayers18 z-score computed on that sample alone carries no information.
