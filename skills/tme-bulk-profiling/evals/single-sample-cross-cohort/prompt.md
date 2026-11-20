---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

We scored our in-house cohort of 40 tumours (RNA-seq, log2(TPM+1), all run through the same
pipeline) with tmelens using the default z-score method. Two follow-ups:

1. A new biopsy from the same pipeline just came back. Can I run tmelens on that single sample, get
   its TIS_Ayers18 score and tell the treating oncologist whether it is above the cut-off from the
   Ayers et al. paper, i.e. whether the patient is likely to respond to pembrolizumab?
2. For our paper I want to show that our cohort has higher TIS_Ayers18 scores than a published
   Affymetrix cohort, which my student scored separately in tmelens with the same z-score method.
   Can I put the two sets of scores in one boxplot and run a Wilcoxon test?
