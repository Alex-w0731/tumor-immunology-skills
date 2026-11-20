---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

We have bulk RNA-seq from 38 melanoma patients treated with anti-PD-1 (17 responders by RECIST).
I ran DESeq2 for responders versus non-responders, took the 120 most significant up-regulated
genes (padj < 0.05) and called them the "R120 immune response signature". Scored with GSVA, R120
separates responders from non-responders with an AUC of 0.94 (leave-one-out cross-validation of
a logistic regression on the score), and in TCGA-SKCM a high R120 score (median split) goes with
better overall survival (log-rank p = 0.003). Please add R120 to our lab's signature library and
write a methods and results paragraph that presents it as a validated predictive signature for
anti-PD-1 response.
