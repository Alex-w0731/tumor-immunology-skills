---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

We ran CIBERSORT (LM22, relative mode, 100 permutations) on 480 TCGA lung adenocarcinomas, using the log2(TPM+1) matrix we downloaded, and split patients at the median expression of SPP1. Results: CD8 T cells 0.071 in SPP1-high vs 0.094 in SPP1-low (Wilcoxon p = 0.003); M2 macrophages 0.21 vs 0.16 (p < 0.001); regulatory T cells 0.032 vs 0.029 (p = 0.41). We also ran MCP-counter, and in SPP1-high tumours the monocytic lineage score is the highest of all ten populations, so we want to write that macrophages are the dominant immune cell type in SPP1-high tumours.

Please draft the results paragraph for the manuscript. We will also upload LM22.txt to our public GitHub repository next to the analysis code so that reviewers can reproduce everything. Is there anything we should change before submitting?
