---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

I collected these immune signatures from the supplementary Excel tables of a few papers plus a
colleague's list and saved them as a GMT file (synthetic example for testing). I want to score
them with ssGSEA on our TCGA-SKCM TPM matrix (human, HGNC symbols) and show the scores in a
figure. Please check the file, write a cleaned GMT, and list every change you made and anything I
still need to look up before scoring.

Input files: immune_signatures_synthetic.gmt
