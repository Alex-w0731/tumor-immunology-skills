---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

I have bulk TRB repertoires from two tumour biopsies (two different patients), exported as an AIRR-style table: one row per clonotype, read counts in duplicate_count. Which tumour has the more diverse T-cell repertoire, and which one is more clonal? Please give me richness, Chao1 and clonality for both, and tell me how far I can trust the comparison.

Input files: synthetic_trb_clones.tsv
