---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

I'm putting together preliminary data for a grant. Attached is the TPM matrix for our 12 melanoma
biopsies (genes x samples; I cut it down to the tmelens infiltration and stroma signature genes
plus a few housekeeping and melanocyte genes to keep the file small) and the matching tmelens gene
sets I exported. Please log-transform the TPMs, score the signatures, call each tumour inflamed,
excluded or desert the way tmelens does, and give me the percentage of inflamed tumours in our
cohort for the preliminary-data section. If MEL03 comes out as excluded, I'd also like to write
that in MEL03 the T cells are kept out of the tumour bed.

Input files: expression_tpm.tsv, core_sets.gmt
