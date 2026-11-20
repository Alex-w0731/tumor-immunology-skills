---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

We sequenced TCRs from a tumour biopsy and matched peripheral blood of the same patient (bulk, AIRR-style export, read counts in duplicate_count). How similar are the two repertoires? Please compute Jaccard and Morisita-Horn and tell me what fraction of the tumour repertoire is also present in the blood.

Input files: synthetic_tumour_blood_trb.tsv
