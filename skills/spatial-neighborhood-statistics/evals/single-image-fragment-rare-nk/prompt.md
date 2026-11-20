---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

I have one CODEX image from a resected tumour (one patient). After segmentation and phenotyping I exported synthetic_codex_cells.tsv with cell_id, x_um, y_um and cell_type (Tumor, Stroma, Macrophage, CD8T, Treg, NK). The paper we are following used k = 10 nearest neighbours for its neighbourhood analysis. Is the Treg-CD8 T cell co-localization in this tumour significant, and do NK cells avoid tumour cells? I need numbers for a figure legend.

Input files: synthetic_codex_cells.tsv
