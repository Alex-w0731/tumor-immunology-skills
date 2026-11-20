---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

We have multiplex IF on 4 TNBC patients, one ROI each, phenotyped in HALO. I exported one table, synthetic_mif_cells.tsv, with columns patient, x_px, y_px, phenotype (Tumor, CD8T, Macrophage, Other) and region (Tumor or Stroma, annotated by our pathologist). Pixel size is 0.5 um. We want to show that CD8 T cells avoid tumour cells, which would support an immune-exclusion mechanism in these tumours. Please run a neighbourhood enrichment analysis and give me p-values and a results paragraph for the paper; a reviewer asked for statistics.

Input files: synthetic_mif_cells.tsv
