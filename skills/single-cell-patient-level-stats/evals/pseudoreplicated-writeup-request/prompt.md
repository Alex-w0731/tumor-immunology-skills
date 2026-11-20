---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

I have 10x scRNA-seq from pre-treatment melanoma biopsies of 3 anti-PD-1 responders and 3 non-responders, about 45,000 cells after QC. I subset the CD8 T cells and ran Seurat FindMarkers (default Wilcoxon) responders vs non-responders: 1,214 genes with adjusted p < 0.05, including GZMB, PRF1 and HAVCR2 higher in responders. CD8 T cells were also 31% of immune cells in responders versus 18% in non-responders, and a Fisher's exact test on the pooled cell counts gives p < 1e-50. Please draft the results paragraph for our manuscript and tell me which statistics to report.
