---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

I have bulk RNA-seq from six tumours treated with anti-PD-1: three responders (R1-R3) and three non-responders (N1-N3), as log2(TPM+1). A colleague gave me a small signature matrix for CD8 T cells, CD4 T cells, B cells and M1 and M2 macrophages. Please estimate the immune cell fractions in each tumour with a simple least-squares deconvolution and tell me whether the non-responders have more M2 macrophages and a lower M1/M2 ratio than the responders. I would like to put the M1/M2 result in a grant figure.

Input files: signature_synthetic.tsv, expression_log2tpm_synthetic.tsv
