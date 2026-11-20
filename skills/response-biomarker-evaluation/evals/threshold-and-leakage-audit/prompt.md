---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

We built a 12-gene response signature by taking the genes most differentially expressed between responders and non-responders in our 40-patient NSCLC anti-PD-1 cohort, then scored the same 40 patients. We picked the cut-off that maximised Youden's index (0.31). At that cut-off sensitivity is 90% and specificity 85%, and the AUC is 0.93. Responders are CR/PR, non-responders SD/PD. Could you add PPV and NPV with confidence intervals and draft the results paragraph for the manuscript? We want to say the signature accurately predicts immunotherapy response.

Input files: synthetic_sig12_scores.tsv
