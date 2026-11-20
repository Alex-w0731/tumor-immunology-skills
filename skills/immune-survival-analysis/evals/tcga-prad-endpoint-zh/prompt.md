---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

我在做TCGA-PRAD（前列腺癌，约500例患者）的免疫微环境与预后分析。现在的方案是：用CIBERSORT估计22种免疫细胞的比例，每种细胞分别用survminer的surv_cutpoint找最佳截断值分成高、低两组，以总生存期（OS）为终点画KM曲线、做log-rank检验，最后把p<0.05的细胞类型写进论文。请帮我审一下这个方案：哪里有问题、应该怎么改？请给出具体可执行的分析步骤（R或Python都可以）。
