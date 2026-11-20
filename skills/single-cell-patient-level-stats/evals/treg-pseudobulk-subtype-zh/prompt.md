---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

我们做了 9 个乳腺癌肿瘤样本的单细胞测序（10x），TNBC 4 个、HR+ 5 个。我把 Treg 细胞的原始 UMI 计数导出成了 cells.tsv（每行一个细胞，后面几列是基因），想比较 TNBC 和 HR+ 之间 Treg 的基因表达差异，比如 FOXP3、CTLA4、CCR8。请帮我做伪 bulk 聚合，告诉我现在的数据能不能做差异分析、应该怎么做，结果可以怎么写。

Input files: cells.tsv
