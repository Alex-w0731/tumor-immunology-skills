---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

我们有一位黑色素瘤患者抗PD-1治疗前（pre）和治疗后（post）肿瘤组织的单细胞V(D)J数据，导出成了contig表（每行一个contig，包含barcode、sample、chain、cdr3、productive、umis这几列）。请帮我计算治疗前后TCR库的克隆性（clonality）和Shannon多样性，看看治疗后有没有出现克隆扩增。结果请用中文写。

Input files: synthetic_scvdj_contigs.csv
