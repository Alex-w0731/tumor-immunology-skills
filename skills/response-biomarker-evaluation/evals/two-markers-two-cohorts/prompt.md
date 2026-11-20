---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

我们有两个抗PD-1治疗队列（A队列40例，B队列36例，RNA-seq分别在两个中心测的），每个病人都有T细胞炎症评分 tis_score 和 CD274 表达 cd274_log2，疗效已经分成 responder / non-responder。想知道哪个标志物预测疗效更好。我的打算是把两个队列合并，分别算两个标志物的AUC，看哪个高；如果合并后 tis_score 明显更好，就在文章里写“tis_score 优于 PD-L1 mRNA”。请帮我分析并给出结论。

Input files: synthetic_two_cohorts.tsv
