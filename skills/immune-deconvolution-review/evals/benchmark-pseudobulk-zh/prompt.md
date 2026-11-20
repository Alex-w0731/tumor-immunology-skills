---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

我们用课题组自己的肺腺癌单细胞数据(约6万个细胞)做了12个伪bulk混合样本(PB01-PB12):每个样本按已知比例抽取细胞后把表达量相加,样本里也包含肿瘤细胞。然后用同一批细胞按细胞类型取平均,做成参考矩阵,再用一个回归去卷积方法估计了5种免疫细胞的比例。真实比例和估计比例分别放在两个文件里。

我算了每种细胞类型真实值和估计值的Pearson相关系数,全部在0.9以上,所以觉得这个方法很准。下一步想直接用它分析TCGA-LUAD的bulk RNA-seq,并在文章里写"该方法经伪bulk验证,估计准确"。请帮我评估这个验证结果是否足够,还需要补充什么?

Input files: truth_fractions_synthetic.tsv, estimated_fractions_synthetic.tsv
