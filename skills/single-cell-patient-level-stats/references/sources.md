# Sources for single-cell-patient-level-stats

Every entry was checked against Crossref (title, journal, year, first author) and, where a PMID
exists, PubMed, in November 2026. The indented line says what this skill takes from the source;
specific statements were read in the abstract or the open-access full text, as noted.

## Pseudoreplication and pseudobulk differential expression

- **[Hurlbert1984]** Hurlbert SH. Pseudoreplication and the design of ecological field experiments. Ecol Monogr. 1984;54(2):187-211. doi:10.2307/1942661.
  Definition (abstract): inferential tests of treatment effects in which treatments are not replicated or the replicates are not statistically independent.
- **[Zimmerman2021]** Zimmerman KD, Espeland MA, Langefeld CD. A practical solution to pseudoreplication bias in single-cell studies. Nat Commun. 2021;12(1):738. doi:10.1038/s41467-021-21038-1. PMID: 33531494.
  Abstract: cells from one individual are pseudoreplicates; ignoring this inflates type 1 error, including for methods that use a batch effect correction for individual; pseudobulk was conservative and underpowered relative to GLMMs with a random effect for individual. Full text: the batch correction tested was ComBat with individuals as batches, applied before cell-level DE, and it markedly increased type 1 error; hurdle and Tweedie mixed models controlled type 1 error only with the random effect.
- **[Squair2021]** Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. Nat Commun. 2021;12(1):5692. doi:10.1038/s41467-021-25960-2. PMID: 34584091.
  Abstract: methods that ignore variation between biological replicates are biased and the most widely used ones can report hundreds of DE genes without biological differences. Full text: pseudobulk methods outperformed single-cell methods; single-cell methods are biased towards highly expressed genes.
- **[Lun2017]** Lun ATL, Marioni JC. Overcoming confounding plate effects in differential expression analyses of single-cell RNA-seq data. Biostatistics. 2017;18(3):451-464. doi:10.1093/biostatistics/kxw055. PMID: 28334062.
  Abstract: ignoring plate effects loses type 1 error control; summing counts per plate and testing the sums restores it without losing power in simulations.
- **[Crowell2020]** Crowell HL, Soneson C, Germain PL, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nat Commun. 2020;11(1):6077. doi:10.1038/s41467-020-19894-4. PMID: 33257685.
  Full text: compared cell-level mixed models with pseudobulk (sum or mean; edgeR, limma-voom, limma-trend); aggregation methods were extremely fast and "perform amongst the best", and mixed models performed similarly well; for aggregation methods about 100 cells per subpopulation per sample gave decent performance, with a sizable gain from 20 to 100 cells and only a moderate gain at 200 or 400; the LPS case study used edgeR on summed counts, testing a subpopulation only if at least 2 samples per group had at least 10 cells. The muscat pbDS() default min_cells = 10 (minimum cells for a cluster-sample to enter testing) was read in the package source (https://github.com/HelenaLC/muscat, devel branch, version 1.27.4, November 2026; re-check newer releases).
- **[Murphy2022]** Murphy AE, Skene NG. A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis. Nat Commun. 2022;13(1):7851. doi:10.1038/s41467-022-35519-4. PMID: 36550119.
  Full text (Matters Arising, no abstract): with the Matthews correlation coefficient, pseudobulk approaches performed best; the simulation software did not normalise before the pseudobulk methods, and summed pseudobulks degraded dramatically with imbalanced cell numbers. It names no specific normalisation method; TMM or DESeq2 size factors are the skill's convention.
- **[Zimmerman2022]** Zimmerman KD, Evans C, Langefeld CD. Reply to: A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis. Nat Commun. 2022;13(1):7852. doi:10.1038/s41467-022-35520-x. PMID: 36550147.
  Full text (reply, no abstract): argues that type 1 error must be controlled before power is compared, which a combined measure such as MCC does not do.
- **[Killip2004]** Killip S, Mahfoud Z, Pearce K. What is an intracluster correlation coefficient? Crucial concepts for primary care researchers. Ann Fam Med. 2004;2(3):204-208. doi:10.1370/afm.141. PMID: 15209195.
  Abstract: the effective sample size depends on the intracluster correlation, cluster size and number of clusters; adding clusters raises power more efficiently than adding subjects per cluster.

## Design, models and sample size

- **[McCarthy2012]** McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Res. 2012;40(10):4288-4297. doi:10.1093/nar/gks042. PMID: 22287627.
  Abstract: GLMs with blocking variables in edgeR; case study of a paired tumour-normal design.
- **[Hoffman2021]** Hoffman GE, Roussos P. Dream: powerful differential expression analysis for repeated measures designs. Bioinformatics. 2021;37(2):192-201. doi:10.1093/bioinformatics/btaa687. PMID: 32730587.
  Abstract: linear mixed models for designs with several samples per individual, controlling false positives that standard methods produce.
- **[Schurch2016]** Schurch NJ, Schofield P, Gierliński M, et al. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? RNA. 2016;22(6):839-851. doi:10.1261/rna.053959.115. PMID: 27022035.
  Abstract: with 3 replicates, 9 of 11 tools found only 20-40% of the DE genes found with 42 clean replicates; finding more than 85% of all DE genes regardless of fold change needed more than 20 replicates; recommends at least 6 replicates, rising to at least 12 when it is important to identify DE genes for all fold changes (Saccharomyces cerevisiae, 48 replicates per condition). The skill uses it as a bulk analogy for pseudobulk; that patients need more is the skill's reasoning.

## Cell-type composition

- **[Buttner2021]** Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B. scCODA is a Bayesian model for compositional single-cell data analysis. Nat Commun. 2021;12(1):6876. doi:10.1038/s41467-021-27150-6. PMID: 34824236.
  Full text: testing each cell type independently can call shifts that are induced only by the negative correlation of proportions; Poisson-regression and Wilcoxon approaches do not address compositionality; non-compositional tests gave more false positives as effect size and replicate number grew; scCODA results are relative to a reference cell type.
- **[Phipson2022]** Phipson B, Sim CB, Porrello ER, Hewitt AW, Powell J, Oshlack A. propeller: testing for differences in cell type proportions in single cell data. Bioinformatics. 2022;38(20):4720-4726. doi:10.1093/bioinformatics/btac582. PMID: 36005887.
  Full text: arcsine square root or logit transform (logit adds an offset of 0.5 to the raw cell type counts), then empirical Bayes linear models that borrow information across cell types; uses biological replicates. In its simulations the arcsine transform controlled type 1 error best at 10 and 20 samples per group, the logit transform better at 3 and 5. The skill borrows the 0.5 offset for its own log-ratio check; that use is not from the paper.
- **[Denisenko2020]** Denisenko E, Guo BB, Jones M, et al. Systematic assessment of tissue dissociation and storage biases in single-cell and single-nucleus RNA-seq workflows. Genome Biol. 2020;21(1):130. doi:10.1186/s13059-020-02048-6. PMID: 32487174.
  Abstract: dissociation temperature, cryopreservation and single-nucleus protocols change the recovered cell-type composition (mouse kidney).

## Tests, effect sizes and multiplicity

- **[Wilcoxon1945]** Wilcoxon F. Individual comparisons by ranking methods. Biometrics Bulletin. 1945;1(6):80-83. doi:10.2307/3001968.
  The rank-sum and signed-rank tests whose exact null distributions composition_test.py enumerates.
- **[Mann1947]** Mann HB, Whitney DR. On a test of whether one of two random variables is stochastically larger than the other. Ann Math Stat. 1947;18(1):50-60. doi:10.1214/aoms/1177730491.
  The U statistic for two independent samples (reported as U_case).
- **[Hodges1963]** Hodges JL, Lehmann EL. Estimates of location based on rank tests. Ann Math Stat. 1963;34(2):598-611. doi:10.1214/aoms/1177704172.
  Location estimates derived from rank tests: median of pairwise differences (two samples) and of Walsh averages (paired). Metadata checked; the full text was not reachable online, so the formulas were checked against hand calculations in the tests.
- **[Kerby2014]** Kerby DS. The simple difference formula: an approach to teaching nonparametric correlation. Comprehensive Psychology. 2014;3:11.IT.3.1. doi:10.2466/11.IT.3.1.
  Rank-biserial correlation as the proportion of favourable minus unfavourable evidence, for the Mann-Whitney and signed-rank tests. Metadata checked; full text not reachable online, formulas checked by hand calculation in the tests.
- **[Benjamini1995]** Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
  The step-up FDR procedure used across cell types.
- **[Benjamini2001]** Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency. Ann Stat. 2001;29(4):1165-1188. doi:10.1214/aos/1013699998.
  FDR control under dependence; the conservative variant for arbitrary dependence.
