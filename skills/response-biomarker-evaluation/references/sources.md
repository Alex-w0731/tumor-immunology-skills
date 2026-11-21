# Sources

Every entry was checked against Crossref and PubMed / Europe PMC in November 2026. The
abstract of each source was read; where a number below comes from the full text, the note
says so. The indented line says what this skill takes from the source.

- **[Altman1994]** Altman DG, Lausen B, Sauerbrei W, Schumacher M. Dangers of using "optimal" cutpoints in the evaluation of prognostic factors. J Natl Cancer Inst. 1994;86(11):829-835. doi:10.1093/jnci/86.11.829. PMID: 8182763.
  No abstract in PubMed; cited only for the warning in its title against data-derived "optimal" cut-points (tmelens cites it for offering median or tertile splits only). The numbers on cut-point optimism come from [Leeflang2008].
- **[Ambroise2002]** Ambroise C, McLachlan GJ. Selection bias in gene extraction on the basis of microarray gene-expression data. Proc Natl Acad Sci U S A. 2002;99(10):6562-6566. doi:10.1073/pnas.102102699. PMID: 11983868.
  Gene selection done outside the cross-validation loop gives near-zero apparent error; selection must be repeated inside each training split; 10-fold preferred to leave-one-out.
- **[Ballman2015]** Ballman KV. Biomarker: predictive or prognostic? J Clin Oncol. 2015;33(33):3968-3971. doi:10.1200/JCO.2015.63.3651. PMID: 26392104.
  A predictive claim needs marker data on treated and untreated patients, preferably randomized, and a significant treatment-by-biomarker interaction.
- **[Benjamini1995]** Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
  FDR control across many screened markers (tmelens applies it across signatures).
- **[Bossuyt2015]** Bossuyt PM, Reitsma JB, Bruns DE, et al. STARD 2015: an updated list of essential items for reporting diagnostic accuracy studies. BMJ. 2015;351:h5527. doi:10.1136/bmj.h5527. PMID: 26511519.
  Reporting checklist when a marker or threshold is evaluated as a classifier of response; item 19 asks for the flow of participants, using a diagram (full text).
- **[Clopper1934]** Clopper CJ, Pearson ES. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika. 1934;26(4):404-413. doi:10.1093/biomet/26.4.404.
  Exact binomial interval used for sensitivity, specificity, PPV and NPV.
- **[Collins2015]** Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. BMJ. 2015;350:g7594. doi:10.1136/bmj.g7594. PMID: 25569120.
  Reporting checklist for developing or validating a multivariable prediction model.
- **[Collins2024]** Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ. 2024;385:e078378. doi:10.1136/bmj-2023-078378. PMID: 38626948.
  Updated TRIPOD guidance covering regression and machine-learning models.
- **[DeLong1988]** DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics. 1988;44(3):837-845. doi:10.2307/2531595. PMID: 3203132.
  Covariance of AUCs measured on the same patients (generalized U-statistics); standard error and paired test in roc_auc.py.
- **[Demler2012]** Demler OV, Pencina MJ, D'Agostino RB. Misuse of DeLong test to compare AUCs for nested models. Stat Med. 2012;31(23):2577-2587. doi:10.1002/sim.5328. PMID: 22415937.
  For nested models developed and tested on the same data, the DeLong test's normal null does not hold (the test becomes conservative); test the added predictor's association instead and, if significant, estimate the AUC change with a CI (abstract).
- **[Efron1986]** Efron B, Tibshirani R. Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy. Stat Sci. 1986;1(1):54-75. doi:10.1214/ss/1177013815.
  Bootstrap percentile interval.
- **[Eisenhauer2009]** Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). Eur J Cancer. 2009;45(2):228-247. doi:10.1016/j.ejca.2008.10.026. PMID: 19097774.
  RECIST 1.1 response categories; confirmation of response required when response is the primary endpoint, not in randomized trials; progression needs a 20% and a 5 mm absolute increase in the sum (abstract).
- **[Hanley1982]** Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology. 1982;143(1):29-36. doi:10.1148/radiology.143.1.7063747. PMID: 7063747.
  AUC equals the probability that a randomly chosen case ranks above a randomly chosen non-case and equals the Wilcoxon (Mann-Whitney) statistic; closed-form standard error used for planning.
- **[Janes2008]** Janes H, Pepe MS. Adjusting for covariates in studies of diagnostic, screening, or prognostic markers: an old concept in a new setting. Am J Epidemiol. 2008;168(1):89-97. doi:10.1093/aje/kwn099. PMID: 18477651.
  Studies of classification accuracy need covariate adjustment just as association studies do (abstract). This skill applies the principle by treating cohort as the covariate (per-cohort and within-cohort AUCs); its pair-weighted within-cohort summary is the skill's own choice, not an estimator from this paper.
- **[Leeflang2008]** Leeflang MM, Moons KG, Reitsma JB, Zwinderman AH. Bias in sensitivity and specificity caused by data-driven selection of optimal cutoff values: mechanisms, magnitude, and solutions. Clin Chem. 2008;54(4):729-737. doi:10.1373/clinchem.2007.096032. PMID: 18258670.
  Youden-optimized cut-offs overstate accuracy, more so in small studies: true sensitivity and specificity of 84% are estimated near 90% with n = 40 and 86% with n = 200 (abstract).
- **[Lijmer1999]** Lijmer JG, Mol BW, Heisterkamp S, et al. Empirical evidence of design-related bias in studies of diagnostic tests. JAMA. 1999;282(11):1061-1066. doi:10.1001/jama.282.11.1061. PMID: 10493205.
  Studies that sample diseased patients and a separate control group overestimate diagnostic performance (relative diagnostic odds ratio 3.0, 95% CI 2.0-4.5; abstract). The paper did not study response labels; this skill applies it by analogy to analyses that drop stable disease.
- **[McShane2005]** McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). J Natl Cancer Inst. 2005;97(16):1180-1184. doi:10.1093/jnci/dji237. PMID: 16106022.
  Reporting checklist for tumour marker prognostic studies.
- **[Obuchowski1998]** Obuchowski NA, Lieber ML. Confidence intervals for the receiver operating characteristic area in studies with small samples. Acad Radiol. 1998;5(8):561-571. doi:10.1016/S1076-6332(98)80208-0. PMID: 9702267.
  For a single AUC, asymptotic intervals under-cover in small samples (high-accuracy AUCs need more than 200 patients); asymptotic intervals for the difference of two AUCs were adequate even with 20 patients; no single best small-sample method (abstract).
- **[Phipson2010]** Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Stat Appl Genet Mol Biol. 2010;9:Article39. doi:10.2202/1544-6115.1585. PMID: 21044043.
  Permutation p-values should never be zero: the naive b/m understates p by about 1/m (abstract). roc_auc.py therefore reports the conservative (b + 1) / (m + 1).
- **[Rizvi2015]** Rizvi NA, Hellmann MD, Snyder A, et al. Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer. Science. 2015;348(6230):124-128. doi:10.1126/science.aaa1348. PMID: 25765070.
  Example of a clinical-benefit label: durable clinical benefit = partial or stable response lasting more than 6 months, with patients followed for less than 6 months kept apart (full text).
- **[Seymour2017]** Seymour L, Bogaerts J, Perrone A, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. Lancet Oncol. 2017;18(3):e143-e152. doi:10.1016/S1470-2045(17)30074-8. PMID: 28271869.
  iUPD must be confirmed (iCPD) at the next assessment 4-8 weeks later; iUPD does not override a later iSD, iPR or iCR in the best overall response (full text).
- **[Varma2006]** Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics. 2006;7:91. doi:10.1186/1471-2105-7-91. PMID: 16504092.
  Tuning a classifier by cross-validation and reporting that same CV error is optimistic: on data with no class difference the estimated error was below 30% in 18.5% (shrunken centroids) and 38% (SVM) of datasets; nested cross-validation with every step repeated in each loop is nearly unbiased (abstract).
