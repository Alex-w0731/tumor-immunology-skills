# Sources for immune-survival-analysis

Every entry was checked against Crossref and PubMed (or OpenAlex where no DOI exists) in
November 2026. The indented line says what the skill takes from the source and how much of
it was read: "abstract" means the claim was checked against the abstract, "full text" against
the article, "metadata" means only the bibliographic record was confirmed.

- **[Altman1994]** Altman DG, Lausen B, Sauerbrei W, Schumacher M. Dangers of using "optimal" cutpoints in the evaluation of prognostic factors. J Natl Cancer Inst. 1994;86(11):829-835. doi:10.1093/jnci/86.11.829. PMID: 8182763.
  The standard warning that choosing the cut-point with the smallest p-value inflates false-positive findings. Metadata only (PubMed has no abstract), so the skill quotes no number from it; cutpoint_inflation.py produces the numbers.
- **[Altman2012]** Altman DG, McShane LM, Sauerbrei W, Taube SE. Reporting Recommendations for Tumor Marker Prognostic Studies (REMARK): explanation and elaboration. PLoS Med. 2012;9(5):e1001216. doi:10.1371/journal.pmed.1001216. PMID: 22675273.
  Full text (PMC3362085): wording of the 20-item REMARK checklist; the items the skill cites (1, 7, 8, 11, 12 and 15 to 18) were checked against the checklist table.
- **[Anderson1983]** Anderson JR, Cain KC, Gelber RD. Analysis of survival by tumor response. J Clin Oncol. 1983;1(11):710-719. doi:10.1200/JCO.1983.1.11.710. PMID: 6668489.
  Abstract: comparing the survival of responders and non-responders is biased in favour of responders, should not be used to infer treatment effect, and valid alternatives exist.
- **[Aran2015]** Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. Nat Commun. 2015;6:8971. doi:10.1038/ncomms9971. PMID: 26634437.
  Abstract: tumour purity in more than 10,000 TCGA samples confounds correlation and clustering of tumours from transcriptomic data.
- **[Benjamini1995]** Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
  The Benjamini-Hochberg step-up procedure implemented in fdr_adjust.py (metadata; the procedure is the paper's subject).
- **[GiobbieHurder2013]** Giobbie-Hurder A, Gelber RD, Regan MM. Challenges of guarantee-time bias. J Clin Oncol. 2013;31(23):2963-2969. doi:10.1200/JCO.2013.49.5283. PMID: 23835712.
  Abstract: guarantee-time (immortal-time) bias arises when groups are defined by an event during follow-up such as response or toxicity; conditional landmark analysis, the extended (time-dependent) Cox model and inverse probability weighting remove it.
- **[Grambsch1994]** Grambsch PM, Therneau TM. Proportional hazards tests and diagnostics based on weighted residuals. Biometrika. 1994;81(3):515-526. doi:10.1093/biomet/81.3.515.
  The residual-based test and plot of the proportional-hazards assumption (R cox.zph). Metadata; also cited for this purpose in [Liu2018].
- **[Holm1979]** Holm S. A simple sequentially rejective multiple test procedure. Scand J Stat. 1979;6:65-70. https://www.jstor.org/stable/4615733.
  The Holm step-down procedure implemented in fdr_adjust.py. No resolvable DOI exists: the record was confirmed in OpenAlex (W2121044470), and the JSTOR DOI that OpenAlex lists returns nothing at Crossref or doi.org (checked November 2026).
- **[Kaplan1958]** Kaplan EL, Meier P. Nonparametric estimation from incomplete observations. J Am Stat Assoc. 1958;53(282):457-481. doi:10.1080/01621459.1958.10501452.
  The product-limit estimator implemented in km_logrank.py (metadata).
- **[Lausen1992]** Lausen B, Schumacher M. Maximally selected rank statistics. Biometrics. 1992;48(1):73. doi:10.2307/2532740.
  Named, by its title, as the analytic alternative to the permutation test for a searched cut-point (implemented in the R package maxstat). Metadata only: the skill takes no number or formula from it.
- **[Liu2018]** Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. Cell. 2018;173(2):400-416.e11. doi:10.1016/j.cell.2018.02.052. PMID: 29625055.
  Full text (PMC6066282), re-read in review: Table 3 endpoint recommendations and event counts per cancer type (all 33 rows of tcga-endpoints.md match it); PFI recommended without reservation for 29 of 33 types in the Results and Table 3 (the Discussion says 27); the text, unlike Table 3, does not recommend OS for PRAD (10 OS events in 500 cases); PFI and DFI record more events than OS and DSS given the short follow-up; DSS mostly approximated; DFI excludes de novo stage IV cases and follow-up under 90 days; stage III-IV versus I-II HR significantly above 1 in most of 14 types, not in MESO, PAAD and UVM; outcomes compared between the largest tissue source sites and the others, with the advice that site may need to be taken into account; the barcode begins with tissue source site and patient.
- **[Mantel1966]** Mantel N. Evaluation of survival data and two new rank order statistics arising in its consideration. Cancer Chemother Rep. 1966;50(3):163-170. PMID: 5910392.
  Origin of the log-rank test (metadata).
- **[McShane2005]** McShane LM, Altman DG, Sauerbrei W, Taube SE, Gion M, Clark GM. Reporting recommendations for tumor marker prognostic studies (REMARK). J Natl Cancer Inst. 2005;97(16):1180-1184. doi:10.1093/jnci/dji237. PMID: 16106022.
  Abstract: the REMARK reporting guideline for prognostic marker studies.
- **[Peduzzi1995]** Peduzzi P, Concato J, Feinstein AR, Holford TR. Importance of events per independent variable in proportional hazards regression analysis. II. Accuracy and precision of regression estimates. J Clin Epidemiol. 1995;48(12):1503-1510. doi:10.1016/0895-4356(95)00048-8. PMID: 8543964.
  Abstract: Monte Carlo simulations at 2 to 25 events per variable (EPV); as EPV fell, coefficients became more biased, 90% confidence limits lost coverage and Z (Wald) statistics lost validity under the null; EPV = 10 "seems most prudent". The paper counts variables; counting parameters is the skill's convention.
- **[Peto1972]** Peto R, Peto J. Asymptotically efficient rank invariant test procedures. J R Stat Soc Ser A. 1972;135(2):185. doi:10.2307/2344317.
  Log-rank scores and test (metadata).
- **[Royston2006]** Royston P, Altman DG, Sauerbrei W. Dichotomizing continuous predictors in multiple regression: a bad idea. Stat Med. 2006;25(1):127-141. doi:10.1002/sim.2331. PMID: 16217841.
  Abstract: dichotomising a continuous predictor loses considerable power and leaves residual confounding; a data-derived optimal cut-point causes serious bias.
- **[Schemper1996]** Schemper M, Smith TL. A note on quantifying follow-up in studies of failure time. Control Clin Trials. 1996;17(4):343-346. doi:10.1016/0197-2456(96)00075-x. PMID: 8889347.
  The reverse Kaplan-Meier estimate of median follow-up computed by km_logrank.py (metadata only; PubMed has no abstract). Its DOI is left out for the same reason as for [Peduzzi1995].
- **[Vittinghoff2007]** Vittinghoff E, McCulloch CE. Relaxing the rule of ten events per variable in logistic and Cox regression. Am J Epidemiol. 2007;165(6):710-718. doi:10.1093/aje/kwk052. PMID: 17182981.
  Abstract: coverage and bias were acceptable below 10 EPV in a range of circumstances; the rule can be relaxed, in particular for sensitivity analyses that show control of confounding.
