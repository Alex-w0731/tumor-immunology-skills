# Sources for spatial-neighborhood-statistics

Every entry was checked against Crossref (title, journal, year, first author) and, where a PMID
exists, PubMed or Europe PMC, in November 2026. The indented line says what this skill takes from
the source and where that was read (abstract or full text). Statements labelled `convention`,
`heuristic`, `derivation` or "this skill's reasoning" are not attributed to any source.

## Neighbourhood permutation tests

- **[Schapiro2017]** Schapiro D, Jackson HW, Raghuraman S, et al. histoCAT: analysis of cell phenotypes and interactions in multiplex image cytometry data. Nat Methods. 2017;14(9):873-876. doi:10.1038/nmeth.4391. PMID: 28783155.
  Full text (PMC5617107): neighbours are found by step-wise pixel expansion of each segmented cell (1 to 6 pixels in most cases; 6 pixels for the pairwise analysis), a neighbour being a cell within the selected pixel distance; interactions in each image are compared with cell-type-randomized controls of the same image, which control for that tissue's connectivity and cell types; one one-tailed test each for interaction and avoidance, significance at p < 0.01; avoidance among rare cell populations will rarely be detected; Equation 1 as printed divides the count of extreme permutations by the number of permutations plus one, without a +1 in the numerator; the cohort display showed a connection only if it was significant in at least 30% of the grouped samples and both phenotypes were present in at least 90%.
- **[Palla2022]** Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. Nat Methods. 2022;19(2):171-178. doi:10.1038/s41592-021-01358-2. PMID: 35102346.
  Full text (PMC8828470): neighbourhood enrichment counts label pairs that are adjacent in the connectivity graph, scrambles the labels while keeping the connectivities (1,000 times by default), and reports z = (observed - permutation mean) / permutation SD, an approach it attributes to Schapiro et al.; graphs can be built from a fixed number of adjacent spots on a grid (hexagonal rings for Visium, n_neigh = 6), k nearest neighbours, a fixed radius in the units of the spatial coordinates, or Delaunay triangulation; it also implements Ripley's L, F and G.
- **[Keren2018]** Keren L, Bosse M, Marquez D, et al. A structured tumor-immune microenvironment in triple negative breast cancer revealed by multiplexed ion beam imaging. Cell. 2018;174(6):1373-1387.e19. doi:10.1016/j.cell.2018.08.039. PMID: 30193111.
  Full text (author manuscript, PMC6132072): 41 TNBC patients, one 800 x 800 um region per patient imaged as 2048 x 2048 pixels "at 500nm resolution", while 100 pixels are stated as 39 um (800/2048 = 0.39 um per pixel; the methods list 500 nm as the probe size); X-positive cells within 100 pixels of Y-positive cells were counted, the locations of Y-positive cells were randomized 500 times with their number kept constant, and a z-score was computed; results were robust to the distance (100 to 500 pixels) and the number of randomizations; a context-dependent version randomized only within the tumour compartment or within lineage categories; a mixing score (immune-tumour interactions divided by immune-immune interactions among cells in direct contact) separated cold, mixed and compartmentalized tumours.

## Null hypotheses, density and edges

- **[Goreaud2003]** Goreaud F, Pélissier R. Avoiding misinterpretation of biotic interactions with the intertype K12-function: population independence vs. random labelling hypotheses. J Veg Sci. 2003;14(5):681-692. doi:10.1111/j.1654-1103.2003.tb02200.x.
  Abstract: the intertype K12 (L12) function is used to study interactions through the spatial structure of bivariate patterns; independence and random labelling are different null hypotheses with different confidence intervals, and choosing the inappropriate one risks misinterpreting interactions, sometimes with opposite conclusions.
- **[Baddeley2000]** Baddeley AJ, Møller J, Waagepetersen R. Non- and semi-parametric estimation of interaction in inhomogeneous point patterns. Stat Neerl. 2000;54(3):329-350. doi:10.1111/1467-9574.00144.
  Abstract (Crossref): methods, including an analogue of the K-function, for analysing interaction between points when the pattern is spatially inhomogeneous.
- **[Wilson2022]** Wilson C, Soupir AC, Thapa R, et al. Tumor immune cell clustering and its association with survival in African American women with ovarian cancer. PLoS Comput Biol. 2022;18(3):e1009900. doi:10.1371/journal.pcbi.1009900. PMID: 35235563.
  Full text (PMC8920290): the edge effect is the loss of neighbours that lie outside the sampled region for cells at the periphery (isotropic and translation corrections are common); TMA cores have more holes and uneven cell density than ROIs; the theoretical value of K under complete spatial randomness is biased, increasingly so with more missing cells or holes, so a permutation-based estimate over the observed cells is used instead; the approach allows analysis within the tumour or stroma compartment separately, and the TMA analyses were restricted to the tumour compartment of each core.

## Replication across samples and patients

- **[Hurlbert1984]** Hurlbert SH. Pseudoreplication and the design of ecological field experiments. Ecol Monogr. 1984;54(2):187-211. doi:10.2307/1942661.
  Definition (abstract): inferential tests of treatment effects in which treatments are not replicated or the replicates are not statistically independent.
- **[Wilson2021]** Wilson CM, Ospina OE, Townsend MK, et al. Challenges and opportunities in the statistical analysis of multiplex immunofluorescence data. Cancers (Basel). 2021;13(12):3031. doi:10.3390/cancers13123031. PMID: 34204319.
  Abstract and full text (PMC8233801): multiple TMA cores or ROIs from the same tumour are repeated measurements, usually handled with mixed or random-effect models; using each subject's mean instead may underestimate variability (inflating type I error) unless the repeated values are very consistent.
- **[Canete2022]** Canete NP, Iyengar SS, Ormerod JT, et al. spicyR: spatial analysis of in situ cytometry data in R. Bioinformatics. 2022;38(11):3099-3105. doi:10.1093/bioinformatics/btac268. PMID: 35438129.
  Abstract and full text (PMC9326848): inference on changes in cell-type co-localization across groups of images; one co-localization score per image, a linear model, or a linear mixed-effects model with a random intercept per subject when subjects contribute several images, with weights that account for images with very different cell numbers; applied to imaging mass cytometry with patient as the random effect.
- **[Bost2023]** Bost P, Schulz D, Engler S, et al. Optimizing multiplexed imaging experimental design through tissue spatial segregation estimation. Nat Methods. 2023;20(3):418-423. doi:10.1038/s41592-022-01692-z. PMID: 36585456.
  Abstract: the number and area of fields of view needed to capture all cell phenotypes of a tissue depend on a measure of tissue spatial segregation.

## Test statistics and multiple testing

- **[Phipson2010]** Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Stat Appl Genet Mol Biol. 2010;9:Article39. doi:10.2202/1544-6115.1585. PMID: 21044043.
  Abstract: permutation p-values are often understated by about 1/m (m permutations), which matters under multiple testing. The skill uses (b + 1) / (m + 1).
- **[Dixon1946]** Dixon WJ, Mood AM. The statistical sign test. J Am Stat Assoc. 1946;41(236):557-566. doi:10.1080/01621459.1946.10501898. PMID: 20279351.
  The exact binomial sign test used for the across-patient direction (metadata; the test is the paper's subject).
- **[Benjamini1995]** Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
  The step-up FDR procedure applied across cell-type pairs (the procedure is the paper's subject).
- **[Benjamini2001]** Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency. Ann Stat. 2001;29(4):1165-1188. doi:10.1214/aos/1013699998.
  FDR control under arbitrary dependence (the conservative variant offered for strongly dependent pair counts).
