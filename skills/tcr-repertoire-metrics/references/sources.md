# Sources for tcr-repertoire-metrics

Every entry was checked against Crossref and, where indexed, PubMed in November 2026. The
indented line says what this skill takes from the source; claims marked "Methods" or "full
text" were read in the PubMed Central full text, the others in the abstract or the formula's
original definition.

## Diversity, evenness and clonality

- **[Shannon1948]** Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948;27:379-423. doi:10.1002/j.1538-7305.1948.tb01338.x.
  Entropy H = -sum p log p, the basis of the Shannon index.
- **[Pielou1966]** Pielou EC. The measurement of diversity in different types of biological collections. J Theor Biol. 1966;13:131-144. doi:10.1016/0022-5193(66)90013-0.
  Evenness J = H / ln S.
- **[Tumeh2014]** Tumeh PC, Harview CL, Yearley JH, et al. PD-1 blockade induces responses by inhibiting adaptive immune resistance. Nature. 2014;515(7528):568-571. doi:10.1038/nature13954. PMID: 25428505.
  Methods (PMC4246418): Shannon entropy on all productive TCR sequences, divided by the log of the number of unique productive sequences, then inverted (1 - normalised entropy) as "clonality"; the Extended Data Fig. 6 legend calls it 1 - Pielou's evenness.
- **[Coghill2025]** Coghill AE, Van Bibber N, Yoder S, Mokhtari S, Forsyth S, Blanck G, Egan KM. T-cell receptor sequencing for detection of Epstein-Barr and cytomegalovirus-specific immune responses in glioma patients: an exploratory study. Neurooncol Adv. 2025;7(1):vdaf216. doi:10.1093/noajnl/vdaf216. PMID: 41624423.
  Methods (PMC12859700): "productive Simpson clonality" from the immunoSEQ Analyzer, defined as the square root of Simpson's index over all productive rearrangements, from 0 (polyclonal) to 1 (one or a few dominant rearrangements). Used here only as an example of a different "clonality".
- **[Simpson1949]** Simpson EH. Measurement of diversity. Nature. 1949;163:688. doi:10.1038/163688a0.
  The Simpson concentration index (sum of squared proportions).
- **[Hill1973]** Hill MO. Diversity and evenness: a unifying notation and its consequences. Ecology. 1973;54:427-432. doi:10.2307/1934352.
  Abstract: species number, Shannon entropy and Simpson's index are members of one family of diversities that estimate an effective number of species and differ in how much weight they give rare species; in Hill's notation richness, exp(Shannon) and inverse Simpson are the orders 0, 1 and 2.

## Richness, rarefaction and sequencing depth

- **[Hurlbert1971]** Hurlbert SH. The nonconcept of species diversity: a critique and alternative parameters. Ecology. 1971;52(4):577-586. doi:10.2307/1934145. PMID: 28973811.
  Abstract: expected number of species in a random sample of n individuals, E(S_n) = sum_i [1 - C(N - N_i, n) / C(N, n)], and the probability of an interspecific encounter [N/(N - 1)][1 - sum (N_i/N)^2].
- **[Gotelli2001]** Gotelli NJ, Colwell RK. Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecol Lett. 2001;4:379-391. doi:10.1046/j.1461-0248.2001.00230.x.
  Abstract: richness depends on abundance and sampling effort; rarefaction allows standardised comparison; asymptotic estimators give lower bounds when richness does not reach an asymptote.
- **[Chiu2014]** Chiu CH, Wang YT, Walther BA, Chao A. An improved nonparametric lower bound of species richness via a modified Good-Turing frequency formula. Biometrics. 2014;70(3):671-682. doi:10.1111/biom.12200. PMID: 24945937.
  Abstract: Chao's (1984, Scand J Stat 11:265-270) estimator is a lower bound built only from singletons and doubletons. Chao 1984 itself has no Crossref or PubMed record and is not cited directly.
- **[Oksanen2026]** Oksanen J, et al. vegan: Community Ecology Package. Reference page "specpool / estimateR", documentation of version 2.8-0, accessed November 2026. https://vegandevs.github.io/vegan/reference/specpool.html
  estimateR uses the bias-corrected Chao1, S + f1(f1 - 1) / (2(f2 + 1)), without the small-sample factor (N - 1)/N, which the page says estimateR does not use; the R source of estimateR.default in 2.8-0 sets SSC <- 1 with the comment "do NOT use small-sample correction". The bias correction reduces but does not remove the downward bias. vegan itself was not run for this skill. Re-check if a later vegan version changes the formula.
- **[Laydon2015]** Laydon DJ, Bangham CR, Asquith B. Estimating T-cell repertoire diversity: limitations of classical estimators and a new approach. Philos Trans R Soc Lond B Biol Sci. 2015;370(1675):20140291. doi:10.1098/rstb.2014.0291. PMID: 26150657.
  Abstract and full text (PMC4528489): total TCR diversity cannot be captured in a blood sample; classical estimators, including bias-corrected Chao1, frequently underestimate TCR diversity, and all estimators tested rose systematically with sample size.
- **[Chaara2018]** Chaara W, Gonzalez-Tort A, Florez LM, Klatzmann D, Mariotti-Ferrandiz E, Six A. RepSeq data representativeness and robustness assessment by Shannon entropy. Front Immunol. 2018;9:1038. doi:10.3389/fimmu.2018.01038. PMID: 29868003.
  Abstract and full text (PMC5962720): in RNA aliquots equivalent to 500 to 50,000 sorted mouse CD4+ Foxp3- T cells, there were about 20-fold (500-cell) to 2-fold (50,000-cell) more clonotypes than cells, and in small samples the excess remained after removing singletons; Pielou evenness of unfiltered data fell strongly as the clonotype-to-cell ratio rose; over-sequencing small samples adds artefactual clonotypes; datasets were normalised by random rarefaction without replacement (vegan rrarefy).
- **[Shugay2015]** Shugay M, et al. VDJtools: unifying post-analysis of T cell receptor repertoires. PLoS Comput Biol. 2015;11(11):e1004503. doi:10.1371/journal.pcbi.1004503. PMID: 26606115.
  Full text (PMC4659587): diversity estimates are normalised by down-sampling datasets to the size of the smallest dataset (mean of three re-samples in that tool).

## Overlap

- **[Jaccard1912]** Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11:37-50. doi:10.1111/j.1469-8137.1912.tb05611.x.
  The coefficient of community: shared species over all species of two samples.
- **[Horn1966]** Horn HS. Measurement of "overlap" in comparative ecological studies. Am Nat. 1966;100:419-424. doi:10.1086/282436.
  The Morisita-Horn overlap index on relative abundances.
- **[Rempala2013]** Rempala GA, Seweryn M. Methods for diversity and overlap analysis in T-cell receptor populations. J Math Biol. 2013;67(6-7):1339-1368. doi:10.1007/s00285-012-0589-7. PMID: 23007599.
  Full text (PMC3543521): MH = 2 sum p1 p2 / (sum p1^2 + sum p2^2) for TCR populations; MH is overly sensitive to the high-abundance components (the most frequent clonotypes), so the Bhattacharyya coefficient suits populations dominated by rare clonotypes better.

## Data standards

- **[AIRR2026]** AIRR Community. AIRR Rearrangement schema, field definitions (documentation, stable release, accessed November 2026). https://docs.airr-community.org/en/stable/datarep/rearrangements.html
  Field names junction_aa, v_call, j_call, locus, productive, duplicate_count, consensus_count, umi_count, cell_id and clone_id. Re-check when the schema version changes.
