# Sources

Every key cited in `SKILL.md`, `method-notes.md` and `validation-design.md`. Each DOI and PMID was
checked against Crossref and PubMed or Europe PMC in November 2026, and the statement the skill
takes from each source was checked against its abstract or full text.

## Methods

- **[Newman2015]** Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods. 2015;12(5):453-457. doi:10.1038/nmeth.3337. PMID: 25822800.
  CIBERSORT; LM22 (547 genes, 22 human hematopoietic phenotypes, Supplementary Table 1) chosen by the lowest 2-norm condition number (11.4, computed with R's `kappa`); works best in non-log linear space; an empirically defined global P value per deconvolution and no P values for individual cell types; CIBERSORT was not explicitly tested on RNA-seq data (full text, Online Methods and Discussion).
- **[Newman2019]** Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol. 2019;37(7):773-782. doi:10.1038/s41587-019-0114-2. PMID: 31061481.
  CIBERSORTx infers cell-type-specific expression and minimises platform-specific variation so that single-cell references can be used on bulk data.
- **[Becht2016]** Becht E, Giraldo NA, Lacroix L, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biol. 2016;17(1):218. doi:10.1186/s13059-016-1070-5. PMID: 27765066.
  MCP-counter was designed to compare the abundance of one population across samples; scores are arbitrary units that depend on the expression platform, so only samples from the same platform can be compared; conceptually close to IHC densities; validated with in vitro mRNA mixtures and IHC.
- **[Aran2017]** Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol. 2017;18(1):220. doi:10.1186/s13059-017-1349-1. PMID: 29141660.
  64 cell types from 1822 pure transcriptomes; raw scores (averaged ssGSEA) are not linear in abundance and do not allow comparison across cell types; the transformed scores resembled known fractions in simulations, which the authors say enables comparison across cell types, yet their Discussion states the scores are enrichment scores that cannot be reliably used as proportions.
- **[Racle2017]** Racle J, de Jonge K, Baumgaertner P, Speiser DE, Gfeller D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife. 2017;6:e26476. doi:10.7554/eLife.26476. PMID: 29130882.
  EPIC: TPM input, reference profiles from circulating and tumour-infiltrating cells, mRNA-content renormalisation, an uncharacterised (cancer) cell fraction.
- **[Finotello2019]** Finotello F, Mayer C, Plattner C, et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome Med. 2019;11(1):34. doi:10.1186/s13073-019-0638-6. PMID: 31126321.
  quanTIseq: ten immune cell types, TPM input, removal of signature genes highly expressed in cancer cell lines, mRNA-content scaling, an "other" fraction that is not tumour purity.
- **[Monaco2019]** Monaco G, Lee B, Xu W, et al. RNA-Seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. Cell Rep. 2019;26(6):1627-1640.e7. doi:10.1016/j.celrep.2019.01.041. PMID: 30726743.
  Cell types differ in mRNA abundance; signatures normalised by it allow absolute deconvolution (PBMC).
- **[Abbas2009]** Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z, Clark HF. Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus. PLoS One. 2009;4(7):e6098. doi:10.1371/journal.pone.0006098. PMID: 19568420.
  Goodness of fit tracked the condition number of the basis matrix, and basis matrices were built partly by minimising it.
- **[Lawson1995]** Lawson CL, Hanson RJ. Solving Least Squares Problems. Philadelphia: Society for Industrial and Applied Mathematics; 1995. ISBN 978-0-89871-356-5. doi:10.1137/1.9781611971217.
  The NNLS active-set algorithm implemented in `scripts/nnls_deconvolve.py`.

## Benchmarks

- **[Sturm2019]** Sturm G, Finotello F, Petitprez F, et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics. 2019;35(14):i436-i445. doi:10.1093/bioinformatics/btz363. PMID: 31510660.
  Table 1 classes each method's output (inter-sample, intra-sample or both; CIBERSORT relative mode intra-sample only, TIMER not across cancer types); deconvolution-based methods mostly predicted absent immune cells at a low background level; minimal detection fraction; spillover between CD8+ and CD4+ T cells, NK and CD8+ T cells, and from DCs to B cells for all methods in simulations, the first two confirmed in validation data.
- **[AvilaCobos2020]** Avila Cobos F, Alquicira-Hernandez J, Powell JE, Mestdagh P, De Preter K. Benchmarking of cell type deconvolution pipelines for transcriptomics data. Nat Commun. 2020;11(1):5650. doi:10.1038/s41467-020-19015-1. PMID: 33159064.
  Linear scale best (log or VST: two- to four-fold higher median RMSE); quantile normalisation always sub-optimal for bulk methods; advises against row scaling, column min-max and column z-score normalisation, each applied to mixtures and reference separately; TPM gave EPIC its best results; omitting a cell type present in the mixture makes results substantially worse; 50:50 training and test split of cells, and keeping each individual's cells in one split gave slightly higher RMSE with the same conclusions.

## Tumour purity and compositional data

- **[Aran2015]** Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. Nat Commun. 2015;6:8971. doi:10.1038/ncomms9971. PMID: 26634437.
  Tumour purity confounds correlation and clustering of tumour transcriptomes.
- **[Yoshihara2013]** Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun. 2013;4:2612. doi:10.1038/ncomms3612. PMID: 24113773.
  ESTIMATE infers stromal and immune admixture from expression; its scores correlate with copy-number-based purity.
- **[Aitchison1982]** Aitchison J. The statistical analysis of compositional data. J R Stat Soc Series B Stat Methodol. 1982;44(2):139-160. doi:10.1111/j.2517-6161.1982.tb01195.x.
  Compositional data (proportions of a whole) need their own methods; the abstract introduces transformed-normal distributions on the simplex, the logistic-normal family defined through log ratios of the parts. Cited for the log-ratio approach.

## Licence terms (date-sensitive)

- **[Stanford2026]** Stanford University. CIBERSORTx Non-Commercial Terms of Use. https://cibersortx.stanford.edu/ (accessed November 2026).
  The service is provided free of charge for non-commercial use only, and use by any commercial entity for any purpose, including research, is prohibited; the terms grant a revocable, non-exclusive, non-transferable right to use the service and reserve all rights not expressly granted. They do not mention redistribution of downloaded files such as LM22, and the LM22-specific notice shown after login was not read. Re-check before relying on it.
