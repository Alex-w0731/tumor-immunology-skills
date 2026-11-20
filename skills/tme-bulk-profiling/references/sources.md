# Sources for tme-bulk-profiling

Every source key used in `SKILL.md`, `references/scoring-methods.md` and
`references/purity-and-phenotypes.md` resolves here. DOIs and PMIDs were checked against Crossref
and PubMed (November 2026). The indented line under each entry says what this skill takes from the
source and where that statement was checked (abstract or full text). Software and database entries
were read directly.

## Primary literature

- **[Aran2015]** Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. Nat Commun. 2015;6:8971. doi:10.1038/ncomms9971. PMID: 26634437.
  Tumor purity confounds correlation and clustering of tumor transcriptomes (abstract); the consensus purity estimate (CPE) is the median of four methods after normalising them to equal means and SDs (full text).
- **[Ayers2017]** Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127(8):2930-2940. doi:10.1172/JCI91190. PMID: 28650338.
  The T cell-inflamed gene expression profile was defined in pembrolizumab-treated patients and developed into a clinical-grade assay (abstract). tmelens copies gene lists from this paper (IFNG_Ayers6, TIS_Ayers18), not the assay.
- **[Barbie2009]** Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462(7269):108-112. doi:10.1038/nature08460. PMID: 19847166.
  Introduced the "single sample" extension of GSEA (the paper does not use the abbreviation ssGSEA; [Hanzelmann2013] attributes ssGSEA to it): genes are ranked by expression within one sample, and the enrichment score is the sum over the ranked list of the difference between a rank-weighted ECDF of the set genes (exponent alpha = 1/4) and the ECDF of the remaining genes (PMC author manuscript, Methods).
- **[Benjamini1995]** Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
  False discovery rate control across the signatures tested in one analysis.
- **[Carter2012]** Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA alterations in human cancer. Nat Biotechnol. 2012;30(5):413-421. doi:10.1038/nbt.2203. PMID: 22544022.
  ABSOLUTE infers tumor purity and ploidy from somatic DNA alterations, independently of expression (abstract).
- **[Chen2017]** Chen DS, Mellman I. Elements of cancer immunity and the cancer-immune set point. Nature. 2017;541(7637):321-330. doi:10.1038/nature21349. PMID: 28102259.
  The cancer-immune set point framework in which immune profiles are discussed (abstract); tmelens cites it for its phenotype framework. The phenotype definitions used in this skill were checked in [Mariathasan2018].
- **[Danaher2017]** Danaher P, Warren S, Dennis L, et al. Gene expression markers of Tumor Infiltrating Leukocytes. J Immunother Cancer. 2017;5:18. doi:10.1186/s40425-017-0215-8. PMID: 28239471.
  Cell-type scores are the simple average log expression of the marker genes; that average equals the log abundance of the cell type plus an unknown constant, so one log2 unit is a doubling within a score (full text). Marker lists used by the tmelens infiltration axis.
- **[Hanzelmann2013]** Hanzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7. doi:10.1186/1471-2105-14-7. PMID: 23323831.
  GSVA evaluates each gene in the context of the sample population (abstract and full text); ssGSEA uses within-sample expression ranks and normalises by the range of values across all gene sets and samples; the combined z-score standardises each gene and then combines the z-scores per sample (full text).
- **[Lee2008]** Lee E, Chuang HY, Kim JW, Ideker T, Lee D. Inferring pathway activity toward precise disease classification. PLoS Comput Biol. 2008;4(11):e1000217. doi:10.1371/journal.pcbi.1000217. PMID: 18989396.
  The combined z-score divides the sum of member-gene z-scores by the square root of the number of genes (full text).
- **[Mariathasan2018]** Mariathasan S, Turley SJ, Nickles D, et al. TGFbeta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018;554(7693):544-548. doi:10.1038/nature25501. PMID: 29443960.
  Lack of response associated with a fibroblast TGF-beta signature, particularly where CD8 T cells were excluded from the tumor parenchyma and found in the peritumoral stroma (abstract). Most solid tumors show an inflamed, excluded or desert phenotype; the metastatic urothelial cohort was 47% excluded, 27% desert and 26% inflamed (text, Fig. 2a,b and Extended Data Fig. 2f; Fig. 2a shows the phenotypes on a combined CD8 immunohistochemistry and trichrome stain); mean CD8 T-effector signature rose from desert to excluded to inflamed; the Pan-F-TBRS signature was generated here (full text).
- **[Rooney2015]** Rooney MS, Shukla SA, Wu CJ, Getz G, Hacohen N. Molecular and genetic properties of tumors associated with local immune cytolytic activity. Cell. 2015;160(1-2):48-61. doi:10.1016/j.cell.2014.12.033. PMID: 25594174.
  Cytolytic activity (CYT) was calculated as the geometric mean of GZMA and PRF1 expressed in TPM with a 0.01 offset (full text, PMC author manuscript).
- **[Yoshihara2013]** Yoshihara K, Shahmoradgoli M, Martinez E, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun. 2013;4:2612. doi:10.1038/ncomms3612. PMID: 24113773.
  ESTIMATE computes ssGSEA stromal and immune scores (141 genes each) that form the basis of the ESTIMATE score used to infer purity (full text); ESTIMATE scores correlate with copy-number-based purity (abstract).

## Software and databases

- **[Wang2026]** Wang A. tmelens: a transparent toolkit to profile the tumour immune microenvironment from bulk transcriptomes. Software, version 0.1.0, released 2026-10-17. MIT licence. https://github.com/Alex-w0731/tmelens
  Read directly (November 2026): `docs/methods.md`, `docs/signatures.md`, `src/tmelens/preprocess.py`, `scoring.py`, `phenotype.py`, `signatures.py`, `cli.py`. Source of every statement about tmelens defaults and behaviour in this skill, including the command line: `cli.main` maps `--log auto|yes|no` to `"auto"`, `True` and `False` before `prepare_expression` is called (read from the code; re-check in later versions).
- **[HGNC2026]** HUGO Gene Nomenclature Committee. HGNC database of approved human gene symbols. https://www.genenames.org (queried through https://rest.genenames.org, November 2026).
  All 48 anchor symbols that `check_expression_matrix.py` uses to recognise the gene axis were "Approved" when queried; CCN2, RFLNB and FAM30A list CTGF, FAM101B and KIAA0125 as previous symbols; MARCHF1 and MTARC1 list MARCH1 and MARC1.
