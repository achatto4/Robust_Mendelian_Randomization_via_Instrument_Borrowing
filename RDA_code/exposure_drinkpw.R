# Exposure-specific multimode MR pipeline.
# Loads one exposure, evaluates multiple outcomes, and writes result tables.

library(TwoSampleMR)
library(ggplot2)
library(MendelianRandomization)
library(dplyr)
library(MRMix)

RSID_FILE="./data/GWAS/MrcML_sumstats/plink2_clump_results/selected_snp_drinkpw.txt"
snp_vector <- readLines(RSID_FILE)[-1]

exposure_dat_drinkpw <- read_exposure_data(
  filename = "./data/GWAS/MrcML_sumstats/summary_txt/drinkpw.sumstats.txt",
  sep = "\t",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

exposure_dat_drinkpw <- exposure_dat_drinkpw %>%
  filter(SNP %in% snp_vector)

bmi_file <- "./data/GWAS/bmi_GWAS_ukbb/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"

outcome_dat_bmi <- read_outcome_data(
  snps = snp_vector,
  filename = bmi_file,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  eaf_col = NA,
  pval_col = "P",
  units_col = NA,
  gene_col = NA,
  samplesize_col = "N"
)

outcome_dat_t2d <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/bmi_GWAS_ukbb/DIAMANTE-EUR.sumstat.txt",
  sep = " ",
  snp_col = "rsID",
  beta_col = "Fixed-effects_beta",
  se_col = "Fixed-effects_SE",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "Fixed-effects_p-value"
)

outcome_dat_cad2 <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/CAD_GWAS_ukbb/cad.add.160614.website.txt",
  sep = "\t",
  snp_col = "markername",
  beta_col = "beta",
  se_col = "se_dgc",
  effect_allele_col = "effect_allele",
  other_allele_col = "noneffect_allele",
  eaf_col = "effect_allele_freq",
  pval_col = "p_dgc"
)

outcome_dat_str <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/Stroke_GWAS/metastroke.all.chr.bp",
  sep = "\t",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value"
)

outcome_dat_LDL <- read_outcome_data(
  filename = "./data/GWAS/lipid_GWAS/filtered_LDL_results_with_header_base_drinkpw.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "EFFECT_SIZE",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "POOLED_ALT_AF",
  pval_col = "pvalue_GC",
  samplesize_col = "N"
)

outcome_dat_HDL <- read_outcome_data(
  filename = "./data/GWAS/lipid_GWAS/filtered_HDL_results_with_header_base_drinkpw.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "EFFECT_SIZE",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "POOLED_ALT_AF",
  pval_col = "pvalue",
  samplesize_col = "N"
)

outcome_dat_logTG <- read_outcome_data(
  filename = "./data/GWAS/lipid_GWAS/filtered_logTG_results_with_header_base_drinkpw.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "EFFECT_SIZE",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "POOLED_ALT_AF",
  pval_col = "pvalue_GC",
  samplesize_col = "N"
)

outcome_dat_nonHDL <- read_outcome_data(
  filename = "./data/GWAS/lipid_GWAS/filtered_nonHDL_results_with_header_base_drinkpw.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "EFFECT_SIZE",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "POOLED_ALT_AF",
  pval_col = "pvalue_GC"
)

outcome_dat_TC <- read_outcome_data(
  filename = "./data/GWAS/lipid_GWAS/filtered_TC_results_with_header_base_drinkpw.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "EFFECT_SIZE",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "POOLED_ALT_AF",
  pval_col = "pvalue_GC"
)

outcome_dat_ckd <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/kidney_GWAS/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
  sep = " ",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
  samplesize_col = "n_total_sum"
)

outcome_dat_bc <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/functionalPRS/data/preprocessed_summary_statistics/EUR_all.UKB_aligned.ma",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p"
)

outcome_dat_pc <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/cancer_GWAS_ukbb/PancreaticC_GWAS_ukbb/GCST90043932_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

outcome_dat_sbp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/BP_GWAS/GCST90310294.h.tsv.gz",
  sep = "\t",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  samplesize_col = "n"
)

outcome_dat_dbp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/BP_GWAS/GCST90310295.h.tsv.gz",
  sep = "\t",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  samplesize_col = "n"
)

outcome_dat_pp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/BP_GWAS/GCST90310296.h.tsv.gz",
  sep = "\t",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  samplesize_col = "n"
)

outcome_dat_BC_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CaBrst.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size"
)

outcome_dat_LC_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CaLung.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size"
)

outcome_dat_PrC_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CaPros.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size"
)

outcome_dat_CC_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CaColon.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size"
)

outcome_dat_CAD_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CircCAD.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_T2D_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_250_2.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_stk_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CircStrk.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_TIA_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CircTIA.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_t1d_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_250_1.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_HLD_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_272_1.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_OSA_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_327_3.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_HTN_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =  "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_401.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_MI_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =  "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_411_2.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_AP_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_411_3.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_PHD_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_415_11.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_PE_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_415.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_CMG_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_416.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_CM_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_425.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_CCD_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_426.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_CD_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_427.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_CHF_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_428.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_AS_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =   "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_440.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_PVD_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =  "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Phe_443.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

outcome_dat_A1C_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.A1C_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_BMI_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.BMI_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_BNP_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.BNP_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_BUN_BSP_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.BUN_BSP_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_CKMB_Abs_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CKMB_Abs_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_Creat_BSP_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Creat_BSP_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_CRP_dL_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =  "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.CRP_dL_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_Diastolic_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Diastolic_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_eGFR_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename =  "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.eGFR_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_Glucose_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Glucose_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_HDLC_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.HDLC_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_LDLC_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.LDLC_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_Systolic_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Systolic_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_Trig_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Trig_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_TroponinI_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.TroponinI_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_Troponin_Mean_INT_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.Troponin_Mean_INT.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "num_samples")

outcome_dat_astm_mvp <- read_outcome_data(
  snps = snp_vector,
  filename = "./data/GWAS/MVP_data/MVP_base_traits/lifted_files/MVP_R4.1000G_AGR.DoAsth.EUR.GIA.dbGaP.lifted_with_header.bed.gz",
  sep = "\t",
  snp_col = "SNP_ID",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  samplesize_col = "eff_sample_size")

n_exp <- max(unique(exposure_dat_drinkpw$samplesize.exposure))

outcome_traits <- list(
  astm_mvp = list(
    data = outcome_dat_astm_mvp,
    aux = list(
      HDL   = outcome_dat_HDL,
      BNP_Mean_INT_mvp = outcome_dat_BNP_Mean_INT_mvp,
      CHF_mvp = outcome_dat_CHF_mvp
    )
  ),
  t2d_mvp = list(
    data = outcome_dat_T2D_mvp,
    aux = list(
      t2d   = outcome_dat_t2d,
      bmi   = outcome_dat_bmi,
      BMI_Mean_INT_mvp = outcome_dat_BMI_Mean_INT_mvp
    )
  ),
  cad_mvp = list(
    data = outcome_dat_CAD_mvp,
    aux = list(
      nHDL  = outcome_dat_nonHDL,
      LDLC_Mean_INT_mvp = outcome_dat_LDLC_Mean_INT_mvp,
      TC    = outcome_dat_TC
    )
  ),
  stk_mvp = list(
    data = outcome_dat_stk_mvp,
    aux = list(
      HDL   = outcome_dat_HDL,
      BUN_BSP_Mean_INT_mvp = outcome_dat_BUN_BSP_Mean_INT_mvp,
      nHDL  = outcome_dat_nonHDL
    )
  )
)

final_results <- data.frame()
# Run multimode MR across trait pairs
# Resolve the shared pipeline relative to this script's location (works from any cwd)
.this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
if (is.null(.this_dir) || !nzchar(.this_dir)) .this_dir <- getwd()
.pipeline <- file.path(.this_dir, "ib_mr_methods_pipeline.R")
if (!file.exists(.pipeline)) .pipeline <- "ib_mr_methods_pipeline.R"  # fallback: cwd = RDA_code/
source(.pipeline)

for (trait1_name in names(outcome_traits)) {
  trait1_entry <- outcome_traits[[trait1_name]]
  trait1_data  <- trait1_entry$data
  aux_traits   <- trait1_entry$aux

  if ("samplesize.outcome" %in% names(trait1_data)) {
    n1_values <- unique(trait1_data$samplesize.outcome)
    if (length(n1_values) > 1) {
      warning(sprintf("Trait '%s': multiple distinct N values found (%s); using max()",
                      trait1_name, paste(n1_values, collapse = ", ")))
    }
    n1 <- max(n1_values, na.rm = TRUE)
  } else {
    warning("no sample size data")
  }

  if (length(aux_traits) == 0) next
  for (count in seq_along(aux_traits)) {
    trait2_name <- names(aux_traits)[count]
    trait2_data <- aux_traits[[trait2_name]]

    n_CML <- min(n_exp, n1, na.rm = TRUE)
    if (!is.finite(n_CML) || is.na(n_CML)) n_CML <- 200000

    result <- harmonize_and_evaluate(
      exposure_data  = exposure_dat_drinkpw,
      outcome_data1  = trait1_data,
      outcome_data2  = trait2_data,
      phi            = 0.5,
      n_boot         = 100,
      alpha          = 0.05,
      count          = count,
      n_CML          = n_CML
    )

    result$trait1 <- trait1_name
    result$trait2 <- trait2_name
    final_results <- rbind(final_results, result)
  }
}

library(dplyr)

# NOTE: GWAS summary statistics must be downloaded separately and placed in ./data/GWAS/
# Output files are written to ./results/ (create this directory before running).


final_results <- final_results %>%
  arrange(trait1, method) %>%
  group_by(trait1) %>%
  filter(
    method %in% c("new_method_weighted", "IB-MR-PRESSO") |
      !duplicated(method)
  ) %>%
  mutate(
    trait2 = ifelse(method %in% c("new_method_weighted", "IB-MR-PRESSO"), trait2, NA_character_)
  ) %>%
  ungroup()

print(final_results, n = nrow(final_results))

write.csv(final_results , file = "../../results/RDA_results/final_results_drinkpw.csv", row.names = FALSE)

quit(save = "no", status = 0)
