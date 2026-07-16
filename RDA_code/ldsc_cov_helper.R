## ===========================================================================
## LDSC cross-trait intercept helper for the IB-Mode bivariate bootstrap (C4).
##
## Provides the cross-trait LDSC intercept I12 for a (primary, auxiliary) outcome
## pair via load_ldsc_intercepts() and ldsc_I12(). This scalar is passed to
## IBMR::IBMODE(ldsc_intercept = I12), which builds the per-SNP cross-trait
## covariance of the two Wald ratios,
##     sigma_{12,k} = I12 * sey1 * sey2 / bx^2 + (by1 * by2) * vbx / bx^4,
## and redraws each instrument's outcome pair from a bivariate normal that
## accounts for outcome-GWAS sample overlap. With no overlap (I12 ~ 0) it reduces
## toward the independent redrawing.
## ===========================================================================

# Saved LDSC intercept matrix = GenomicSEM ldsc()$I.
# Not distributed with this repo: place LDSC_intercepts.csv in RDA_code/ to enable the
# C4 correction. Resolved relative to the working directory (RDA_code/ when the exposure
# scripts run); if absent, the cross-trait covariance is treated as 0 (independent bootstrap).
LDSC_INTERCEPT_PATH <- "LDSC_intercepts.csv"

# Row/column order of that matrix == the trait.names passed to ldsc() when
# LDSC_intercepts.csv was written.
LDSC_TRAIT_ORDER <- c(
  "astm","t2d","bmi","cad","str","HDL","LDL","logTG","nHDL","TC",
  "ckd","bc","pc","sbp","dbp","pp","t2d_mvp","cad_mvp","stk_mvp","TIA_mvp",
  "t1d_mvp","HLD_mvp","OSA_mvp","HTN_mvp","MI_mvp","AP_mvp","PHD_mvp","PE_mvp",
  "CMG_mvp","CM_mvp","CCD_mvp","CD_mvp","CHF_mvp","AS_mvp","PVD_mvp",
  "A1C_Mean_INT_mvp","BMI_Mean_INT_mvp","BNP_Mean_INT_mvp","BUN_BSP_Mean_INT_mvp",
  "CKMB_Abs_Mean_INT_mvp","Creat_BSP_Mean_INT_mvp","CRP_dL_Mean_INT_mvp",
  "Diastolic_Mean_INT_mvp","eGFR_Mean_INT_mvp","Glucose_Mean_INT_mvp",
  "HDLC_Mean_INT_mvp","LDLC_Mean_INT_mvp","Systolic_Mean_INT_mvp",
  "Trig_Mean_INT_mvp","TroponinI_Mean_INT_mvp"
)

# RDA outcome short-name -> name in LDSC_TRAIT_ORDER.
#   astm_mvp : confirmed to be the position-1 'astm' row (same MVP asthma GWAS).
#   cad2/nonHDL : consortium-named outcomes mapping to cad / nHDL.
LDSC_NAME_ALIAS <- c(
  astm_mvp = "astm",
  cad2     = "cad",
  nonHDL   = "nHDL"
)

# Load the intercept matrix once; returns a named matrix or NULL if absent.
load_ldsc_intercepts <- function(path = LDSC_INTERCEPT_PATH) {
  # Try the given path, then in-repo fallbacks (never escapes the repo via ../).
  candidates <- unique(c(path, "LDSC_intercepts.csv", "RDA_code/LDSC_intercepts.csv"))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    warning(sprintf("LDSC intercept file not found (tried: %s); cross-trait covariance treated as 0.",
                    paste(candidates, collapse = ", ")))
    return(NULL)
  }
  path <- hit[1]
  # write.table() saved a leading row-index column. read.csv may either auto-use
  # it as row names (-> 50 cols) or keep it as a data column (-> 51 cols); handle
  # both, then require an exact n x n matrix.
  M0 <- read.csv(path, header = TRUE, sep = "\t", check.names = FALSE)
  if (ncol(M0) == length(LDSC_TRAIT_ORDER) + 1) M0 <- M0[, -1, drop = FALSE]
  M <- as.matrix(M0)
  if (nrow(M) != length(LDSC_TRAIT_ORDER) || ncol(M) != length(LDSC_TRAIT_ORDER)) {
    stop(sprintf("LDSC intercept matrix is %dx%d but LDSC_TRAIT_ORDER has %d traits.",
                 nrow(M), ncol(M), length(LDSC_TRAIT_ORDER)))
  }
  storage.mode(M) <- "double"
  rownames(M) <- colnames(M) <- LDSC_TRAIT_ORDER
  M
}

# Resolve an RDA outcome name to a matrix trait name (alias + case-insensitive).
.ldsc_resolve_name <- function(nm) {
  if (is.null(nm) || is.na(nm)) return(NA_character_)
  if (nm %in% LDSC_TRAIT_ORDER)        return(nm)
  if (nm %in% names(LDSC_NAME_ALIAS))  return(unname(LDSC_NAME_ALIAS[nm]))
  hit <- LDSC_TRAIT_ORDER[tolower(LDSC_TRAIT_ORDER) == tolower(nm)]
  if (length(hit) == 1) return(hit)
  NA_character_
}

# Off-diagonal cross-trait LDSC intercept I12 for an outcome pair.
# STOPs (rather than silently using 0) if a name cannot be mapped, so a naming
# mismatch can never quietly corrupt the standard errors.
ldsc_I12 <- function(name1, name2, M) {
  if (is.null(M)) return(0)
  r1 <- .ldsc_resolve_name(name1); r2 <- .ldsc_resolve_name(name2)
  if (is.na(r1) || is.na(r2)) {
    stop(sprintf("ldsc_I12: outcome name(s) not in LDSC trait list: '%s'->%s, '%s'->%s",
                 name1, r1, name2, r2))
  }
  M[r1, r2]
}
