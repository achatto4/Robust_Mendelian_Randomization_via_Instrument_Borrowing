# IB-Mode and IB-PRESSO are provided by the IBMR package.
# Install with: remotes::install_github("achatto4/IBMR")
if (!requireNamespace("IBMR", quietly = TRUE)) {
  stop("Package 'IBMR' is required. Install with: remotes::install_github('achatto4/IBMR')")
}

# Shared helper: reads MR-PRESSO's outlier-corrected row (repository root).
.pc_cand <- c("presso_corrected.R", "../presso_corrected.R")
.pc_hit  <- .pc_cand[file.exists(.pc_cand)]
if (!length(.pc_hit)) stop("Cannot find presso_corrected.R at the repository root.")
source(.pc_hit[1])

# Cross-trait LDSC intercept helper for the IB-Mode C4 (outcome-sample-overlap) correction.
# Resolved relative to common run locations; if absent, .LDSC_I stays NULL (independent bootstrap).
.ldsc_helper_cand <- c("ldsc_cov_helper.R", "RDA_code/ldsc_cov_helper.R")
.ldsc_helper_hit  <- .ldsc_helper_cand[file.exists(.ldsc_helper_cand)]
if (length(.ldsc_helper_hit)) source(.ldsc_helper_hit[1])
# Resolve LDSC_intercepts.csv next to the helper (in-repo only); NULL if not shipped.
.LDSC_I <- if (exists("load_ldsc_intercepts") && length(.ldsc_helper_hit))
  load_ldsc_intercepts(file.path(dirname(.ldsc_helper_hit[1]), "LDSC_intercepts.csv")) else NULL

empty_result_row <- function(method, nsnp = NA) {
  data.frame(method = method, nsnp = nsnp, b = NA, se = NA, pval = NA)
}

safe_try <- function(expr) {
  out <- try(expr, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

append_outcome_labels <- function(result_df, outcome_data1, outcome_data2, add_col) {
  if (!add_col) return(result_df)

  result_df$outcome1 <- gsub("outcome_dat_", "", deparse(substitute(outcome_data1)))
  result_df$outcome2 <- gsub("outcome_dat_", "", deparse(substitute(outcome_data2)))
  result_df
}

prepare_common_instruments <- function(exposure_data, outcome_data1, outcome_data2, pval_threshold) {
  h1 <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data1)
  h2 <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data2)

  f11 <- h1 %>% dplyr::filter(pval.exposure <= pval_threshold & mr_keep == TRUE)
  f21 <- h2 %>% dplyr::filter(pval.exposure <= pval_threshold & mr_keep == TRUE)

  snp_common <- intersect(f11$SNP, f21$SNP)
  f1 <- f11 %>% dplyr::filter(SNP %in% snp_common)
  f2 <- f21 %>% dplyr::filter(SNP %in% snp_common)

  list(filtered_data1 = f1, filtered_data2 = f2, filtered_data11 = f11, filtered_data21 = f21)
}

make_presso_subset <- function(filtered_data1, filtered_data11, filtered_data21, presso_nmax, pressed_seed) {
  n_presso <- nrow(filtered_data1)
  idx_presso <- seq_len(n_presso)

  if (n_presso > presso_nmax) {
    if (!is.null(pressed_seed)) {
      snp_key <- paste(sort(as.character(filtered_data1$SNP)), collapse = "|")
      hash <- sum(as.integer(charToRaw(snp_key)))
      set.seed(abs(hash) + as.integer(pressed_seed))
    }
    idx_presso <- sort(sample(n_presso, presso_nmax))
  }

  sampled <- filtered_data1[idx_presso, , drop = FALSE]
  f11_presso <- filtered_data11[filtered_data11$SNP %in% sampled$SNP, , drop = FALSE]
  f21_presso <- filtered_data21[filtered_data21$SNP %in% sampled$SNP, , drop = FALSE]

  list(filtered_data11_presso = f11_presso, filtered_data21_presso = f21_presso)
}

run_ibmode <- function(filtered_data1, filtered_data2, phi, n_boot, alpha, ldsc_intercept = NULL) {
  # ldsc_intercept: cross-trait LD-score intercept for the (primary, auxiliary)
  # outcome pair. NULL -> independent bootstrap; supplying it makes IBMODE draw the
  # overlap-aware bivariate bootstrap used in the paper.
  fit <- IBMR::IBMODE(
    BetaXG = filtered_data1$beta.exposure,
    BetaYG_matrix = cbind(filtered_data1$beta.outcome, filtered_data2$beta.outcome),
    seBetaXG = filtered_data1$se.exposure,
    seBetaYG_matrix = cbind(filtered_data1$se.outcome, filtered_data2$se.outcome),
    phi = phi,
    n_boot = n_boot,
    alpha = alpha,
    ldsc_intercept = ldsc_intercept
  )

  # IBMODE returns one row per phi value; take the first row (primary phi).
  # Columns are named Estimate_Outcome1, SE_Outcome1, P_Outcome1 for the primary outcome.
  data.frame(
    method = "IB-Mode",
    nsnp = nrow(filtered_data1),
    b = fit$Estimate_Outcome1[1],
    se = fit$SE_Outcome1[1],
    pval = fit$P_Outcome1[1]
  )
}

run_ibpresso <- function(filtered_data11_presso, filtered_data21_presso, pressed_seed) {
  if (nrow(filtered_data11_presso) != nrow(filtered_data21_presso) || nrow(filtered_data11_presso) <= 3) {
    return(empty_result_row("IB-MR-PRESSO", nsnp = NA))
  }

  ibpresso_df <- data.frame(
    BetaOutcome = filtered_data11_presso$beta.outcome,
    BetaExposure = filtered_data11_presso$beta.exposure,
    BetaAux = filtered_data21_presso$beta.outcome,
    SdOutcome = filtered_data11_presso$se.outcome,
    SdExposure = filtered_data11_presso$se.exposure,
    SdAux = filtered_data21_presso$se.outcome
  )

  nsnp <- nrow(ibpresso_df)
  fit <- safe_try(
    IBMR::IBPRESSO(
      BetaOutcome = "BetaOutcome",
      BetaExposure = "BetaExposure",
      BetaAux = "BetaAux",
      SdOutcome = "SdOutcome",
      SdExposure = "SdExposure",
      SdAux = "SdAux",
      data = ibpresso_df,
      OUTLIERtest = TRUE,
      DISTORTIONtest = FALSE,
      NbDistribution = 10000,
      seed = pressed_seed,
      SignifThreshold = 0.05
    )
  )

  if (is.null(fit)) {
    return(empty_result_row("IB-MR-PRESSO", nsnp = nsnp))
  }

  b_val <- if (!is.null(fit$corrected_beta) && !is.na(fit$corrected_beta)) fit$corrected_beta else fit$raw_beta
  se_val <- if (!is.null(fit$corrected_se) && !is.na(fit$corrected_se)) fit$corrected_se else fit$raw_se

  n_outliers <- if (!is.null(fit$outlier_idx) && !all(is.na(fit$outlier_idx))) length(fit$outlier_idx) else 0
  n_post <- nsnp - n_outliers

  if (
    !is.null(fit$corrected_beta) && !is.na(fit$corrected_beta) &&
    !is.null(fit$corrected_se) && !is.na(fit$corrected_se) &&
    n_outliers > 0 && n_post > 1
  ) {
    pval <- 2 * pt(-abs(b_val / se_val), df = n_post - 1)
  } else if (nsnp > 1) {
    pval <- 2 * pt(-abs(b_val / se_val), df = nsnp - 1)
  } else {
    pval <- NA
  }

  data.frame(method = "IB-MR-PRESSO", nsnp = nsnp, b = b_val, se = se_val, pval = pval)
}

run_mrpresso <- function(filtered_data11_presso, count, pressed_seed) {
  if (count != 1 || !requireNamespace("MRPRESSO", quietly = TRUE) || nrow(filtered_data11_presso) <= 3) {
    return(empty_result_row("MR-PRESSO", nsnp = NA))
  }

  presso_df <- data.frame(
    BetaOutcome = filtered_data11_presso$beta.outcome,
    BetaExposure = filtered_data11_presso$beta.exposure,
    SdOutcome = filtered_data11_presso$se.outcome,
    SdExposure = filtered_data11_presso$se.exposure
  )

  nsnp <- nrow(presso_df)
  fit <- safe_try(
    MRPRESSO::mr_presso(
      BetaOutcome = "BetaOutcome",
      BetaExposure = "BetaExposure",
      SdOutcome = "SdOutcome",
      SdExposure = "SdExposure",
      data = presso_df,
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      NbDistribution = 10000,
      seed = pressed_seed,
      SignifThreshold = 0.05
    )
  )

  if (is.null(fit)) {
    return(empty_result_row("MR-PRESSO", nsnp = nsnp))
  }

  ## Outlier-corrected estimate; the raw fit when no instrument was removed.
  mrp <- .presso_corrected(fit)
  n_post <- nsnp - mrp$n_outliers
  df_p   <- if (mrp$corrected && n_post > 1) n_post - 1 else
            if (nsnp > 1) nsnp - 1 else NA
  data.frame(
    method = "MR-PRESSO",
    nsnp = nsnp,
    b = mrp$b,
    se = mrp$se,
    pval = if (!is.na(df_p)) 2 * stats::pt(-abs(mrp$b / mrp$se), df = df_p) else NA
  )
}

run_standard_mr_methods <- function(filtered_data11, n_CML, phi) {
  n <- nrow(filtered_data11)

  mr_obj <- safe_try(
    mr_input(
      bx = filtered_data11$beta.exposure,
      bxse = filtered_data11$se.exposure,
      by = filtered_data11$beta.outcome,
      byse = filtered_data11$se.outcome
    )
  )

  result_conmix <- if (!is.null(mr_obj)) safe_try(mr_conmix(mr_obj)) else NULL
  result_cml <- if (!is.null(mr_obj)) safe_try(mr_cML(mr_obj, n = n_CML)) else NULL

  result_mix <- safe_try(
    MRMix(
      filtered_data11$beta.exposure,
      filtered_data11$beta.outcome,
      sx = filtered_data11$se.exposure,
      sy = filtered_data11$se.outcome,
      theta_temp_vec = seq(-1, 1, by = 0.01),
      pi_init = 0.6,
      sigma_init = 1e-5
    )
  )

  param <- default_parameters()
  param$phi <- phi
  result_base <- safe_try(mr(filtered_data11, parameters = param)[, 5:9])

  if (is.null(result_base)) {
    result_base <- data.frame(
      method = c("Inverse variance weighted", "MR Egger", "Weighted median", "Simple mode", "Weighted mode"),
      nsnp = n,
      b = NA,
      se = NA,
      pval = NA
    )
  }
  if (!"nsnp" %in% names(result_base)) result_base$nsnp <- n

  row_conmix <- if (is.null(result_conmix)) {
    empty_result_row("MR-ConMix", nsnp = n)
  } else {
    ci_len <- result_conmix$CIUpper - result_conmix$CILower
    se <- if (is.null(ci_len) || any(is.na(ci_len))) NA else sum(ci_len) / 1.96 / 2
    data.frame(method = "MR-ConMix", nsnp = n, b = result_conmix$Estimate, se = se, pval = result_conmix$Pvalue)
  }

  row_cml <- if (is.null(result_cml)) {
    empty_result_row("MR-cML", nsnp = n)
  } else {
    data.frame(method = "MR-cML", nsnp = n, b = result_cml$Estimate, se = result_cml$StdError, pval = result_cml$Pvalue)
  }

  row_mix <- if (is.null(result_mix)) {
    empty_result_row("MR-Mix", nsnp = n)
  } else {
    data.frame(method = "MR-Mix", nsnp = n, b = result_mix$theta, se = result_mix$SE_theta, pval = result_mix$pvalue_theta)
  }

  do.call(rbind, list(result_base, row_conmix, row_cml, row_mix))
}

harmonize_and_evaluate <- function(
    exposure_data,
    outcome_data1,
    outcome_data2,
    phi = 1,
    n_boot = 100,
    alpha = 0.05,
    pval_threshold = 5e-8,
    add_col = FALSE,
    n_CML,
    count = NA,
    pressed_seed = 2025,
    ## Computational cap on the instruments passed to the two PRESSO arms,
    ## whose cost is O(NbDistribution * K); the other methods use the full
    ## instrument set.  At K <= 200 and NbDistribution = 10000 the bootstrap
    ## resolution K/NbDistribution = 0.02 stays below the 0.05 per-SNP
    ## threshold.
    presso_nmax = 200,
    outcome1_name = NULL, outcome2_name = NULL
) {
  # Step 1: harmonize and keep common SNPs.
  prepared <- prepare_common_instruments(exposure_data, outcome_data1, outcome_data2, pval_threshold)
  filtered_data1 <- prepared$filtered_data1
  filtered_data2 <- prepared$filtered_data2
  filtered_data11 <- prepared$filtered_data11
  filtered_data21 <- prepared$filtered_data21

  # Not enough instruments: return method table with NA values.
  if (nrow(filtered_data1) < 4 || nrow(filtered_data2) < 4) {
    methods <- c(
      "Inverse variance weighted", "MR Egger", "Weighted median", "Simple mode",
      "Weighted mode", "MR-ConMix", "MR-Mix", "MR-cML", "IB-Mode",
      "MR-PRESSO", "IB-MR-PRESSO"
    )
    out <- data.frame(method = methods, nsnp = NA, b = NA, se = NA, pval = NA)
    out <- append_outcome_labels(out, outcome_data1, outcome_data2, add_col)
    warning("Too few SNPs after harmonisation/intersection; returning NA rows for all methods.")
    return(out)
  }

  # Step 2: build PRESSO subset and run dual-outcome methods.
  presso_subset <- make_presso_subset(filtered_data1, filtered_data11, filtered_data21, presso_nmax, pressed_seed)
  filtered_data11_presso <- presso_subset$filtered_data11_presso
  filtered_data21_presso <- presso_subset$filtered_data21_presso

  # C4: cross-trait LD-score intercept for this (primary, auxiliary) outcome pair.
  ldsc_intercept <- if (!is.null(outcome1_name) && !is.null(outcome2_name) && !is.null(.LDSC_I))
    ldsc_I12(outcome1_name, outcome2_name, .LDSC_I) else NULL
  row_ibmode <- run_ibmode(filtered_data1, filtered_data2, phi, n_boot, alpha, ldsc_intercept = ldsc_intercept)
  row_ibpresso <- run_ibpresso(filtered_data11_presso, filtered_data21_presso, pressed_seed)

  # Step 3: optionally run single-outcome methods for the primary block.
  if (count == 1) {
    rows_standard <- run_standard_mr_methods(filtered_data11, n_CML, phi)
    row_presso <- run_mrpresso(filtered_data11_presso, count, pressed_seed)
    out <- do.call(rbind, list(rows_standard, row_ibmode, row_presso, row_ibpresso))
  } else {
    out <- do.call(rbind, list(row_ibmode, row_ibpresso))
  }

  out <- append_outcome_labels(out, outcome_data1, outcome_data2, add_col)
  print(out)
  out
}
