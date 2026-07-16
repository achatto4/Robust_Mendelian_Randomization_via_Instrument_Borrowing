# ============================================================
# CONDITIONAL INFERENCE FOR SELECTED INSTRUMENTS
#
# TARGET: Coheterogeneity among SELECTED SNPs -- the conditional estimand
#         rho_CH^(n), the precision-weighted correlation of the ratio-scale
#         pleiotropic deviations across the K SNPs that passed selection.
#
# DGP: confounder model (paper eq. 7). Invalid instruments load on a shared
#      confounder U0 (a fraction D_ov of them) or on outcome-specific
#      confounders U1 / U2. D_ov is the invalid-instrument overlap and is
#      the only source of cross-trait pleiotropic correlation, so the
#      coheterogeneity increases monotonically with D_ov (Lemma 1).
#
# CONDITIONAL INFERENCE:
#   - Oracle: TRUE coheterogeneity among selected SNPs (conditional target)
#   - Estimator: bias-corrected weighted estimate of that correlation
#   - Coverage: Does the CI cover the oracle? (Should be ~95% for large K)
# ============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(IBMR)
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

select <- dplyr::select
filter <- dplyr::filter
eps <- 1e-12

# ============================================================
# ARCHITECTURE
# ============================================================

# Confounder-overlap (D_ov) DGP -- paper eq. (7).
# Invalid instruments draw their pleiotropy from a shared confounder U0
# (a fraction D_ov of them) or from outcome-specific confounders U1 / U2.
# Only the shared confounder U0 creates cross-trait pleiotropic correlation,
# so the coheterogeneity grows monotonically with the overlap D_ov (Lemma 1).
generate_architecture <- function(M, D_ov,
                                  sigma2x = 2e-4,
                                  pi_assoc = 0.01,
                                  prop_invalid = 0.7,
                                  theta_ux = 0.3,
                                  theta_uy1 = 0.5,
                                  theta_uy2 = 0.5,
                                  sigma2u = 5e-4,
                                  sigma2d = 1e-4,
                                  theta_causal = c(0.2, 0.3),
                                  r1 = 0.5, r2 = 0.5,
                                  seed = 999) {
  set.seed(seed)

  # Exposure-variance attenuation: invalid instruments draw part of their
  # X-effect from a confounder, so the residual direct effect is shrunk to
  # keep Var(beta_x) = sigma2x for valid and invalid instruments alike.
  sigma2x_tilde <- sigma2x - theta_ux^2 * sigma2u
  if (sigma2x_tilde < 0) {
    stop("sigma2x - theta_ux^2 * sigma2u < 0: reduce sigma2u or theta_ux.")
  }

  # SNPs associated with the exposure form the candidate instrument pool.
  n_assoc <- round(M * pi_assoc)
  assoc_snps <- sample.int(M, n_assoc)

  n_invalid <- round(n_assoc * prop_invalid)
  n_valid <- n_assoc - n_invalid
  valid_snps <- if (n_valid > 0) assoc_snps[1:n_valid] else integer(0)
  invalid_snps <- if (n_invalid > 0) assoc_snps[(n_valid + 1):n_assoc] else integer(0)

  beta_x <- alpha1 <- alpha2 <- rep(0, M)

  # Valid instruments: exposure effect only, no pleiotropy.
  if (n_valid > 0) {
    beta_x[valid_snps] <- rnorm(n_valid, 0, sqrt(sigma2x))
  }

  # Invalid instruments split by confounder pathway:
  #   U0 = shared confounder      (fraction D_ov of the invalid instruments)
  #   U1 = Y1-specific confounder, U2 = Y2-specific confounder.
  U0_snps <- U1_snps <- U2_snps <- integer(0)
  if (n_invalid > 0) {
    probs <- c(D_ov, (1 - D_ov) * r1, (1 - D_ov) * r2)
    probs <- probs / sum(probs)
    counts <- as.integer(rmultinom(1, n_invalid, probs))
    shuffled <- sample(invalid_snps)
    cuts <- cumsum(counts)
    if (counts[1] > 0) U0_snps <- shuffled[1:cuts[1]]
    if (counts[2] > 0) U1_snps <- shuffled[(cuts[1] + 1):cuts[2]]
    if (counts[3] > 0) U2_snps <- shuffled[(cuts[2] + 1):cuts[3]]

    groups <- list(list(idx = U0_snps, on1 = TRUE,  on2 = TRUE),
                   list(idx = U1_snps, on1 = TRUE,  on2 = FALSE),
                   list(idx = U2_snps, on1 = FALSE, on2 = TRUE))

    for (g in groups) {
      n <- length(g$idx)
      if (n == 0) next
      gamma_g <- rnorm(n, 0, sqrt(sigma2x_tilde))  # residual direct X-effect
      phi_g   <- rnorm(n, 0, sqrt(sigma2u))        # confounder loading
      d1_g    <- rnorm(n, 0, sqrt(sigma2d))        # intrinsic direct effect on Y1
      d2_g    <- rnorm(n, 0, sqrt(sigma2d))        # intrinsic direct effect on Y2
      beta_x[g$idx] <- gamma_g + theta_ux * phi_g
      alpha1[g$idx] <- (if (g$on1) theta_uy1 * phi_g else 0) + d1_g
      alpha2[g$idx] <- (if (g$on2) theta_uy2 * phi_g else 0) + d2_g
    }
  }

  betay1_true <- theta_causal[1] * beta_x + alpha1
  betay2_true <- theta_causal[2] * beta_x + alpha2

  list(M = M,
       betax_true = beta_x,
       betay1_true = betay1_true,
       betay2_true = betay2_true,
       alpha1 = alpha1,
       alpha2 = alpha2,
       assoc_snps = assoc_snps,
       valid_snps = valid_snps,
       invalid_snps = invalid_snps,
       U0_snps = U0_snps,
       U1_snps = U1_snps,
       U2_snps = U2_snps,
       n_valid = n_valid,
       n_invalid = n_invalid,
       D_ov = D_ov,
       theta_causal = theta_causal,
       pi_assoc = pi_assoc)
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================

gwas_cov_by_overlap <- function(N_y1, N_y2, overlap, rho_e) {
  m <- overlap * min(N_y1, N_y2)
  m <- max(0, min(m, min(N_y1, N_y2)))
  rho_e * m / (N_y1 * N_y2)
}

compute_sigma_terms <- function(bx, by1, by2,
                                var_bx, var_by1, var_by2, cov_by12,
                                eps = 1e-12) {
  bx2 <- bx^2 + eps
  bx4 <- bx2^2
  s1 <- (var_by1 / bx2) + (by1^2 * var_bx / bx4)
  s2 <- (var_by2 / bx2) + (by2^2 * var_bx / bx4)
  s12 <- (cov_by12 / bx2) + (by1 * by2 * var_bx / bx4)
  s1[!is.finite(s1)] <- NA
  s2[!is.finite(s2)] <- NA
  s12[!is.finite(s12)] <- NA
  list(s1 = s1, s2 = s2, s12 = s12)
}

# ============================================================
# ESTIMATOR
# ============================================================

.cohet_pkg <- function(bx, by1, by2, N_x, N_y1, N_y2, overlap = 0, rho_e = 0) {
  seX  <- rep(sqrt(1 / N_x),  length(bx))
  seY1 <- rep(sqrt(1 / N_y1), length(bx))
  seY2 <- rep(sqrt(1 / N_y2), length(bx))
  I12  <- rho_e * overlap * min(N_y1, N_y2) / sqrt(N_y1 * N_y2)
  IBMR::coheterogeneity_Q(
    BetaXG = bx,
    BetaYG_matrix = cbind(by1, by2),
    seBetaXG = seX,
    seBetaYG_matrix = cbind(seY1, seY2),
    ldsc_intercepts = matrix(c(1, I12, I12, 1), 2, 2),
    use_ldsc = TRUE,
    min_K_pair = 3
  )
}

rho_hat_T <- function(bx, by1, by2, N_x, N_y1, N_y2, overlap = 0, rho_e = 0, eps = 1e-12)
  .cohet_pkg(bx, by1, by2, N_x, N_y1, N_y2, overlap, rho_e)$rho[1, 2]

# ============================================================
# INSTRUMENT SELECTION
# ============================================================

select_instruments <- function(arch, N_x, F_threshold=10, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  M <- arch$M
  bx_true <- arch$betax_true
  
  bx_hat <- bx_true + rnorm(M, 0, sqrt(1/N_x))
  z_scores <- sqrt(N_x) * bx_hat
  pval <- 2 * pnorm(-abs(z_scores))
  
  pthr <- 0.05 / M
  
  idx <- which(pval < pthr)
  if (length(idx) == 0) return(list(idx=NULL, K=0, bx_hat=bx_hat))
  
  F_stat <- N_x * bx_hat[idx]^2
  idx_strong <- idx[F_stat > F_threshold]
  if (length(idx_strong) == 0) return(list(idx=NULL, K=0, bx_hat=bx_hat))
  
  idx_sorted <- idx_strong[order(pval[idx_strong])]
  
  list(idx=idx_sorted, K=length(idx_sorted), bx_hat=bx_hat)
}

# ============================================================
# ORACLE TARGET: TRUE correlation among SELECTED SNPs
# ============================================================

oracle_target <- function(arch, idx, N_x, N_y1, N_y2, overlap, rho_e) {
  if (is.null(idx) || length(idx) < 3) return(NA_real_)
  
  bx_true  <- arch$betax_true[idx]
  by1_true <- arch$betay1_true[idx]
  by2_true <- arch$betay2_true[idx]
  
  valid <- abs(bx_true) > 1e-10
  if (sum(valid) < 3) return(NA_real_)
  
  bx_true  <- bx_true[valid]
  by1_true <- by1_true[valid]
  by2_true <- by2_true[valid]
  
  theta1_true <- by1_true / bx_true
  theta2_true <- by2_true / bx_true
  
  var_bx  <- 1 / N_x
  var_by1 <- 1 / N_y1
  var_by2 <- 1 / N_y2
  cov_by12 <- gwas_cov_by_overlap(N_y1, N_y2, overlap, rho_e)
  
  sig <- compute_sigma_terms(bx_true, by1_true, by2_true,
                             var_bx, var_by1, var_by2, cov_by12, eps)
  
  s1 <- sig$s1
  s2 <- sig$s2
  s12 <- sig$s12
  
  ok <- is.finite(s1) & is.finite(s2) & is.finite(s12) & s1 > 0 & s2 > 0
  if (sum(ok) < 3) return(NA_real_)
  
  theta1_true <- theta1_true[ok]
  theta2_true <- theta2_true[ok]
  s1 <- s1[ok]
  s2 <- s2[ok]
  s12 <- s12[ok]
  
  w_raw <- 1 / sqrt((s1 + eps) * (s2 + eps))
  w <- w_raw / sum(w_raw)
  
  theta1_bar <- sum(w * theta1_true)
  theta2_bar <- sum(w * theta2_true)
  
  d1 <- theta1_true - theta1_bar
  d2 <- theta2_true - theta2_bar
  
  S11 <- sum(w * d1^2)
  S22 <- sum(w * d2^2)
  S12 <- sum(w * d1 * d2)
  
  if (S11 <= eps || S22 <= eps) return(NA_real_)
  
  rho <- S12 / sqrt(S11 * S22)
  max(-1, min(1, rho))
}

# ============================================================
# SIMULATE OUTCOMES
# ============================================================

simulate_outcomes <- function(arch, idx, N_y1, N_y2, overlap, rho_e) {
  K <- length(idx)
  if (K == 0) return(NULL)
  
  by1_true <- arch$betay1_true[idx]
  by2_true <- arch$betay2_true[idx]
  
  cov_by12 <- gwas_cov_by_overlap(N_y1, N_y2, overlap, rho_e)
  SigmaY <- matrix(c(1/N_y1, cov_by12, cov_by12, 1/N_y2), 2, 2)
  
  E <- mvrnorm(K, mu=c(0,0), Sigma=SigmaY)
  
  list(by1=by1_true + E[,1], by2=by2_true + E[,2])
}

estimate_rho_with_se <- function(bx, by1, by2,
                                 N_x, N_y1, N_y2,
                                 overlap=0, rho_e=0,
                                 alpha=0.05) {
  K <- length(bx)
  res <- .cohet_pkg(bx, by1, by2, N_x, N_y1, N_y2, overlap, rho_e)
  rho <- res$rho[1, 2]
  se <- res$se[1, 2]
  z <- qnorm(1 - alpha / 2)
  if (!is.finite(rho) || !is.finite(se)) {
    return(list(rho = rho, se = NA, ci_lower = NA, ci_upper = NA, pvalue = NA, K = K))
  }
  list(
    rho = rho,
    se = se,
    ci_lower = max(-1, rho - z * se),
    ci_upper = min(1, rho + z * se),
    pvalue = if (se > 0) 2 * pnorm(-abs(rho / se)) else NA,
    var = se^2,
    K = K
  )
}

# ============================================================
# SINGLE REPLICATION
# ============================================================

simulate_one_rep <- function(arch, N_x, N_y1, N_y2, overlap, rho_e,
                             F_threshold=10, alpha=0.05, seed=NULL) {
  
  sel <- select_instruments(arch, N_x, F_threshold, seed)
  
  if (is.null(sel$idx) || sel$K < 5) {
    return(list(rho=NA, se=NA, ci_lower=NA, ci_upper=NA, pvalue=NA,
                K=sel$K, rho_target=NA, error=NA))
  }
  
  bx_hat <- sel$bx_hat[sel$idx]
  out <- simulate_outcomes(arch, sel$idx, N_y1, N_y2, overlap, rho_e)
  
  rho_target <- oracle_target(arch, sel$idx, N_x, N_y1, N_y2, overlap, rho_e)
  est <- estimate_rho_with_se(bx_hat, out$by1, out$by2, N_x, N_y1, N_y2, overlap, rho_e, alpha)
  
  error <- if(is.finite(est$rho) && is.finite(rho_target)) est$rho - rho_target else NA
  
  list(rho=est$rho, se=est$se, ci_lower=est$ci_lower, ci_upper=est$ci_upper,
       pvalue=est$pvalue, K=sel$K, rho_target=rho_target, error=error)
}

# ============================================================
# RUN FOR SINGLE SAMPLE SIZE
# ============================================================

run_for_sample_size <- function(arch, N_x, N_y_ratio=1, N_rep=200,
                                overlap=0.5, rho_e=0.2,
                                F_threshold=10, alpha=0.05, verbose=TRUE) {
  
  N_y <- N_x * N_y_ratio
  
  results <- data.frame(
    rep = 1:N_rep,
    N_x = N_x,
    M = arch$M,
    K = NA_integer_,
    rho = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    pvalue = NA_real_,
    rho_target = NA_real_,
    error = NA_real_,
    cover = NA
  )
  
  for (r in 1:N_rep) {
    res <- simulate_one_rep(arch, N_x, N_y, N_y, overlap, rho_e,
                            F_threshold, alpha, seed = 10000 + r)
    
    results$K[r] <- res$K
    results$rho[r] <- res$rho
    results$se[r] <- res$se
    results$ci_lower[r] <- res$ci_lower
    results$ci_upper[r] <- res$ci_upper
    results$pvalue[r] <- res$pvalue
    results$rho_target[r] <- res$rho_target
    results$error[r] <- res$error
    
    results$cover[r] <- is.finite(res$rho_target) && is.finite(res$ci_lower) &&
      res$rho_target >= res$ci_lower && res$rho_target <= res$ci_upper
    
    if (verbose && r %% 50 == 0) {
      cat(sprintf("  N=%s, M=%s: rep %d/%d, K=%d\n", 
                  format(N_x, big.mark=","), format(arch$M, big.mark=","),
                  r, N_rep, res$K))
    }
  }
  
  results
}

# ============================================================
# MAIN STUDY
# ============================================================

run_study <- function(N_vec = c(5e4, 8e4, 1e5, 1.25e5, 1.625e5, 2e5),
                      M_func = function(N) ceiling(N^0.6),
                      N_rep = 200, alpha = 0.05,
                      D_ov = 0.75) {
  
  cat("\n")
  cat(rep("=", 100), "\n", sep="")
  cat("CONDITIONAL INFERENCE FOR SELECTED INSTRUMENTS\n")
  cat("Target: Correlation among selected SNPs (rep-specific oracle)\n")
  cat(rep("=", 100), "\n", sep="")
  
  all_results <- list()
  
  for (i in seq_along(N_vec)) {
    N_x <- N_vec[i]
    M <- M_func(N_x)
    
    cat(sprintf("\n[%d/%d] N = %s, M = %s\n", 
                i, length(N_vec), 
                format(N_x, big.mark=","),
                format(M, big.mark=",")))
    
    arch <- generate_architecture(M = M,
                                  D_ov = D_ov,
                                  sigma2x = 2e-4,
                                  pi_assoc = 0.01,
                                  prop_invalid = 0.7,
                                  theta_ux = 0.3,
                                  theta_uy1 = 0.5,
                                  theta_uy2 = 0.5,
                                  sigma2u = 5e-4,
                                  sigma2d = 1e-4,
                                  theta_causal = c(0.2, 0.3),
                                  seed = 22345 + i)
    
    results <- run_for_sample_size(arch, N_x, N_y_ratio = 1, N_rep = N_rep,
                                   overlap = 0.5, rho_e = 0.2,
                                   F_threshold = 10, alpha = alpha, verbose = TRUE)
    all_results[[i]] <- results
  }
  
  all_data <- do.call(rbind, all_results)
  
  # Summary statistics
  summary_stats <- all_data %>%
    filter(is.finite(rho) & is.finite(se) & K >= 5) %>%
    group_by(N_x, M) %>%
    summarise(
      n_valid = n(),
      K_mean = mean(K),
      K_sd = sd(K),
      rho_target_mean = mean(rho_target, na.rm=TRUE),
      rho_mean = mean(rho),
      rho_sd = sd(rho),
      bias = mean(error, na.rm=TRUE),
      rmse = sqrt(mean(error^2, na.rm=TRUE)),
      se_mean = mean(se),
      se_empirical = sd(rho),
      se_ratio = mean(se) / sd(rho),
      coverage = mean(cover, na.rm=TRUE),
      .groups = "drop"
    )
  
  # Print summary
  cat("\n")
  cat(rep("=", 100), "\n", sep="")
  cat("SUMMARY STATISTICS\n")
  cat(rep("=", 100), "\n", sep="")
  
  print_df <- data.frame(
    N = format(summary_stats$N_x, big.mark=","),
    K = sprintf("%.0f", summary_stats$K_mean),
    `Oracle(ρ)` = sprintf("%.3f", summary_stats$rho_target_mean),
    RMSE = sprintf("%.4f", summary_stats$rmse),
    Coverage = sprintf("%.1f%%", 100*summary_stats$coverage),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  print(print_df, row.names = FALSE)
  cat(rep("=", 100), "\n", sep="")
  
  cat(sprintf("\nMean coverage: %.1f%%\n", 100*mean(summary_stats$coverage)))
  cat(sprintf("Mean oracle ρ: %.3f (pleiotropy correlation)\n", mean(summary_stats$rho_target_mean)))
  
  list(data = all_data, summary = summary_stats)
}

# ============================================================
# CREATE THE TWO PLOTS (rounded N labels WITHOUT changing axis order)
# ============================================================

# ============================================================
# CREATE TWO READABLE PLOTS (coverage + ridge distribution)
# - Rounded N labels (e.g., 62,345 -> 62K)
# - Preserve correct axis ordering by numeric N_x
# - Coverage plot: fewer x labels + rotated + label only outliers
# - Ridge plot: fewer y labels + less overlap + zoomed x-range
# ============================================================
nice_round <- function(x) {
  ifelse(
    x < 1e5,
    round(x / 5000) * 5000,   # 5-digit: nearest 5,000
    round(x / 50000) * 50000  # 6-digit: nearest 50,000
  )
}

# Pretty formatter: 62345 -> 62K, 1000000 -> 1M
N_pretty <- function(x) {
  scales::label_number(
    scale_cut = scales::cut_short_scale(),
    accuracy  = 1
  )(nice_round(x))
}

create_two_plots <- function(results,
                             # readability knobs
                             x_label_every = 2,         # show every 2nd x tick label
                             y_label_every = 2,         # show every 2nd y tick label
                             outlier_low  = 0.90,       # label coverage < 90%
                             outlier_high = 0.975,      # label coverage > 97.5%
                             ridge_scale  = 1.8,        # lower => less overlap (was ~2.5)
                             ridge_xlim   = c(-0.2, 0.9) # zoom ridge x-range for readability
) {
  
  all_data <- results$data
  summary_stats <- results$summary
  
  #-------------------------------
  # Prep ridge data (rep-level)
  #-------------------------------
  ridge_data <- all_data %>%
    dplyr::filter(is.finite(rho) & is.finite(se) & K >= 5) %>%
    dplyr::mutate(
      N_label_pretty = paste0("n=", N_pretty(N_x))
    ) %>%
    dplyr::arrange(N_x) %>%
    dplyr::mutate(
      # preserve ordering by numeric N_x
      N_label_pretty = factor(N_label_pretty, levels = unique(N_label_pretty))
    )
  
  mean_oracle <- mean(ridge_data$rho_target, na.rm = TRUE)
  
  #====================================================================
  # PLOT 1: Coverage Probability
  #====================================================================
  coverage_data <- summary_stats %>%
    dplyr::mutate(
      N_label_pretty = paste0("n=", N_pretty(N_x)),
      N_K_label = sprintf("%s\nK=%d", N_label_pretty, round(K_mean))
    ) %>%
    dplyr::arrange(N_x) %>%
    dplyr::mutate(
      # preserve ordering by numeric N_x
      N_K_label = factor(N_K_label, levels = unique(N_K_label))
    )
  
  x_breaks <- levels(coverage_data$N_K_label)[
    seq(1, length(levels(coverage_data$N_K_label)), by = x_label_every)
  ]
  
  outlier_points <- coverage_data %>%
    dplyr::filter(is.finite(coverage) & (coverage < outlier_low | coverage > outlier_high))
  
  p1 <- ggplot(coverage_data, aes(x = N_K_label, y = coverage)) +
    geom_hline(yintercept = 0.95, linetype = "solid",
               color = "#E63946", linewidth = 1.1) +
    geom_line(aes(group = 1), color = "#1D3557", linewidth = 1.1) +
    geom_point(size = 3.5, color = "#1D3557") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0.5, 1.0)) +
    scale_x_discrete(breaks = x_breaks) +
    labs(
      title = "Coverage of 95% CI",
      x = "Sample size (n) and number of instruments (K)",
      y = "Coverage Probability"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "#1D3557"),
      plot.subtitle = element_text(color = "gray40", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.title = element_text(face = "bold", size = 12),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  #====================================================================
  # PLOT 2: Ridge distribution (readability improvements)
  #====================================================================
  y_levels <- levels(ridge_data$N_label_pretty)
  y_breaks <- y_levels[seq(1, length(y_levels), by = y_label_every)]
  
  p2 <- ggplot(ridge_data, aes(x = rho, y = N_label_pretty, fill = after_stat(x))) +
    ggridges::stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantile_lines = TRUE,
      quantiles = 2,            # median
      alpha = 0.85,
      scale = ridge_scale,      # less overlap than 2.5
      rel_min_height = 0.01,
      linewidth = 0.25
    ) +
    scale_fill_viridis_c(option = "mako", name = "ρ̂") +
    scale_x_continuous(breaks = seq(-1, 1, 0.25)) +
    scale_y_discrete(breaks = y_breaks) +
    coord_cartesian(xlim = ridge_xlim) +
    labs(
      title = "Distribution of coheterogeneity statistics",
      x = expression(hat(rho)[CH]^"(n)"),
      y = "Sample size (n)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "#1D3557"),
      plot.subtitle = element_text(color = "gray40", size = 11),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  list(coverage = p1, distribution = p2)
  }

# ============================================================
# EXECUTE
# ============================================================

# Load ggridges for ridge plot
if (!require("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
  library(ggridges)
}

# Run study
results <- run_study(
  N_vec = round(exp(seq(log(1e4), log(2e5), length.out = 15))),
  M_func = function(N) 100000,
  N_rep = 200,
  alpha = 0.05,
  D_ov = 0.75
)

# Create plots
plots <- create_two_plots(results)

# Display
cat("\n")
cat(rep("=", 100), "\n", sep="")
cat("DISPLAYING PLOTS\n")
cat(rep("=", 100), "\n", sep="")

print(plots$coverage)
print(plots$distribution)

# Combine into single figure
combined <- (plots$distribution | plots$coverage) +
  plot_annotation(
    title = "",
     theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
    )
  )

print(combined)

# save as EPS
ggsave(
  filename = "combined_plots_coverage.eps",
  plot     = combined,
  width    = 12,     # inches
  height   = 6.5,    # inches
  device   = cairo_ps
)


# ============================================================

# Population reference: the mean selected-set coheterogeneity rho_CH^(n)
# realised across the run. Under the confounder DGP this tracks D_ov
# (Lemma 1) rather than being a fixed constant.
lemma_val <- mean(results$summary$rho_target_mean, na.rm = TRUE)

# Plot the estimate distribution against the mean selected-set target
dist_plot_with_lemma <- results$data %>%
  filter(is.finite(rho) & is.finite(se) & K >= 5) %>%
  mutate(
    N_label = factor(paste0("N=", scales::label_number(scale_cut = scales::cut_short_scale())(N_x)),
                     levels = unique(paste0("N=", scales::label_number(scale_cut = scales::cut_short_scale())(sort(unique(N_x))))))
  ) %>%
  ggplot(aes(x = rho, y = N_label, fill = after_stat(x))) +
  geom_vline(xintercept = lemma_val, linetype = "dashed", color = "#E63946", linewidth = 1) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = 2,
    alpha = 0.7,
    scale = 1.8,
    rel_min_height = 0.01
  ) +
  annotate("text", x = lemma_val + 0.02, y = 14.5, label = sprintf("Mean selected-set rho = %.2f", lemma_val),
           color = "#E63946", fontface = "bold", angle = 90, vjust = 0, size = 3.5) +
  scale_fill_viridis_c(option = "mako", name = expression(hat(rho)[CH])) +
  coord_cartesian(xlim = c(-0.2, 0.9)) +
  labs(
    title = "Co-heterogeneity Distributions vs. Population Target",
    subtitle = "Confounder DGP (paper eq. 7); coheterogeneity set by D_ov",
    x = expression(hat(rho)[CH]^"(n)"),
    y = "Sample size (n)"
  ) +
  theme_minimal()

print(dist_plot_with_lemma)
