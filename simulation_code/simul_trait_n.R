#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#######################################################################
#######################################################################
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
#######################################################################

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

#######################################################################
# ARCHITECTURE
#######################################################################

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

#######################################################################
# HELPER FUNCTIONS
#######################################################################

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

#######################################################################
# ESTIMATOR
#######################################################################

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

#######################################################################
# INSTRUMENT SELECTION
#######################################################################

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

#######################################################################
# ORACLE TARGET: TRUE correlation among SELECTED SNPs
#######################################################################

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

#######################################################################
# SIMULATE OUTCOMES
#######################################################################

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

#######################################################################
# SINGLE REPLICATION
#######################################################################

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

#######################################################################
# RUN FOR SINGLE SAMPLE SIZE
#######################################################################

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

#######################################################################
# MAIN STUDY
#######################################################################

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

#######################################################################
# CREATE THE TWO PLOTS (rounded N labels WITHOUT changing axis order)
#######################################################################

#######################################################################
# CREATE TWO READABLE PLOTS (coverage + ridge distribution)
# - Rounded N labels (e.g., 62,345 -> 62K)
# - Preserve correct axis ordering by numeric N_x
# - Coverage plot: fewer x labels + rotated + label only outliers
# - Ridge plot: fewer y labels + less overlap + zoomed x-range
#######################################################################
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

#######################################################################
# EXECUTE
#######################################################################

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


############################################################

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

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

#######################################################################
#######################################################################
# PUBLICATION SIMULATION FOR WEIGHTED CH CORRELATION ESTIMATOR
#
# Outputs:
#   Plot 2: Weighted estimator and weighted oracle ridge plots.
#   Plot 3: Weighted mean trajectories.
#   Plot 4A: Unrestricted-beta ridge plot for weighted estimator.
#   Plot 4B: Unrestricted-beta ridge plot for weighted oracle.
#   Plot 4 combined: Plot 4A and 4B stacked with shared legend.
#   Plot 5: Oracle narrowness using IQR.
#
# Notes:
#   - Plot headings are removed for manuscript captions.
#   - N labels are rounded to nearest 10k by default.
#   - Set N_round_base <- 25000 for nearest 25k.
#######################################################################
#######################################################################

rm(list = setdiff(
  ls(),
  c("gwas_cov_by_overlap", "compute_sigma_terms", ".cohet_pkg",
    "rho_hat_T", "estimate_rho_with_se")
))

suppressPackageStartupMessages({
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggridges)
  library(patchwork)
})

select <- dplyr::select
filter <- dplyr::filter

eps <- 1e-12

#######################################################################
# GLOBAL CONTROLS
#######################################################################

RUN_SIMULATION <- TRUE

N_rep_main <- 200

M0_main <- 1e5
N0_main <- 1e5
M_gamma_main <- 0.5

# 10000 = nearest 10k; 25000 = nearest 25k.
N_round_base <- 10000

scenario_labels <- c(
  all    = "Unrestricted true betas",
  strong = "Strong instruments"
)

#######################################################################
# HELPERS
#######################################################################

clip_corr <- function(x) {
  ifelse(is.finite(x), pmax(-1, pmin(1, x)), NA_real_)
}

make_M_N <- function(N_x, M0 = 1e5, N0 = 1e5, gamma = 0.5, max_M = Inf) {
  M <- ceiling(M0 * (N_x / N0)^gamma)
  M <- max(10L, as.integer(M))
  M <- min(M, max_M)
  as.integer(M)
}

safe_mean <- function(x) {
  if (all(!is.finite(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(is.finite(x)) < 2) return(NA_real_)
  sd(x, na.rm = TRUE)
}

round_to_nearest <- function(x, base = N_round_base) {
  as.integer(base * round(x / base))
}

make_N_factor <- function(N_x, base = N_round_base) {
  rounded <- round_to_nearest(N_x, base)
  levels_rounded <- sort(unique(rounded))
  
  factor(
    paste0("N=", label_number(scale_cut = cut_short_scale())(rounded)),
    levels = paste0(
      "N=",
      label_number(scale_cut = cut_short_scale())(levels_rounded)
    )
  )
}

rho_legend_labels <- c(
  expression(hat(rho)[CH]),
  expression(rho[CH])
)

#######################################################################
# UNWEIGHTED GENOME-WIDE ORACLE FORMULA
#######################################################################

lemma_formula <- function(D_ov,
                          pi_inv,
                          r1,
                          r2,
                          theta_U1,
                          theta_U2,
                          sigma2_eta1,
                          sigma2_eta2,
                          sigma2_UR) {
  num <- pi_inv * D_ov * theta_U1 * theta_U2 * sigma2_UR
  
  d1 <- sigma2_eta1 +
    pi_inv * (D_ov + (1 - D_ov) * r1) * theta_U1^2 * sigma2_UR
  
  d2 <- sigma2_eta2 +
    pi_inv * (D_ov + (1 - D_ov) * r2) * theta_U2^2 * sigma2_UR
  
  clip_corr(num / sqrt(d1 * d2))
}

#######################################################################
# DGP
#######################################################################

generate_lemma_arch <- function(M,
                                pi_inv = 0.5,
                                D_ov = 0.5,
                                r1 = 0.5,
                                r2 = 0.5,
                                allow_other_invalid = FALSE,
                                theta_U1 = 0.3,
                                theta_U2 = 0.3,
                                sigma2_x = 5e-5,
                                sigma2_U = 1e-4,
                                sigma2_eta1 = 1e-4,
                                sigma2_eta2 = 1e-4,
                                theta_causal = c(0.2, 0.3),
                                scenario = c("all", "strong"),
                                strong_c_mult = 2,
                                seed = 999) {
  scenario <- match.arg(scenario)
  set.seed(seed)
  
  if (!allow_other_invalid && abs(r1 + r2 - 1) > 1e-8) {
    stop("For the default DGP, require r1 + r2 = 1.")
  }
  
  if (D_ov < 0 || D_ov > 1) stop("D_ov must be in [0, 1].")
  if (pi_inv < 0 || pi_inv > 1) stop("pi_inv must be in [0, 1].")
  if (r1 < 0 || r2 < 0) stop("r1 and r2 must be nonnegative.")
  
  M <- as.integer(M)
  
  n_inv <- round(M * pi_inv)
  inv_snps <- if (n_inv > 0) sample.int(M, n_inv) else integer(0)
  
  p_U0 <- D_ov
  p_U1 <- (1 - D_ov) * r1
  p_U2 <- (1 - D_ov) * r2
  p_other <- max(0, 1 - p_U0 - p_U1 - p_U2)
  
  if (!allow_other_invalid && p_other > 1e-8) {
    stop("Conditional invalid-pathway probabilities do not sum to 1.")
  }
  
  probs <- c(U0 = p_U0, U1 = p_U1, U2 = p_U2, other = p_other)
  probs <- probs / sum(probs)
  
  counts <- if (n_inv > 0) {
    as.integer(rmultinom(1, n_inv, probs))
  } else {
    rep(0L, 4)
  }
  names(counts) <- names(probs)
  
  shuf <- if (n_inv > 0) sample(inv_snps) else integer(0)
  cuts <- cumsum(counts)
  
  U0 <- if (counts["U0"] > 0) {
    shuf[1:cuts["U0"]]
  } else {
    integer(0)
  }
  
  U1 <- if (counts["U1"] > 0) {
    shuf[(cuts["U0"] + 1):cuts["U1"]]
  } else {
    integer(0)
  }
  
  U2 <- if (counts["U2"] > 0) {
    shuf[(cuts["U1"] + 1):cuts["U2"]]
  } else {
    integer(0)
  }
  
  U_other <- if (counts["other"] > 0) {
    shuf[(cuts["U2"] + 1):cuts["other"]]
  } else {
    integer(0)
  }
  
  pathway <- rep("valid", M)
  pathway[U0] <- "U0_shared"
  pathway[U1] <- "U1_Y1_specific"
  pathway[U2] <- "U2_Y2_specific"
  pathway[U_other] <- "other_invalid"
  
  if (scenario == "all") {
    beta_x <- rnorm(M, 0, sqrt(sigma2_x))
  }
  
  if (scenario == "strong") {
    c_thresh <- strong_c_mult * sqrt(sigma2_x)
    mags <- c_thresh + abs(rnorm(M, 0, sqrt(sigma2_x)))
    signs <- sample(c(-1, 1), M, replace = TRUE)
    beta_x <- signs * mags
  }
  
  sigma2_UR <- sigma2_U / sigma2_x
  
  psi_0 <- rep(0, M)
  psi_1 <- rep(0, M)
  psi_2 <- rep(0, M)
  
  if (length(U0) > 0) {
    psi_0[U0] <- rnorm(length(U0), 0, sqrt(sigma2_UR))
  }
  
  if (length(U1) > 0) {
    psi_1[U1] <- rnorm(length(U1), 0, sqrt(sigma2_UR))
  }
  
  if (length(U2) > 0) {
    psi_2[U2] <- rnorm(length(U2), 0, sqrt(sigma2_UR))
  }
  
  eta_1 <- rnorm(M, 0, sqrt(sigma2_eta1))
  eta_2 <- rnorm(M, 0, sqrt(sigma2_eta2))
  
  alpha_R_1 <- eta_1 + theta_U1 * psi_0 + theta_U1 * psi_1
  alpha_R_2 <- eta_2 + theta_U2 * psi_0 + theta_U2 * psi_2
  
  betay1_true <- (theta_causal[1] + alpha_R_1) * beta_x
  betay2_true <- (theta_causal[2] + alpha_R_2) * beta_x
  
  list(
    M = M,
    scenario = scenario,
    betax_true = beta_x,
    betay1_true = betay1_true,
    betay2_true = betay2_true,
    alpha_R_1 = alpha_R_1,
    alpha_R_2 = alpha_R_2,
    pathway = pathway,
    U0 = U0,
    U1 = U1,
    U2 = U2,
    U_other = U_other,
    pi_inv = pi_inv,
    D_ov = D_ov,
    r1 = r1,
    r2 = r2,
    theta_U1 = theta_U1,
    theta_U2 = theta_U2,
    sigma2_UR = sigma2_UR,
    sigma2_eta1 = sigma2_eta1,
    sigma2_eta2 = sigma2_eta2,
    theta_causal = theta_causal,
    sigma2_x = sigma2_x,
    sigma2_U = sigma2_U,
    strong_c = if (scenario == "strong") {
      strong_c_mult * sqrt(sigma2_x)
    } else {
      NA_real_
    }
  )
}

#######################################################################
# WEIGHTED ORACLE HELPERS
#######################################################################

compute_oracle_weights <- function(bx, by1, by2,
                                   N_x, N_y1, N_y2,
                                   overlap = 0,
                                   rho_e = 0) {
  var_bx <- 1 / N_x
  var_by1 <- 1 / N_y1
  var_by2 <- 1 / N_y2
  cov_by12 <- gwas_cov_by_overlap(N_y1, N_y2, overlap, rho_e)
  
  sig <- compute_sigma_terms(
    bx = bx,
    by1 = by1,
    by2 = by2,
    var_bx = var_bx,
    var_by1 = var_by1,
    var_by2 = var_by2,
    cov_by12 = cov_by12
  )
  
  ok <- complete.cases(sig$s1, sig$s2, sig$s12) &
    sig$s1 > 0 &
    sig$s2 > 0
  
  if (sum(ok) < 3) return(NULL)
  
  w_raw <- 1 / sqrt((sig$s1[ok] + eps) * (sig$s2[ok] + eps))
  w <- w_raw / sum(w_raw)
  
  list(
    ok = ok,
    w = w,
    s1 = sig$s1[ok],
    s2 = sig$s2[ok],
    s12 = sig$s12[ok]
  )
}

weighted_selected_target <- function(arch, idx,
                                     N_x, N_y1, N_y2,
                                     overlap = 0,
                                     rho_e = 0) {
  if (is.null(idx) || length(idx) < 3) return(NA_real_)
  
  bx_t <- arch$betax_true[idx]
  by1_t <- arch$betay1_true[idx]
  by2_t <- arch$betay2_true[idx]
  
  ok0 <- complete.cases(bx_t, by1_t, by2_t) &
    is.finite(bx_t) &
    abs(bx_t) > eps
  
  if (sum(ok0) < 3) return(NA_real_)
  
  bx_t <- bx_t[ok0]
  by1_t <- by1_t[ok0]
  by2_t <- by2_t[ok0]
  
  th1 <- by1_t / bx_t
  th2 <- by2_t / bx_t
  
  wobj <- compute_oracle_weights(
    bx = bx_t,
    by1 = by1_t,
    by2 = by2_t,
    N_x = N_x,
    N_y1 = N_y1,
    N_y2 = N_y2,
    overlap = overlap,
    rho_e = rho_e
  )
  
  if (is.null(wobj)) return(NA_real_)
  
  ok <- wobj$ok
  
  th1 <- th1[ok]
  th2 <- th2[ok]
  w <- wobj$w
  
  th1b <- sum(w * th1)
  th2b <- sum(w * th2)
  
  d1 <- th1 - th1b
  d2 <- th2 - th2b
  
  S11 <- sum(w * d1^2)
  S22 <- sum(w * d2^2)
  S12 <- sum(w * d1 * d2)
  
  if (!is.finite(S11) || !is.finite(S22) || S11 <= eps || S22 <= eps) {
    return(NA_real_)
  }
  
  clip_corr(S12 / sqrt(S11 * S22))
}

#######################################################################
# SELECTION AND OUTCOME SIMULATION
#######################################################################

select_instruments <- function(arch,
                               N_x,
                               p_threshold = 5e-8,
                               F_threshold = 10,
                               use_p = TRUE,
                               use_F = TRUE,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  M <- arch$M
  
  bx_hat <- arch$betax_true + rnorm(M, 0, sqrt(1 / N_x))
  z <- sqrt(N_x) * bx_hat
  pval <- 2 * pnorm(-abs(z))
  F_obs <- N_x * bx_hat^2
  
  keep <- rep(TRUE, M)
  
  if (use_p) keep <- keep & (pval < p_threshold)
  if (use_F) keep <- keep & (F_obs > F_threshold)
  
  idx <- which(keep)
  
  if (length(idx) == 0) {
    return(list(
      idx = NULL,
      K = 0L,
      bx_hat = bx_hat,
      pval = pval,
      F_obs = F_obs,
      z = z
    ))
  }
  
  idx <- idx[order(pval[idx])]
  
  list(
    idx = idx,
    K = length(idx),
    bx_hat = bx_hat,
    pval = pval,
    F_obs = F_obs,
    z = z
  )
}

simulate_outcomes <- function(arch, idx,
                              N_y1, N_y2,
                              overlap = 0,
                              rho_e = 0) {
  K <- length(idx)
  
  if (K == 0) return(NULL)
  
  by1_t <- arch$betay1_true[idx]
  by2_t <- arch$betay2_true[idx]
  
  cov12 <- gwas_cov_by_overlap(N_y1, N_y2, overlap, rho_e)
  Sigma <- matrix(c(1 / N_y1, cov12, cov12, 1 / N_y2), 2, 2)
  
  E <- mvrnorm(K, c(0, 0), Sigma)
  
  list(
    by1 = by1_t + E[, 1],
    by2 = by2_t + E[, 2]
  )
}

simulate_one_rep <- function(arch,
                             N_x,
                             N_y1,
                             N_y2,
                             overlap = 0,
                             rho_e = 0,
                             p_threshold = 5e-8,
                             F_threshold = 10,
                             use_p = TRUE,
                             use_F = TRUE,
                             min_K = 5,
                             seed = NULL) {
  sel <- select_instruments(
    arch = arch,
    N_x = N_x,
    p_threshold = p_threshold,
    F_threshold = F_threshold,
    use_p = use_p,
    use_F = use_F,
    seed = seed
  )
  
  rho_hat <- NA_real_
  rho_weighted_selected <- NA_real_
  
  if (sel$K >= min_K) {
    bx_hat_sel <- sel$bx_hat[sel$idx]
    
    out <- simulate_outcomes(
      arch = arch,
      idx = sel$idx,
      N_y1 = N_y1,
      N_y2 = N_y2,
      overlap = overlap,
      rho_e = rho_e
    )
    
    rho_hat <- rho_hat_T(
      bx = bx_hat_sel,
      by1 = out$by1,
      by2 = out$by2,
      N_x = N_x,
      N_y1 = N_y1,
      N_y2 = N_y2,
      overlap = overlap,
      rho_e = rho_e
    )
    
    rho_weighted_selected <- weighted_selected_target(
      arch = arch,
      idx = sel$idx,
      N_x = N_x,
      N_y1 = N_y1,
      N_y2 = N_y2,
      overlap = overlap,
      rho_e = rho_e
    )
  }
  
  data.frame(
    rho_hat = rho_hat,
    rho_weighted_selected = rho_weighted_selected,
    K = sel$K
  )
}

#######################################################################
# RUNNER
#######################################################################

run_vary_N <- function(scenario,
                       N_grid = round(exp(seq(log(2e4), log(5e5), length.out = 10))),
                       D_ov = 0.5,
                       M0 = 1e5,
                       N0 = 1e5,
                       M_gamma = 0.5,
                       max_M = Inf,
                       N_rep = 200,
                       pi_inv = 0.5,
                       r1 = 0.5,
                       r2 = 0.5,
                       theta_U1 = 0.3,
                       theta_U2 = 0.3,
                       sigma2_x = 5e-5,
                       sigma2_U = 1e-4,
                       sigma2_eta1 = 1e-4,
                       sigma2_eta2 = 1e-4,
                       overlap_y = 0.5,
                       rho_e = 0.2,
                       strong_c_mult = 2,
                       p_threshold = 5e-8,
                       F_threshold = 10,
                       use_p = TRUE,
                       use_F = TRUE,
                       base_seed = 123400) {
  out <- list()
  
  for (i in seq_along(N_grid)) {
    N_x <- N_grid[i]
    M_N <- make_M_N(
      N_x,
      M0 = M0,
      N0 = N0,
      gamma = M_gamma,
      max_M = max_M
    )
    
    cat(sprintf(
      "[%s | N %d/%d] N = %s, M_N = %s\n",
      scenario,
      i,
      length(N_grid),
      format(N_x, big.mark = ","),
      format(M_N, big.mark = ",")
    ))
    
    arch <- generate_lemma_arch(
      M = M_N,
      pi_inv = pi_inv,
      D_ov = D_ov,
      r1 = r1,
      r2 = r2,
      theta_U1 = theta_U1,
      theta_U2 = theta_U2,
      sigma2_x = sigma2_x,
      sigma2_U = sigma2_U,
      sigma2_eta1 = sigma2_eta1,
      sigma2_eta2 = sigma2_eta2,
      scenario = scenario,
      strong_c_mult = strong_c_mult,
      seed = base_seed + i
    )
    
    pop_formula <- lemma_formula(
      D_ov = D_ov,
      pi_inv = pi_inv,
      r1 = r1,
      r2 = r2,
      theta_U1 = theta_U1,
      theta_U2 = theta_U2,
      sigma2_eta1 = sigma2_eta1,
      sigma2_eta2 = sigma2_eta2,
      sigma2_UR = arch$sigma2_UR
    )
    
    reps <- vector("list", N_rep)
    
    for (r in seq_len(N_rep)) {
      reps[[r]] <- simulate_one_rep(
        arch = arch,
        N_x = N_x,
        N_y1 = N_x,
        N_y2 = N_x,
        overlap = overlap_y,
        rho_e = rho_e,
        p_threshold = p_threshold,
        F_threshold = F_threshold,
        use_p = use_p,
        use_F = use_F,
        seed = 100000 * i + r
      )
    }
    
    df <- bind_rows(reps) %>%
      mutate(
        rep = seq_len(N_rep),
        scenario = scenario,
        D_ov = D_ov,
        N_x = N_x,
        M = M_N,
        M_gamma = M_gamma,
        pop_formula = pop_formula
      )
    
    out[[i]] <- df
  }
  
  bind_rows(out)
}

#######################################################################
# PUBLICATION PLOTS
#######################################################################

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title = element_text(face = "bold")
    )
}

make_pub_plot2_weighted_ridges <- function(d, N_round_base = 10000) {
  long <- d %>%
    pivot_longer(
      cols = c(rho_hat, rho_weighted_selected),
      names_to = "kind",
      values_to = "rho"
    ) %>%
    mutate(
      kind = recode(
        kind,
        rho_hat = "hat(rho)[CH]",
        rho_weighted_selected = "rho[CH]"
      ),
      scenario_label = scenario_labels[scenario],
      N_label = make_N_factor(N_x, N_round_base)
    ) %>%
    filter(is.finite(rho))
  
  formula_df <- d %>%
    distinct(scenario, pop_formula) %>%
    mutate(scenario_label = scenario_labels[scenario]) %>%
    filter(is.finite(pop_formula))
  
  ggplot(long, aes(x = rho, y = N_label, fill = kind)) +
    geom_density_ridges(
      alpha = 0.65,
      scale = 1.35,
      rel_min_height = 0.01,
      color = "white",
      linewidth = 0.25
    ) +
    geom_vline(
      data = formula_df,
      aes(xintercept = pop_formula, linetype = "Genome-wide oracle"),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.85
    ) +
    facet_grid(
      kind ~ scenario_label,
      labeller = labeller(kind = label_parsed)
    ) +
    scale_linetype_manual(values = c("Genome-wide oracle" = "dashed")) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = expression(rho),
      y = "Sample size (n)",
      linetype = NULL
    ) +
    theme_pub(base_size = 13)
}

make_pub_plot3_weighted_trajectories <- function(d, N_round_base = 10000) {
  summ <- d %>%
    filter(
      is.finite(rho_hat),
      is.finite(rho_weighted_selected)
    ) %>%
    group_by(scenario, N_x, M) %>%
    summarise(
      mean_hat = safe_mean(rho_hat),
      sd_hat = safe_sd(rho_hat),
      mean_weighted = safe_mean(rho_weighted_selected),
      sd_weighted = safe_sd(rho_weighted_selected),
      pop_formula = first(pop_formula),
      .groups = "drop"
    ) %>%
    mutate(
      scenario_label = scenario_labels[scenario],
      N_round = round_to_nearest(N_x, N_round_base)
    )
  
  long_mean <- bind_rows(
    summ %>%
      transmute(
        scenario,
        scenario_label,
        N_x,
        N_round,
        M,
        target = "Estimator",
        mean = mean_hat,
        sd = sd_hat
      ),
    summ %>%
      transmute(
        scenario,
        scenario_label,
        N_x,
        N_round,
        M,
        target = "Oracle",
        mean = mean_weighted,
        sd = sd_weighted
      )
  )
  
  formula_df <- summ %>%
    distinct(scenario, scenario_label, pop_formula) %>%
    filter(is.finite(pop_formula))
  
  ggplot(long_mean, aes(x = N_round, y = mean, color = target, fill = target)) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd),
      alpha = 0.16,
      color = NA
    ) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 2.25) +
    geom_hline(
      data = formula_df,
      aes(yintercept = pop_formula, linetype = "Genome-wide oracle"),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.85
    ) +
    facet_wrap(~ scenario_label) +
    scale_x_log10(
      breaks = sort(unique(long_mean$N_round)),
      labels = label_number(scale_cut = cut_short_scale())
    ) +
    scale_color_discrete(
      breaks = c("Estimator", "Oracle"),
      labels = rho_legend_labels
    ) +
    scale_fill_discrete(
      breaks = c("Estimator", "Oracle"),
      labels = rho_legend_labels
    ) +
    scale_linetype_manual(values = c("Genome-wide oracle" = "dashed")) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Sample Size (N)",
      y = expression(rho),
      color = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    theme_pub(base_size = 13)
}

make_single_unrestricted_gradient_ridge <- function(
    d,
    quantity = c("rho_hat", "rho_weighted_selected"),
    N_round_base = 10000
) {
  quantity <- match.arg(quantity)
  
  x_lab <- if (quantity == "rho_hat") {
    expression(hat(rho)[CH])
  } else {
    expression(rho[CH])
  }
  
  d_all <- d %>%
    filter(scenario == "all") %>%
    filter(is.finite(.data[[quantity]])) %>%
    mutate(
      N_round = round_to_nearest(N_x, N_round_base),
      N_label = make_N_factor(N_x, N_round_base)
    )
  
  formula_df <- d_all %>%
    distinct(pop_formula) %>%
    filter(is.finite(pop_formula))
  
  ggplot(
    d_all,
    aes(
      x = .data[[quantity]],
      y = N_label,
      fill = after_stat(x)
    )
  ) +
    geom_density_ridges_gradient(
      scale = 1.45,
      rel_min_height = 0.01,
      color = "#17324D",
      linewidth = 0.35,
      alpha = 0.95
    ) +
    geom_vline(
      data = formula_df,
      aes(xintercept = pop_formula, linetype = "Genome-wide oracle"),
      inherit.aes = FALSE,
      color = "#B22222",
      linewidth = 0.85
    ) +
    scale_fill_viridis_c(
      option = "D",
      limits = c(-1, 1),
      name = expression(rho)
    ) +
    scale_linetype_manual(
      values = c("Genome-wide oracle" = "dashed")
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      linetype = guide_legend(order = 2)
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = x_lab,
      y = "Sample size (n)",
      linetype = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.35),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

make_pub_plot4_unrestricted_two_ridges <- function(d, N_round_base = 10000) {
  p4_est <- make_single_unrestricted_gradient_ridge(
    d = d,
    quantity = "rho_hat",
    N_round_base = N_round_base
  )
  
  p4_oracle <- make_single_unrestricted_gradient_ridge(
    d = d,
    quantity = "rho_weighted_selected",
    N_round_base = N_round_base
  )
  
  p4_combined <- (p4_est / p4_oracle) +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )
    )
  
  list(
    estimator = p4_est,
    oracle = p4_oracle,
    combined = p4_combined
  )
}

make_pub_plot5_oracle_narrowness <- function(d, N_round_base = 10000) {
  width_df <- d %>%
    filter(scenario == "all") %>%
    mutate(
      N_round = round_to_nearest(N_x, N_round_base)
    ) %>%
    pivot_longer(
      cols = c(rho_hat, rho_weighted_selected),
      names_to = "quantity",
      values_to = "rho"
    ) %>%
    filter(is.finite(rho)) %>%
    mutate(
      quantity = recode(
        quantity,
        rho_hat = "Estimator",
        rho_weighted_selected = "Oracle"
      )
    ) %>%
    group_by(N_round, quantity) %>%
    summarise(
      width_iqr = IQR(rho, na.rm = TRUE),
      width_sd = sd(rho, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      width_iqr_plot = pmax(width_iqr, eps)
    )
  
  ggplot(width_df, aes(x = N_round, y = width_iqr_plot, color = quantity)) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 2.3) +
    scale_x_log10(
      breaks = sort(unique(width_df$N_round)),
      labels = label_number(scale_cut = cut_short_scale())
    ) +
    scale_y_log10() +
    scale_color_discrete(
      breaks = c("Estimator", "Oracle"),
      labels = rho_legend_labels
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Sample Size (N)",
      y = "Distribution width: IQR",
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.title = element_text(face = "bold")
    )
}

#######################################################################
# MAIN RUN
#######################################################################

if (RUN_SIMULATION) {
  cat("\n############################################################\n")
  cat("# PUBLICATION RUN: WEIGHTED PLOTS ONLY\n")
  cat("############################################################\n\n")
  
  scenarios_to_run <- c("all", "strong")
  
  N_grid <- round(exp(seq(log(2e4), log(5e5), length.out = 10)))
  
  cat("\n=== Vary N with M_N growing slower than N ===\n")
  cat(sprintf(
    "Using M_N = M0 * (N / N0)^gamma with gamma = %.2f\n",
    M_gamma_main
  ))
  
  d2 <- bind_rows(lapply(scenarios_to_run, function(sc) {
    run_vary_N(
      scenario = sc,
      N_grid = N_grid,
      D_ov = 0.5,
      M0 = M0_main,
      N0 = N0_main,
      M_gamma = M_gamma_main,
      N_rep = N_rep_main,
      pi_inv = 0.5,
      r1 = 0.5,
      r2 = 0.5,
      overlap_y = 0.5,
      rho_e = 0.2,
      strong_c_mult = 2,
      p_threshold = 5e-8,
      F_threshold = 10,
      use_p = TRUE,
      use_F = TRUE
    )
  }))
  
  p2_pub <- make_pub_plot2_weighted_ridges(
    d2,
    N_round_base = N_round_base
  )
  
  p3_pub <- make_pub_plot3_weighted_trajectories(
    d2,
    N_round_base = N_round_base
  )
  
  p4_list <- make_pub_plot4_unrestricted_two_ridges(
    d2,
    N_round_base = N_round_base
  )
  
  p4_estimator <- p4_list$estimator
  p4_oracle <- p4_list$oracle
  p4_combined <- p4_list$combined
  
  p5_narrowness <- make_pub_plot5_oracle_narrowness(
    d2,
    N_round_base = N_round_base
  )
  
  print(p2_pub)
  print(p3_pub)
  print(p4_estimator)
  print(p4_oracle)
  print(p4_combined)
  print(p5_narrowness)

  cat("\n=== Done: saved publication plots ===\n")
}
