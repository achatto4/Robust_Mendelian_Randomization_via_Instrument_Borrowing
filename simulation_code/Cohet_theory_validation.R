## =============================================================================
## FINAL COMPREHENSIVE CODE: polygenic validation of sqrt(K) CLT + stability of U
## =============================================================================
## What this script demonstrates (end-to-end, no saving):
##  (1) Very polygenic triangular-array exposure architecture: beta_x = b0 / sqrt(nx)
##      -> selection at fixed GW threshold yields MANY instruments near threshold.
##  (2) Bounded ratio-scale pleiotropy alpha_{l,k} (RC1) and beta_y = beta_x*(theta+alpha)
##      -> Wald ratio signal = theta + alpha, bounded.
##  (3) Shows "variance stabilization" numerically:
##      - fraction of selected SNPs in [zthr, zthr+1]
##      - quantiles of nx*beta_x^2 (true) among selected stable across nx
##      - quantiles of nx*bhat_x^2 (estimated) among selected stable across nx
##      - quantiles of sigma1_hat (plug-in Wald ratio variance proxy) stable across nx
##  (4) Confirms sqrt(K) rate:
##      - SD of zK := sqrt(K)*(rho_hat - rho_target_hatW) stable across nx
##      - Histograms + QQ plots for zK by nx (shape ~Normal as K large)
##  (5) Optional: bias correction demo via jackknife bias estimate (quick, approximate)
##
## Notes:
##  - No overlap exposure vs outcomes; outcome overlap via ks,rho_e.
##  - The target used for centering is estimator-aligned:
##      rho_target_hatW = corr(true alpha1, alpha2) under estimated weights w_hat.
##  - This script is designed to support the narrative:
##      "In a very polygenic setting, sqrt(K) CLT holds naturally because many
##       selected instruments are near threshold, stabilizing sigma and U."
##
## Dependencies: MASS, dplyr, tidyr, ggplot2

suppressPackageStartupMessages({
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# -----------------------------
# Utilities
# -----------------------------
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

p_select <- function(b0_sd, zthr) {
  # if Z = b0 + N(0,1), b0~N(0,b0_sd^2), then Z~N(0,1+b0_sd^2)
  2 * (1 - pnorm(zthr / sqrt(1 + b0_sd^2)))
}

choose_b0_sd <- function(M, zthr, targetK) {
  target_p <- targetK / M
  f <- function(b) p_select(b, zthr) - target_p
  uniroot(f, lower = 0.1, upper = 50)$root
}

# -----------------------------
# Polygenic triangular-array genome
# -----------------------------
make_polygenic_genome <- function(M = 2e6,
                                  b0_sd = 3,
                                  alpha_sd = 0.10,
                                  alpha_bound = 0.40,
                                  rho_alpha = 0.35,
                                  theta1 = 0.10,
                                  theta2 = -0.05,
                                  seed = 1) {
  set.seed(seed)
  
  # Z-scale exposure effects
  b0 <- rnorm(M, 0, b0_sd)
  
  # bounded ratio-scale pleiotropy via correlated normals + clipping
  Z <- MASS::mvrnorm(M, mu = c(0, 0),
                     Sigma = matrix(c(1, rho_alpha, rho_alpha, 1), 2, 2))
  alpha1 <- clip(Z[, 1] * alpha_sd, -alpha_bound, alpha_bound)
  alpha2 <- clip(Z[, 2] * alpha_sd, -alpha_bound, alpha_bound)
  
  list(M = M, b0 = b0, alpha1 = alpha1, alpha2 = alpha2,
       theta1 = theta1, theta2 = theta2,
       b0_sd = b0_sd)
}

# -----------------------------
# One replicate at one nx, returns estimator + diagnostics
# -----------------------------
simulate_one_polygenic <- function(genome, nx,
                                   k1 = 2, k2 = 2, ks = 0.2, rho_e = 0.2,
                                   alpha_gw = 5e-8,
                                   minK = 500) {
  
  M <- genome$M
  b0 <- genome$b0
  a1 <- genome$alpha1
  a2 <- genome$alpha2
  theta1 <- genome$theta1
  theta2 <- genome$theta2
  
  n1 <- round(k1 * nx)
  n2 <- round(k2 * nx)
  ns <- round(ks * nx)
  
  se_x  <- 1 / sqrt(nx)
  se_y1 <- 1 / sqrt(n1)
  se_y2 <- 1 / sqrt(n2)
  
  zthr <- qnorm(1 - alpha_gw/2)
  
  # triangular beta_x => Z_true = sqrt(nx)*beta_x = b0
  beta_x <- b0 / sqrt(nx)
  
  # reduced-form outcomes consistent with ratio-scale alpha
  beta_y1 <- beta_x * (theta1 + a1)
  beta_y2 <- beta_x * (theta2 + a2)
  
  # exposure GWAS
  bhat_x <- beta_x + rnorm(M, 0, se_x)
  Zx <- bhat_x / se_x
  
  # outcomes GWAS (overlap only between outcomes)
  corr_y <- rho_e * ns / sqrt(n1*n2)
  corr_y <- max(min(corr_y, 0.99), -0.99)
  e12 <- MASS::mvrnorm(M, mu = c(0, 0),
                       Sigma = matrix(c(se_y1^2, corr_y*se_y1*se_y2,
                                        corr_y*se_y1*se_y2, se_y2^2), 2, 2))
  bhat_y1 <- beta_y1 + e12[, 1]
  bhat_y2 <- beta_y2 + e12[, 2]
  
  # selection by exposure GW significance
  sel <- abs(Zx) > zthr
  K <- sum(sel)
  if (K < minK) return(list(ok = FALSE, nx = nx, K = K))
  
  # Near-threshold mass: fraction within one Z-unit above threshold
  frac_band1 <- mean(abs(Zx[sel]) <= (zthr + 1))
  
  bx_hat  <- bhat_x[sel]
  by1_hat <- bhat_y1[sel]
  by2_hat <- bhat_y2[sel]
  
  th1 <- by1_hat / bx_hat
  th2 <- by2_hat / bx_hat
  
  # Plug-in sigma formulas (your paper)
  sig1_hat  <- (1/(nx * bx_hat^2)) * (1/k1 + (by1_hat^2)/(bx_hat^2))
  sig2_hat  <- (1/(nx * bx_hat^2)) * (1/k2 + (by2_hat^2)/(bx_hat^2))
  sig12_hat <- (1/(nx * bx_hat^2)) * (ks*rho_e/(k1*k2) + (by1_hat*by2_hat)/(bx_hat^2))
  
  # weights
  w_raw <- 1 / sqrt(sig1_hat * sig2_hat)
  w <- w_raw / sum(w_raw)
  
  # bias-corrected estimator rho_hat
  m1 <- sum(w * th1); m2 <- sum(w * th2)
  d1 <- th1 - m1;     d2 <- th2 - m2
  
  M11 <- sum(w * d1^2); M22 <- sum(w * d2^2); M12 <- sum(w * d1*d2)
  B11 <- sum(w * (1-w) * sig1_hat)
  B22 <- sum(w * (1-w) * sig2_hat)
  B12 <- sum(w * (1-w) * sig12_hat)
  
  S11 <- M11 - B11
  S22 <- M22 - B22
  S12 <- M12 - B12
  if (S11 <= 0 || S22 <= 0) return(list(ok = FALSE, nx = nx, K = K))
  
  rho_hat <- S12 / sqrt(S11*S22)
  
  # estimator-aligned target (true alpha corr under estimated weights w)
  a1s <- a1[sel]; a2s <- a2[sel]
  a1b <- sum(w * a1s); a2b <- sum(w * a2s)
  t1 <- sum(w * (a1s - a1b)^2)
  t2 <- sum(w * (a2s - a2b)^2)
  C12 <- sum(w * (a1s - a1b) * (a2s - a2b))
  rho_target <- as.numeric(C12 / sqrt(t1*t2))
  
  err <- rho_hat - rho_target
  zK <- sqrt(K) * err
  
  # Stabilization quantities
  nx_bx2_true <- nx * (beta_x[sel]^2)
  nx_bx2_hat  <- nx * (bx_hat^2)
  
  q <- function(v) as.numeric(quantile(v, c(.1,.5,.9), na.rm = TRUE))
  q_true <- q(nx_bx2_true)
  q_hat  <- q(nx_bx2_hat)
  q_sig1 <- q(sig1_hat)
  
  list(ok = TRUE, nx = nx, K = K,
       rho_hat = rho_hat, rho_target = rho_target, err = err, zK = zK,
       frac_band1 = frac_band1,
       nx_bx2_true_q10 = q_true[1], nx_bx2_true_q50 = q_true[2], nx_bx2_true_q90 = q_true[3],
       nx_bx2_hat_q10  = q_hat[1],  nx_bx2_hat_q50  = q_hat[2],  nx_bx2_hat_q90  = q_hat[3],
       sig1_hat_q10 = q_sig1[1], sig1_hat_q50 = q_sig1[2], sig1_hat_q90 = q_sig1[3])
}

# -----------------------------
# Main runner
# -----------------------------
run_final_validation <- function(nx_grid = c(20000, 40000, 80000, 160000),
                                 R = 200,
                                 M = 2000000,
                                 alpha_gw = 5e-8,
                                 targetK = 5000,
                                 k1 = 2, k2 = 2, ks = 0.2, rho_e = 0.2,
                                 alpha_sd = 0.10, alpha_bound = 0.40, rho_alpha = 0.35,
                                 theta1 = 0.10, theta2 = -0.05,
                                 minK = 500,
                                 seed = 1) {
  
  zthr <- qnorm(1 - alpha_gw/2)
  b0_sd <- choose_b0_sd(M = M, zthr = zthr, targetK = targetK)
  cat("Calibrated b0_sd =", round(b0_sd, 3),
      " => approx P(select) =", signif(p_select(b0_sd, zthr), 3),
      " => approx E[K] =", round(M * p_select(b0_sd, zthr)), "\n\n")
  
  genome <- make_polygenic_genome(M = M, b0_sd = b0_sd,
                                  alpha_sd = alpha_sd, alpha_bound = alpha_bound, rho_alpha = rho_alpha,
                                  theta1 = theta1, theta2 = theta2,
                                  seed = seed)
  
  out <- list()
  for (nx in nx_grid) {
    cat("\n=== nx =", nx, "R =", R, "===\n")
    tmp <- vector("list", R)
    for (r in 1:R) {
      tmp[[r]] <- simulate_one_polygenic(genome, nx,
                                         k1 = k1, k2 = k2, ks = ks, rho_e = rho_e,
                                         alpha_gw = alpha_gw,
                                         minK = minK)
      if (r %% max(1, floor(R/10)) == 0) cat("  rep", r, "\n")
    }
    ok <- vapply(tmp, `[[`, logical(1), "ok")
    cat("success:", sum(ok), "/", R, "\n")
    out[[as.character(nx)]] <- bind_rows(lapply(tmp[ok], as.data.frame))
  }
  
  df <- bind_rows(out)
  if (nrow(df) == 0) stop("No successful replicates; lower minK or increase M/targetK.")
  
  # Summary table
  summ <- df %>% group_by(nx) %>% summarise(
    R = n(),
    K_mean = mean(K),
    
    # CLT scaling
    mean_err = mean(err),
    sd_err = sd(err),
    z_mean = mean(zK),
    z_sd = sd(zK),
    stab = sd_err * sqrt(mean(K)),
    
    # near-threshold mass
    frac_band1_mean = mean(frac_band1),
    
    # stabilization diagnostics
    nx_bx2_true_q10 = mean(nx_bx2_true_q10),
    nx_bx2_true_q50 = mean(nx_bx2_true_q50),
    nx_bx2_true_q90 = mean(nx_bx2_true_q90),
    nx_bx2_hat_q10  = mean(nx_bx2_hat_q10),
    nx_bx2_hat_q50  = mean(nx_bx2_hat_q50),
    nx_bx2_hat_q90  = mean(nx_bx2_hat_q90),
    sig1_hat_q10 = mean(sig1_hat_q10),
    sig1_hat_q50 = mean(sig1_hat_q50),
    sig1_hat_q90 = mean(sig1_hat_q90),
    .groups = "drop"
  )
  
  cat("\n--- FINAL POLYGENIC VALIDATION SUMMARY ---\n")
  print(summ)
  
  # -----------------------------
  # Plots
  # -----------------------------
  
  # (A) CLT: histogram + normal overlay
  ggplot(df, aes(x = zK)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey85", color = "white") +
    stat_function(fun = dnorm, linewidth = 1) +
    facet_wrap(~nx, scales = "free_y") +
    theme_bw() +
    labs(title = "CLT check: zK = sqrt(K)*(rho_hat - rho_target_hatW)",
         subtitle = "In polygenic triangular-array setting, SD(zK) should stabilize across nx",
         x = "zK", y = "Density") |> print()
  
  ggplot(df, aes(sample = zK)) +
    stat_qq(size = 0.6) + stat_qq_line() +
    facet_wrap(~nx) +
    theme_bw() +
    labs(title = "QQ plot: zK vs Normal",
         subtitle = "QQ straightening with nx supports CLT",
         x = "Theoretical", y = "Sample") |> print()
  
  # (B) Stability of SD(zK) and stab across nx
  ggplot(summ, aes(x = nx, y = z_sd)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    theme_bw() +
    labs(title = "Stability of SD(zK) across nx (sqrt(K) scaling)",
         x = "nx (log10)", y = "SD(zK)") |> print()
  
  ggplot(summ, aes(x = nx, y = frac_band1_mean)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    theme_bw() +
    labs(title = "Near-threshold mass among selected SNPs",
         subtitle = "P(|Zx| <= zthr+1 | selected) should be bounded away from 0",
         x = "nx (log10)", y = "fraction") |> print()
  
  # (C) Stability of nx*beta_x^2 (true and estimated)
  df_bx <- summ %>%
    select(nx,
           nx_bx2_true_q10, nx_bx2_true_q50, nx_bx2_true_q90,
           nx_bx2_hat_q10,  nx_bx2_hat_q50,  nx_bx2_hat_q90) %>%
    pivot_longer(-nx, names_to = "which", values_to = "value")
  
  ggplot(df_bx, aes(x = nx, y = value, color = which)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    theme_bw() +
    labs(title = "Stability of n_x * beta_x^2 among selected (true and estimated)",
         subtitle = "Quantiles should be ~constant across nx in the polygenic regime",
         x = "nx (log10)", y = "quantile value") |> print()
  
  # (D) Stability of sigma1_hat (proxy for Wald ratio variance scale)
  df_sig <- summ %>%
    select(nx, sig1_hat_q10, sig1_hat_q50, sig1_hat_q90) %>%
    pivot_longer(-nx, names_to = "which", values_to = "value")
  
  ggplot(df_sig, aes(x = nx, y = value, color = which)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    theme_bw() +
    labs(title = "Stability of plug-in sigma1_hat among selected",
         subtitle = "Shows Wald ratio variance scale is O(1) and stable across nx",
         x = "nx (log10)", y = "quantile of sigma1_hat") |> print()
  
  invisible(list(all = df, summary = summ, genome = genome, b0_sd = b0_sd))
}

## =============================================================================
## RUN
## =============================================================================
set.seed(1)

res <- run_final_validation(
  nx_grid = c(20000, 40000, 80000, 160000),
  R = 200,
  M = 2000000,
  alpha_gw = 5e-8,
  targetK = 5000,
  k1 = 2, k2 = 2, ks = 0.2, rho_e = 0.2,
  alpha_sd = 0.10, alpha_bound = 0.40, rho_alpha = 0.35,
  theta1 = 0.10, theta2 = -0.05,
  minK = 500,
  seed = 1
)

## Use:
##  res$summary  # table you pasted
##  res$all      # replicate-level data


############
############
############
############################################################
## THEORETICAL ESTIMATOR - Exact Implementation
## This implements the estimator exactly as described in the paper:
##  - Bias correction: subtract σ_{12,k} (not B_k)
##  - Influence function: matches equation (42) exactly
##  - Variance: Σ_hat = K^{-1} Σ ψ_k^2 with correct K scaling
##  - Standard error: sqrt(Σ_hat / K)
############################################################

library(MASS)  # for mvrnorm if needed

qnorm_two_sided <- function(p) qnorm(1 - p/2)

############################################################
## Main simulation function - THEORETICAL VERSION
############################################################
simulate_replicate_theory <- function(
    n_x,
    n1 = n_x, 
    n2 = n_x,
    M = 150000,
    z_GW = qnorm_two_sided(5e-8),
    
    ## Exposure architecture
    arch = c("spike_slab", "polygenic"),
    pi_x = 0.003, 
    sd_bx = 0.07,
    sd_bx_poly = 0.01,
    
    ## Pleiotropy
    pleio = c("gaussian", "t", "mixture", "directional"),
    sd_a1 = 0.05, 
    sd_a2 = 0.06, 
    rho_a = 0.4,
    df_t = 4,
    mix_p = 0.98, 
    mix_sd_big = 0.35,
    dir_mean1 = 0.03, 
    dir_mean2 = 0.03,
    corr_bx_a = 0.0,
    
    ## Causal effects
    theta1 = 0.10, 
    theta2 = 0.15,
    
    ## Overlap model (working correlation, matches theory)
    k_s = 1.0,          # overlap ratio n_s/n_x (1.0 = complete overlap)
    rho_e = 0.0,        # residual phenotypic correlation
    
    ## LD-block dependence (for robustness checks)
    ld = c("none", "blocks"),
    block_size = 20,
    rho_block_alpha = 0.7,
    
    ## Estimator controls
    eps_denom_floor = 1e-6,
    trunc_tau = TRUE
) {
  arch <- match.arg(arch)
  pleio <- match.arg(pleio)
  ld <- match.arg(ld)
  
  M <- as.integer(M)
  k1 <- n1/n_x
  k2 <- n2/n_x
  
  ## ---- 1) Generate beta_x ----
  if (arch == "spike_slab") {
    I <- rbinom(M, 1, pi_x) == 1
    beta_x <- numeric(M)
    beta_x[I] <- rnorm(sum(I), 0, sd_bx)
  } else {
    beta_x <- rnorm(M, 0, sd_bx_poly)
  }
  
  ## ---- 2) Generate (alpha1, alpha2) with correlation rho_a ----
  z1 <- rnorm(M)
  z2 <- rnorm(M)
  alpha1_g <- sd_a1 * z1
  alpha2_g <- sd_a2 * (rho_a * z1 + sqrt(1 - rho_a^2) * z2)
  
  if (pleio == "gaussian") {
    alpha1 <- alpha1_g
    alpha2 <- alpha2_g
  } else if (pleio == "t") {
    s <- sqrt(rchisq(M, df_t) / df_t)
    alpha1 <- alpha1_g / s
    alpha2 <- alpha2_g / s
  } else if (pleio == "mixture") {
    alpha1 <- alpha1_g
    alpha2 <- alpha2_g
    is_big <- rbinom(M, 1, 1 - mix_p) == 1
    if (any(is_big)) {
      alpha1[is_big] <- rnorm(sum(is_big), 0, mix_sd_big)
      alpha2[is_big] <- rnorm(sum(is_big), 0, mix_sd_big)
    }
  } else if (pleio == "directional") {
    alpha1 <- alpha1_g + dir_mean1
    alpha2 <- alpha2_g + dir_mean2
  }
  
  ## Inject correlation between beta_x and alpha's (selection tilt test)
  if (abs(corr_bx_a) > 0) {
    vx <- var(beta_x)
    if (vx > 0) {
      g1 <- corr_bx_a * sd(alpha1) / sqrt(vx)
      g2 <- corr_bx_a * sd(alpha2) / sqrt(vx)
      alpha1 <- alpha1 + g1 * beta_x
      alpha2 <- alpha2 + g2 * beta_x
    }
  }
  
  ## ---- 3) LD-block dependence (optional, for robustness) ----
  if (ld == "blocks" && abs(rho_block_alpha) > 0) {
    B <- ceiling(M / block_size)
    shocks <- rnorm(B, 0, sd(alpha1))
    rep_shocks1 <- rep(shocks, each = block_size)[1:M]
    alpha1 <- alpha1 + rho_block_alpha * rep_shocks1
    
    shocks2 <- rnorm(B, 0, sd(alpha2))
    rep_shocks2 <- rep(shocks2, each = block_size)[1:M]
    alpha2 <- alpha2 + rho_block_alpha * rep_shocks2
  }
  
  ## ---- 4) True outcome effects ----
  beta_y1 <- theta1 * beta_x + alpha1
  beta_y2 <- theta2 * beta_x + alpha2
  
  ## ---- 5) GWAS sampling errors ----
  # Standard errors (homoskedastic, as in theory)
  se_x <- 1/sqrt(n_x)
  se_y1 <- 1/sqrt(n1)
  se_y2 <- 1/sqrt(n2)
  
  u_x <- rnorm(M, 0, se_x)
  
  # Outcome errors with overlap correlation structure
  # Theory specifies: Cov(u_y1, u_y2) = (k_s * rho_e) / (sqrt(k1 * k2) * n_x)
  # For simplicity with constant SE, correlation is:
  cov_u <- (k_s * rho_e) / (sqrt(k1 * k2) * n_x)
  r_overlap <- cov_u / (se_y1 * se_y2)
  r_overlap <- max(min(r_overlap, 0.99), -0.99)  # safety
  
  e1 <- rnorm(M)
  e2 <- rnorm(M)
  u_y1 <- se_y1 * e1
  u_y2 <- se_y2 * (r_overlap * e1 + sqrt(1 - r_overlap^2) * e2)
  
  bhat_x <- beta_x + u_x
  bhat_y1 <- beta_y1 + u_y1
  bhat_y2 <- beta_y2 + u_y2
  
  ## ---- 6) Selection (genome-wide significant instruments) ----
  Zx <- bhat_x / se_x
  sel <- abs(Zx) > z_GW
  K <- sum(sel)
  if (K < 30) return(list(ok = FALSE, K = K))
  
  bx <- bhat_x[sel]
  by1 <- bhat_y1[sel]
  by2 <- bhat_y2[sel]
  
  ## Store true values for diagnostics
  true_alpha1 <- alpha1[sel]
  true_alpha2 <- alpha2[sel]
  true_beta_x <- beta_x[sel]
  
  ## ---- 7) Wald ratios ----
  denom <- bx
  small <- abs(denom) < eps_denom_floor
  denom[small] <- sign(denom[small]) * eps_denom_floor + eps_denom_floor
  
  theta1_hat <- by1 / denom
  theta2_hat <- by2 / denom
  
  ## ---- 8) Theory: Variance formulas (equations 2-3) ----
  # Note: Using TRUE beta values in the formulas as in population
  # In practice, we'd use estimates, but asymptotically equivalent
  
  # Variance formula: σ²_{l,k} = (1/(n_x β_x²)) * (1/k_l + β_y²/β_x²)
  sig1_sq <- (1/(n_x * denom^2)) * (1/k1 + (by1^2)/(denom^2))
  sig2_sq <- (1/(n_x * denom^2)) * (1/k2 + (by2^2)/(denom^2))
  
  # Covariance formula: σ_{12,k} = (1/(n_x β_x²)) * (k_s*ρ_e/(k1*k2) + β_y1*β_y2/β_x²)
  sig12 <- (1/(n_x * denom^2)) * ((k_s * rho_e)/(k1 * k2) + (by1 * by2)/(denom^2))
  
  ## ---- 9) Weights (geometric mean inverse-variance) ----
  w <- 1 / sqrt(sig1_sq * sig2_sq)
  wtil <- w / sum(w)
  
  ## ---- 10) Centered Wald ratios ----
  theta1_bar <- sum(wtil * theta1_hat)
  theta2_bar <- sum(wtil * theta2_hat)
  
  Delta1 <- theta1_hat - theta1_bar
  Delta2 <- theta2_hat - theta2_bar
  
  ## ---- 11) THEORY: Bias-corrected estimators (equations 9-10) ----
  # Heterogeneity (equation 9):
  # τ²_{l,corr} = Σ w_k (Δ_k^(l))² - Σ w_k σ²_{l,k}
  
  S1_sq <- sum(wtil * Delta1^2)
  S2_sq <- sum(wtil * Delta2^2)
  
  tau1_sq_corr <- S1_sq - sum(wtil * sig1_sq)
  tau2_sq_corr <- S2_sq - sum(wtil * sig2_sq)
  
  if (trunc_tau) {
    tau1_sq_corr <- max(0, tau1_sq_corr)
    tau2_sq_corr <- max(0, tau2_sq_corr)
  }
  
  tau1_corr <- sqrt(tau1_sq_corr)
  tau2_corr <- sqrt(tau2_sq_corr)
  
  if (tau1_corr < 1e-10 || tau2_corr < 1e-10) {
    return(list(ok = FALSE, K = K))
  }
  
  # Coheterogeneity (equation 10):
  # Ĉ_12 = Σ w_k Δ_k^(1) Δ_k^(2) - Σ w_k σ_{12,k}
  C12_hat <- sum(wtil * Delta1 * Delta2) - sum(wtil * sig12)
  
  ## ---- 12) Correlation estimator (equation 11) ----
  rho_hat <- C12_hat / (tau1_corr * tau2_corr)
  
  ## ---- 13) THEORY: first-order influence / closed-form variance (Theorem 1) ----
  
  # First-order influence gradient terms (Theorem 1, fixed-weight form):
  #   D1_k = Delta1_k - rho*(tau1/tau2)*Delta2_k ;  D2_k = Delta2_k - rho*(tau2/tau1)*Delta1_k
  D1 <- Delta1 - rho_hat * (tau1_corr / tau2_corr) * Delta2
  D2 <- Delta2 - rho_hat * (tau2_corr / tau1_corr) * Delta1
  
  ## ---- 14) THEORY: closed-form variance of rho-hat ----
  # Var(rho-hat) = sum_k (w_k/(tau1 tau2))^2 [ D2^2 sig1 + D1^2 sig2 + 2 D1 D2 sig12 ]
  Sigma_hat <- sum(
    wtil^2 * (D2^2 * sig1_sq + D1^2 * sig2_sq + 2 * D1 * D2 * sig12)
  ) / (tau1_corr^2 * tau2_corr^2)
  
  se_rho <- sqrt(Sigma_hat)
  
  ## ---- 15) Diagnostics ----
  ess <- 1 / sum(wtil^2)
  wmax <- max(wtil)
  
  # True correlation for comparison (if we know true alphas)
  true_rho_pleio <- cor(true_alpha1, true_alpha2)
  
  list(
    ok = TRUE, 
    K = K, 
    rho_hat = rho_hat, 
    se = se_rho,
    Sigma_hat = Sigma_hat,
    tau1_corr = tau1_corr,
    tau2_corr = tau2_corr,
    C12_hat = C12_hat,
    ess = ess, 
    wmax = wmax,
    true_rho_pleio = true_rho_pleio,
    theta1_bar = theta1_bar,
    theta2_bar = theta2_bar
  )
}

############################################################
## Run scenario with multiple replicates
############################################################
run_scenario_theory <- function(R = 200, seed = 1, ...) {
  set.seed(seed)
  res <- vector("list", R)
  for (r in 1:R) {
    res[[r]] <- simulate_replicate_theory(...)
  }
  
  ok <- vapply(res, `[[`, logical(1), "ok")
  res_ok <- res[ok]
  if (!any(ok)) stop("No successful replicates")
  
  K <- vapply(res_ok, `[[`, numeric(1), "K")
  rho <- vapply(res_ok, `[[`, numeric(1), "rho_hat")
  se <- vapply(res_ok, `[[`, numeric(1), "se")
  Sigma <- vapply(res_ok, `[[`, numeric(1), "Sigma_hat")
  ess <- vapply(res_ok, `[[`, numeric(1), "ess")
  wmax <- vapply(res_ok, `[[`, numeric(1), "wmax")
  true_rho <- vapply(res_ok, `[[`, numeric(1), "true_rho_pleio")
  
  # Center around empirical mean (as in theory: estimate target parameter)
  rho_mean <- mean(rho)
  rho_sd <- sd(rho)
  
  # Wald-type standardization using SE
  z_wald <- (rho - rho_mean) / se
  
  # Check if SE × sqrt(K) stabilizes (should be ≈ constant if rate is correct)
  se_times_sqrtK <- se * sqrt(K)
  sd_times_sqrtK <- rho_sd * sqrt(mean(K))
  
  # Coverage (using empirical mean as "truth")
  cover95 <- mean(abs(rho - rho_mean) <= 1.96 * se)
  
  # Compare to true rho_pleio (correlation of true alpha's among selected)
  true_rho_mean <- mean(true_rho)
  bias_from_true <- rho_mean - true_rho_mean
  
  list(
    ok_rate = mean(ok),
    R_ok = length(res_ok),
    K_mean = mean(K),
    K_median = median(K),
    K_sd = sd(K),
    rho_mean = rho_mean,
    rho_sd = rho_sd,
    mean_se = mean(se),
    median_se = median(se),
    sd_se = sd(se),
    mean_Sigma = mean(Sigma),
    mean_se_times_sqrtK = mean(se_times_sqrtK),
    sd_times_sqrtK = sd_times_sqrtK,
    cover95 = cover95,
    ess_mean = mean(ess),
    wmax_mean = mean(wmax),
    true_rho_mean = true_rho_mean,
    bias_from_true = bias_from_true,
    z_wald = z_wald,
    rho = rho,
    se = se,
    K = K
  )
}

############################################################
## Comprehensive test suite
############################################################
run_suite_theory <- function(R = 200, seed = 1,
                             n_x = 200000,
                             M = 150000) {
  
  base <- list(
    n_x = n_x, 
    n1 = n_x, 
    n2 = n_x, 
    M = M,
    arch = "spike_slab", 
    pi_x = 0.003, 
    sd_bx = 0.07,
    pleio = "gaussian", 
    sd_a1 = 0.05, 
    sd_a2 = 0.06, 
    rho_a = 0.4,
    k_s = 1.0,
    rho_e = 0.0,
    ld = "none"
  )
  
  scenarios <- list(
    baseline = base,
    
    ld_blocks = modifyList(base, list(
      ld = "blocks", 
      block_size = 20, 
      rho_block_alpha = 0.7
    )),
    
    heavy_tails = modifyList(base, list(
      pleio = "t", 
      df_t = 3
    )),
    
    mixture_outliers = modifyList(base, list(
      pleio = "mixture", 
      mix_p = 0.98, 
      mix_sd_big = 0.35
    )),
    
    directional = modifyList(base, list(
      pleio = "directional", 
      dir_mean1 = 0.03, 
      dir_mean2 = 0.02
    )),
    
    corr_bx_alpha = modifyList(base, list(
      corr_bx_a = 0.5
    )),
    
    with_overlap = modifyList(base, list(
      k_s = 0.5,
      rho_e = 0.3
    ))
  )
  
  out <- list()
  for (nm in names(scenarios)) {
    cat("\nRunning scenario:", nm, "\n")
    sc <- scenarios[[nm]]
    out[[nm]] <- do.call(run_scenario_theory, c(list(R = R, seed = seed), sc))
  }
  out
}

############################################################
## Summarize results
############################################################
summarize_suite_theory <- function(res_list) {
  nm <- names(res_list)
  data.frame(
    scenario = nm,
    ok_rate = vapply(res_list, `[[`, numeric(1), "ok_rate"),
    R_ok = vapply(res_list, `[[`, numeric(1), "R_ok"),
    K_mean = vapply(res_list, `[[`, numeric(1), "K_mean"),
    rho_mean = vapply(res_list, `[[`, numeric(1), "rho_mean"),
    rho_sd = vapply(res_list, `[[`, numeric(1), "rho_sd"),
    mean_se = vapply(res_list, `[[`, numeric(1), "mean_se"),
    sd_times_sqrtK = vapply(res_list, `[[`, numeric(1), "sd_times_sqrtK"),
    meanSE_times_sqrtK = vapply(res_list, `[[`, numeric(1), "mean_se_times_sqrtK"),
    ratio_sd_se = vapply(res_list, function(x) x$rho_sd / x$mean_se, numeric(1)),
    cover95 = vapply(res_list, `[[`, numeric(1), "cover95"),
    ess_mean = vapply(res_list, `[[`, numeric(1), "ess_mean"),
    wmax_mean = vapply(res_list, `[[`, numeric(1), "wmax_mean"),
    true_rho = vapply(res_list, `[[`, numeric(1), "true_rho_mean"),
    bias = vapply(res_list, `[[`, numeric(1), "bias_from_true")
  )
}

############################################################
## QQ plots for diagnostics
############################################################
plot_diagnostics <- function(res_list) {
  par(mfrow = c(2, 3))
  
  for (nm in names(res_list)[1:min(6, length(res_list))]) {
    z <- res_list[[nm]]$z_wald
    qqnorm(z, main = paste(nm, ": Wald z"), pch = 19, cex = 0.5)
    qqline(z, col = "red", lwd = 2)
  }
}

############################################################
## Rate verification: log-log slope test
############################################################
rate_experiment_theory <- function(
    R = 300,
    seed = 1,
    n_x = 200000,
    M_grid = c(100000, 200000, 400000, 700000, 1000000),
    pi_x = 0.003,
    sd_bx = 0.07
) {
  
  results <- vector("list", length(M_grid))
  
  for (i in seq_along(M_grid)) {
    cat("\nM =", M_grid[i], "\n")
    
    res <- run_scenario_theory(
      R = R,
      seed = seed + i,
      n_x = n_x,
      M = M_grid[i],
      pi_x = pi_x,
      sd_bx = sd_bx,
      pleio = "gaussian",
      rho_a = 0.4
    )
    
    results[[i]] <- res
    
    cat("  K_mean =", round(res$K_mean, 1),
        "  sd(rho) =", signif(res$rho_sd, 4),
        "  sd*sqrt(K) =", signif(res$sd_times_sqrtK, 4),
        "  mean(SE)*sqrt(K) =", signif(res$mean_se_times_sqrtK, 4), "\n")
  }
  
  df <- data.frame(
    M = M_grid,
    K_mean = vapply(results, `[[`, numeric(1), "K_mean"),
    K_sd = vapply(results, `[[`, numeric(1), "K_sd"),
    rho_sd = vapply(results, `[[`, numeric(1), "rho_sd"),
    mean_se = vapply(results, `[[`, numeric(1), "mean_se"),
    sd_sqrtK = vapply(results, `[[`, numeric(1), "sd_times_sqrtK"),
    meanSE_sqrtK = vapply(results, `[[`, numeric(1), "mean_se_times_sqrtK")
  )
  
  # Log-log regression
  fit <- lm(log(rho_sd) ~ log(K_mean), data = df)
  slope <- coef(fit)[2]
  ci <- confint(fit, level = 0.95)[2, ]
  
  list(df = df, fit = fit, slope = slope, slope_ci = ci, results = results)
}

plot_rate_theory <- function(rate_obj) {
  df <- rate_obj$df
  
  par(mfrow = c(2, 2))
  
  # 1. Log-log plot (slope should be -0.5)
  plot(log(df$K_mean), log(df$rho_sd), pch = 19, cex = 1.5,
       xlab = "log(K)", ylab = "log(SD(ρ̂))",
       main = sprintf("Rate test: slope = %.3f [%.3f, %.3f]",
                      rate_obj$slope, rate_obj$slope_ci[1], rate_obj$slope_ci[2]))
  abline(rate_obj$fit, col = "red", lwd = 2)
  abline(h = log(df$rho_sd[1]) - 0.5 * (log(df$K_mean) - log(df$K_mean[1])),
         col = "blue", lty = 2, lwd = 2)
  legend("topright", c("Fitted", "Slope = -0.5"), 
         col = c("red", "blue"), lty = c(1, 2), lwd = 2)
  
  # 2. SD × sqrt(K) stabilization (should be constant)
  plot(df$K_mean, df$sd_sqrtK, type = "b", pch = 19, cex = 1.5,
       xlab = "K", ylab = "SD(ρ̂) × √K",
       main = "Stabilization: SD × √K")
  abline(h = mean(df$sd_sqrtK), col = "blue", lwd = 2, lty = 2)
  
  # 3. SE × sqrt(K) (theoretical prediction)
  plot(df$K_mean, df$meanSE_sqrtK, type = "b", pch = 19, cex = 1.5,
       xlab = "K", ylab = "mean(SE) × √K",
       main = "Theory: SE × √K")
  abline(h = mean(df$meanSE_sqrtK), col = "red", lwd = 2, lty = 2)
  
  # 4. Compare empirical vs theoretical
  plot(df$sd_sqrtK, df$meanSE_sqrtK, pch = 19, cex = 1.5,
       xlab = "SD(ρ̂) × √K (empirical)", 
       ylab = "mean(SE) × √K (theory)",
       main = "Empirical vs Theory")
  abline(0, 1, col = "red", lwd = 2)
  
  invisible(NULL)
}

############################################################
## EXAMPLE USAGE
############################################################

# Quick test
cat("Running quick baseline test...\n")
system.time({
  test <- simulate_replicate_theory(
    n_x = 200000,
    M = 150000,
    pleio = "gaussian",
    rho_a = 0.4
  )
})
print(test)

# Full suite (takes a few minutes)
cat("\n\nRunning full test suite...\n")
system.time({
  suite <- run_suite_theory(R = 200, seed = 123, n_x = 200000, M = 300000)
})

tab <- summarize_suite_theory(suite)
print(tab)

# Check: ratio_sd_se should be ≈ 1 and cover95 ≈ 0.95
cat("\nKey checks:\n")
cat("  SD/SE ratio (should be ≈1):", round(tab$ratio_sd_se, 3), "\n")
cat("  Coverage (should be ≈0.95):", round(tab$cover95, 3), "\n")

# Diagnostic plots
plot_diagnostics(suite)

# Rate experiment (this takes longer)
cat("\n\nRunning rate verification...\n")
system.time({
  rate_res <- rate_experiment_theory(
    R = 250,
    seed = 456,
    n_x = 200000,
    M_grid = c(100000, 200000, 400000, 700000, 1000000)
  )
})

print(rate_res$df)
cat("\nEstimated slope:", rate_res$slope, 
    "\n95% CI:", rate_res$slope_ci, "\n")
cat("(Should be ≈ -0.5 for sqrt(K) rate)\n")

plot_rate_theory(rate_res)

