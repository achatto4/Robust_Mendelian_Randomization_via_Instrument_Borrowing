#!/usr/bin/env Rscript

# ============================================================================
# MR2 Calibration: Type I Error and Power Analysis
# ============================================================================

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript mr2_calibration.R <N> <output_dir>")
}

N_INPUT <- as.numeric(args[1])
OUTPUT_DIR <- args[2]

cat(sprintf("\n==========================================================\n"))
cat(sprintf("Starting MR2 Calibration for N = %d\n", N_INPUT))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("==========================================================\n\n"))

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(MR2)
  library(dplyr)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

TARGET_ALPHA <- 0.05

# Simulation parameters
N <- N_INPUT
nx <- N
ny <- N / 2
ny_alt <- N / 1
pthr <- 5e-8
M <- 2e5

# Pleiotropy parameters
thetaU <- 0.3
thetaU_alt <- 0.3
theta_alt <- 0.3

# IV parameters
prop_invalid <- 0.3
pi1 <- 0.02 * (1 - prop_invalid)
pi2 <- 0.02 * prop_invalid
pi3 <- 0.01

# Effect size parameters
sigma2x <- sigma2y <- 5e-5
sigma2u <- 1e-4
sigma2x_td <- sigma2y_td <- (5e-5) - thetaU * thetaU * sigma2u

# Overlap between traits
overlap <- 0.5
overlap_adj <- (2 * overlap) / (overlap + 1)

# Number of replications
N_REP_NULL <- 200
N_REP_ALT <- 200

# MR2 parameters
MR2_NITER <- 1500
MR2_BURNIN <- 500
MR2_THIN <- 5

cat(sprintf("Configuration:\n"))
cat(sprintf("  N = %d (nx=%d, ny=%d, ny_alt=%d)\n", N, nx, ny, ny_alt))
cat(sprintf("  Proportion invalid IVs: %.2f\n", prop_invalid))
cat(sprintf("  Replications: Null=%d, Alt=%d\n", N_REP_NULL, N_REP_ALT))
cat(sprintf("  MCMC: niter=%d, burnin=%d, thin=%d\n", MR2_NITER, MR2_BURNIN, MR2_THIN))
cat(sprintf("\n"))

# ============================================================================
# FUNCTION: Generate Data and Run MR2
# ============================================================================

run_mr2_simulation <- function(theta, seed_offset = 0) {
  
  set.seed(seed_offset)
  
  # Generate SNP indices
  ind1 <- sample(M, round(M * pi1))
  ind2 <- sample(setdiff(1:M, ind1), round(M * pi2))
  ind3 <- sample(setdiff(1:M, c(ind1, ind2)), round(M * pi3))
  
  # Overlapping indices for secondary outcome
  shared_ind1 <- sample(ind1, round(length(ind1) * overlap_adj))
  ind1_alt <- c(shared_ind1, sample(setdiff(1:M, shared_ind1), 
                                    round(M * pi1) - length(shared_ind1)))
  
  shared_ind2 <- sample(ind2, round(length(ind2) * overlap_adj))
  ind2_alt <- c(shared_ind2, sample(setdiff(1:M, c(ind1_alt, shared_ind2)), 
                                    round(M * pi2) - length(shared_ind2)))
  
  shared_ind3 <- sample(ind3, round(length(ind3) * overlap_adj))
  ind3_alt <- c(shared_ind3, sample(setdiff(1:M, c(ind1_alt, ind2_alt, shared_ind3)), 
                                    round(M * pi3) - length(shared_ind3)))
  
  # Initialize effect vectors
  gamma <- phi <- phi_alt <- alpha <- alpha_alt <- rep(0, M)
  
  # Generate true effects
  gamma[union(ind1, ind1_alt)] <- rnorm(length(union(ind1, ind1_alt)), 
                                        mean = 0, sd = sqrt(sigma2x))
  gamma[union(ind2, ind2_alt)] <- rnorm(length(union(ind2, ind2_alt)), 
                                        mean = 0, sd = sqrt(sigma2x_td))
  
  alpha[ind2] <- rnorm(length(ind2), mean = 0.005, sd = sqrt(sigma2y_td))
  alpha[ind3] <- rnorm(length(ind3), mean = 0.005, sd = sqrt(sigma2y))
  alpha_alt[ind2_alt] <- rnorm(length(ind2_alt), mean = 0.003, sd = sqrt(sigma2y_td))
  alpha_alt[ind3_alt] <- rnorm(length(ind3_alt), mean = 0.003, sd = sqrt(sigma2y))
  
  phi[union(ind2, ind2_alt)] <- rnorm(length(union(ind2, ind2_alt)), 
                                      mean = 0, sd = sqrt(sigma2u))
  phi_alt[ind2_alt] <- phi[ind2_alt]
  phi[-ind2] <- 0
  
  # Calculate true effects
  betax <- gamma + thetaU * phi
  betay <- alpha + theta * betax + thetaU * phi
  betay_alt <- alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
  
  # Generate observed effects with noise
  betahat_x <- betax + rnorm(M, mean = 0, sd = sqrt(1 / nx))
  betahat_y <- betay + rnorm(M, mean = 0, sd = sqrt(1 / ny))
  betahat_y_alt <- betay_alt + rnorm(M, mean = 0, sd = sqrt(1 / ny_alt))
  
  # Filter significant SNPs
  ind_filter <- which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
  numIV <- length(ind_filter)
  
  # Need at least 3 IVs
  if (numIV < 3) {
    return(list(
      pip = NA,
      estimate = NA,
      se = NA,
      numIV = numIV,
      success = FALSE
    ))
  }
  
  betahat_x.flt <- betahat_x[ind_filter]
  betahat_y.flt <- betahat_y[ind_filter]
  betahat_y_alt.flt <- betahat_y_alt[ind_filter]
  
  # Run MR2
  result <- tryCatch({
    betaHat_Y_matrix <- cbind(betahat_y.flt, betahat_y_alt.flt)
    betaHat_X_matrix <- matrix(betahat_x.flt, ncol = 1)
    
    # Run MR2
    MR2_fit <- MR2(
      betaHat_Y = betaHat_Y_matrix,
      betaHat_X = betaHat_X_matrix,
      EVgamma = 0.5,
      niter = MR2_NITER,
      burnin = MR2_BURNIN,
      thin = MR2_THIN,
      monitor = 0  # Suppress output
    )
    
    # Post-processing
    PostProc_res <- PostProc(MR2_fit, betaHat_Y_matrix, betaHat_X_matrix, 
                             alpha = 0.05)
    
    # Extract results for PRIMARY outcome (column 1)
    pip <- as.numeric(MR2_fit$postMean$gamma[1, 1])
    estimate <- as.numeric(PostProc_res$thetaPost[1, 1])
    
    # SE from posterior samples
    theta_samples <- MR2_fit$theta[, 1]
    se_samples <- sd(theta_samples)
    
    # Clean up
    rm(MR2_fit, PostProc_res, betaHat_Y_matrix, betaHat_X_matrix)
    gc(verbose = FALSE)
    
    # Return results
    list(
      pip = pip,
      estimate = estimate,
      se = se_samples,
      numIV = numIV,
      success = TRUE
    )
    
  }, error = function(e) {
    list(
      pip = NA,
      estimate = NA,
      se = NA,
      numIV = numIV,
      success = FALSE
    )
  })
  
  return(result)
}

# ============================================================================
# STEP 1: CALIBRATION - Find PIP threshold for Type I Error = 0.05
# ============================================================================

OVERALL_START <- Sys.time()

cat("\n============================================================\n")
cat("STEP 1: CALIBRATION (Under Null: theta = 0)\n")
cat("============================================================\n")

# Storage for null simulations
null_pips <- numeric(N_REP_NULL)
null_estimates <- numeric(N_REP_NULL)
null_ses <- numeric(N_REP_NULL)
null_numIVs <- numeric(N_REP_NULL)
null_success <- logical(N_REP_NULL)

# Run simulations under null
start_time <- Sys.time()
for (i in 1:N_REP_NULL) {
  if (i %% 20 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    estimated_total <- elapsed * N_REP_NULL / i
    remaining <- estimated_total - elapsed
    cat(sprintf("  Progress: %d/%d (%.1f%%) | Elapsed: %.1f min | Remaining: %.1f min\n", 
                i, N_REP_NULL, 100*i/N_REP_NULL, elapsed, remaining))
  }
  
  result <- run_mr2_simulation(theta = 0, seed_offset = 10000 + i)
  
  null_pips[i] <- result$pip
  null_estimates[i] <- result$estimate
  null_ses[i] <- result$se
  null_numIVs[i] <- result$numIV
  null_success[i] <- result$success
}

# Remove failed runs
null_pips_valid <- null_pips[null_success & !is.na(null_pips)]
n_valid_null <- length(null_pips_valid)
success_rate_null <- n_valid_null / N_REP_NULL

cat(sprintf("\nCompleted! Valid runs: %d/%d (%.1f%%)\n", 
            n_valid_null, N_REP_NULL, 100*success_rate_null))

# Calculate Type I error rate for various PIP thresholds
pip_thresholds <- seq(0.1, 0.9, by = 0.01)
type1_error_rates <- sapply(pip_thresholds, function(threshold) {
  mean(null_pips_valid > threshold)
})

# Find threshold closest to target alpha
optimal_idx <- which.min(abs(type1_error_rates - TARGET_ALPHA))
optimal_pip_threshold <- pip_thresholds[optimal_idx]
achieved_type1_error <- type1_error_rates[optimal_idx]

cat("\n--- CALIBRATION RESULTS ---\n")
cat(sprintf("Target Type I error rate: %.3f\n", TARGET_ALPHA))
cat(sprintf("Optimal PIP threshold: %.3f\n", optimal_pip_threshold))
cat(sprintf("Achieved Type I error: %.3f\n", achieved_type1_error))
cat(sprintf("Difference: %.4f\n", abs(achieved_type1_error - TARGET_ALPHA)))

# Summary statistics of PIPs under null
cat("\n--- PIP Distribution Under Null ---\n")
cat(sprintf("Mean PIP: %.3f\n", mean(null_pips_valid)))
cat(sprintf("Median PIP: %.3f\n", median(null_pips_valid)))
cat(sprintf("SD PIP: %.3f\n", sd(null_pips_valid)))
cat(sprintf("Range: [%.3f, %.3f]\n", min(null_pips_valid), max(null_pips_valid)))
cat(sprintf("Q25, Q50, Q75: %.3f, %.3f, %.3f\n", 
            quantile(null_pips_valid, 0.25),
            quantile(null_pips_valid, 0.50),
            quantile(null_pips_valid, 0.75)))

# Save calibration curve data
calibration_curve <- data.frame(
  pip_threshold = pip_thresholds,
  type1_error = type1_error_rates
)

# ============================================================================
# STEP 2: POWER CALCULATION - Apply calibrated threshold under alternative
# ============================================================================

cat("\n============================================================\n")
cat("STEP 2: POWER CALCULATION (Under Alternative)\n")
cat("============================================================\n")

# Test multiple alternative values of theta
theta_alternatives <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)

power_results <- data.frame(
  N = integer(),
  theta_true = numeric(),
  n_valid = integer(),
  success_rate = numeric(),
  power = numeric(),
  mean_pip = numeric(),
  sd_pip = numeric(),
  mean_estimate = numeric(),
  sd_estimate = numeric(),
  bias = numeric(),
  rmse = numeric(),
  mean_se = numeric(),
  mean_numIV = numeric()
)

# Detailed results storage
detailed_results <- data.frame()

for (theta_alt_val in theta_alternatives) {
  cat(sprintf("\n--- Testing theta = %.2f ---\n", theta_alt_val))
  
  # Storage
  alt_pips <- numeric(N_REP_ALT)
  alt_estimates <- numeric(N_REP_ALT)
  alt_ses <- numeric(N_REP_ALT)
  alt_numIVs <- numeric(N_REP_ALT)
  alt_success <- logical(N_REP_ALT)
  
  # Run simulations
  start_time <- Sys.time()
  for (i in 1:N_REP_ALT) {
    if (i %% 20 == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      estimated_total <- elapsed * N_REP_ALT / i
      remaining <- estimated_total - elapsed
      cat(sprintf("  Progress: %d/%d (%.1f%%) | Elapsed: %.1f min | Remaining: %.1f min\n", 
                  i, N_REP_ALT, 100*i/N_REP_ALT, elapsed, remaining))
    }
    
    result <- run_mr2_simulation(theta = theta_alt_val, seed_offset = 20000 + i)
    
    alt_pips[i] <- result$pip
    alt_estimates[i] <- result$estimate
    alt_ses[i] <- result$se
    alt_numIVs[i] <- result$numIV
    alt_success[i] <- result$success
  }
  
  # Remove failed runs
  valid_idx <- alt_success & !is.na(alt_pips)
  alt_pips_valid <- alt_pips[valid_idx]
  alt_estimates_valid <- alt_estimates[valid_idx]
  alt_ses_valid <- alt_ses[valid_idx]
  alt_numIVs_valid <- alt_numIVs[valid_idx]
  n_valid_alt <- length(alt_pips_valid)
  
  # Calculate metrics
  power <- mean(alt_pips_valid > optimal_pip_threshold)
  mean_pip <- mean(alt_pips_valid)
  sd_pip <- sd(alt_pips_valid)
  mean_estimate <- mean(alt_estimates_valid)
  sd_estimate <- sd(alt_estimates_valid)
  bias <- mean_estimate - theta_alt_val
  rmse <- sqrt(mean((alt_estimates_valid - theta_alt_val)^2))
  mean_se <- mean(alt_ses_valid)
  mean_numIV <- mean(alt_numIVs_valid)
  
  cat(sprintf("\nResults for theta = %.2f:\n", theta_alt_val))
  cat(sprintf("  Valid runs: %d/%d (%.1f%%)\n", n_valid_alt, N_REP_ALT, 100*n_valid_alt/N_REP_ALT))
  cat(sprintf("  Power (PIP > %.3f): %.3f\n", optimal_pip_threshold, power))
  cat(sprintf("  Mean PIP: %.3f (SD: %.3f)\n", mean_pip, sd_pip))
  cat(sprintf("  Mean estimate: %.4f (SD: %.4f, true = %.2f)\n", mean_estimate, sd_estimate, theta_alt_val))
  cat(sprintf("  Bias: %.4f\n", bias))
  cat(sprintf("  RMSE: %.4f\n", rmse))
  cat(sprintf("  Mean SE: %.4f\n", mean_se))
  cat(sprintf("  Mean # IVs: %.1f\n", mean_numIV))
  
  # Store summary results
  power_results <- rbind(power_results, data.frame(
    N = N,
    theta_true = theta_alt_val,
    n_valid = n_valid_alt,
    success_rate = n_valid_alt / N_REP_ALT,
    power = power,
    mean_pip = mean_pip,
    sd_pip = sd_pip,
    mean_estimate = mean_estimate,
    sd_estimate = sd_estimate,
    bias = bias,
    rmse = rmse,
    mean_se = mean_se,
    mean_numIV = mean_numIV
  ))
  
  # Store detailed results for this theta
  detailed_results <- rbind(detailed_results, data.frame(
    N = N,
    theta_true = theta_alt_val,
    rep = 1:N_REP_ALT,
    pip = alt_pips,
    estimate = alt_estimates,
    se = alt_ses,
    numIV = alt_numIVs,
    success = alt_success,
    reject = alt_pips > optimal_pip_threshold
  ))
}

OVERALL_END <- Sys.time()
total_runtime <- difftime(OVERALL_END, OVERALL_START, units = "mins")

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n============================================================\n")
cat("SAVING RESULTS\n")
cat("============================================================\n")

# Summary table
summary_table <- data.frame(
  Parameter = c(
    "N", "nx", "ny", "ny_alt", 
    "prop_invalid", "thetaU", "overlap",
    "target_alpha", "optimal_pip_threshold", "achieved_type1_error",
    "n_rep_null", "n_rep_alt", "success_rate_null", "runtime_minutes"
  ),
  Value = c(
    N, nx, ny, ny_alt,
    prop_invalid, thetaU, overlap,
    TARGET_ALPHA, optimal_pip_threshold, achieved_type1_error,
    N_REP_NULL, N_REP_ALT, success_rate_null, as.numeric(total_runtime)
  )
)

# Save all results
output_prefix <- file.path(OUTPUT_DIR, sprintf("MR2_N%d", N))

write.csv(summary_table, 
          paste0(output_prefix, "_summary.csv"), 
          row.names = FALSE)

write.csv(power_results, 
          paste0(output_prefix, "_power.csv"), 
          row.names = FALSE)

write.csv(calibration_curve,
          paste0(output_prefix, "_calibration_curve.csv"),
          row.names = FALSE)

write.csv(detailed_results,
          paste0(output_prefix, "_detailed_results.csv"),
          row.names = FALSE)

# Save null distribution
null_distribution <- data.frame(
  N = N,
  rep = 1:N_REP_NULL,
  pip = null_pips,
  estimate = null_estimates,
  se = null_ses,
  numIV = null_numIVs,
  success = null_success
)

write.csv(null_distribution,
          paste0(output_prefix, "_null_distribution.csv"),
          row.names = FALSE)

# Save RData with all objects
save(
  N, nx, ny, ny_alt,
  prop_invalid, thetaU, overlap,
  TARGET_ALPHA, optimal_pip_threshold, achieved_type1_error,
  power_results, calibration_curve, detailed_results,
  null_distribution, summary_table,
  file = paste0(output_prefix, "_results.rda")
)

cat(sprintf("\nResults saved with prefix: %s\n", output_prefix))

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("CALIBRATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("Sample size: N = %d\n", N))
cat(sprintf("Optimal PIP threshold: %.3f\n", optimal_pip_threshold))
cat(sprintf("Achieved Type I error: %.3f (target: %.3f)\n", 
            achieved_type1_error, TARGET_ALPHA))
cat(sprintf("Total runtime: %.1f minutes (%.2f hours)\n", total_runtime, total_runtime/60))
cat("\n--- Power Results ---\n")
print(power_results[, c("theta_true", "power", "mean_pip", "bias", "rmse")])
cat("\n============================================================\n")

cat(sprintf("\nAll results saved to: %s\n", OUTPUT_DIR))
cat("Done!\n")