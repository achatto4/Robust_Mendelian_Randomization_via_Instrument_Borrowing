rm(list = ls())

# Load required libraries
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(penalized)
library(ks)
library(MRPRESSO)

# Check if MR2 is installed, if not, install it
if (!requireNamespace("MR2", quietly = TRUE)) {
  message("MR2 package not found. Installing from GitHub...")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("lb664/MR2")
}
library(MR2)

# Set working directory and source files
code_dir <- "."  # set to simulation_code/ directory
setwd(code_dir)

# Source required functions
library(IBMR)  # coheterogeneity_Q, IBMODE, IBPRESSO all available via IBMR
# MBE_multi_old_wt1 -> use IBMR::IBMODE
# mr_presso_ib -> use IBMR::IBPRESSO

# ============================================================================
# CONFIGURATION: Methods to silence (skip in analysis)
# ============================================================================
# Set to TRUE to skip a method, FALSE to run it
SILENCE_METHODS <- list(
  MRMode = FALSE,
  ContMix = FALSE,
  MRMix = FALSE,
  mode_new_phi1 = FALSE,
  MRPRESSO = TRUE,      # SILENCED
  MRcML = FALSE,
  IVW = FALSE,
  median = FALSE,
  Egger = FALSE,
  IBPRESSO = TRUE,      # SILENCED
  MR2 = FALSE           # ACTIVE
)
# ============================================================================

# Parameter vectors
thetavec = c(0.2, 0, -0.2, 0.1, -0.1)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5, 1.5e5, 2e5, 5e5, 1e6) # 1:7
prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7)
overlap_vec = c(0.5, 0.75, 1)

# Get command line arguments or use defaults
temp = as.integer(commandArgs(trailingOnly = TRUE))
# Default parameter cell when run without command-line arguments
# (indices into thetavec / thetaUvec / Nvec / prop_invalid_vec / overlap_vec)
if (length(temp) == 0) temp = c(4, 1, 4, 2, 2)

# Set defaults if no command line arguments provided
if (length(temp) == 0) {
  message("No command line arguments provided. Using defaults:")
  temp <- c(1, 1, 3, 1, 1)  # Default: theta=0.2, thetaU=0.3, N=1e5, prop_invalid=0.1, overlap=0.5
  message(paste("  theta index =", temp[1], "(theta =", thetavec[temp[1]], ")"))
  message(paste("  thetaU index =", temp[2], "(thetaU =", thetaUvec[temp[2]], ")"))
  message(paste("  N index =", temp[3], "(N =", Nvec[temp[3]], ")"))
  message(paste("  prop_invalid index =", temp[4], "(prop_invalid =", prop_invalid_vec[temp[4]], ")"))
  message(paste("  overlap index =", temp[5], "(overlap =", overlap_vec[temp[5]], ")"))
}

# Extract parameters
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
theta_alt = -0.3 
thetaU_alt = -0.3
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs
overlap = overlap_vec[temp[5]]

# Other parameters
pthr = 5e-8 # p-value threshold for instrument selection
NxNy_ratio = 2 # Ratio of sample sizes for X and Y
NxNy_alt_ratio = 1 # Ratio of sample sizes for X and Y
M = 2e5 # Total number of independent SNPs representing the common variants in the genome

# Model parameters for effect size distribution
pi1 = 0.02 * (1 - prop_invalid)
pi3 = 0.01
pi2 = 0.02 * prop_invalid
sigma2x = sigma2y = 5e-5
sigma2u = 1e-4
sigma2x_td = sigma2y_td = (5e-5) - thetaU * thetaUx * sigma2u

# Print configuration
cat("\n========================================\n")
cat("SIMULATION CONFIGURATION\n")
cat("========================================\n")
print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, 
            "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio, "overlap", overlap))

cat("\nMethods Status:\n")
active_methods <- c()
for (method_name in names(SILENCE_METHODS)) {
  status <- ifelse(SILENCE_METHODS[[method_name]], "SILENCED", "ACTIVE")
  cat(sprintf("  %-15s: %s\n", method_name, status))
  if (!SILENCE_METHODS[[method_name]]) {
    active_methods <- c(active_methods, method_name)
  }
}
cat("========================================\n\n")

overlap_adj = (2 * overlap) / (overlap + 1) # to make sure that intersection / union is equal to overlap value.
N_rep = 100
nx = N
ny = N / NxNy_ratio
ny_alt = N / NxNy_alt_ratio

# Create column names based on active methods only
all_methods_map = c(
  MRMode = "MRMode",
  ContMix = "Cont-Mix",
  MRMix = "MR-Mix",
  mode_new_phi1 = "mode_new_phi1",
  MRPRESSO = "MR-PRESSO",
  MRcML = "MR-cML",
  IVW = "IVW",
  median = "median",
  Egger = "Egger",
  IBPRESSO = "IB-PRESSO",
  MR2 = "MR2"
)

# Filter to only active methods
active_method_names <- all_methods_map[active_methods]

# Initialize results matrix with dynamic columns
# MODIFIED: Add extra column for MR2_PIP
n_methods <- length(active_method_names)
has_MR2 <- "MR2" %in% active_method_names

# Base columns: numIV, varX_expl, varY_expl
# Method columns: estimates, SEs, times
# Additional: MR2_PIP (if MR2 is active)
n_cols <- 3 + n_methods * 3 + ifelse(has_MR2, 1, 0)

est = matrix(NA, nrow = N_rep, ncol = n_cols)

# Build column names
base_cols <- c("numIV", "varX_expl", "varY_expl")
estimate_cols <- active_method_names
se_cols <- paste0(active_method_names, "_se")
time_cols <- paste0(active_method_names, "_time")

if (has_MR2) {
  colnames(est) <- c(base_cols, estimate_cols, se_cols, time_cols, "MR2_PIP")
} else {
  colnames(est) <- c(base_cols, estimate_cols, se_cols, time_cols)
}

# Create lookup for column indices
get_col_idx <- function(method_display_name, type = "estimate") {
  base_cols <- 3
  method_idx <- which(active_method_names == method_display_name)
  if (length(method_idx) == 0) return(NULL)
  
  if (type == "estimate") {
    return(base_cols + method_idx)
  } else if (type == "se") {
    return(base_cols + n_methods + method_idx)
  } else if (type == "time") {
    return(base_cols + 2 * n_methods + method_idx)
  } else if (type == "pip" && method_display_name == "MR2" && has_MR2) {
    return(ncol(est))  # Last column
  }
  return(NULL)
}

boot_num = 100

# Main simulation loop
for (repind in 1:N_rep) {
  set.seed(6765 * repind)
  cat("\n========================================\n")
  cat(sprintf("Starting Replication %d/%d\n", repind, N_rep))
  cat("========================================\n")
  
  # Generate SNP indices
  ind1 = sample(M, round(M * pi1))
  ind2 = sample(setdiff(1:M, ind1), round(M * pi2))
  ind3 = sample(setdiff(1:M, c(ind1, ind2)), round(M * pi3))
  
  # Create overlapping indices for second outcome
  shared_ind1 = sample(ind1, round(length(ind1) * overlap_adj))
  ind1_1 = c(shared_ind1, sample(setdiff(1:M, shared_ind1), round(M * pi1) - length(shared_ind1)))
  shared_ind2 = sample(ind2, round(length(ind2) * overlap_adj))
  ind2_1 = c(shared_ind2, sample(setdiff(1:M, c(ind1_1, shared_ind2)), round(M * pi2) - length(shared_ind2)))
  shared_ind3 = sample(ind3, round(length(ind3) * overlap_adj))
  ind3_1 = c(shared_ind3, sample(setdiff(1:M, c(ind1_1, ind2_1, shared_ind3)), round(M * pi3) - length(shared_ind3)))
  
  # Initialize effect vectors
  gamma = phi = phi_alt = alpha = alpha_alt = rep(0, M)
  
  # Generate true effects
  gamma[union(ind1, ind1_1)] = rnorm(length(union(ind1, ind1_1)), mean = 0, sd = sqrt(sigma2x))
  gamma[union(ind2, ind2_1)] = rnorm(length(union(ind2, ind2_1)), mean = 0, sd = sqrt(sigma2x_td))
  alpha[ind2] = rnorm(length(ind2), mean = 0.005, sd = sqrt(sigma2y_td))
  alpha[ind3] = rnorm(length(ind3), mean = 0.005, sd = sqrt(sigma2y))
  alpha_alt[ind2_1] = rnorm(length(ind2_1), mean = 0.003, sd = sqrt(sigma2y_td))
  alpha_alt[ind3_1] = rnorm(length(ind3_1), mean = 0.003, sd = sqrt(sigma2y))
  phi[union(ind2, ind2_1)] = rnorm(length(union(ind2, ind2_1)), mean = 0, sd = sqrt(sigma2u))
  phi_full = phi
  phi_alt[ind2_1] = phi[ind2_1]
  phi[-ind2] = 0
  
  # Calculate true effects
  betax = gamma + thetaUx * phi
  betay = alpha + theta * betax + thetaU * phi
  betay_alt = alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
  
  # Generate observed effects with noise
  betahat_x = betax + rnorm(M, mean = 0, sd = sqrt(1 / nx))
  betahat_y = betay + rnorm(M, mean = 0, sd = sqrt(1 / ny))
  betahat_y_alt = betay_alt + rnorm(M, mean = 0, sd = sqrt(1 / ny_alt))
  
  # Filter the SNPs that reach genome-wide significance
  ind_filter = which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
  numIV = length(ind_filter)
  est[repind, 1] = numIV
  est[repind, 2] = sum(betax[ind_filter]^2)
  est[repind, 3] = sum(betay[ind_filter]^2)
  
  cat(sprintf("Number of IVs: %d\n\n", numIV))
  
  # MR analysis with all methods
  if (numIV > 2) {
    betahat_x.flt = betahat_x[ind_filter]
    betahat_y.flt = betahat_y[ind_filter]
    betahat_y_alt.flt = betahat_y_alt[ind_filter]
    
    # Compute CoHeterogeneity Q statistic for first rep
    if (repind == 1) {
      CoHetQ = coheterogeneity_Q(
        BetaXG = betahat_x.flt,
        BetaYG_matrix = cbind(betahat_y.flt, betahat_y_alt.flt),
        seBetaXG = rep(1 / sqrt(nx), length(betahat_x.flt)),
        seBetaYG_matrix = cbind(
          rep(1 / sqrt(ny), length(betahat_y.flt)), 
          rep(1 / sqrt(ny_alt), length(betahat_y_alt.flt))
        )
      )
      cat("DN: Absolute CoHeterogeneity Q value:", abs(CoHetQ$Q_corr_matrix)[1, 2], "\n\n")
    }
    
    # Create MR input object for standard methods
    mr.obj = mr_input(bx = betahat_x.flt, 
                      bxse = rep(1 / sqrt(nx), length(betahat_x.flt)),
                      by = betahat_y.flt, 
                      byse = rep(1 / sqrt(ny), length(betahat_y.flt)))
    
    # 1. MR Mode
    if (!SILENCE_METHODS$MRMode) {
      cat("Running MRMode... ")
      T0 = proc.time()[3]
      res = mr_mbe(mr.obj, weighting = "weighted")
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("MRMode", "estimate")] = res$Estimate
      est[repind, get_col_idx("MRMode", "se")] = res$StdError
      est[repind, get_col_idx("MRMode", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("MRMode... [SILENCED]\n")
    }
    
    # 2. Contamination mixture
    if (!SILENCE_METHODS$ContMix) {
      cat("Running Cont-Mix... ")
      T0 = proc.time()[3]
      res = mr_conmix(mr.obj)
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("Cont-Mix", "estimate")] = res$Estimate
      CIlength = res$CIUpper - res$CILower
      if (length(CIlength) > 1) cat(" [multimodal] ")
      est[repind, get_col_idx("Cont-Mix", "se")] = sum(CIlength) / 1.96 / 2
      est[repind, get_col_idx("Cont-Mix", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("Cont-Mix... [SILENCED]\n")
    }
    
    # 3. MRMix
    if (!SILENCE_METHODS$MRMix) {
      cat("Running MR-Mix... ")
      theta_temp_vec = seq(-0.5, 0.5, by = 0.01)
      T0 = proc.time()[3]
      res = MRMix(betahat_x.flt, betahat_y.flt, 
                  sx = 1 / sqrt(nx), sy = 1 / sqrt(ny), 
                  theta_temp_vec, pi_init = 0.6, sigma_init = 1e-5)
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("MR-Mix", "estimate")] = res$theta
      est[repind, get_col_idx("MR-Mix", "se")] = res$SE_theta
      est[repind, get_col_idx("MR-Mix", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("MR-Mix... [SILENCED]\n")
    }
    
    # 4. mode_new_phi1
    if (!SILENCE_METHODS$mode_new_phi1) {
      cat("Running mode_new_phi1... ")
      T0 = proc.time()[3]
      res = IBMR::IBMODE(
        BetaXG = betahat_x.flt,
        BetaYG_matrix = cbind(betahat_y.flt, betahat_y_alt.flt),
        seBetaXG = rep(1 / sqrt(nx), length(betahat_x.flt)),
        seBetaYG_matrix = cbind(
          rep(1 / sqrt(ny), length(betahat_y.flt)), 
          rep(1 / sqrt(ny_alt), length(betahat_y_alt.flt))
        ),
        phi = 1,
        n_boot = boot_num,
        alpha = 0.05
      )
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("mode_new_phi1", "estimate")] = res$Estimate.1
      est[repind, get_col_idx("mode_new_phi1", "se")] = res$SE.1
      est[repind, get_col_idx("mode_new_phi1", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("mode_new_phi1... [SILENCED]\n")
    }
    
    # 5. MR-PRESSO
    if (!SILENCE_METHODS$MRPRESSO) {
      cat("Running MR-PRESSO... ")
      T0 = proc.time()[3]
      res <- tryCatch({
        presso <- MRPRESSO::mr_presso(
          BetaOutcome = "BetaOutcome",
          BetaExposure = "BetaExposure",
          SdOutcome = "SdOutcome",
          SdExposure = "SdExposure",
          data = data.frame(
            BetaOutcome = betahat_y.flt,
            BetaExposure = betahat_x.flt,
            SdOutcome = rep(1 / sqrt(ny), length(betahat_y.flt)),
            SdExposure = rep(1 / sqrt(nx), length(betahat_x.flt))
          ),
          OUTLIERtest = TRUE,
          DISTORTIONtest = FALSE,
          NbDistribution = 5000,
          seed = repind,
          SignifThreshold = 0.05
        )
        results <- presso[["Main MR results"]]
        corrected <- results[results$`MR Analysis` == "Outlier-corrected", ]
        if (is.na(corrected$`Causal Estimate`)) {
          corrected <- results[results$`MR Analysis` == "Raw", ]
        }
        list(b = corrected$`Causal Estimate`, se = corrected$Sd)
      }, error = function(e) {
        cat("[ERROR] ")
        list(b = NA, se = NA)
      })
      
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("MR-PRESSO", "estimate")] = res$b
      est[repind, get_col_idx("MR-PRESSO", "se")] = res$se
      est[repind, get_col_idx("MR-PRESSO", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("MR-PRESSO... [SILENCED]\n")
    }
    
    # 6. MR-cML
    if (!SILENCE_METHODS$MRcML) {
      cat("Running MR-cML... ")
      T0 = proc.time()[3]
      res = mr_cML(mr.obj, n = ny)
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("MR-cML", "estimate")] = res$Estimate
      est[repind, get_col_idx("MR-cML", "se")] = res$StdError
      est[repind, get_col_idx("MR-cML", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("MR-cML... [SILENCED]\n")
    }
    
    # 7. IVW
    if (!SILENCE_METHODS$IVW) {
      cat("Running IVW... ")
      T0 = proc.time()[3]
      res = mr_ivw(mr.obj)
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("IVW", "estimate")] = res$Estimate
      est[repind, get_col_idx("IVW", "se")] = res$StdError
      est[repind, get_col_idx("IVW", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("IVW... [SILENCED]\n")
    }
    
    # 8. Median
    if (!SILENCE_METHODS$median) {
      cat("Running median... ")
      T0 = proc.time()[3]
      res = mr_median(mr.obj)
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("median", "estimate")] = res$Estimate
      est[repind, get_col_idx("median", "se")] = res$StdError
      est[repind, get_col_idx("median", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("median... [SILENCED]\n")
    }
    
    # 9. Egger
    if (!SILENCE_METHODS$Egger) {
      cat("Running Egger... ")
      T0 = proc.time()[3]
      res = mr_egger(mr.obj)
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("Egger", "estimate")] = res$Estimate
      est[repind, get_col_idx("Egger", "se")] = res$StdError.Est
      est[repind, get_col_idx("Egger", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("Egger... [SILENCED]\n")
    }
    
    # 10. IB-PRESSO
    if (!SILENCE_METHODS$IBPRESSO) {
      cat("Running IB-PRESSO... ")
      T0 = proc.time()[3]
      res <- tryCatch({
        ibpresso <- IBMR::IBPRESSO(
          BetaOutcome = "BetaOutcome",
          BetaExposure = "BetaExposure",
          BetaAux = "BetaAux",
          SdOutcome = "SdOutcome",
          SdExposure = "SdExposure",
          SdAux = "SdAux",
          data = data.frame(
            BetaOutcome = betahat_y.flt,
            BetaExposure = betahat_x.flt,
            BetaAux = betahat_y_alt.flt,
            SdOutcome = rep(1 / sqrt(ny), length(betahat_y.flt)),
            SdExposure = rep(1 / sqrt(nx), length(betahat_x.flt)),
            SdAux = rep(1 / sqrt(ny_alt), length(betahat_y_alt.flt))
          ),
          OUTLIERtest = TRUE,
          DISTORTIONtest = FALSE,
          NbDistribution = 5000,
          seed = repind,
          SignifThreshold = 0.05
        )
        b_val <- if (!is.null(ibpresso$corrected_beta)) ibpresso$corrected_beta else ibpresso$raw_beta
        se_val <- if (!is.null(ibpresso$corrected_se)) ibpresso$corrected_se else ibpresso$raw_se
        list(b = b_val, se = se_val)
      }, error = function(e) {
        cat("[ERROR] ")
        list(b = NA, se = NA)
      })
      
      T1 = proc.time()[3]
      runtime = T1 - T0
      est[repind, get_col_idx("IB-PRESSO", "estimate")] = res$b
      est[repind, get_col_idx("IB-PRESSO", "se")] = res$se
      est[repind, get_col_idx("IB-PRESSO", "time")] = runtime
      cat(sprintf("Done! (%.3f sec)\n", runtime))
      rm(res)
    } else {
      cat("IB-PRESSO... [SILENCED]\n")
    }
    
    # 11. MR2 (WITH PIP EXTRACTION - CORRECTED)
    if (!SILENCE_METHODS$MR2) {
      cat("Running MR2... ")
      T0 = proc.time()[3]
      res <- tryCatch({
        # Prepare matrices
        betaHat_Y_matrix <- cbind(betahat_y.flt, betahat_y_alt.flt)
        betaHat_X_matrix <- matrix(betahat_x.flt, ncol = 1)
        
        # Validation checks
        if (any(is.na(betaHat_Y_matrix)) || any(is.na(betaHat_X_matrix))) {
          stop("NA values in input")
        }
        if (any(!is.finite(betaHat_Y_matrix)) || any(!is.finite(betaHat_X_matrix))) {
          stop("Non-finite values in input")
        }
        if (nrow(betaHat_Y_matrix) < 3) {
          stop("Insufficient IVs")
        }
        
        # Run MR2
        MR2_fit <- MR2(
          betaHat_Y = betaHat_Y_matrix, 
          betaHat_X = betaHat_X_matrix, 
          EVgamma = 0.5,
          niter = 7500,
          burnin = 2500,
          thin = 5,
          monitor = 1000
        )
        
        # Post-processing
        PostProc_res <- PostProc(MR2_fit, betaHat_Y_matrix, betaHat_X_matrix, alpha = 0.05)
        
        # Extract point estimate for outcome 1 (first column)
        b_val <- as.numeric(PostProc_res$thetaPost[1, 1])
        
        # Calculate SE from 95% credible interval
        ci_lower <- as.numeric(PostProc_res$thetaPost_CI[1, 1, 1])
        ci_upper <- as.numeric(PostProc_res$thetaPost_CI[1, 1, 2])
        se_val <- (ci_upper - ci_lower) / (2 * qnorm(0.975))
        
        # Extract PIP from postMean$gamma (posterior inclusion probability)
        # This is a (n_exposures × n_outcomes) matrix
        # For exposure 1 → outcome 1, use [1, 1]
        pip_val <- as.numeric(MR2_fit$postMean$gamma[1, 1])
        
        # Validate outputs
        if (any(is.na(c(b_val, se_val, pip_val))) ||
            any(!is.finite(c(b_val, se_val, pip_val))) ||
            any(c(length(b_val), length(se_val), length(pip_val)) != 1)) {
          stop("Non-finite or non-scalar output")
        }
        
        # Clean up
        rm(MR2_fit, PostProc_res, betaHat_Y_matrix, betaHat_X_matrix)
        
        list(b = b_val, se = se_val, pip = pip_val, success = TRUE)
        
      }, error = function(e) {
        err_msg <- substr(e$message, 1, 40)
        cat(sprintf("\n[FAILED: %s, rep=%d, IVs=%d, seed=%d]\n", 
                    err_msg, repind, numIV, 6765 * repind))
        list(b = NA, se = NA, pip = NA, success = FALSE)
      })
      
      T1 = proc.time()[3]
      est[repind, get_col_idx("MR2", "estimate")] = res$b
      est[repind, get_col_idx("MR2", "se")] = res$se
      est[repind, get_col_idx("MR2", "time")] = T1 - T0
      est[repind, get_col_idx("MR2", "pip")] = res$pip  # Store PIP
      
      if (res$success) {
        cat(sprintf("Done! (%.2f sec, PIP=%.3f)\n", T1 - T0, res$pip))
      }
      rm(res)
    }
    
    # Print summary of runtimes for this replication
    cat("\n--- Runtime Summary ---\n")
    for (method_name in active_method_names) {
      method_time <- est[repind, get_col_idx(method_name, "time")]
      if (!is.na(method_time)) {
        if (method_name == "MR2" && has_MR2) {
          pip_val <- est[repind, get_col_idx("MR2", "pip")]
          cat(sprintf("%-15s: %8.3f sec (PIP=%.3f)\n", method_name, method_time, pip_val))
        } else {
          cat(sprintf("%-15s: %8.3f sec\n", method_name, method_time))
        }
      }
    }
    cat("----------------------\n")
  }
  
  # Save intermediate results every 5 iterations
  if (repind %% 5 == 0) {
    cat(sprintf("\n==> Saving intermediate results at rep %d\n", repind))
    save(est, file = file.path(
      code_dir,
      paste0("DN_est_theta", theta,
             "_thetaU", thetaU,
             "_N", format(N, scientific = FALSE),
             "_prop_invalid", prop_invalid,
             "_overlap", overlap,
             ".rda")
    ))
    cat("\n=== CUMULATIVE STATISTICS (Reps 1 to", repind, ") ===\n")
    for (method_name in active_method_names) {
      est_col <- get_col_idx(method_name, "estimate")
      se_col <- get_col_idx(method_name, "se")
      
      # Get non-NA values up to current replication
      estimates <- est[1:repind, est_col]
      ses <- est[1:repind, se_col]
      
      valid_est <- estimates[!is.na(estimates)]
      valid_se <- ses[!is.na(ses)]
      
      if (length(valid_est) > 0) {
        cum_mean_est <- mean(valid_est)
        cum_sd_est <- sd(valid_est)
        cum_mean_se <- mean(valid_se)
        success_rate <- (length(valid_est) / repind) * 100
        
        # Add PIP statistics for MR2
        if (method_name == "MR2" && has_MR2) {
          pip_col <- get_col_idx("MR2", "pip")
          pips <- est[1:repind, pip_col]
          valid_pip <- pips[!is.na(pips)]
          if (length(valid_pip) > 0) {
            mean_pip <- mean(valid_pip)
            cat(sprintf("%-15s: Mean=%.4f (SD=%.4f), SE=%.4f, PIP=%.3f, Success=%.1f%%\n", 
                        method_name, cum_mean_est, cum_sd_est, cum_mean_se, mean_pip, success_rate))
          } else {
            cat(sprintf("%-15s: Mean=%.4f (SD=%.4f), SE=%.4f, Success=%.1f%%\n", 
                        method_name, cum_mean_est, cum_sd_est, cum_mean_se, success_rate))
          }
        } else {
          cat(sprintf("%-15s: Mean=%.4f (SD=%.4f), SE=%.4f, Success=%.1f%%\n", 
                      method_name, cum_mean_est, cum_sd_est, cum_mean_se, success_rate))
        }
      } else {
        cat(sprintf("%-15s: [No successful runs]\n", method_name))
      }
    }
    cat(sprintf("True theta = %.2f\n", theta))
    cat("=============================================\n\n")
  }
}

# Save final results
output_file <- file.path(
  code_dir,
  paste0("DN_est_theta", theta,
         "_thetaU", thetaU,
         "_N", format(N, scientific = FALSE),
         "_prop_invalid", prop_invalid,
         "_overlap", overlap,
         ".rda")
)

save(est, file = output_file)
cat("\n========================================\n")
cat("SIMULATION COMPLETE!\n")
cat("========================================\n")
message(paste("Results saved to:", output_file))

# Print overall runtime statistics
cat("\n=== OVERALL RUNTIME STATISTICS ===\n")
for (method_name in active_method_names) {
  method_times <- est[, get_col_idx(method_name, "time")]
  method_times <- method_times[!is.na(method_times)]
  if (length(method_times) > 0) {
    cat(sprintf("%-15s: Mean=%.3f sec, Median=%.3f sec, Total=%.2f min\n", 
                method_name, 
                mean(method_times), 
                median(method_times),
                sum(method_times)/60))
  }
}

# Print MR2 PIP summary if available
if (has_MR2) {
  cat("\n=== MR2 PIP SUMMARY ===\n")
  pip_col <- get_col_idx("MR2", "pip")
  all_pips <- est[, pip_col]
  valid_pips <- all_pips[!is.na(all_pips)]
  
  if (length(valid_pips) > 0) {
    cat(sprintf("Mean PIP: %.3f\n", mean(valid_pips)))
    cat(sprintf("Median PIP: %.3f\n", median(valid_pips)))
    cat(sprintf("SD PIP: %.3f\n", sd(valid_pips)))
    cat(sprintf("Min PIP: %.3f\n", min(valid_pips)))
    cat(sprintf("Max PIP: %.3f\n", max(valid_pips)))
    cat(sprintf("Proportion PIP > 0.5: %.1f%%\n", mean(valid_pips > 0.5) * 100))
    cat(sprintf("Proportion PIP > 0.8: %.1f%%\n", mean(valid_pips > 0.8) * 100))
  }
}

cat("===================================\n")