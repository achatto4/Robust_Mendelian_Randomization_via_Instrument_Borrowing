# Summary level simulations with directional pleiotropy but InSIDE assumption not satisfied (NULL setting)
rm(list = ls())

library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(penalized)
library(ks)

library(IBMR)  # coheterogeneity_Q is available as IBMR::coheterogeneity_Q
library(IBMR)  # IBMODE is available as IBMR::IBMODE

# Null setting for type-I error evaluation
theta <- 0
phi_val <- 1

# Secondary trait effect grid (varied in null setting)
theta_alt_vec <- c(-0.2, -0.1, 0, 0.1, 0.2)
thetaU_alt <- 0.3

thetaUvec <- c(0.3)
Nvec <- c(5e4, 8e4, 1e5, 2e5)
prop_invalid_vec <- c(0.1, 0.3, 0.5, 0.7)
overlap_vec <- c(0, 0.25, 0.5, 0.75, 1)
alpha_levels <- c(0.005, 0.01, 0.05, 1)

# Backward compatible argument format: theta_alt_idx thetaU_idx N_idx prop_invalid_idx phi_idx overlap_idx
args <- as.integer(commandArgs(trailingOnly = TRUE))
if (length(args) < 6 || any(is.na(args))) {
  stop("Expected 6 integer args: theta_alt_idx thetaU_idx N_idx prop_invalid_idx phi_idx overlap_idx")
}

theta_alt <- theta_alt_vec[args[1]]
thetaU <- thetaUx <- thetaUvec[1]
N <- Nvec[args[3]]
prop_invalid <- prop_invalid_vec[args[4]]
overlap <- overlap_vec[args[6]]

pthr <- 5e-8
NxNy_ratio <- 2
NxNy_alt_ratio <- 1
M <- 2e5
N_rep <- 1000
boot_num <- 100

pi1 <- 0.02 * (1 - prop_invalid)
pi3 <- 0.01
pi2 <- 0.02 * prop_invalid
sigma2x <- 5e-5
sigma2y <- 5e-5

sigma2u <- 1e-4
sigma2x_td <- sigma2y_td <- (5e-5) - thetaU * thetaUx * sigma2u

print(paste(
  "Scenario", "DN",
  "theta", theta,
  "phi", phi_val,
  "N", N,
  "thetaU", thetaU,
  "theta_alt", theta_alt,
  "prop_invalid", prop_invalid,
  "overlap", overlap
))

overlap_adj <- (2 * overlap) / (overlap + 1)
nx <- N
ny <- N / NxNy_ratio
ny_alt <- N / NxNy_alt_ratio

est <- matrix(NA, nrow = N_rep, ncol = 3 + 2 * 3)
mr_methods <- c("MRMode", "mode_new")
colnames(est) <- c(
  "numIV", "varX_expl", "varY_expl",
  mr_methods,
  paste0(mr_methods, "_se"),
  paste0(mr_methods, "_time")
)

for (repind in seq_len(N_rep)) {
  set.seed(6765 * repind)

  ind1 <- sample(M, round(M * pi1))
  ind2 <- sample(setdiff(1:M, ind1), round(M * pi2))
  ind3 <- sample(setdiff(1:M, c(ind1, ind2)), round(M * pi3))

  shared_ind1 <- sample(ind1, round(length(ind1) * overlap_adj))
  ind1_1 <- c(shared_ind1, sample(setdiff(1:M, shared_ind1), round(M * pi1) - length(shared_ind1)))
  shared_ind2 <- sample(ind2, round(length(ind2) * overlap_adj))
  ind2_1 <- c(shared_ind2, sample(setdiff(1:M, c(ind1_1, shared_ind2)), round(M * pi2) - length(shared_ind2)))
  shared_ind3 <- sample(ind3, round(length(ind3) * overlap_adj))
  ind3_1 <- c(shared_ind3, sample(setdiff(1:M, c(ind1_1, ind2_1, shared_ind3)), round(M * pi3) - length(shared_ind3)))

  gamma <- phi <- phi_alt <- alpha <- alpha_alt <- rep(0, M)

  gamma[union(ind1, ind1_1)] <- rnorm(length(union(ind1, ind1_1)), mean = 0, sd = sqrt(sigma2x))
  gamma[union(ind2, ind2_1)] <- rnorm(length(union(ind2, ind2_1)), mean = 0, sd = sqrt(sigma2x_td))

  if (TRUE) {
    alpha[ind2] <- rnorm(length(ind2), mean = 0.005, sd = sqrt(sigma2y_td))
    alpha[ind3] <- rnorm(length(ind3), mean = 0.005, sd = sqrt(sigma2y))
    alpha_alt[ind2_1] <- rnorm(length(ind2_1), mean = 0.003, sd = sqrt(sigma2y_td))
    alpha_alt[ind3_1] <- rnorm(length(ind3_1), mean = 0.003, sd = sqrt(sigma2y))
  } else {
    alpha[ind2] <- rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
    alpha[ind3] <- rnorm(length(ind3), mean = 0, sd = sqrt(sigma2y))
    alpha_alt[ind2_1] <- rnorm(length(ind2_1), mean = 0, sd = sqrt(sigma2y_td))
    alpha_alt[ind3_1] <- rnorm(length(ind3_1), mean = 0, sd = sqrt(sigma2y))
  }

  if (!FALSE) {
    phi[union(ind2, ind2_1)] <- rnorm(length(union(ind2, ind2_1)), mean = 0, sd = sqrt(sigma2u))
    phi_alt[ind2_1] <- phi[ind2_1]
    phi[-ind2] <- 0
  }

  betax <- gamma + thetaUx * phi
  betay <- alpha + theta * betax + thetaU * phi
  betay_alt <- alpha_alt + theta_alt * betax + thetaU_alt * phi_alt

  betahat_x <- betax + rnorm(M, mean = 0, sd = sqrt(1 / nx))
  betahat_y <- betay + rnorm(M, mean = 0, sd = sqrt(1 / ny))
  betahat_y_alt <- betay_alt + rnorm(M, mean = 0, sd = sqrt(1 / ny_alt))

  ind_filter <- which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
  numIV <- length(ind_filter)

  est[repind, 1] <- numIV
  est[repind, 2] <- sum(betax[ind_filter]^2)
  est[repind, 3] <- sum(betay[ind_filter]^2)

  if (numIV > 2) {
    bx <- betahat_x[ind_filter]
    by <- betahat_y[ind_filter]
    by_alt <- betahat_y_alt[ind_filter]

    if (repind == 1) {
      cohet <- coheterogeneity_Q(
        BetaXG = bx,
        BetaYG_matrix = cbind(by, by_alt),
        seBetaXG = rep(1 / sqrt(nx), length(bx)),
        seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(by)), rep(1 / sqrt(ny_alt), length(by_alt)))
      )
      cat("DN: Absolute CoHeterogeneity Q value:", abs(cohet$Q_corr_matrix)[1, 2], "\n")
    }

    mr.obj <- mr_input(
      bx = bx,
      bxse = rep(1 / sqrt(nx), length(bx)),
      by = by,
      byse = rep(1 / sqrt(ny), length(by))
    )

    t0 <- proc.time()[3]
    res <- mr_mbe(mr.obj, weighting = "weighted", phi = phi_val)
    t1 <- proc.time()[3]
    est[repind, 4] <- res$Estimate
    est[repind, 6] <- res$StdError
    est[repind, 8] <- t1 - t0

    t0 <- proc.time()[3]
    res <- IBMR::IBMODE(
      BetaXG = bx,
      BetaYG_matrix = cbind(by, by_alt),
      seBetaXG = rep(1 / sqrt(nx), length(bx)),
      seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(by)), rep(1 / sqrt(ny_alt), length(by_alt))),
      phi = phi_val,
      n_boot = boot_num,
      alpha = 0.05
    )
    t1 <- proc.time()[3]
    est[repind, 5] <- res$Estimate.1
    est[repind, 7] <- res$SE.1
    est[repind, 9] <- t1 - t0
  }
}

mode_p <- rep(NA_real_, N_rep)
ib_mode_p <- rep(NA_real_, N_rep)

valid_mode <- !is.na(est[, 4]) & !is.na(est[, 6]) & est[, 6] > 0
valid_ib <- !is.na(est[, 5]) & !is.na(est[, 7]) & est[, 7] > 0

mode_p[valid_mode] <- 2 * pnorm(-abs(est[valid_mode, 4] / est[valid_mode, 6]))
ib_mode_p[valid_ib] <- 2 * pnorm(-abs(est[valid_ib, 5] / est[valid_ib, 7]))

type1_summary <- do.call(rbind, lapply(alpha_levels, function(a) {
  data.frame(
    scenario = "DN",
    theta = theta,
    theta_alt = theta_alt,
    phi = phi_val,
    alpha = a,
    n_rep = N_rep,
    n_valid_mode = sum(!is.na(mode_p)),
    n_valid_ib_mode = sum(!is.na(ib_mode_p)),
    type1_mrmode = mean(mode_p < a, na.rm = TRUE),
    type1_ib_mode = mean(ib_mode_p < a, na.rm = TRUE),
    thetaU = thetaU,
    N = N,
    prop_invalid = prop_invalid,
    overlap = overlap
  )
}))

base_name <- paste0(
  "DN_NULL_theta0_thetaAlt", theta_alt, "_phi1_thetaU", thetaU,
  "_N", format(N, scientific = FALSE),
  "_prop_invalid", prop_invalid,
  "_overlap", overlap
)

save(
  est,
  type1_summary,
  file = file.path("../../../results/simulation_results/mode_comp", paste0(base_name, ".rda"))
)

write.csv(
  type1_summary,
  file = file.path("../../../results/simulation_results/mode_comp", paste0(base_name, "_type1_summary.csv")),
  row.names = FALSE
)
