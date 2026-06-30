#————————————————————————————————————————————#
#  Simulation comparing univariate MR‐PRESSO
#             vs. IB‐MR‐PRESSO (beta, p‐values, distortion)
#————————————————————————————————————————————#

rm(list=ls())
library(data.table); library(dplyr)
library(MASS)
library(MendelianRandomization)
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
if (!requireNamespace("MRPRESSO", quietly=TRUE)) devtools::install_github("rondolab/MR-PRESSO", force=TRUE)
library(MRPRESSO)
library(IBMR)  # IBPRESSO is available as IBMR::IBPRESSO
#source("IBPresso.R")
#————————————————————————————————————————————#
#  Fixed parameters
#————————————————————————————————————————————#
thetavec         <- c(0.2, 0, -0.2)
thetaUvec        <- c(0.3, 0.5)
Nvec             <- c(5e4,8e4,1e5,1.5e5,2e5,5e5,1e6)
prop_invalid_vec <- c(0.1,0.3,0.5)
overlap_vec      <- c(0,0.25,0.5,0.75,1)

args         <- as.integer(commandArgs(trailingOnly=TRUE))
theta        <- thetavec[args[1]]
thetaU    = thetaUx    <- thetaUvec[args[2]]
N            <- Nvec[args[3]]
prop_invalid <- prop_invalid_vec[args[4]]
overlap      <- overlap_vec[args[5]]

pthr        <- 5e-8
NxNy_ratio  <- 2
M           <- 2e5
pi1 <- 0.02*(1-prop_invalid)
pi2 <- 0.02*prop_invalid
pi3 <- 0.01
sigma2x    <- sigma2y    <- 5e-5
sigma2x_td <- sigma2y_td <- 5e-5
overlap_adj <- (2*overlap)/(1+overlap)

N_rep  <- 100
nx     <- N
ny     <- N/NxNy_ratio
ny_alt <- N

# Results matrix now includes SEs
est <- matrix(NA, nrow=N_rep, ncol=11)
colnames(est) <- c(
  "numIV","varX_expl","varY_expl",
  "uni_beta","uni_se","uni_p","uni_time",
  "ib_beta","ib_se","ib_p","ib_time"
)

for(repind in 1:N_rep) {
  set.seed(1234 + repind)
  
  # simulate SNP‐sets and true effects
  ind1 <- sample(M, round(M*pi1))
  ind2 <- sample(setdiff(1:M,ind1), round(M*pi2))
  ind3 <- sample(setdiff(1:M,c(ind1,ind2)), round(M*pi3))
  shared1 <- sample(ind1, round(length(ind1)*overlap_adj))
  ind1_1 <- c(shared1,
              sample(setdiff(1:M,shared1), length(ind1)-length(shared1)))
  shared2 <- sample(ind2, round(length(ind2)*overlap_adj))
  ind2_1 <- c(shared2,
              sample(setdiff(1:M,c(ind1_1,shared2)), length(ind2)-length(shared2)))
  shared3 <- sample(ind3, round(length(ind3)*overlap_adj))
  ind3_1 <- c(shared3,
              sample(setdiff(1:M,c(ind1_1,ind2_1,shared3)), length(ind3)-length(shared3)))
  
  phi <- rnorm(M)*0
  gamma <- alpha <- alpha_alt <- rep(0,M)
  gamma[ union(ind1,ind1_1) ] <- rnorm(length(union(ind1,ind1_1)),0,sqrt(sigma2x))
  gamma[ union(ind2,ind2_1) ] <- rnorm(length(union(ind2,ind2_1)),0,sqrt(sigma2x_td))
  alpha[ind2]       <- rnorm(length(ind2), 0, sqrt(sigma2y_td))
  alpha[ind3]       <- rnorm(length(ind3), 0, sqrt(sigma2y))
  alpha_alt[ind2_1] <- rnorm(length(ind2_1),0, sqrt(sigma2y_td))
  alpha_alt[ind3_1] <- rnorm(length(ind3_1),0, sqrt(sigma2y))
  
  betax     <- gamma + thetaU * phi
  betay     <- alpha + theta  * betax + thetaU * phi
  betay_alt <- alpha_alt + 0.3 * betax + 0.3 * phi
  
  betahat_x     <- betax     + rnorm(M,0,sqrt(1/nx))
  betahat_y     <- betay     + rnorm(M,0,sqrt(1/ny))
  betahat_y_alt <- betay_alt + rnorm(M,0,sqrt(1/ny_alt))
  
  sel   <- which(2*pnorm(-sqrt(nx)*abs(betahat_x)) < pthr)
  numIV <- length(sel)
  est[repind,"numIV"]     <- numIV
  est[repind,"varX_expl"] <- sum(betax[sel]^2)
  est[repind,"varY_expl"] <- sum(betay[sel]^2)
  
  if(numIV > 2) {
    df <- data.frame(
      BetaExposure = betahat_x[sel],
      SdExposure   = rep(1/sqrt(nx), numIV),
      BetaOutcome  = betahat_y[sel],
      SdOutcome    = rep(1/sqrt(ny), numIV)
    )
    df_ib <- df %>% mutate(
      BetaAux = betahat_y_alt[sel],
      SdAux   = rep(1/sqrt(ny_alt), numIV)
    )
    
    # Univariate MR‐PRESSO (with outlier correction)
    t0 <- proc.time()[3]
    uni <- mr_presso(
      BetaOutcome    = "BetaOutcome",
      BetaExposure   = "BetaExposure",
      SdOutcome      = "SdOutcome",
      SdExposure     = "SdExposure",
      data           = df,
      OUTLIERtest    = TRUE,
      DISTORTIONtest = FALSE,
      NbDistribution = 5000,
      seed           = 2025
    )
    t1 <- proc.time()[3]
    
    mr_results <- uni[["Main MR results"]]
    corrected <- mr_results[mr_results$`MR Analysis` == "Outlier-corrected", ]
    
    if (is.na(corrected$`Causal Estimate`)) {
      corrected <- mr_results[mr_results$`MR Analysis` == "Raw", ]
    }
    
    est[repind,"uni_beta"] <- corrected$`Causal Estimate`
    est[repind,"uni_se"  ] <- corrected$Sd
    est[repind,"uni_p"   ] <- uni[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    est[repind,"uni_time"] <- t1 - t0
    
    # IB‐MR‐PRESSO
    t0 <- proc.time()[3]
    ib <- mr_presso_ib(
      BetaOutcome    = "BetaOutcome",
      BetaExposure   = "BetaExposure",
      BetaAux        = "BetaAux",
      SdOutcome      = "SdOutcome",
      SdExposure     = "SdExposure",
      SdAux          = "SdAux",
      data           = df_ib,
      OUTLIERtest    = TRUE,
      DISTORTIONtest = FALSE,
      NbDistribution = 5000,
      seed           = 2025
    )
    t1 <- proc.time()[3]
    
    est[repind,"ib_beta"] <- if (!is.na(ib$corrected_beta)) ib$corrected_beta else ib$raw_beta
    est[repind,"ib_se"  ] <- if (!is.na(ib$corrected_se))    ib$corrected_se    else ib$raw_se
    est[repind,"ib_p"   ] <- ib$p_value
    est[repind,"ib_time"] <- t1 - t0
  }
  
  if(repind %% 1 == 0) message("Done replicate ", repind)
  if(repind %% 1 == 0){
    print(est[1:repind,])
  }
  if(repind %% 5 == 0) save(est,
                            file = sprintf(
                              "../../../results/simulation_results/presso_comp/sim_MRpressos_BI_beta_p_se_theta%s_U%s_N%s_invalid%s_overlap%s.rda",
                              theta, thetaU, N, prop_invalid, overlap
                            ))
}

save(est,
     file = sprintf(
       "../../../results/simulation_results/presso_comp/sim_MRpressos_BI_beta_p_se_theta%s_U%s_N%s_invalid%s_overlap%s.rda",
       theta, thetaU, N, prop_invalid, overlap
     ))
