# Summary level simulations with balanced pleiotropy and InSIDE assumption satisfied
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(penalized)
library(ks)
library(MRPRESSO)
library(IBMR)  # coheterogeneity_Q is available as IBMR::coheterogeneity_Q
library(IBMR)  # IBMODE is available as IBMR::IBMODE
library(IBMR)  # IBPRESSO is available as IBMR::IBPRESSO

thetavec = c(0.2, 0, -0.2, 0.1, -0.1)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5, 1.5e5, 2e5, 5e5, 1e6) # 1:7
prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7)
overlap_vec = c(0.5, 0.75, 1)

temp = as.integer(commandArgs(trailingOnly = TRUE))
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
theta_alt = 0.3
thetaU_alt = 0.3
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs
overlap = overlap_vec[temp[5]]

pthr = 5e-8 # p-value threshold for instrument selection
NxNy_ratio = 2 # Ratio of sample sizes for X and Y
NxNy_alt_ratio = 1 # Ratio of sample sizes for X and Y
M = 2e5 # Total number of independent SNPs representing the common variants in the genome

# Model parameters for effect size distribution
pi1=0.02*(1-prop_invalid); pi3=0.01
pi2=0.02*prop_invalid;
sigma2x = sigma2y = 5e-5; sigma2u = 1e-4; sigma2x_td = sigma2y_td = (5e-5)-thetaU*thetaUx*sigma2u

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio, "overlap", overlap))
overlap_adj  = (2*overlap)/(overlap+1) # to make sure that intersection / union is equal to overlap value.
N_rep = 100
nx = N; ny = N/NxNy_ratio; ny_alt = N/NxNy_alt_ratio
est = matrix(NA, nrow = N_rep, ncol = 3+10*3)
mr_methods = c("MRMode","Cont-Mix","MR-Mix","mode_new_phi1","MR-PRESSO",
               "MR-cML","IVW","median","Egger","IB-PRESSO")
colnames(est) = c("numIV", "varX_expl","varY_expl", mr_methods, paste0(mr_methods,"_se"), paste0(mr_methods,"_time"))
# varX_expl: variance of X explained by IVs; varY_expl: variance of Y explained by IVs
boot_num = 100
for (repind in 1:N_rep){
  set.seed(6765 * repind)

  ind1 = sample(M, round(M * pi1))
  ind2 = sample(setdiff(1:M, ind1), round(M * pi2))
  ind3 = sample(setdiff(1:M, c(ind1, ind2)), round(M * pi3))

  shared_ind1 = sample(ind1, round(length(ind1) * overlap_adj))
  ind1_1 = c(shared_ind1, sample(setdiff(1:M, shared_ind1), round(M * pi1) - length(shared_ind1)))
  shared_ind2 = sample(ind2, round(length(ind2) * overlap_adj))
  ind2_1 = c(shared_ind2, sample(setdiff(1:M, c(ind1_1, shared_ind2)), round(M * pi2) - length(shared_ind2)))
  shared_ind3 = sample(ind3, round(length(ind3) * overlap_adj))
  ind3_1 = c(shared_ind3, sample(setdiff(1:M, c(ind1_1, ind2_1, shared_ind3)), round(M * pi3) - length(shared_ind3)))

  gamma = phi = phi_alt = alpha = alpha_alt = rep(0, M)

  gamma[union(ind1,ind1_1)] = rnorm(length(union(ind1,ind1_1)), mean = 0, sd = sqrt(sigma2x))
  gamma[union(ind2,ind2_1)] = rnorm(length(union(ind2,ind2_1)), mean = 0, sd = sqrt(sigma2x_td))
  alpha[ind2] = rnorm(length(ind2), mean = 0.005, sd = sqrt(sigma2y_td))
  alpha[ind3] = rnorm(length(ind3), mean = 0.005, sd = sqrt(sigma2y))
  alpha_alt[ind2_1] = rnorm(length(ind2_1), mean = 0.003, sd = sqrt(sigma2y_td))
  alpha_alt[ind3_1] = rnorm(length(ind3_1), mean = 0.003, sd = sqrt(sigma2y))
  phi[union(ind2,ind2_1)] = rnorm(length(union(ind2,ind2_1)), mean = 0, sd = sqrt(sigma2u))
  phi_full = phi
  phi_alt[ind2_1] = phi[ind2_1]
  phi[-ind2] = 0
  betax = gamma + thetaUx * phi
  betay = alpha + theta * betax + thetaU * phi
  betay_alt = alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
  betahat_x = betax + rnorm(M, mean = 0, sd = sqrt(1/nx))
  betahat_y = betay + rnorm(M, mean = 0, sd = sqrt(1/ny))
  betahat_y_alt = betay_alt + rnorm(M, mean = 0, sd = sqrt(1/ny_alt))

  # Filter the SNPs that reach genome-wide significance in the study associated with X.
  ind_filter = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<pthr)
  numIV = length(ind_filter)
  est[repind,1] = numIV
  est[repind,2] = sum(betax[ind_filter]^2)
  est[repind,3] = sum(betay[ind_filter]^2)

  # MR analysis with 2 methods
  if (numIV>2){
    betahat_x.flt = betahat_x[ind_filter]
    betahat_y.flt = betahat_y[ind_filter]
    betahat_y_alt.flt = betahat_y_alt[ind_filter]

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

      cat("DN: Absolute CoHeterogeneity Q value:", abs(CoHetQ$Q_corr_matrix)[1, 2], "\n")
    }

    mr.obj = mr_input(bx = betahat_x.flt, bxse = rep(1/sqrt(nx), length(betahat_x.flt)),
                      by = betahat_y.flt, byse = rep(1/sqrt(ny), length(betahat_y.flt)))

    # 1. mode
    T0 = proc.time()[3]
    res = mr_mbe(mr.obj, weighting = "weighted")
    T1 = proc.time()[3]
    est[repind,3+1] = res$Estimate
    est[repind,3+11] = res$StdError
    est[repind,3+21] = T1-T0
    rm(res)

    # 2. contamination mixture
    T0 = proc.time()[3]
    res = mr_conmix(mr.obj)
    T1 = proc.time()[3]
    est[repind,3+2] = res$Estimate
    CIlength = res$CIUpper-res$CILower
    if (length(CIlength)>1) print(paste("Repind",repind,"conmix multimodal"))
    est[repind,3+12] = sum(CIlength)/1.96/2 ## Caution: this may be problematic
    est[repind,3+22] = T1-T0
    rm(res)

    # 3. MRMix
    theta_temp_vec = seq(-0.5,0.5,by=0.01)
    T0 = proc.time()[3]
    res = MRMix(betahat_x.flt, betahat_y.flt, sx=1/sqrt(nx), sy=1/sqrt(ny), theta_temp_vec, pi_init = 0.6, sigma_init = 1e-5)
    T1 = proc.time()[3]
    est[repind,3+3] = res$theta
    est[repind,3+13] = res$SE_theta
    est[repind,3+23] = T1-T0
    rm(res)

    # 4. mode_new_phi1
    T0 = proc.time()[3]
    res = IBMR::IBMODE(
      BetaXG = betahat_x.flt,
      BetaYG_matrix = cbind(betahat_y.flt, betahat_y_alt.flt),
      seBetaXG = rep(1 / sqrt(nx), length(betahat_x.flt)),
      seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(betahat_y.flt)), rep(
        1 / sqrt(ny_alt), length(betahat_y_alt.flt)
      )),
      phi = 1,
      n_boot = boot_num,
      alpha = 0.05
    )
    T1 = proc.time()[3]
    est[repind,3+4] = res$Estimate.1
    est[repind,3+14] = res$SE.1
    est[repind,3+24] = T1-T0
    rm(res)

    # 5. MR-PRESSO
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
          SdOutcome = rep(1/sqrt(ny), length(betahat_y.flt)),
          SdExposure = rep(1/sqrt(nx), length(betahat_x.flt))
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
    }, error=function(e) list(b=NA,se=NA))

    T1 = proc.time()[3]
    est[repind,3+5] = res$b
    est[repind,3+15] = res$se
    est[repind,3+25] = T1-T0
    rm(res)

    # 6. MR-cML
    T0 = proc.time()[3]
    res = mr_cML(mr.obj, n = ny)
    T1 = proc.time()[3]
    est[repind,3+6] = res$Estimate
    est[repind,3+16] = res$StdError
    est[repind,3+26] = T1-T0
    rm(res)

    # 7. IVW
    T0 = proc.time()[3]
    res = mr_ivw(mr.obj)
    T1 = proc.time()[3]
    est[repind,3+7] = res$Estimate
    est[repind,3+17] = res$StdError
    est[repind,3+27] = T1-T0
    rm(res)

    # 8. median
    T0 = proc.time()[3]
    res = mr_median(mr.obj)
    T1 = proc.time()[3]
    est[repind,3+8] = res$Estimate
    est[repind,3+18] = res$StdError
    est[repind,3+28] = T1-T0
    rm(res)

    # 9. Egger
    T0 = proc.time()[3]
    res = mr_egger(mr.obj)
    T1 = proc.time()[3]
    est[repind,3+9] = res$Estimate
    est[repind,3+19] = res$StdError.Est
    est[repind,3+29] = T1-T0
    rm(res)

    # 10. IB-PRESSO
    T0 = proc.time()[3]
    res <- tryCatch({
      ibpresso <- mr_presso_ib(
        BetaOutcome = "BetaOutcome",
        BetaExposure = "BetaExposure",
        BetaAux     = "BetaAux",
        SdOutcome   = "SdOutcome",
        SdExposure  = "SdExposure",
        SdAux       = "SdAux",
        data = data.frame(
          BetaOutcome = betahat_y.flt,
          BetaExposure = betahat_x.flt,
          BetaAux     = betahat_y_alt.flt,
          SdOutcome = rep(1/sqrt(ny), length(betahat_y.flt)),
          SdExposure = rep(1/sqrt(nx), length(betahat_x.flt)),
          SdAux     = rep(1/sqrt(ny_alt), length(betahat_y_alt.flt))
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
    }, error=function(e) list(b=NA,se=NA))

    T1 = proc.time()[3]
    est[repind,3+10] = res$b
    est[repind,3+20] = res$se
    est[repind,3+30] = T1-T0
    rm(res)

  }
  if (repind%%5==0){
    print(paste("Rep",repind,"numIV",numIV))
    save(est, file = paste0("../../results/simulation_results/DN_est_theta",
                            theta,
                            "_thetaU",
                            thetaU,
                            "_N",
                            format(N, scientific = FALSE),
                            "_prop_invalid",
                            prop_invalid,
                            "_overlap",
                            overlap,
                            ".rda"))
  }}
save(est, file = paste0("../../results/simulation_results/DN_est_theta",
                        theta,
                        "_thetaU",
                        thetaU,
                        "_N",
                        format(N, scientific = FALSE),
                        "_prop_invalid",
                        prop_invalid,
                        "_overlap",
                        overlap,
                        ".rda"))
#warnings()
