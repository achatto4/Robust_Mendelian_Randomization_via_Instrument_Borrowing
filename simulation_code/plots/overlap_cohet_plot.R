# # Summary level simulations with balanced pleiotropy and InSIDE assumption satisfied
# rm(list = ls())
# library(data.table)
# library(dplyr)
# library(MASS)
# library(MendelianRandomization)
# library(MRMix)
# library(penalized)
# library(ks)
# library(IBMR)  # coheterogeneity_Q is available as IBMR::coheterogeneity_Q
# 
# # Parameters
# N = 1e5
# theta = 50 / sqrt(N)
# thetaU = thetaUx = 0.3
# theta_alt = 0.3
# thetaU_alt = 0.3
# prop_invalid = 0.5
# pthr = 5e-8
# M = 2e5
# 
# pi1 = 0.02 * (1 - prop_invalid)
# pi2 = 0.02 * prop_invalid
# pi3 = 0.01
# 
# sigma2x = sigma2y = 5e-5
# sigma2u = 1e-4
# sigma2x_td = sigma2y_td = sigma2x - thetaU * thetaUx * sigma2u
# 
# nx = N
# ny = N / 2
# ny_alt = N
# 
# # Grid
# overlap_adj_seq = seq(0, 1, by = 0.1)
# n_rep = 100
# Q_corr_mat = matrix(NA, nrow = n_rep, ncol = length(overlap_adj_seq))
# sig_count_vec = numeric(length(overlap_adj_seq))
# 
# set.seed(123)
# 
# # Fixed SNPs
# ind1_master = sample(M, round(M * pi1))
# ind2_master = sample(setdiff(1:M, ind1_master), round(M * pi2))
# ind3_master = sample(setdiff(1:M, c(ind1_master, ind2_master)), round(M * pi3))
# 
# for (i in seq_along(overlap_adj_seq)) {
#   overlap_adj = overlap_adj_seq[i]
#   
#   for (rep in 1:n_rep) {
#     set.seed(10000 + rep)
#     
#     shared_ind1 = sample(ind1_master, round(length(ind1_master) * overlap_adj))
#     ind1_1 = c(shared_ind1, sample(setdiff(1:M, shared_ind1), round(M * pi1) - length(shared_ind1)))
#     
#     shared_ind2 = sample(ind2_master, round(length(ind2_master) * overlap_adj))
#     ind2_1 = c(shared_ind2, sample(setdiff(1:M, c(ind1_1, shared_ind2)), round(M * pi2) - length(shared_ind2)))
#     
#     shared_ind3 = sample(ind3_master, round(length(ind3_master) * overlap_adj))
#     ind3_1 = c(shared_ind3, sample(setdiff(1:M, c(ind1_1, ind2_1, shared_ind3)), round(M * pi3) - length(shared_ind3)))
#     
#     gamma = phi = phi_alt = alpha = alpha_alt = rep(0, M)
#     
#     gamma[union(ind1_master, ind1_1)] = rnorm(length(union(ind1_master, ind1_1)), 0, sqrt(sigma2x))
#     gamma[union(ind2_master, ind2_1)] = rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2x_td))
#     
#     alpha[ind2_master] = rnorm(length(ind2_master), 0.005, sqrt(sigma2y_td))
#     alpha[ind3_master] = rnorm(length(ind3_master), 0.005, sqrt(sigma2y))
#     alpha_alt[ind2_1] = rnorm(length(ind2_1), 0.003, sqrt(sigma2y_td))
#     alpha_alt[ind3_1] = rnorm(length(ind3_1), 0.003, sqrt(sigma2y))
#     
#     phi[union(ind2_master, ind2_1)] = rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2u))
#     phi[-ind2_master] = 0
#     phi_alt[ind2_1] = phi[ind2_1]
#     
#     betax = gamma + thetaUx * phi
#     betay = alpha + theta * betax + thetaU * phi
#     betay_alt = alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
#     
#     betahat_x = betax + rnorm(M, 0, sqrt(1 / nx))
#     betahat_y = betay + rnorm(M, 0, sqrt(1 / ny))
#     betahat_y_alt = betay_alt + rnorm(M, 0, sqrt(1 / ny_alt))
#     
#     ind_filter = which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
#     if (length(ind_filter) < 3) next
#     
#     CoHetQ = coheterogeneity_Q(
#       BetaXG = betahat_x[ind_filter],
#       BetaYG_matrix = cbind(betahat_y[ind_filter], betahat_y_alt[ind_filter]),
#       seBetaXG = rep(1 / sqrt(nx), length(ind_filter)),
#       seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(ind_filter)),
#                               rep(1 / sqrt(ny_alt), length(ind_filter)))
#     )
#     
#     r_val = abs(CoHetQ$Q_corr_matrix[1, 2])
#     Q_corr_mat[rep, i] = r_val
#     
#     # t-statistic for correlation significance
#     df = length(ind_filter) - 2
#     if (df > 1 && r_val < 1) {
#       t_stat = r_val * sqrt(df) / sqrt(1 - r_val^2)
#       p_val = 2 * pt(-abs(t_stat), df)
#       if (!is.na(p_val) && p_val < 0.05) sig_count_vec[i] = sig_count_vec[i] + 1
#     }
#   }
#   cat("Done overlap_adj =", overlap_adj, "\n")
# }
# 
# mean_Q = colMeans(Q_corr_mat, na.rm = TRUE)
# sd_Q = apply(Q_corr_mat, 2, sd, na.rm = TRUE)
# sig_frac = sig_count_vec / n_rep  # Between 0 and 1
# 
# # Color palette
# cols = colorRampPalette(c("white", "blue"))(100)
# col_idx = as.integer(1 + (sig_frac * 99))
# point_colors = cols[col_idx]
# 
# # Plot
# plot(overlap_adj_seq, mean_Q, type = "n", ylim = c(0, 1),
#      xlab = "Adjusted Overlap", ylab = "Mean CoHeterogeneity Correlation",
#      main = "CoHeterogeneity vs Overlap\n(Colored by % Significant)")
# 
# grid()
# 
# polygon(c(overlap_adj_seq, rev(overlap_adj_seq)), 
#         c(mean_Q + sd_Q, rev(mean_Q - sd_Q)), 
#         col = adjustcolor("firebrick", alpha.f = 0.2), border = NA)
# 
# lines(overlap_adj_seq, mean_Q, lwd = 2, col = "firebrick")
# 
# points(overlap_adj_seq, mean_Q, col = point_colors, pch = 19)
# 
# legend("topleft", legend = c("0%", "50%", "100%"),
#        fill = colorRampPalette(c("white", "blue"))(3), 
#        title = "% reps with p < 0.05", border = NA)


#---------------------------------------------
# Extended Summary‐Level Simulations for CoHeterogeneity
#---------------------------------------------

# 1. Setup -------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(penalized)
library(ks)
library(ggplot2)

# Source your CoHeterogeneity function
library(IBMR)  # coheterogeneity_Q is available as IBMR::coheterogeneity_Q

# 2. Fixed parameters --------------------------------------------------------
N            <- 1e5
theta        <- 50 / sqrt(N)
thetaU       <- thetaUx <- 0.3
theta_alt    <- 0.3
thetaU_alt   <- 0.3

# SNP proportions (will split valid vs invalid dynamically)
pi_total1    <- 0.02
pi_total2    <- 0.02
pi3          <- 0.01

sigma2x      <- sigma2y      <- 5e-5
sigma2u      <- 1e-4
sigma2x_td   <- sigma2y_td   <- sigma2x - thetaU * thetaUx * sigma2u

nx           <- N
ny           <- N / 2
ny_alt       <- N

# Simulation grid
overlap_seq      <- seq(0, 1, by = 0.1)
prop_invalid_seq <- c(0.25, 0.5, 0.75)    # example settings
n_rep            <- 100
pthr             <- 5e-8
M                <- 2e5

set.seed(123)

# 3. Prepare result container ------------------------------------------------
results <- CJ(prop_invalid = prop_invalid_seq,
              overlap     = overlap_seq,
              rep         = 1:n_rep)[
                , .(Q_corr = NA_real_, sig      = NA_integer_), 
                by = .(prop_invalid, overlap, rep)
              ]

# 4. Main simulation loop ----------------------------------------------------
for (pi_inv in prop_invalid_seq) {
  
  # derive valid vs invalid proportions
  pi1 <- pi_total1 * (1 - pi_inv)
  pi2 <- pi_total2 * pi_inv
  
  # pick master SNP sets once per invalid proportion
  ind1_master <- sample(M, round(M * pi1))
  ind2_master <- sample(setdiff(1:M, ind1_master), round(M * pi2))
  ind3_master <- sample(setdiff(1:M, c(ind1_master, ind2_master)), round(M * pi3))
  
  for (ov in overlap_seq) {
    for (rep in seq_len(n_rep)) {
      set.seed(10000 + rep)
      
      # create overlapping subsets
      shared1   <- sample(ind1_master, round(length(ind1_master) * ov))
      ind1_1    <- c(shared1,
                     sample(setdiff(1:M, shared1), round(M*pi1) - length(shared1)))
      
      shared2   <- sample(ind2_master, round(length(ind2_master) * ov))
      ind2_1    <- c(shared2,
                     sample(setdiff(1:M, c(ind1_1, shared2)), round(M*pi2) - length(shared2)))
      
      shared3   <- sample(ind3_master, round(length(ind3_master) * ov))
      ind3_1    <- c(shared3,
                     sample(setdiff(1:M, c(ind1_1, ind2_1, shared3)), round(M*pi3) - length(shared3)))
      
      # initialize effect vectors
      gamma        <- phi <- phi_alt <- alpha <- alpha_alt <- rep(0, M)
      
      # genetic effects on exposure
      gamma[union(ind1_master, ind1_1)]            <- 
        rnorm(length(union(ind1_master, ind1_1)), 0, sqrt(sigma2x))
      gamma[union(ind2_master, ind2_1)]            <- 
        rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2x_td))
      
      # pleiotropic effects on outcome
      alpha[ind2_master]                           <- 
        rnorm(length(ind2_master), 0.005, sqrt(sigma2y_td))
      alpha[ind3_master]                           <- 
        rnorm(length(ind3_master), 0.005, sqrt(sigma2y))
      alpha_alt[ind2_1]                            <- 
        rnorm(length(ind2_1), 0.003, sqrt(sigma2y_td))
      alpha_alt[ind3_1]                            <- 
        rnorm(length(ind3_1), 0.003, sqrt(sigma2y))
      
      # hidden confounder effects
      phi[union(ind2_master, ind2_1)]             <- 
        rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2u))
      phi <- 0
      phi_alt[ind2_1]                             <- phi[ind2_1]
      
      # true betaX & betaY
      betax     <- gamma + thetaUx * phi
      betay     <- alpha + theta * betax + thetaU * phi
      betay_alt <- alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
      
      # observed estimates w/ sampling error
      betahat_x       <- betax     + rnorm(M, 0, sqrt(1 / nx))
      betahat_y       <- betay     + rnorm(M, 0, sqrt(1 / ny))
      betahat_y_alt   <- betay_alt + rnorm(M, 0, sqrt(1 / ny_alt))
      
      # instrument selection
      sel_idx <- which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
      if (length(sel_idx) < 3) next
      
      # compute CoHeterogeneity
      CoHetQ <- coheterogeneity_Q(
        BetaXG          = betahat_x[sel_idx],
        BetaYG_matrix   = cbind(betahat_y[sel_idx], betahat_y_alt[sel_idx]),
        seBetaXG        = rep(1 / sqrt(nx), length(sel_idx)),
        seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(sel_idx)),
                                rep(1 / sqrt(ny_alt), length(sel_idx)))
      )
      
      r_val <- abs(CoHetQ$Q_corr_matrix[1,2])
      df    <- length(sel_idx) - 2
      p_val <- if (df > 1 && r_val < 1) {
        t_stat <- r_val * sqrt(df) / sqrt(1 - r_val^2)
        2 * pt(-abs(t_stat), df)
      } else NA_real_
      
      # store
      idx_row <- which(results$prop_invalid == pi_inv &
                         results$overlap     == ov      &
                         results$rep         == rep)
      results[idx_row, `:=`(Q_corr = r_val,
                            sig    = as.integer(!is.na(p_val) && p_val < 0.05))]
    }
  }
  message("Finished prop_invalid = ", pi_inv)
}

# 5. Summarize ---------------------------------------------------------------
summary_dt <- results[, .(
  mean_Q   = mean(Q_corr, na.rm = TRUE),
  sd_Q     =   sd(Q_corr, na.rm = TRUE),
  sig_frac =   mean(sig, na.rm = TRUE)
), by = .(prop_invalid, overlap)]

# 6. Plot --------------------------------------------------------------------
ggplot(summary_dt, aes(x = overlap, y = mean_Q, group = 1)) +
  facet_wrap(~ prop_invalid,
             labeller = labeller(prop_invalid = function(x) paste0("Prop. invalid = ", x))) +
  geom_ribbon(aes(ymin = mean_Q - sd_Q, ymax = mean_Q + sd_Q),
              fill = "firebrick", alpha = 0.2) +
  geom_line(colour = "firebrick", size = 1) +
  geom_point(aes(colour = sig_frac, size = sig_frac)) +
  scale_colour_gradient(low = "white", high = "blue",
                        name = "Sig. prop.") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(limits = c(0,0.5)) +
  labs(x    = "Adjusted Overlap",
       y    = "Mean |CoHeterogeneity correlation|",
       title= "CoHeterogeneity captures shared‐instrument overlap",
       subtitle = "Ribbon = ±1 SD; point colour ∝ fraction of reps p<0.05") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")


###################################################
###################################################
###################################################

# 1. Setup -------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(penalized)
library(ks)
library(ggplot2)

# Source your CoHeterogeneity function
library(IBMR)  # coheterogeneity_Q is available as IBMR::coheterogeneity_Q

# 2. Fixed parameters --------------------------------------------------------
#N_seq         <- c(7.5e4, 1e5, 1.5e5)     # Different sample sizes
N_seq         <- c(1e5)     # Different sample sizes
thetaU        <- thetaUx <- 0.3
theta_alt     <- 0.3
thetaU_alt    <- 0.3

# SNP proportions
pi_total1     <- 0.02
pi_total2     <- 0.02
pi3           <- 0.01

sigma2x       <- sigma2y      <- 5e-5
sigma2u       <- 1e-4
sigma2x_td    <- sigma2y_td   <- sigma2x - thetaU * thetaUx * sigma2u

# Simulation grid
#prop_invalid_seq <- c(0.25, 0.5, 0.75)
prop_invalid_seq <- c(0.5)
overlap_seq      <- seq(0, 1, by = 0.1)
n_rep            <- 100
pthr             <- 5e-8
M                <- 2e5

set.seed(123)

# 3. Prepare result container ------------------------------------------------
results <- CJ(N            = N_seq,
              prop_invalid = prop_invalid_seq,
              overlap      = overlap_seq,
              rep          = 1:n_rep)[
                , .(Q_corr = NA_real_, sig = NA_integer_), 
                by = .(N, prop_invalid, overlap, rep)
              ]

# 4. Main simulation loop ----------------------------------------------------
for (N in N_seq) {
  
  theta   <- 50 / sqrt(N)
  nx      <- N
  ny      <- N / 2
  ny_alt  <- N
  
  for (pi_inv in prop_invalid_seq) {
    
    # derive valid vs invalid proportions
    pi1 <- pi_total1 * (1 - pi_inv)
    pi2 <- pi_total2 * pi_inv
    
    # pick master SNP sets
    ind1_master <- sample(M, round(M * pi1))
    ind2_master <- sample(setdiff(1:M, ind1_master), round(M * pi2))
    ind3_master <- sample(setdiff(1:M, c(ind1_master, ind2_master)), round(M * pi3))
    
    for (ov in overlap_seq) {
      for (rep in seq_len(n_rep)) {
        set.seed(10000 + rep)
        
        # create overlapping subsets
        shared1   <- sample(ind1_master, round(length(ind1_master) * ov))
        ind1_1    <- c(shared1,
                       sample(setdiff(1:M, shared1), round(M*pi1) - length(shared1)))
        
        shared2   <- sample(ind2_master, round(length(ind2_master) * ov))
        ind2_1    <- c(shared2,
                       sample(setdiff(1:M, c(ind1_1, shared2)), round(M*pi2) - length(shared2)))
        
        shared3   <- sample(ind3_master, round(length(ind3_master) * ov))
        ind3_1    <- c(shared3,
                       sample(setdiff(1:M, c(ind1_1, ind2_1, shared3)), round(M*pi3) - length(shared3)))
        
        # initialize effect vectors
        gamma        <- phi <- phi_alt <- alpha <- alpha_alt <- rep(0, M)
        
        # genetic effects on exposure
        gamma[union(ind1_master, ind1_1)] <- rnorm(length(union(ind1_master, ind1_1)), 0, sqrt(sigma2x))
        gamma[union(ind2_master, ind2_1)] <- rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2x_td))
        
        # pleiotropic effects
        alpha[ind2_master]    <- rnorm(length(ind2_master), 0.008, sqrt(sigma2y_td))
        alpha[ind3_master]    <- rnorm(length(ind3_master), 0.008, sqrt(sigma2y))
        alpha_alt[ind2_1]     <- rnorm(length(ind2_1), 0.003, sqrt(sigma2y_td))
        alpha_alt[ind3_1]     <- rnorm(length(ind3_1), 0.003, sqrt(sigma2y))
        
        # confounder
        phi[union(ind2_master, ind2_1)] <- rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2u))
        phi <- 0
        phi_alt[ind2_1]                 <- phi[ind2_1]
        
        # true effects
        betax     <- gamma + thetaUx * phi
        betay     <- alpha + theta * betax + thetaU * phi
        betay_alt <- alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
        
        # observed estimates
        betahat_x     <- betax + rnorm(M, 0, sqrt(1 / nx))
        betahat_y     <- betay + rnorm(M, 0, sqrt(1 / ny))
        betahat_y_alt <- betay_alt + rnorm(M, 0, sqrt(1 / ny_alt))
        
        # instrument selection
        sel_idx <- which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
        if (length(sel_idx) < 3) next
        
        # CoHeterogeneity calculation
        CoHetQ <- coheterogeneity_Q(
          BetaXG          = betahat_x[sel_idx],
          BetaYG_matrix   = cbind(betahat_y[sel_idx], betahat_y_alt[sel_idx]),
          seBetaXG        = rep(1 / sqrt(nx), length(sel_idx)),
          seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(sel_idx)),
                                  rep(1 / sqrt(ny_alt), length(sel_idx)))
        )
        
        r_val <- abs(CoHetQ$Q_corr_matrix[1, 2])
        df    <- length(sel_idx) - 2
        p_val <- if (df > 1 && r_val < 1) {
          t_stat <- r_val * sqrt(df) / sqrt(1 - r_val^2)
          2 * pt(-abs(t_stat), df)
        } else NA_real_
        
        # store result
        idx_row <- which(results$N == N &
                           results$prop_invalid == pi_inv &
                           results$overlap == ov &
                           results$rep == rep)
        results[idx_row, `:=`(Q_corr = r_val,
                              sig    = as.integer(!is.na(p_val) && p_val < 0.05))]
      }
    }
    message("Finished N = ", N, " | prop_invalid = ", pi_inv)
  }
}

# 5. Summarize ---------------------------------------------------------------
summary_dt <- results[, .(
  mean_Q   = mean(Q_corr, na.rm = TRUE),
  sd_Q     = sd(Q_corr, na.rm = TRUE),
  sig_frac = mean(sig, na.rm = TRUE)
), by = .(N, prop_invalid, overlap)]

# Convert N and prop_invalid to factors with labels
summary_dt[, N_factor := factor(N, levels = N_seq,
                                labels = paste0("N = ", format(N_seq, scientific = FALSE)))]

summary_dt[, prop_invalid_factor := factor(prop_invalid,
                                           levels = prop_invalid_seq,
                                           labels = paste0("Prop. invalid = ", prop_invalid_seq))]
library(grid)  # for unit()

p <- ggplot(summary_dt, aes(x = overlap, y = mean_Q, group = 1)) +
  facet_grid(N_factor ~ prop_invalid_factor) +
  
  # Ribbon and line
  geom_ribbon(aes(ymin = mean_Q - sd_Q, ymax = mean_Q + sd_Q),
              fill = "firebrick", alpha = 0.4) +
  geom_line(colour = "firebrick", size = 1.2) +
  geom_point(aes(colour = sig_frac, size = sig_frac)) +
  
  # Scales
  scale_colour_gradient(low = "white", high = "blue", name = "Proportion Signif.") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 0.3)) +
  
  #, title = "CoHeterogeneity captures shared-instrument overlap across N and invalidity"
  # Labels
  labs(x = expression(D[ov]),
       y = "CoHeterogeneity correlation") +
  
  # Custom theme
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#F0F0F0", colour = NA),  # light gray
    plot.background  = element_rect(fill = "#FFFFFF", colour = NA),  # white outer background
    strip.background = element_rect(fill = "#D9D9D9", colour = NA),  # facet label background
    strip.text       = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13),
    legend.key.size  = unit(2, "cm")
  )
p
# Save as EPS
ggsave("CoHeterogeneity_plot_DI.eps", plot = p, device = "eps",
       width = 10, height = 8, units = "in")

#############################################################################

# 1. Setup -------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(dplyr)
library(MASS)
library(MendelianRandomization)
library(MRMix)
library(penalized)
library(ks)
library(ggplot2)

# Source your CoHeterogeneity function
library(IBMR)  # coheterogeneity_Q is available as IBMR::coheterogeneity_Q

# 2. Fixed parameters --------------------------------------------------------
N_seq         <- c(1e5)     # Different sample sizes
thetaU        <- thetaUx <- 0.3
theta_alt     <- 0.3
thetaU_alt    <- 0.3

# SNP proportions
pi_total1     <- 0.02
pi_total2     <- 0.02
pi3           <- 0.01

sigma2x       <- sigma2y      <- 5e-5
sigma2u       <- 1e-4
sigma2x_td    <- sigma2y_td   <- sigma2x - thetaU * thetaUx * sigma2u

# Simulation grid
prop_invalid_seq <- c(0.5)
overlap_seq      <- seq(0, 1, by = 0.1)
n_rep            <- 100
pthr             <- 5e-8
M                <- 2e5

set.seed(123)

# 3. Prepare result container ------------------------------------------------
results <- CJ(N            = N_seq,
              prop_invalid = prop_invalid_seq,
              overlap      = overlap_seq,
              rep          = 1:n_rep)[
                , .(Q_corr = NA_real_, sig = NA_integer_), 
                by = .(N, prop_invalid, overlap, rep)
              ]

# 4. Main simulation loop ----------------------------------------------------
for (N in N_seq) {
  
  theta   <- 50 / sqrt(N)
  nx      <- N
  ny      <- N / 2
  ny_alt  <- N
  
  for (pi_inv in prop_invalid_seq) {
    
    # derive valid vs invalid proportions
    pi1 <- pi_total1 * (1 - pi_inv)
    pi2 <- pi_total2 * pi_inv
    
    # pick master SNP sets
    ind1_master <- sample(M, round(M * pi1))
    ind2_master <- sample(setdiff(1:M, ind1_master), round(M * pi2))
    ind3_master <- sample(setdiff(1:M, c(ind1_master, ind2_master)), round(M * pi3))
    
    for (ov in overlap_seq) {
      for (rep in seq_len(n_rep)) {
        set.seed(10000 + rep)
        
        # create overlapping subsets
        shared1   <- sample(ind1_master, round(length(ind1_master) * ov))
        ind1_1    <- c(shared1,
                       sample(setdiff(1:M, shared1), round(M*pi1) - length(shared1)))
        
        shared2   <- sample(ind2_master, round(length(ind2_master) * ov))
        ind2_1    <- c(shared2,
                       sample(setdiff(1:M, c(ind1_1, shared2)), round(M*pi2) - length(shared2)))
        
        shared3   <- sample(ind3_master, round(length(ind3_master) * ov))
        ind3_1    <- c(shared3,
                       sample(setdiff(1:M, c(ind1_1, ind2_1, shared3)), round(M*pi3) - length(shared3)))
        
        # initialize effect vectors
        gamma        <- phi <- phi_alt <- alpha <- alpha_alt <- rep(0, M)
        
        # genetic effects on exposure
        gamma[union(ind1_master, ind1_1)] <- rnorm(length(union(ind1_master, ind1_1)), 0, sqrt(sigma2x))
        gamma[union(ind2_master, ind2_1)] <- rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2x_td))
        
        # pleiotropic effects
        alpha[ind2_master]    <- rnorm(length(ind2_master), 0.008, sqrt(sigma2y_td))
        alpha[ind3_master]    <- rnorm(length(ind3_master), 0.008, sqrt(sigma2y))
        alpha_alt[ind2_1]     <- rnorm(length(ind2_1), 0.008, sqrt(sigma2y_td))
        alpha_alt[ind3_1]     <- rnorm(length(ind3_1), 0.008, sqrt(sigma2y))
        
        # confounder
        phi[union(ind2_master, ind2_1)] <- rnorm(length(union(ind2_master, ind2_1)), 0, sqrt(sigma2u))
        phi_alt[ind2_1]                 <- phi[ind2_1]
        
        # true effects
        betax     <- gamma + thetaUx * phi
        betay     <- alpha + theta * betax + thetaU * phi
        betay_alt <- alpha_alt + theta_alt * betax + thetaU_alt * phi_alt
        
        # observed estimates
        betahat_x     <- betax + rnorm(M, 0, sqrt(1 / nx))
        betahat_y     <- betay + rnorm(M, 0, sqrt(1 / ny))
        betahat_y_alt <- betay_alt + rnorm(M, 0, sqrt(1 / ny_alt))
        
        # instrument selection
        sel_idx <- which(2 * pnorm(-sqrt(nx) * abs(betahat_x)) < pthr)
        if (length(sel_idx) < 3) next
        
        # CoHeterogeneity calculation
        CoHetQ <- coheterogeneity_Q(
          BetaXG          = betahat_x[sel_idx],
          BetaYG_matrix   = cbind(betahat_y[sel_idx], betahat_y_alt[sel_idx]),
          seBetaXG        = rep(1 / sqrt(nx), length(sel_idx)),
          seBetaYG_matrix = cbind(rep(1 / sqrt(ny), length(sel_idx)),
                                  rep(1 / sqrt(ny_alt), length(sel_idx)))
        )
        
        r_val <- abs(CoHetQ$Q_corr_matrix[1, 2])
        df    <- length(sel_idx) - 2
        p_val <- if (df > 1 && r_val < 1) {
          t_stat <- r_val * sqrt(df) / sqrt(1 - r_val^2)
          2 * pt(-abs(t_stat), df)
        } else NA_real_
        
        # store result
        idx_row <- which(results$N == N &
                           results$prop_invalid == pi_inv &
                           results$overlap == ov &
                           results$rep == rep)
        results[idx_row, `:=`(Q_corr = r_val,
                              sig    = as.integer(!is.na(p_val) && p_val < 0.05))]
      }
    }
    message("Finished N = ", N, " | prop_invalid = ", pi_inv)
  }
}

# 5. Summarize ---------------------------------------------------------------
summary_dt <- results[, .(
  mean_Q   = mean(Q_corr, na.rm = TRUE),
  sd_Q     = sd(Q_corr, na.rm = TRUE)
), by = .(N, prop_invalid, overlap)]

# Convert N and prop_invalid to factors with labels
summary_dt[, N_factor := factor(N, levels = N_seq,
                                labels = paste0("N = ", format(N_seq, scientific = FALSE)))]

summary_dt[, prop_invalid_factor := factor(prop_invalid,
                                           levels = prop_invalid_seq,
                                           labels = paste0("Prop. invalid = ", prop_invalid_seq))]
library(grid)  # for unit()

# Points to highlight
special_overlaps <- seq(0.2, 1, by = 0.2)  # 0.2, 0.4, ..., 1.0
highlight_df <- summary_dt %>% filter(overlap %in% special_overlaps)

p <- ggplot(summary_dt, aes(x = overlap, y = mean_Q, group = 1)) +
  # Uncertainty ribbon, light blue with alpha
  geom_ribbon(aes(ymin = mean_Q - sd_Q, ymax = mean_Q + sd_Q),
              fill = "#E0ECF8", alpha = 0.75) +
  # Thick, bold blue line
  geom_line(colour = "#1565C0", size = 8, lineend = "round") +
  # Regular points (for context, smaller)
  geom_point(shape = 21, fill = "#1976D2", color = "#FFFFFF", size = 4, stroke = 1.5, alpha = 0.9) +
  # Highlighted points at overlaps 0.2, 0.4, ...
  geom_point(
    data = highlight_df,
    aes(x = overlap, y = mean_Q),
    shape = 21, fill = "#1565C0", color = "#FFFFFF", size = 12, stroke = 3, alpha = 1
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(D[ov]), y = expression("Coheterogeneity statistics ("*rho[H]*")")) +
  theme_minimal(base_size = 32) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF", colour = NA),
    plot.background  = element_rect(fill = "#FFFFFF", colour = NA),
    strip.background = element_rect(fill = "#FFFFFF", colour = NA),
    strip.text       = element_text(face = "bold", size = 34),
    axis.text        = element_text(size = 30, colour = "black"),
    axis.title       = element_text(size = 38, face = "bold"),
    panel.grid.major = element_line(colour = "grey85"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

ggsave("CoHeterogeneity_plot_grant.eps", plot = p, device = "eps",
       width = 12, height = 10, units = "in")
