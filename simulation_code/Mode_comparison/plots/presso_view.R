

analyse_new_scenarios_eff <- function(power_dir = 1,
                                      thetaUvec = c(0.3, 0.5),
                                      Nvec = c(5e4,8e4,1e5,1.5e5,2e5),
                                      prop_invalid_vec = c(0.1, 0.3, 0.5),
                                      overlap_vec = c(0,0.25,0.5,0.75,1),
                                      scenario_vec = c("BI", "DI"),
                                      color_var = "overlap", plot_type = 1) {
  library(ggplot2)
  library(dplyr)
  library(splines)
  
  opt_phi_values <- data.frame()
  power_results <- data.frame()
  
  for (scenario in scenario_vec) {
    for (thetaU in thetaUvec) {
      for (N in Nvec) {
        for (prop_invalid in prop_invalid_vec) {
          for (overlap in overlap_vec) {
            
            file_path <- sprintf(
              "./results/presso_comp/sim_MRpressos_BI_beta_p_distort_theta%s_U%s_N%s_invalid%s_overlap%s.rda",
              theta, thetaU, N, prop_invalid, overlap
            )
            
            error_data <- data.frame()
            
            for (p in phi_range) {
              load(sprintf(file_path, p))
              if (!exists("est") || nrow(est) == 0) {
                next
              }
              est <- na.omit(est)
              methods <- c("MRMode", "mode_new")
              if (!all(c(methods, paste0(methods, "_se")) %in% colnames(est))) {
                next
              }
              type1_error <- colMeans(abs(est[, methods]) > 1.96 * est[, paste0(methods, "_se")]) * 100
              error_data <- rbind(error_data, data.frame(Method = methods, Phi = p, Type1Error = type1_error))
            }
            
            # Smooth and interpolate Type 1 Error curve
            opt_phi <- data.frame()
            for (method in unique(error_data$Method)) {
              df <- error_data %>% filter(Method == method)
              spline_fit <- smooth.spline(log(df$Phi), df$Type1Error, spar = 0.7)
              interp_phi <- exp(seq(log(0.1), log(50), length.out = 100))
              interp_error <- predict(spline_fit, log(interp_phi))$y
              interp_error[interp_error < 0] <- 0  # Ensuring non-negativity
              
              best_phi <- ifelse(any(interp_error <= 5), max(interp_phi[interp_error <= 5]), 50)
              opt_phi <- rbind(opt_phi, data.frame(Method = method, OptPhi = best_phi))
            }
            opt_phi_values <- rbind(opt_phi_values, cbind(opt_phi, scenario, thetaU, N, prop_invalid, overlap))
            
            # Compute Power at optimal phi using interpolation
            power_data <- data.frame()
            power_values <- data.frame()
            
            file_path_power <- if (power_dir == 1) {
              sprintf(
                "../../../results/simulation_results/mode_comp/%s_est_theta0.2_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
                scenario, thetaU, format(N, scientific = TRUE), prop_invalid, "%s", overlap)
            } else {
              sprintf(
                "../../../results/simulation_results/mode_comp/%s_est_theta-0.2_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
                scenario, thetaU, format(N, scientific = TRUE), prop_invalid, "%s", overlap)
            }
            
            for (p in phi_range) {
              load(sprintf(file_path_power, p))
              if (!exists("est") || nrow(est) == 0) {
                next
              }
              est <- na.omit(est)
              methods <- c("MRMode", "mode_new")
              if (!all(c(methods, paste0(methods, "_se")) %in% colnames(est))) {
                next
              }
              type2_error <- colMeans(abs(est[, methods]) <= 1.96 * est[, paste0(methods, "_se")]) * 100
              power_values <- rbind(power_values, data.frame(Method = methods, Phi = p, Power = pmin(100, 100 - type2_error)))
            }
            
            for (method in unique(opt_phi$Method)) {
              df <- power_values %>% filter(Method == method)
              spline_fit <- smooth.spline(log(df$Phi), df$Power, spar = 0.7)
              interp_phi <- exp(seq(log(0.1), log(50), length.out = 100))
              interp_power <- predict(spline_fit, log(interp_phi))$y
              interp_power <- pmin(100, pmax(0, interp_power))  # Ensuring power is within [0,100]
              
              opt_power <- interp_power[which.min(abs(interp_phi - opt_phi$OptPhi[opt_phi$Method == method]))]
              power_data <- rbind(power_data, data.frame(Method = method, Power = opt_power))
            }
            
            power_results <- rbind(power_results, cbind(power_data, scenario, thetaU, N, prop_invalid, overlap))
          }
        }
      }
    }
  }
  
  library(RColorBrewer)
  
  
  # Reshape data for plotting
  power_results_wide <- power_results %>%
    select(Method, Power, scenario, thetaU, N, prop_invalid, overlap) %>%
    spread(key = Method, value = Power)
  
  # Get the number of unique values for N to ensure sufficient shapes
  unique_N <- length(unique(power_results_wide$N))
  
  # Generate plot comparing optimal power values with differentiation by prop_invalid, N, and faceting by overlap
  if(plot_type == 1){
    # Plot comparing optimal power values with differentiation by prop_invalid and N
    power_plot <- ggplot(power_results_wide, aes(x = MRMode, y = mode_new, 
                                                 color = prop_invalid,  # treat overlap as continuous
                                                 shape = factor(N))) +
      geom_point(size = 3) +  # Points to represent the data
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = "Optimal Power of MR-MODE", 
           y = "Optimal Power of IB-MODE", 
           title = "Optimal Power Comparison") +
      scale_color_viridis_c(option = "turbo", name = "Invalid %") +
      scale_shape_manual(name = "Sample Size (N)", 
                         values = seq(1, unique_N)) +  # Dynamic shape scale
      theme_minimal() +
      theme(legend.position = "bottom",   # Improved legend positioning
            legend.title = element_text(size = 10),  # Improve legend title size
            legend.text = element_text(size = 8),  # Improve legend text size
            plot.title = element_text(hjust = 0.5))  # Center title
  } else {
    power_results_wide$overlap_label <- factor(
      paste0("overlap = ", power_results_wide$overlap)
    )
    
    power_plot <- ggplot(power_results_wide, aes(x = MRMode, y = mode_new, 
                                                 color = prop_invalid,
                                                 shape = factor(N))) +
      geom_point(size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = "Optimal Power of MR-MODE", 
           y = "Optimal Power of IB-MODE", 
           title = "Optimal Power Comparison") +
      scale_color_viridis_c(option = "turbo", name = "Invalid %") +
      scale_shape_manual(name = "Sample Size (N)", 
                         values = seq(1, length(unique(power_results_wide$N)))) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ overlap_label, scales = "free")
  }
  
  # Construct base path
  base_path <- "./plots"
  filename_comp <- sprintf("%s/overlap_power_comparison.eps",
                           base_path)
  ggsave(filename_comp, plot = power_plot, device = "eps", width = 10, height = 6) }

###################################################################################
###################################################################################
###################################################################################
###################################################################################


compare_and_plot_pvals <- function(
    presso_dir,
    theta_vec      = c(0, 0.3, 0.5),
    thetaU_vec     = c(0.2, 0.5),
    N_vec          = c(5e4, 1e5, 2e5),
    prop_invalid_vec = c(0.1, 0.3, 0.5),
    overlap_vec    = c(0, 0.25, 0.5, 0.75, 1)
) {
  library(dplyr)
  library(ggplot2)
  
  all_res <- expand.grid(
    theta        = theta_vec,
    thetaU       = thetaU_vec,
    N            = N_vec,
    prop_invalid = prop_invalid_vec,
    overlap      = overlap_vec,
    stringsAsFactors = FALSE
  ) %>% rowwise() %>% do({
    TH   <- .$theta
    TU   <- .$thetaU
    Nn   <- .$N
    PI   <- .$prop_invalid
    Ov   <- .$overlap
    
    # adjust this pattern to exactly match your filenames:
    pat <- sprintf(
      "sim_MRpressos_BI_beta_p_distort_theta%s_U%s_N%s_invalid%s_overlap%s\\.rda$",
      gsub("\\.", "_", as.character(TH)),
      TU, Nn, PI, Ov
    )
    files <- list.files(presso_dir, pat, full.names = TRUE)
    
    # For each file, extract the columns we need:
    df_list <- lapply(files, function(f) {
      load(f)   # loads `est`
      if (!exists("est") || nrow(est)==0) return(NULL)
      e <- na.omit(est)
      
      if (TH == 0) {
        req0 <- c("numIV","varX_expl","varY_expl",
                  "uni_beta","uni_p","uni_time")
        if (!all(req0 %in% names(e))) return(NULL)
        out <- e[1, req0]
        out$ib_beta <- NA; out$ib_p <- NA; out$ib_time <- NA
        
      } else {
        req1 <- c("numIV","varX_expl","varY_expl",
                  "uni_beta","uni_p","uni_time",
                  "ib_beta","ib_p","ib_time")
        if (!all(req1 %in% names(e))) return(NULL)
        out <- e[1, req1]
      }
      # add metadata
      out$theta        <- TH
      out$thetaU       <- TU
      out$N            <- Nn
      out$prop_invalid <- PI
      out$overlap      <- Ov
      out
    })
    bind_rows(df_list)
  }) %>% ungroup()
  
  # ensure numeric
  all_res <- all_res %>%
    mutate(across(c(theta, thetaU, N, prop_invalid, overlap,
                    numIV, varX_expl, varY_expl,
                    uni_beta, uni_p, uni_time,
                    ib_beta, ib_p, ib_time),
                  as.numeric))
  
  # --- Scatter plot for theta ≠ 0 ---
  df_nz <- filter(all_res, theta != 0, !is.na(ib_p))
  df_nz <- df_nz %>%
    mutate(
      facet_label = paste0("θ=", theta, ", overlap=", overlap)
    )
  
  p_scatter <- ggplot(df_nz,
                      aes(x = uni_p, y = ib_p,
                          color = prop_invalid,
                          shape = factor(N))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ facet_label, scales = "free") +
    labs(
      x = expression(uni[p]),
      y = expression(ib[p]),
      color = "Invalid %",
      shape = "Sample size\n(N)",
      title = "uni_p vs ib_p across all scenarios"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8)
    )
  
  # --- Histogram for theta = 0 (optional) ---
  df0 <- filter(all_res, theta == 0)
  p_hist0 <- ggplot(df0, aes(x = uni_p)) +
    geom_histogram(bins = 30) +
    facet_grid(prop_invalid ~ overlap,
               labeller = label_both) +
    labs(
      x = expression(uni[p]~"(θ=0)"),
      y = "Count",
      title = "Distribution of uni_p (θ=0) by invalid% and overlap"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 7)
    )
  
  return(list(
    data        = all_res,
    scatter     = p_scatter,
    hist_theta0 = p_hist0
  ))
}



res <- compare_and_plot_pvals(
  presso_dir      = "./results/presso_comp",
  theta_vec       = c(0.2,-0.2),
  thetaU_vec      = c(0.3, 0.5),
  N_vec           = c(5e4,8e4,1e5,1.5e5),
  prop_invalid_vec = c(0.1, 0.3, 0.5),
  overlap_vec     = c(0,0.25,0.5,0.75,1)
)

# your big scatter of uni_p vs ib_p:
print(res$scatter)

# optional: histogram for θ = 0
print(res$hist_theta0)

# combined table
head(res$data)
