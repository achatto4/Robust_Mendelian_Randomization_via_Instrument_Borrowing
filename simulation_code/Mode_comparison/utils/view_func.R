library(ggplot2)
library(tidyr)

# file_path <- sprintf("./results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                      scenario, est_theta, thetaU, format(N, scientific = sctfc), prop_invalid,  overlap)

#setwd(".")  # set working directory as needed

formatted_phi <- formatC(exp(seq(log(0.1), log(50), length.out = 30)), format = "g", digits = 14, drop0trailing = TRUE)

library(ggplot2)

analyze_all_scenarios_coverage <- function(est_theta,
                                           thetaU,
                                           N,
                                           phi = exp(seq(log(0.1), log(50), length.out = 30)),
                                           prop_invalid,
                                           overlap,
                                           scenario,
                                           sctfc = FALSE) {
  df_all <- data.frame()
  power_df <- data.frame()
  
  for (p in phi) {
    file_path <- sprintf(
      "../../../results/simulation_results/mode_comp/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      p,
      overlap
    )
    
    load(file_path)
    est <- na.omit(est)
    
    methods <- c("MRMode", "mode_new")
    within_95 <- colMeans((est[, methods] - 1.96 * est[, paste0(methods, "_se")] <= est_theta) &
                            (est[, methods] + 1.96 * est[, paste0(methods, "_se")] >= est_theta)) * 100
    
    power_df_scenario <- data.frame(Method = methods,
                                    Power = within_95,
                                    Scenario = scenario,
                                    Phi = p)
    power_df <- rbind(power_df, power_df_scenario)
  }
  
  plot_title <- sprintf("Coverage vs Phi\n est_theta=%g, thetaU=%g, N=%s, prop_invalid=%g, overlap=%g, Scenario=%s",
                        est_theta, thetaU, format(N, scientific = sctfc), prop_invalid, overlap, scenario)
  
  p_power <- ggplot(power_df, aes(x = Phi, y = Power, color = Method, group = Method)) +
    geom_line() +
    geom_point() +
    labs(title = plot_title, y = "Coverage (%)", x = "Phi") +
    theme_minimal()
  
  print(p_power)
}



analyze_all_scenarios_error <- function(est_theta,
                                        thetaU,
                                        N,
                                        phi = exp(seq(log(0.1), log(50), length.out = 30)),
                                        prop_invalid,
                                        overlap,
                                        scenario,
                                        sctfc = TRUE) {
  library(ggplot2)
  library(dplyr)
  
  error_df <- data.frame()
  
  for (p in phi) {
    file_path <- sprintf(
      "../../../results/simulation_results/mode_comp/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      p,
      overlap
    )
    
    load(file_path)
    est <- na.omit(est)
    
    methods <- c("MRMode", "mode_new")
    
    if (est_theta == 0) {
      error_rate <- colMeans(abs(est[, methods]) > 1.96 * est[, paste0(methods, "_se")]) * 100
      error_label <- "Type 1 Error"
    } else {
      type2_error_rate <- colMeans(abs(est[, methods]) <= 1.96 * est[, paste0(methods, "_se")]) * 100
      error_rate <- 100 - type2_error_rate
      error_label <- "Power (1 - Type 2 Error)"
    }
    
    error_df_scenario <- data.frame(Method = methods,
                                    ErrorRate = error_rate,
                                    Scenario = scenario,
                                    Phi = p)
    error_df <- rbind(error_df, error_df_scenario)
  }
  
  plot_title <- sprintf("%s vs Phi\n est_theta=%g, thetaU=%g, N=%s, prop_invalid=%g, overlap=%g, Scenario=%s",
                        error_label, est_theta, thetaU, format(N, scientific = sctfc), prop_invalid, overlap, scenario)
  
  p_error <- ggplot(error_df, aes(x = Phi, y = ErrorRate, color = Method, group = Method)) +
    geom_line() +
    geom_point() +
    labs(title = plot_title, y = paste(error_label, "Rate (%)"), x = "Phi") +
    theme_minimal()
  
  print(p_error)
}

analyze_all_scenarios_mse <- function(est_theta,
                                      thetaU,
                                      N,
                                      phi = exp(seq(log(0.1), log(50), length.out = 30)),
                                      prop_invalid,
                                      overlap,
                                      scenario,
                                      sctfc = FALSE) {
  mse_df <- data.frame()    
  
  for (p in phi) {
    file_path <- sprintf(
      "../../../results/simulation_results/mode_comp/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      p,
      overlap
    )
    
    load(file_path)
    est <- na.omit(est)
    methods <- c("MRMode", "mode_new")
    
    mse <- colMeans((est[, methods] - est_theta) ^ 2)
    
    mse_df_scenario <- data.frame(Method = methods,
                                  MSE = mse,
                                  Scenario = scenario,
                                  Phi = p)
    mse_df <- rbind(mse_df, mse_df_scenario)
  }
  
  plot_title <- sprintf("MSE vs Phi\n est_theta=%g, thetaU=%g, N=%s, prop_invalid=%g, overlap=%g, Scenario=%s",
                        est_theta, thetaU, format(N, scientific = sctfc), prop_invalid, overlap, scenario)
  
  p_mse <- ggplot(mse_df, aes(x = Phi, y = MSE, color = Method, group = Method)) +
    geom_line() +
    geom_point() +
    labs(title = plot_title, y = "Mean Squared Error", x = "Phi") +
    theme_minimal()
  
  print(p_mse)
}

analyse_new_scenarios_eff <- function(power_dir = 1,
                                      thetaUvec = c(0.3, 0.5),
                                      Nvec = c(5e4, 8e4, 1e5, 2e5),
                                      prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7),
                                      overlap_vec = c(0, 0.1, 0.3, 0.5, 0.75, 1),
                                      phi_range = exp(seq(log(0.1), log(10), length.out = 25)),
                                      scenario_vec = c("BI", "BN", "DI", "DN"),
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
              "../../../results/simulation_results/mode_comp/%s_est_theta0_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
              scenario, thetaU, format(N, scientific = FALSE), prop_invalid, "%s", overlap)
            
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
                scenario, thetaU, format(N, scientific = FALSE), prop_invalid, "%s", overlap)
            } else {
              sprintf(
                "../../../results/simulation_results/mode_comp/%s_est_theta-0.2_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
                scenario, thetaU, format(N, scientific = FALSE), prop_invalid, "%s", overlap)
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
  
  library(ggplot2)
  library(dplyr)
  
  power_results_wide_clean <- power_results_wide %>% 
    filter(abs(MRMode - mode_new) <= 60)
  
  point_size <- 6
  legend_point_size <- 6
  abline_width <- 2
  
  colorbar_breaks <- range(power_results_wide_clean$prop_invalid)
  
  if (plot_type == 1) {
    power_plot <- ggplot(
      power_results_wide_clean,
      aes(
        x = MRMode, 
        y = mode_new, 
        color = prop_invalid,
        shape = factor(N)
      )
    ) +
      geom_point(size = point_size, stroke = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = abline_width) +
      labs(
        x = "Power of MR-MODE", 
        y = "Power of IB-MODE", 
        title = "Optimal Power Comparison"
      ) +
      scale_color_viridis_c(
        option = "turbo",
        name = "% of invalid instruments",
        breaks = colorbar_breaks,         
        guide = guide_colorbar(
          title.position = "top",
          barwidth = 32,   # *** Here is the extra-long scale! ***
          barheight = 1.5,
          ticks = TRUE
        )
      ) +
      scale_shape_manual(
        name = "Sample Size (N)",
        values = seq(1, length(unique(power_results_wide_clean$N)))
      ) +
      theme_minimal(base_size = 21) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 20)
      ) +
      guides(
        color = guide_colorbar(order = 1, title.vjust = 1),
        shape = guide_legend(order = 2, override.aes = list(size = legend_point_size))
      )
  } else {
    power_results_wide_clean$overlap_label <- factor(
      paste0("overlap = ", power_results_wide_clean$overlap)
    )
    power_plot <- ggplot(
      power_results_wide_clean,
      aes(
        x = MRMode, 
        y = mode_new, 
        color = prop_invalid,
        shape = factor(N)
      )
    ) +
      geom_point(size = point_size, stroke = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = abline_width) +
      labs(
        x = "Optimal Power of MR-MODE", 
        y = "Optimal Power of IB-MODE", 
        title = "Optimal Power Comparison"
      ) +
      scale_color_viridis_c(
        option = "turbo",
        name = "Invalid %",
        breaks = colorbar_breaks,
        guide = guide_colorbar(
          title.position = "top",
          barwidth = unit(1.5, "cm"),  # *** make it very long ***
          barheight = 1.5,
          ticks = TRUE
        )
      ) +
      scale_shape_manual(
        name = "Sample Size (N)",
        values = seq(1, length(unique(power_results_wide_clean$N)))
      ) +
      theme_minimal(base_size = 21) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 20)
      ) +
      guides(
        color = guide_colorbar(order = 1, title.vjust = 1),
        shape = guide_legend(order = 2, override.aes = list(size = legend_point_size))
      ) +
      facet_wrap(~ overlap_label, scales = "free", nrow = 2, ncol = 3)
  }
  
  base_path <- "./plots"
  filename_comp <- sprintf("%s/overlap_power_comparison.eps", base_path)
  ggsave(filename_comp, plot = power_plot, device = "eps", width = 10, height = 8)
  # Generate plot comparing optimal power values with differentiation by prop_invalid, N, and faceting by overlap
  if(plot_type == 1){
    # Plot comparing optimal power values with differentiation by prop_invalid and N
    power_plot <- ggplot(power_results_wide, aes(x = MRMode, y = mode_new,
                                                 color = prop_invalid,  # treat overlap as continuous
                                                 shape = factor(N))) +
      geom_point(size = 3) +  # Points to represent the data
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = "Power of MR-MODE",
           y = "Power of IB-MODE",
           title = "Optimal Power Comparison") +
      scale_color_viridis_c(option = "turbo", name = "% of invalid instruments") +
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
  ggsave(filename_comp, plot = power_plot, device = "eps", width = 10, height = 6)

  
   }

analyse_new_scenarios_phi_hist <- function(thetaUvec = c(0.3, 0.5),
                                           Nvec = c(5e4, 8e4, 1e5, 2e5),
                                           prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7),
                                           overlap_vec = c(0, 0.1, 0.3, 0.5, 0.75, 1),
                                           phi_range = exp(seq(log(0.1), log(50), length.out =30)),
                                           scenario_vec = c("BI", "BN", "DI", "DN")) {
  library(ggplot2)
  library(dplyr)
  library(splines)
  
  opt_phi_values <- data.frame()
  
  for (scenario in scenario_vec) {
    for (thetaU in thetaUvec) {
      for (N in Nvec) {
        for (prop_invalid in prop_invalid_vec) {
          for (overlap in overlap_vec) {
            
            file_path <- sprintf(
              "../../../results/simulation_results/mode_comp/%s_est_theta0_thetaU%g_N%s_prop_invalid%g_phi%s_overlap%g.rda",
              scenario, thetaU, format(N, scientific = FALSE), prop_invalid, "%s", overlap)
            
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
              
              best_phi <- ifelse(any(interp_error <= 5), max(interp_phi[interp_error <= 5]), NA)
              opt_phi <- rbind(opt_phi, data.frame(Method = method, OptPhi = best_phi))
            }
            opt_phi_values <- rbind(opt_phi_values, cbind(opt_phi, scenario, thetaU, N, prop_invalid, overlap))
          }
        }
      }
    }
  }
  
  # Generate line-based histograms (density plots) for each method
  phi_hist_plot <- ggplot(opt_phi_values, aes(x = OptPhi)) +
    geom_density(data = subset(opt_phi_values, Method == "MRMode"), aes(color = "MRMode"), fill = "blue", alpha = 0.4) +
    geom_density(data = subset(opt_phi_values, Method == "mode_new"), aes(color = "mode_new"), fill = "red", alpha = 0.4) +
    labs(x = "Optimal Phi Values", y = "Density", title = "Density Plot of Optimal Phi Values Across Methods") +
    scale_color_manual(values = c("MRMode" = "blue", "mode_new" = "red")) +
    theme_minimal()
  
  print(phi_hist_plot)
  
  base_path <- "./plots"
  filename_comp <- sprintf("phiDepen_comparison.eps",
                           base_path)
  ggsave(filename_comp, plot = power_plot, device = "eps", width = 10, height = 6)
}

analyze_all_scenarios_error_alt <- function(est_theta,
                                            theta_alt_vec = seq(-1, 1, by = 0.1),
                                            NxNy_alt_ratio_vec = c(1:5)) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(knitr)
  
  error_df <- data.frame()
  methods <- c("MRMode", "mode_new")
  
  for (theta_alt in theta_alt_vec) {
    for (NxNy_alt_ratio in NxNy_alt_ratio_vec) {
      file_path <- sprintf(
        "../../../results/simulation_results/mode_comp/ThetaU_Nalt_est_theta%g_thetaU0.3_N80000_prop_invalid0.5_phi1_overlap0.3_theta_alt%g_NxNy_alt_ratio%g.rda",
        est_theta, theta_alt, NxNy_alt_ratio
      )
      
      if (!file.exists(file_path)) {
        warning(paste("File does not exist:", file_path))
        next
      }
      
      load(file_path)
      est <- na.omit(est)
      
      if (!all(c(methods, paste0(methods, "_se")) %in% colnames(est))) {
        warning("Missing required columns in 'est'")
        next
      }
      
      if (est_theta == 0) {
        error_rate <- colMeans(abs(est[, methods]) > 1.96 * est[, paste0(methods, "_se")]) * 100
        error_label <- "Type 1 Error"
      } else {
        type2_error_rate <- colMeans(abs(est[, methods]) <= 1.96 * est[, paste0(methods, "_se")]) * 100
        error_rate <- 100 - type2_error_rate
        error_label <- "Power (1 - Type 2 Error)"
      }
      
      error_df_scenario <- data.frame(
        Method = methods,
        ErrorRate = error_rate,
        theta_alt = theta_alt,
        NxNy_alt_ratio = NxNy_alt_ratio,
        error_label = error_label
      )
      
      error_df <- rbind(error_df, error_df_scenario)
    }
  }
  
  if (nrow(error_df) == 0) {
    warning("No data available to plot.")
    return(NULL)
  }
  
  if (est_theta == 0) {
    # ---- Plot Type 1 Error heatmap ----
    p_t1 <- ggplot(error_df, aes(x = NxNy_alt_ratio, y = theta_alt, fill = ErrorRate)) +
      geom_tile(color = "white") +
      facet_wrap(~ Method) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 5,
                           name = "Type 1 Error (%)") +
      labs(
        title = "Type 1 Error across Scenarios",
        x = "NxNy_alt_ratio",
        y = expression(theta[alt])
      ) +
      theme_minimal(base_size = 14)
    
    print(p_t1)
    
  } else {
    # ---- Plot Power difference heatmap ----
    power_wide <- error_df %>%
      filter(error_label == "Power (1 - Type 2 Error)") %>%
      pivot_wider(names_from = Method, values_from = ErrorRate) %>%
      mutate(PowerDiff = mode_new - MRMode)
    
    p_diff <- ggplot(power_wide, aes(x = NxNy_alt_ratio, y = theta_alt, fill = PowerDiff)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           name = "mode_new - MRMode\nPower Diff") +
      labs(
        title = "Power Difference: mode_new vs MRMode",
        x = "NxNy_alt_ratio",
        y = expression(theta[alt])
      ) +
      theme_minimal(base_size = 14)
    
    print(p_diff)
    
    # ---- Summary Tables ----
    power_df <- error_df %>%
      filter(error_label == "Power (1 - Type 2 Error)")
    
    # 1. Average Power by Method
    avg_power_method <- power_df %>%
      group_by(Method) %>%
      summarise(AveragePower = mean(ErrorRate)) %>%
      arrange(desc(AveragePower))
    
    cat("\n\n### Average Power Across All Scenarios (by Method):\n")
    print(knitr::kable(avg_power_method, digits = 2))
    
    # 2. Average Power by theta_alt
    avg_power_theta <- power_df %>%
      group_by(Method, theta_alt) %>%
      summarise(AveragePower = mean(ErrorRate)) %>%
      arrange(Method, theta_alt)
    
    cat("\n\n### Average Power by Method and theta_alt:\n")
    print(knitr::kable(avg_power_theta, digits = 2))
    
    # 3. Average Power by NxNy_alt_ratio
    avg_power_ratio <- power_df %>%
      group_by(Method, NxNy_alt_ratio) %>%
      summarise(AveragePower = mean(ErrorRate)) %>%
      arrange(Method, NxNy_alt_ratio)
    
    cat("\n\n### Average Power by Method and NxNy_alt_ratio:\n")
    print(knitr::kable(avg_power_ratio, digits = 2))
  }
}

#analyze_all_scenarios_error_alt(est_theta = 0.2)


# analyse_new_scenarios_eff(power_dir = -1, thetaUvec = c(0.3, 0.5),
#                           Nvec = c("50000", "80000", "100000", "200000"),
#                           prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7),
#                           overlap_vec = c(0, 0.1, 0.3, 0.5, 0.75, 1),
#                           phi_range = exp(seq(log(0.1), log(10), length.out = 10)),
#                           scenario_vec = c("BI", "BN", "DI", "DN"), color = "overlap",plot_type = 0)
# # 
# analyse_new_scenarios_phi_hist(thetaUvec = c(0.3, 0.5),
#                                Nvec = c("80000", "100000", "200000"),
#                                prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7),
#                                overlap_vec = c(0, 0.1, 0.3, 0.5, 0.75, 1),
#                                phi_range = exp(seq(log(0.1), log(10), length.out = 10)),
#                                scenario_vec = c("BI", "BN", "DI", "DN"))
# analyze_all_scenarios_coverage(est_theta = 0.2,
#                                            thetaU = 0.5,
#                                            N = 50000,
#                                            prop_invalid = 0.5,
#                                            overlap = 0.75,
#                                            scenario = "BI",
#                                            sctfc = FALSE)
# 
# analyze_all_scenarios_error(est_theta = 0.2,
#                             thetaU = 0.5,
#                             N = 80000,
#                             prop_invalid = 0.5,
#                             overlap = 0.75,
#                             scenario = "DI",
#                             sctfc = FALSE)

# analyze_all_scenarios_mse(est_theta = 0.2,
#                             thetaU = 0.5,
#                             N = 80000,
#                             prop_invalid = 0.5,
#                             overlap = 0.75,
#                             scenario = "DI",
#                             sctfc = FALSE)


#####################################################IB-Mix

# analyse_mrmix_ibmix_power <- function(power_dir = 1,
#                                       thetaUvec = c(0.3,0.5),
#                                       Nvec = c("50000", "80000", "100000", "150000", "200000"),
#                                       prop_invalid_vec = c(0.3, 0.5, 0.7),
#                                       overlap_vec = c(0.5, 0.75, 1),
#                                       scenario_vec = c("BI", "BN", "DI", "DN"),
#                                       color_var = "overlap", plot_type = 1) {
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(RColorBrewer)
# 
#   power_results <- data.frame()
# 
#   for (scenario in scenario_vec) {
#     for (thetaU in thetaUvec) {
#       for (N in Nvec) {
#         for (prop_invalid in prop_invalid_vec) {
#           for (overlap in overlap_vec) {
# 
#             file_path_power <- if (power_dir == 1) {
#               sprintf("./results/%s_est_theta0.2_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                       scenario, thetaU, format(N, scientific = FALSE), prop_invalid, overlap)
#             } else if (power_dir == 2) {
#               sprintf("./results/%s_est_theta-0.2_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                       scenario, thetaU, format(N, scientific = FALSE), prop_invalid, overlap)
#             } else if (power_dir == 3) {
#               sprintf("./results/%s_est_theta0.1_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                       scenario, thetaU, format(N, scientific = FALSE), prop_invalid, overlap)
#             } else if (power_dir == 4) {
#               sprintf("./results/%s_est_theta-0.1_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                       scenario, thetaU, format(N, scientific = FALSE), prop_invalid, overlap)
#             }
# 
#             if (!file.exists(file_path_power)) next
# 
#             load(file_path_power)
#             if (!exists("est") || nrow(est) == 0) next
# 
#             est <- na.omit(est)
#             methods <- c("MR-Mix", "IB-Mix")
#             se_methods <- paste0(methods, "_se")
# 
#             if (!all(c(methods, se_methods) %in% colnames(est))) next
# 
#             type2_error <- colMeans(abs(est[, methods]) <= 1.96 * est[, se_methods]) * 100
#             power_vals <- pmin(100, 100 - type2_error)
# 
#             power_data <- data.frame(Method = methods,
#                                      Power = power_vals,
#                                      scenario = scenario,
#                                      thetaU = thetaU,
#                                      N = N,
#                                      prop_invalid = prop_invalid,
#                                      overlap = overlap)
# 
#             power_results <- rbind(power_results, power_data)
#           }
#         }
#       }
#     }
#   }
# 
#   # after pivoting:
#   power_results_wide <- power_results %>%
#     pivot_wider(names_from = Method, values_from = Power) %>%
#     mutate(overlap_label = factor(paste0("overlap = ", overlap)))
# 
#   unique_N <- length(unique(power_results_wide$N))
# 
#   power_plot <- ggplot(power_results_wide,
#                        aes(x = `MR-Mix`,
#                            y = `IB-Mix`,
#                            color = prop_invalid,
#                            shape = factor(N))) +
#     geom_point(size = 3) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     facet_wrap(~ overlap_label, scales = "free") +    # <-- here’s your 3 panels
#     labs(
#       x = "Power of MR-Mix",
#       y = "Power of IB-Mix",
#       title = "Power Comparison by Overlap",
#       color = "Invalid %",
#       shape = "Sample size\n(N)"
#     ) +
#     scale_color_viridis_c(option = "turbo") +
#     scale_shape_manual(values = seq_len(unique_N)) +
#     theme_minimal() +
#     theme(
#       legend.position = "bottom",
#       plot.title     = element_text(hjust = 0.5)
#     )
# 
#   print(power_plot)
#   base_path <- "./plots"
#   filename    <- file.path(base_path, "mrmix_ibmix_power_comparison.eps")
#   ggsave(filename, plot = power_plot,
#          device = "eps", width = 10, height = 6)
# 
# }


analyse_mrmix_ibmix_power <- function(thetaUvec = c(0.3,0.5),
                                      Nvec = c("50000", "80000", "100000", "150000", "200000"),
                                      prop_invalid_vec = c(0.3, 0.5, 0.7),
                                      overlap_vec = c(0.5, 0.75, 1),
                                      scenario_vec = c("BI", "BN", "DI", "DN"),
                                      color_var = "overlap", plot_type = 1) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  
  power_results <- data.frame()
  
  for (scenario in scenario_vec) {
    for (thetaU in thetaUvec) {
      for (N in Nvec) {
        for (prop_invalid in prop_invalid_vec) {
          for (overlap in overlap_vec) {
            
            file_path_power <-
              sprintf("./results/%s_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
                      scenario, thetaU, format(N, scientific = FALSE), prop_invalid, overlap)
            
            
            if (!file.exists(file_path_power)) next
            
            load(file_path_power)
            if (!exists("est") || nrow(est) == 0) next
            
            est <- na.omit(est)
            methods <- c("MR-Mix", "IB-Mix")
            se_methods <- paste0(methods, "_se")
            
            if (!all(c(methods, se_methods) %in% colnames(est))) next
            
            type2_error <- colMeans(abs(est[, methods]) <= 1.96 * est[, se_methods]) * 100
            power_vals <- pmin(100, 100 - type2_error)
            
            power_data <- data.frame(Method = methods,
                                     Power = power_vals,
                                     scenario = scenario,
                                     thetaU = thetaU,
                                     N = N,
                                     prop_invalid = prop_invalid,
                                     overlap = overlap)
            
            power_results <- rbind(power_results, power_data)
          }
        }
      }
    }
  }
  
  # after pivoting:
  power_results_wide <- power_results %>%
    pivot_wider(names_from = Method, values_from = Power) %>%
    mutate(overlap_label = factor(paste0("overlap = ", overlap)))
  
  unique_N <- length(unique(power_results_wide$N))
  
  power_plot <- ggplot(power_results_wide,
                       aes(x = `MR-Mix`,
                           y = `IB-Mix`,
                           color = prop_invalid,
                           shape = factor(N))) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ overlap_label, scales = "free") +    # <-- here’s your 3 panels
    labs(
      x = "Power of MR-Mix",
      y = "Power of IB-Mix",
      title = "Power Comparison by Overlap",
      color = "Invalid %",
      shape = "Sample size\n(N)"
    ) +
    scale_color_viridis_c(option = "turbo") +
    scale_shape_manual(values = seq_len(unique_N)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title     = element_text(hjust = 0.5)
    )
  
  print(power_plot)
  base_path <- "./plots"
  filename    <- file.path(base_path, "mrmix_ibmix_power_comparison.eps")
  ggsave(filename, plot = power_plot,
         device = "eps", width = 10, height = 6)
  
}

