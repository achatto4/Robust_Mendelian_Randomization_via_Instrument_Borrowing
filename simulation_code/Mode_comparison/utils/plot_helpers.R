library(ggplot2)
library(tidyr)

# file_path <- sprintf("./results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                      scenario, est_theta, thetaU, format(N, scientific = sctfc), prop_invalid,  overlap)

#setwd(".")  # set working directory as needed

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
