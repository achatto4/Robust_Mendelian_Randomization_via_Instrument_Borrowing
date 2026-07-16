library(ggplot2)
if (!require("RColorBrewer")) install.packages("RColorBrewer", dependencies = TRUE)
library(RColorBrewer)

# Define the methods used for extracting values from the est object
METHODS <- c(
  "IVW", "median", "Egger",
  "Cont-Mix", "MR-cML", "MR-Mix",
  "mode_new_phi1", "MR-PRESSO",
  "MRMode", "IB-PRESSO"
)

NEW_ORDER <- c(
  "IVW", "median", "Egger",
  "Cont-Mix", "MR-cML", "MR-Mix",
  "MRMode", "mode_new_phi1", "MR-PRESSO",
  "IB-PRESSO"
)

NEW_ORDER1 <-  c(
  "IVW", "Median", "Egger",
  "Cont-Mix", "MR-cML", "MR-Mix",
   "MR-Mode", "IB-MODE", "MR-PRESSO",
  "IB-PRESSO"
)
##COVERAGE

analyze_all_scenarios_multi_coverage <- function(est_theta,
                                                 scenario,
                                                 thetaUvec,
                                                 Nvec,
                                                 prop_invalid_vec,
                                                 overlap_vec,
                                                 sctfc = FALSE) {
  library(ggplot2)
  library(dplyr)

  facets <- c()
  if (length(thetaUvec) > 1) facets <- c(facets, "thetaU")
  if (length(Nvec) > 1) facets <- c(facets, "N")
  if (length(prop_invalid_vec) > 1) facets <- c(facets, "prop_invalid")
  if (length(overlap_vec) > 1) facets <- c(facets, "overlap")

  if (length(facets) == 0) stop("At least one parameter must be a vector for faceting.")
  if (length(facets) > 2) stop("At most two parameters can be vectors for faceting.")

  df_all <- data.frame()
  # methods <- c("MRMode", "Cont-Mix", "MR-Mix", "mode_new_phi1", "MR-PRESSO",
  #              "MR-cML", "IVW", "median", "Egger", "IB-PRESSO")
  methods <- METHODS

  for (thetaU in thetaUvec) {
    for (N in Nvec) {
      for (prop_invalid in prop_invalid_vec) {
        for (overlap in overlap_vec) {

          N_str <- if (sctfc) format(N, scientific = TRUE) else as.character(N)

          file_path <- sprintf(
            "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
            scenario, est_theta, thetaU, N_str, prop_invalid, overlap
          )

          if (!file.exists(file_path)) {
            warning(paste("File does not exist:", file_path))
            next
          }

          load(file_path)
          est <- na.omit(est)

          coverage <- colMeans(
            (est[, methods] - 1.96 * est[, paste0(methods, "_se")] <= est_theta) &
              (est[, methods] + 1.96 * est[, paste0(methods, "_se")] >= est_theta)
          ) * 100

          df_scenario <- data.frame(
            Method = methods,
            Coverage = coverage,
            Scenario = scenario
          )

          # Add facet values
          if ("thetaU" %in% facets) {
            if (!"Facet1" %in% names(df_scenario)) {
              df_scenario$Facet1 <- thetaU
              facet_labels <- list(Facet1 = "thetaU")
            } else {
              df_scenario$Facet2 <- thetaU
              facet_labels$Facet2 <- "thetaU"
            }
          }
          if ("N" %in% facets) {
            if (!"Facet1" %in% names(df_scenario)) {
              df_scenario$Facet1 <- factor(N, levels = Nvec)
              facet_labels <- list(Facet1 = "N")
            } else {
              df_scenario$Facet2 <- factor(N, levels = Nvec)
              facet_labels$Facet2 <- "N"
            }
          }
          if ("prop_invalid" %in% facets) {
            if (!"Facet1" %in% names(df_scenario)) {
              df_scenario$Facet1 <- prop_invalid
              facet_labels <- list(Facet1 = "prop_invalid")
            } else {
              df_scenario$Facet2 <- prop_invalid
              facet_labels$Facet2 <- "prop_invalid"
            }
          }
          if ("overlap" %in% facets) {
            if (!"Facet1" %in% names(df_scenario)) {
              df_scenario$Facet1 <- overlap
              facet_labels <- list(Facet1 = "overlap")
            } else {
              df_scenario$Facet2 <- overlap
              facet_labels$Facet2 <- "overlap"
            }
          }

          df_all <- rbind(df_all, df_scenario)
        }
      }
    }
  }
  if (nrow(df_all) == 0) stop("No data was loaded into df_all. Check file paths or NA filtering.")

  # 1. Define the new horizontal order, grouping the two “mode” and two “mix” at the end:
  # new_order <- c(
  #   "IVW", "median", "Egger","Cont-Mix", "MR-cML",
  #   "mode_new_phi1", "MR-PRESSO", "MRMode",
  #   "MR-Mix", "IB-PRESSO"
  # )
  new_order <- NEW_ORDER

  # 2. Apply to df_all:
  df_all$Method <- factor(df_all$Method, levels = new_order)

  facet1_levels <- length(unique(df_all$Facet1))
  facet2_levels <- length(unique(df_all$Facet2))

  facet_formula <- if (facet1_levels > 1 & facet2_levels == 1) {
    Facet1 ~ .
  } else if (facet2_levels > 1 & facet1_levels == 1) {
    . ~ Facet2
  } else {
    Facet1 ~ Facet2
  }

  facet_labeller <- labeller(
    Facet1 = setNames(
      paste0(facet_labels$Facet1, " = ", unique(df_all$Facet1)),
      unique(df_all$Facet1)
    ),
    Facet2 = setNames(
      paste0(facet_labels$Facet2, " = ", unique(df_all$Facet2)),
      unique(df_all$Facet2)
    )
  )

  # Compute dynamic y-axis lower bound
  ymin_dynamic <- max(min(df_all$Coverage) - 10, 0)

  p_coverage <- ggplot(df_all, aes(x = Method, y = Coverage, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6, color = "black") +
    geom_vline(xintercept = c(6.5),
               linetype   = "dotted",
               linewidth = 0.3) +
    facet_grid(facet_formula, labeller = facet_labeller) +
    labs(
      title = "Coverage: True Value Within 95% CI",
      y     = "Coverage (%)",
      x     = ""
    ) +
    scale_fill_manual(
      values = brewer.pal(length(new_order), "Set3"),
      breaks = new_order,
      labels = NEW_ORDER1
    ) +
    # instead of scale_y_continuous(limits=...):
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(ymin_dynamic, NA)) +   # <— this just zooms, doesn't drop data
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(size = 8),
      axis.title       = element_text(size = 11),
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 9)
    )

  # Define scenario as a variable
  scenario_str <- scenario  # or set this dynamically if within a function

  # Construct base path
  base_path <- "./plots"

  # Construct parts of filename from parameters
  thetaU_str <- paste(thetaUvec, collapse = "-")
  N_str <- paste(Nvec, collapse = "-")
  prop_invalid_str <- paste(prop_invalid_vec, collapse = "-")
  overlap_str <- paste(overlap_vec, collapse = "-")

 # --- Save Coverage plot ---
  filename_coverage <- sprintf("%s/coverage_plot_est_%s_scenario_%s_thetaU_%s_N_%s_propInvalid_%s_overlap_%s.eps",
                               base_path, est_theta, scenario_str, thetaU_str, N_str, prop_invalid_str, overlap_str)
  ggsave(filename_coverage, plot = p_coverage, device = "eps", width = 10, height = 6)

  #print(p_coverage)
}

analyze_all_scenarios_multi_error <- function(est_theta,
                                              scenario,
                                              thetaUvec,
                                              Nvec,
                                              prop_invalid_vec,
                                              overlap_vec,
                                              sctfc = TRUE) {
  
  library(ggplot2)
  library(dplyr)
  
  facets <- c()
  if (length(thetaUvec) > 1) facets <- c(facets, "thetaU")
  if (length(Nvec) > 1) facets <- c(facets, "N")
  if (length(prop_invalid_vec) > 1) facets <- c(facets, "prop_invalid")
  if (length(overlap_vec) > 1) facets <- c(facets, "overlap")
  
  if (length(facets) == 0) stop("At least one parameter must be a vector for faceting.")
  if (length(facets) > 2) stop("At most two parameters can be vectors for faceting.")
  
  error_df <- data.frame()
  methods <- METHODS
  
  for (thetaU in thetaUvec) {
    for (N in Nvec) {
      for (prop_invalid in prop_invalid_vec) {
        for (overlap in overlap_vec) {
          
          N_str <- if (sctfc) format(N, scientific = TRUE) else as.character(N)
          
          file_path <- sprintf(
            "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
            scenario, est_theta, thetaU, N_str, prop_invalid, overlap
          )
          
          if (!file.exists(file_path)) {
            warning(paste("File does not exist:", file_path))
            next
          }
          
          load(file_path)
          est <- na.omit(est)
          
          if (est_theta == 0) {
            error_rate <- colMeans(abs(est[, methods]) > 1.96 * est[, paste0(methods, "_se")]) * 100
            error_label <- "Type 1 Error"
            
          } else {
            type2_error_rate <- colMeans(abs(est[, methods]) <= 1.96 * est[, paste0(methods, "_se")]) * 100
            error_rate <- 100 - type2_error_rate
            error_label <- "Power (1 - Type 2 Error)"
          }
          
          error_df_scenario <- data.frame(Method = methods, ErrorRate = error_rate, Scenario = scenario)
          
          # Add facets and labels
          if ("thetaU" %in% facets) {
            if (!"Facet1" %in% names(error_df_scenario)) {
              error_df_scenario$Facet1 <- thetaU
              facet_labels <- list(Facet1 = "thetaU")
            } else {
              error_df_scenario$Facet2 <- thetaU
              facet_labels$Facet2 <- "thetaU"
            }
          }
          if ("N" %in% facets) {
            if (!"Facet1" %in% names(error_df_scenario)) {
              error_df_scenario$Facet1 <- factor(N, levels = Nvec)
              facet_labels <- list(Facet1 = "N")
            } else {
              error_df_scenario$Facet2 <- factor(N, levels = Nvec)
              facet_labels$Facet2 <- "N"
            }
          }
          if ("prop_invalid" %in% facets) {
            if (!"Facet1" %in% names(error_df_scenario)) {
              error_df_scenario$Facet1 <- prop_invalid
              facet_labels <- list(Facet1 = "prop_invalid")
            } else {
              error_df_scenario$Facet2 <- prop_invalid
              facet_labels$Facet2 <- "prop_invalid"
            }
          }
          if ("overlap" %in% facets) {
            if (!"Facet1" %in% names(error_df_scenario)) {
              error_df_scenario$Facet1 <- overlap
              facet_labels <- list(Facet1 = "overlap")
            } else {
              error_df_scenario$Facet2 <- overlap
              facet_labels$Facet2 <- "overlap"
            }
          }
          
          error_df <- rbind(error_df, error_df_scenario)
        }
      }
    }
  }
  
  if (nrow(error_df) == 0) stop("No data was loaded into df_all. Check file paths or NA filtering.")
  
  new_order <- NEW_ORDER
  error_df$Method <- factor(error_df$Method, levels = new_order)
  
  facet1_levels <- length(unique(error_df$Facet1))
  facet2_levels <- length(unique(error_df$Facet2))
  
  facet_formula <- if (facet1_levels > 1 & facet2_levels == 1) {
    Facet1 ~ .
  } else if (facet2_levels > 1 & facet1_levels == 1) {
    . ~ Facet2
  } else {
    Facet1 ~ Facet2
  }
  
  facet_labeller <- labeller(
    Facet1 = setNames(
      paste0(facet_labels$Facet1, " = ", unique(error_df$Facet1)),
      unique(error_df$Facet1)
    ),
    Facet2 = setNames(
      paste0(facet_labels$Facet2, " = ", unique(error_df$Facet2)),
      unique(error_df$Facet2)
    )
  )
  
  p_error <- ggplot(error_df, aes(x = Method, y = ErrorRate, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6, color = "black") +
    geom_vline(xintercept = c(6.5),
               linetype   = "dotted",
               linewidth = 0.3) +
    facet_grid(facet_formula,  labeller = facet_labeller) +
    scale_fill_manual(
      values = brewer.pal(length(new_order), "Set3"),
      breaks = new_order,
      labels = NEW_ORDER1
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9)
    ) +
    labs(
      title = paste0(error_label, " by Method"),
      y = paste0(error_label, " (%)"),
      x = ""
    )
  
  if (est_theta == 0) {
    p_error <- p_error +
      geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 0.7)
  }
  
  scenario_str <- scenario
  base_path <- "./plots"
  thetaU_str <- paste(thetaUvec, collapse = "-")
  N_str <- paste(Nvec, collapse = "-")
  prop_invalid_str <- paste(prop_invalid_vec, collapse = "-")
  overlap_str <- paste(overlap_vec, collapse = "-")
  
  filename_error <- sprintf("%s/error_plot_est_%s_scenario_%s_thetaU_%s_N_%s_propInvalid_%s_overlap_%s.eps",
                            base_path, est_theta, scenario_str, thetaU_str, N_str, prop_invalid_str, overlap_str)
  ggsave(filename_error, plot = p_error, device = "eps", width = 10, height = 6)
  # p_error
}

analyze_all_scenarios_multi_mse <- function(
    est_theta,
    scenario,
    thetaUvec,
    Nvec,
    prop_invalid_vec,
    overlap_vec,
    sctfc = FALSE,
    use_log_scale = FALSE,
    truncation_multiplier = 4
) {
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # ============================================================================
  # 1. PARAMETER VALIDATION & FACET SETUP
  # ============================================================================
  
  facets <- c()
  if (length(thetaUvec)        > 1) facets <- c(facets, "thetaU")
  if (length(Nvec)             > 1) facets <- c(facets, "N")
  if (length(prop_invalid_vec) > 1) facets <- c(facets, "prop_invalid")
  if (length(overlap_vec)      > 1) facets <- c(facets, "overlap")
  
  if (length(facets) == 0) stop("At least one parameter must vary for faceting.")
  if (length(facets) > 2)  stop("At most two parameters can vary for faceting.")
  
  methods   <- METHODS
  new_order <- NEW_ORDER
  
  # ============================================================================
  # 2. DATA LOADING
  # ============================================================================
  
  mse_df <- data.frame()
  
  for (thetaU in thetaUvec) {
    for (N in Nvec) {
      for (prop_invalid in prop_invalid_vec) {
        for (overlap in overlap_vec) {
          
          N_str <- if (sctfc) format(N, scientific = TRUE) else as.character(N)
          file_path <- sprintf(
            "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
            scenario, est_theta, thetaU, N_str, prop_invalid, overlap
          )
          
          if (!file.exists(file_path)) {
            warning("Missing file: ", file_path)
            next
          }
          
          load(file_path)
          est <- na.omit(est)
          mse <- colMeans((est[, methods] - est_theta)^2, na.rm = TRUE)
          
          tmp <- data.frame(
            Method       = methods,
            MSE          = mse,
            thetaU       = thetaU,
            N            = N,
            prop_invalid = prop_invalid,
            overlap      = overlap,
            stringsAsFactors = FALSE
          )
          
          mse_df <- bind_rows(mse_df, tmp)
        }
      }
    }
  }
  
  if (nrow(mse_df) == 0) stop("No data loaded; check file paths.")
  
  # ============================================================================
  # 3. DATA PREPARATION
  # ============================================================================
  
  mse_df$Method <- factor(mse_df$Method, levels = new_order)
  
  # Create facet columns
  if (length(facets) == 1) {
    if (facets[1] == "N") {
      mse_df$Facet1 <- factor(mse_df$N, levels = Nvec)
    } else {
      mse_df$Facet1 <- mse_df[[facets[1]]]
    }
    mse_df$Facet2 <- NA
    facet_label1 <- facets[1]
    facet_label2 <- NULL
  } else {
    if (facets[1] == "N") {
      mse_df$Facet1 <- factor(mse_df$N, levels = Nvec)
    } else {
      mse_df$Facet1 <- mse_df[[facets[1]]]
    }
    if (facets[2] == "N") {
      mse_df$Facet2 <- factor(mse_df$N, levels = Nvec)
    } else {
      mse_df$Facet2 <- mse_df[[facets[2]]]
    }
    facet_label1 <- facets[1]
    facet_label2 <- facets[2]
  }
  
  # ============================================================================
  # 4. ROW-AWARE TRUNCATION LOGIC (4× MEDIAN)
  # ============================================================================
  
  # Identify which facet represents rows
  f1_n <- n_distinct(mse_df$Facet1)
  f2_n <- n_distinct(mse_df$Facet2)
  
  # Determine row grouping variable
  if (f1_n > 1 && f2_n == 1) {
    row_var <- "Facet1"
    mse_df$RowGroup <- mse_df$Facet1
  } else if (f2_n > 1 && f1_n == 1) {
    row_var <- "Facet2"
    mse_df$RowGroup <- mse_df$Facet2
  } else {
    row_var <- "Facet1"
    mse_df$RowGroup <- mse_df$Facet1
  }
  
  # Apply truncation WITHIN each row (shared y-axis)
  mse_df <- mse_df %>%
    group_by(RowGroup) %>%
    mutate(
      median_mse_in_row = median(MSE, na.rm = TRUE),
      mad_mse_in_row = mad(MSE, na.rm = TRUE),
      truncate_threshold = truncation_multiplier * median_mse_in_row,
      is_outlier = MSE > truncate_threshold,
      MSE_plot = pmin(MSE, truncate_threshold),
      show_actual = is_outlier & (MSE > truncate_threshold),
      fold_over_threshold = MSE / truncate_threshold
    ) %>%
    ungroup()
  
  # ============================================================================
  # 5. BUILD FACET FORMULA AND LABELLER
  # ============================================================================
  
  facet_formula <- if (f1_n > 1 && f2_n == 1) {
    Facet1 ~ .
  } else if (f2_n > 1 && f1_n == 1) {
    . ~ Facet2
  } else {
    Facet1 ~ Facet2
  }
  
  lab1 <- setNames(
    paste0(facet_label1, " = ", unique(mse_df$Facet1)),
    unique(mse_df$Facet1)
  )
  
  if (!is.null(facet_label2)) {
    lab2 <- setNames(
      paste0(facet_label2, " = ", unique(mse_df$Facet2)),
      unique(mse_df$Facet2)
    )
  } else {
    lab2 <- label_value
  }
  
  # ============================================================================
  # 6. CREATE PLOT (LOG SCALE OR TRUNCATED)
  # ============================================================================
  
  if (use_log_scale) {
    # ========== LOG SCALE PLOT ==========
    
    mse_df$MSE_plot_log <- mse_df$MSE + 1e-10
    
    plot_mse <- ggplot(mse_df, aes(x = Method, y = MSE_plot_log, fill = Method)) +
      geom_bar(
        stat = "identity", 
        position = position_dodge(0.8),
        width = 0.65, 
        color = "black", 
        linewidth = 0.3
      ) +
      geom_vline(xintercept = 6.5, linetype = "dotted", linewidth = 0.4) +
      
      geom_text(
        data = filter(mse_df, is_outlier),
        aes(label = sprintf("%.2e", MSE)),
        size = 2.8, 
        vjust = -0.3, 
        color = "red", 
        fontface = "bold"
      ) +
      
      facet_grid(
        facet_formula,
        scales   = "free_y",
        labeller = labeller(Facet1 = lab1, Facet2 = lab2)
      ) +
      
      scale_fill_manual(
        values = brewer.pal(length(new_order), "Set3"),
        breaks = new_order,
        labels = NEW_ORDER1
      ) +
      
      scale_y_log10(
        expand = expansion(mult = c(0, 0.1)),
        labels = scales::label_scientific()
      ) +
      
      annotation_logticks(sides = "l", size = 0.3, alpha = 0.5) +
      
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_text(size = 9),
        axis.title       = element_text(size = 11, face = "bold"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.caption     = element_text(size = 9, hjust = 0, color = "grey40"),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text       = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
        panel.spacing    = unit(0.8, "lines"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 9)
      ) +
      labs(
        title = "Mean Squared Error by Method (Log Scale)",
        y = "MSE (log₁₀ scale)",
        x = "",
        caption = "Note: Y-axis uses logarithmic scale. Red values indicate methods with exceptionally high MSE."
      )
    
  } else {
    # ========== TRUNCATED PLOT - ONLY TEXT ON TOP ==========
    
    plot_mse <- ggplot(mse_df, aes(x = Method, y = MSE_plot, fill = Method)) +
      geom_bar(
        stat = "identity", 
        position = position_dodge(0.8),
        width = 0.65, 
        color = "black", 
        linewidth = 0.3
      ) +
      
      geom_vline(xintercept = 6.5, linetype = "dotted", linewidth = 0.4) +
      
      # ONLY actual MSE values ON TOP of truncated bars - small red text
      geom_text(
        data = filter(mse_df, show_actual),
        aes(x = Method, y = MSE_plot, label = sprintf("(%.4f)", MSE)),
        color = "red", 
        size = 2.0,              # Small text
        fontface = "bold",
        vjust = -0.3,            # Just above the bar
        inherit.aes = FALSE
      ) +
      
      facet_grid(
        facet_formula,
        scales   = "free_y",
        labeller = labeller(Facet1 = lab1, Facet2 = lab2)
      ) +
      
      scale_fill_manual(
        values = brewer.pal(length(new_order), "Set3"),
        breaks = new_order,
        labels = NEW_ORDER1
      ) +
      
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +  # Extra space for text
      
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_text(size = 9),
        axis.title       = element_text(size = 11, face = "bold"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.caption     = element_text(size = 9, hjust = 0, color = "grey40"),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text       = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(0.8, "lines"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 9)
      ) +
      labs(
        title = "Mean Squared Error by Method",
        y = "Mean Squared Error",
        x = ""
      )
  }
  
  # ============================================================================
  # 7. DIAGNOSTIC OUTPUT
  # ============================================================================
  
  truncation_summary <- mse_df %>%
    filter(show_actual) %>%
    group_by(RowGroup) %>%
    summarise(
      n_truncated = n(),
      methods_truncated = paste(unique(Method), collapse = ", "),
      median_mse_row = first(median_mse_in_row),
      truncate_at = first(truncate_threshold),
      max_actual_MSE = max(MSE),
      max_fold_over = max(fold_over_threshold),
      .groups = "drop"
    )
  
  if (nrow(truncation_summary) > 0) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("TRUNCATION SUMMARY (by row)\n")
    cat(rep("=", 70), "\n", sep = "")
    cat(sprintf("Truncation rule: %d× median MSE\n\n", truncation_multiplier))
    print(truncation_summary, n = Inf)
    cat(rep("=", 70), "\n\n", sep = "")
  } else {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("No methods exceeded truncation threshold (", 
        truncation_multiplier, "× median MSE)\n", sep = "")
    cat(rep("=", 70), "\n\n", sep = "")
  }
  
  method_summary <- mse_df %>%
    group_by(Method) %>%
    summarise(
      mean_MSE = mean(MSE, na.rm = TRUE),
      median_MSE = median(MSE, na.rm = TRUE),
      sd_MSE = sd(MSE, na.rm = TRUE),
      min_MSE = min(MSE, na.rm = TRUE),
      max_MSE = max(MSE, na.rm = TRUE),
      n_times_truncated = sum(show_actual, na.rm = TRUE),
      max_fold_over_threshold = ifelse(
        any(show_actual), 
        max(fold_over_threshold[show_actual], na.rm = TRUE), 
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    arrange(median_MSE)
  
  cat("METHOD PERFORMANCE SUMMARY\n")
  cat(rep("=", 70), "\n", sep = "")
  print(method_summary, n = Inf)
  cat(rep("=", 70), "\n\n", sep = "")
  
  row_stats <- mse_df %>%
    group_by(RowGroup) %>%
    summarise(
      median_mse = first(median_mse_in_row),
      mad_mse = first(mad_mse_in_row),
      truncate_threshold = first(truncate_threshold),
      n_outliers = sum(is_outlier),
      .groups = "drop"
    )
  
  cat("ROW-SPECIFIC STATISTICS\n")
  cat(rep("=", 70), "\n", sep = "")
  print(row_stats, n = Inf)
  cat(rep("=", 70), "\n\n", sep = "")
  
  # ============================================================================
  # 8. SAVE PLOTS
  # ============================================================================
  
  base_path <- "./plots"
  
  thetaU_str       <- paste(thetaUvec, collapse = "-")
  N_str            <- paste(Nvec, collapse = "-")
  prop_invalid_str <- paste(prop_invalid_vec, collapse = "-")
  overlap_str      <- paste(overlap_vec, collapse = "-")
  
  scale_suffix <- if (use_log_scale) "_logscale" else ""
  
  filename_error <- sprintf(
    "%s/mse_plot%s_est_%s_scenario_%s_thetaU_%s_N_%s_propInvalid_%s_overlap_%s.eps",
    base_path, scale_suffix, est_theta, scenario, thetaU_str, N_str, 
    prop_invalid_str, overlap_str
  )
  
  ggsave(filename_error, plot = plot_mse, device = "eps", width = 12, height = 7)
  
  cat("Plot saved to:", filename_error, "\n\n")
}
