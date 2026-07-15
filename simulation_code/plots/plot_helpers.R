library(ggplot2)

# file_path <- sprintf("./results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#                      scenario, est_theta, thetaU, format(N, scientific = sctfc), prop_invalid,  overlap)

#setwd(".")  # set working directory as needed

# Define a function to generalize the analysis and plot all methods and scenarios in a 2x2 layout
analyze_all_scenarios <- function(est_theta,
                                  thetaU,
                                  N,
                                  prop_invalid,
                                  overlap,
                                  sctfc = FALSE) {
  # Define the specific scenario filenames
  scenarios <- c("BI", "DN", "BN", "DI")
  df_all <- data.frame()
  
  for (scenario in scenarios) {
    # Construct the file path based on input parameters for each scenario
    file_path <- sprintf(
      "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      overlap
    )
    # Load the data
    load(file_path)
    est <- na.omit(est)
    # Extract estimates and standard deviations for all methods
    estimates <- apply(est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )], 2, mean)
    se <- apply(est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )], 2, sd)
    methods <- c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )
    
    # Create a data frame for the current scenario and add a scenario column
    df_scenario <- data.frame(
      Method = methods,
      Estimate = estimates,
      SE = se,
      Scenario = scenario
    )
    
    # Bind the scenario data to the main data frame
    df_all <- rbind(df_all, df_scenario)
  }
  
  # Create a 2x2 grid of plots
  p <- ggplot(df_all,
              aes(
                x = Method,
                y = Estimate,
                color = Scenario,
                group = Scenario
              )) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 *
                        SE), width = 0.2) +
    geom_hline(yintercept = est_theta,
               color = "red",
               linetype = "dashed") +
    labs(title = "Comparison of Method Estimates with SD Bands", y = "Estimate", x = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Scenario, ncol = 2)  # Creates 2x2 grid of plots
  
  # Print the plot
  print(p)
}


library(ggplot2)

# Define a function to generalize the analysis, plot all methods and scenarios,
# and add a bar plot for the power corresponding to the true value.
analyze_all_scenarios_coverage <- function(est_theta,
                                           thetaU,
                                           N,
                                           prop_invalid,
                                           overlap,
                                           sctfc = FALSE) {
  # Define the specific scenario filenames
  scenarios <- c("BI", "DN", "BN", "DI")
  df_all <- data.frame()
  power_df <- data.frame()  # Data frame for power calculation
  
  for (scenario in scenarios) {
    # Construct the file path based on input parameters for each scenario
    file_path <- sprintf(
      "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      overlap
    )
    # Load the data
    load(file_path)
    est <- na.omit(est)
    # Extract estimates and standard deviations for all methods
    estimates <- apply(est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )], 2, mean)
    se <- apply(est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )], 2, sd)
    methods <- c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )
    
    # Calculate the percentage of times the true value falls within the 95% CI
    within_95 <- colMeans((est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )] -
      1.96 * est[, c(
        "MRMode_se",
        "Cont-Mix_se",
        "MR-Mix_se",
        "mode_new_phi1_se",
        "MR-PRESSO_se",
        "MR-cML_se",
        "IVW_se",
        "median_se",
        "Egger_se", "IB-PRESSO_se"
      )] <= est_theta) &
      (est[, c(
        "MRMode",
        "Cont-Mix",
        "MR-Mix",
        "mode_new_phi1",
        "MR-PRESSO",
        "MR-cML",
        "IVW",
        "median",
        "Egger", "IB-PRESSO"
      )] +
        1.96 * est[, c(
          "MRMode_se",
          "Cont-Mix_se",
          "MR-Mix_se",
          "mode_new_phi1_se",
          "MR-PRESSO_se",
          "MR-cML_se",
          "IVW_se",
          "median_se",
          "Egger_se", "IB-PRESSO_se"
        )] >= est_theta)) * 100
    
    # Create data frame for power (95% CI coverage) for each method and scenario
    power_df_scenario <- data.frame(Method = methods,
                                    Power = within_95,
                                    Scenario = scenario)
    power_df <- rbind(power_df, power_df_scenario)
    
  }
  
  # Create a bar plot for power (percentage of times true value is within the 95% CI)
  p_power <- ggplot(power_df, aes(x = Method, y = Power, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap( ~ Scenario, ncol = 2) +
    labs(title = "Coverage: True Value Within 95% CI", y = "Coverage (%)", x = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plots
  print(p_power)
}

analyze_all_scenarios_error <- function(est_theta,
                                        thetaU,
                                        N,
                                        prop_invalid,
                                        overlap,
                                        sctfc = TRUE) {
  library(ggplot2)
  library(dplyr)
  # Define the specific scenario filenames
  scenarios <- c("BI", "DN", "BN", "DI")
  df_all <- data.frame()
  error_df <- data.frame()  # Data frame for errors
  
  for (scenario in scenarios) {
    # Construct the file path based on input parameters for each scenario
    file_path <- sprintf(
      "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      overlap
    )
    # Load the data
    load(file_path)
    est <- na.omit(est)
    # Extract estimates and standard deviations for all methods
    estimates <- apply(est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )], 2, mean)
    se <- apply(est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )], 2, sd)
    methods <- c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )
    
    # Calculate Type 1 Error or Power
    if (est_theta == 0) {
      # Type 1 Error: Rejecting null hypothesis when it's true
      error_rate <- colMeans(abs(est[, methods]) >
                               1.96 * est[, paste0(methods, "_se")]) * 100
      error_label <- "Type 1 Error"
    } else {
      # Type 2 Error and calculate Power (1 - Type 2 Error)
      type2_error_rate <- colMeans(abs(est[, methods]) <=
                                     1.96 * est[, paste0(methods, "_se")]) * 100
      error_rate <- 100 - type2_error_rate
      error_label <- "Power (1 - Type 2 Error)"
    }
    
    # Create data frame for errors
    error_df_scenario <- data.frame(Method = methods,
                                    ErrorRate = error_rate,
                                    Scenario = scenario)
    error_df <- rbind(error_df, error_df_scenario)
  }
  
  # Create a bar plot for Type 1 Error or Power
  p_error <- ggplot(error_df, aes(x = Method, y = ErrorRate, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap( ~ Scenario, ncol = 2) +
    labs(
      title = paste(error_label, ": Across Methods and Scenarios"),
      y = paste(error_label, "Rate (%)"),
      x = "Method"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_error)
}

analyze_all_scenarios_mse <- function(est_theta,
                                      thetaU,
                                      N,
                                      prop_invalid,
                                      overlap,
                                      sctfc = FALSE) {
  # Define the specific scenario filenames
  scenarios <- c("BI", "DN", "BN", "DI")
  mse_df <- data.frame()    # Data frame for MSE calculation
  
  for (scenario in scenarios) {
    # Construct the file path based on input parameters for each scenario
    file_path <- sprintf(
      "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      overlap
    )
    # Load the data
    load(file_path)
    est <- na.omit(est)
    # Extract estimates and standard deviations for all methods
    methods <- c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )
    
    # Calculate MSE for each method
    mse <- colMeans((est[, c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )] - est_theta) ^ 2)
    
    # Add MSE to the dataframe
    mse_df_scenario <- data.frame(Method = methods,
                                  MSE = mse,
                                  Scenario = scenario)
    mse_df <- rbind(mse_df, mse_df_scenario)
  }
  
  p_mse <- ggplot(mse_df, aes(x = Method, y = MSE, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Scenario, ncol = 2) +
    labs(title = "MSE Across Methods and Scenarios", y = "Mean Squared Error", x = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Print the MSE plot
  print(p_mse)
}

analyze_all_scenarios_time <- function(est_theta,
                                       thetaU,
                                       N,
                                       prop_invalid,
                                       overlap,
                                       sctfc = FALSE) {
  # Define the specific scenario filenames
  scenarios <- c("BI", "DN", "BN", "DI")
  time_df <- data.frame()  # Data frame for time calculation
  
  for (scenario in scenarios) {
    # Construct the file path based on input parameters for each scenario
    file_path <- sprintf(
      "../../results/simulation_results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
      scenario,
      est_theta,
      thetaU,
      format(N, scientific = sctfc),
      prop_invalid,
      overlap
    )
    # Load the data
    load(file_path)
    est <- na.omit(est)
    
    # Extract time taken for all methods
    methods <- c(
      "MRMode",
      "Cont-Mix",
      "MR-Mix",
      "mode_new_phi1",
      "MR-PRESSO",
      "MR-cML",
      "IVW",
      "median",
      "Egger", "IB-PRESSO"
    )
    
    # Extract time columns dynamically
    time_columns <- paste0(methods, "_time")
    time_values <- colMeans(est[, time_columns], na.rm = TRUE)
    
    # Create dataframe for this scenario
    time_df_scenario <- data.frame(Method = methods,
                                   Time = time_values,
                                   Scenario = scenario)
    
    # Combine results
    time_df <- rbind(time_df, time_df_scenario)
  }
  
  # Compute scenario-wise averages
  time_df_avg <- time_df %>%
    group_by(Method) %>%
    summarize(Average_Time = mean(Time, na.rm = TRUE)) %>%
    ungroup()
  
  # Plot the averaged time across scenarios
  p_time <- ggplot(time_df_avg, aes(x = Method, y = Average_Time, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Average Time Taken Across Methods", y = "Average Time (seconds)", x = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(p_time)
}


# check_se <- function(est_theta,
#                      thetaU,
#                      N,
#                      prop_invalid,
#                      overlap,
#                      sctfc = FALSE) {
#   # Define the specific scenario filenames
#   scenarios <- c("BI", "DN", "BN", "DI")
#   df_all <- data.frame()
#   df_all_summary <- data.frame()
#   
#   for (scenario in scenarios) {
#     # Construct the file path based on input parameters for each scenario
#     file_path <- sprintf(
#       "./results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#       scenario,
#       est_theta,
#       thetaU,
#       format(N, scientific = sctfc),
#       prop_invalid,
#       overlap
#     )
#     # Load the data
#     load(file_path)
#     est <- na.omit(est)
#     
#     # Calculate estimated SE using the standard deviation across iterations
#     estimated_se <- apply(est[, c(
#       "MRMode",
#       "Cont-Mix",
#       "MR-Mix",
#       "mode_new_phi1",
#       "MR-PRESSO",
#       "MR-cML",
#       "IVW",
#       "median",
#       "Egger", "IB-PRESSO"
#     )], 2, sd)
#     
#     # Extract built-in SE values across iterations to create an empirical distribution
#     builtin_se_values <- est[, c(
#       "MRMode_se",
#       "Cont-Mix_se",
#       "MR-Mix_se",
#       "mode_new_phi1_se",
#       "MR-PRESSO_se",
#       "MR-cML_se",
#       "IVW_se",
#       "median_se",
#       "Egger_se", "IB-PRESSO_se"
#     )]
#     
#     # Calculate the mean of built-in SEs
#     builtin_se_means <- apply(builtin_se_values, 2, median, na.rm = TRUE)
# #colMeans(builtin_se_values, na.rm = TRUE)
#     
#     # Calculate empirical p-values for estimated SEs
#     p_values <- sapply(1:ncol(builtin_se_values), function(i) {
#       mean(builtin_se_values[, i] >= estimated_se[i])
#     })
#     
#     methods <- c(
#       "MRMode",
#       "Cont-Mix",
#       "MR-Mix",
#       "mode_new_phi1",
#       "MR-PRESSO",
#       "MR-cML",
#       "IVW",
#       "median",
#       "Egger", "IB-PRESSO"
#     )
#     
#     # Create a data frame for the current scenario with full SE distribution
#     df_scenario <- data.frame(
#       Method = rep(methods, each = nrow(builtin_se_values)),
#       Builtin_SE = as.vector(builtin_se_values),
#       Scenario = scenario
#     )
#     
#     # Create a summary data frame for mean built-in SEs, estimated SEs, and p-values
#     df_scenario_summary <- data.frame(
#       Method = methods,
#       Mean_Builtin_SE = builtin_se_means,
#       Estimated_SE = estimated_se,
#       Empirical_P = p_values,
#       Scenario = scenario
#     )
#     
#     # Combine the data into the main data frames
#     df_all <- rbind(df_all, df_scenario)
#     df_all_summary <- rbind(df_all_summary, df_scenario_summary)
#   }
#   
#   # Plot the comparison of SEs using box plots and overlaying mean and estimated SE points
#   p <- ggplot(df_all, aes(x = Method, y = Builtin_SE, fill = Scenario)) +
#     geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#     geom_point(
#       data = df_all_summary,
#       aes(x = Method, y = Mean_Builtin_SE, color = "Mean Builtin SE"),
#       size = 3,
#       shape = 16
#     ) +
#     geom_point(
#       data = df_all_summary,
#       aes(x = Method, y = Estimated_SE, color = "Estimated SE"),
#       size = 3,
#       shape = 18
#     ) +
#     labs(title = "Comparison of Built-in SE Distribution and Estimated SE", y = "Standard Error (SE)", x = "Method") +
#     scale_color_manual(
#       name = "SE Type",
#       values = c(
#         "Mean Builtin SE" = "blue",
#         "Estimated SE" = "red"
#       )
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     facet_wrap( ~ Scenario, ncol = 2)  # Creates a 2x2 grid of plots
#   
#   # Print the plot and the summary data frame with empirical p-values
#   print(p)
#   print(df_all_summary)
# }

####################################################################(FOR MR-Mix vs IB-Mix)
# check_se <- function(est_theta,
#                      thetaU,
#                      N,
#                      prop_invalid,
#                      overlap,
#                      sctfc = FALSE) {
#   # Define the specific scenario filenames
#   df_all <- data.frame()
#   df_all_summary <- data.frame()
# 
#     # Construct the file path based on input parameters for each scenario
#     # trial_est_theta0_thetaU0.3_N80000_prop_invalid0.5_overlap1.rda
#     file_path <- sprintf(
#       "./results/mode_comp/trial_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#       est_theta,
#       thetaU,
#       format(N, scientific = sctfc),
#       prop_invalid,
#       overlap
#     )
#     # Load the data
#     load(file_path)
#     est <- na.omit(est)
# 
#     # Calculate estimated SE using the standard deviation across iterations
#     estimated_se <- apply(est[, c(
#       "MR-Mix",
#       "IB-Mix")], 2, sd)
# 
#     # Extract built-in SE values across iterations to create an empirical distribution
#     builtin_se_values <- est[, c(
#       "MR-Mix_se",
#       "IB-Mix_se")]
# 
#     # Calculate the mean of built-in SEs
#     builtin_se_means <- apply(builtin_se_values, 2, median, na.rm = TRUE)
#     #colMeans(builtin_se_values, na.rm = TRUE)
# 
#     # Calculate empirical p-values for estimated SEs
#     p_values <- sapply(1:ncol(builtin_se_values), function(i) {
#       mean(builtin_se_values[, i] >= estimated_se[i])
#     })
# 
#     methods <- c(
#       "MR-Mix",
#       "IB-Mix")
# 
#     # Create a data frame for the current scenario with full SE distribution
#     df_scenario <- data.frame(
#       Method = rep(methods, each = nrow(builtin_se_values)),
#       Builtin_SE = as.vector(builtin_se_values)
#     )
# 
#     # Create a summary data frame for mean built-in SEs, estimated SEs, and p-values
#     df_scenario_summary <- data.frame(
#       Method = methods,
#       Mean_Builtin_SE = builtin_se_means,
#       Estimated_SE = estimated_se,
#       Empirical_P = p_values
#     )
# 
#     # Combine the data into the main data frames
#     df_all <- rbind(df_all, df_scenario)
#     df_all_summary <- rbind(df_all_summary, df_scenario_summary)
# 
#   # Plot the comparison of SEs using box plots and overlaying mean and estimated SE points
#   p <- ggplot(df_all, aes(x = Method, y = Builtin_SE)) +
#     geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#     geom_point(
#       data = df_all_summary,
#       aes(x = Method, y = Mean_Builtin_SE, color = "Mean Builtin SE"),
#       size = 3,
#       shape = 16
#     ) +
#     geom_point(
#       data = df_all_summary,
#       aes(x = Method, y = Estimated_SE, color = "Estimated SE"),
#       size = 3,
#       shape = 18
#     ) +
#     labs(title = "Comparison of Built-in SE Distribution and Estimated SE", y = "Standard Error (SE)", x = "Method") +
#     scale_color_manual(
#       name = "SE Type",
#       values = c(
#         "Mean Builtin SE" = "blue",
#         "Estimated SE" = "red"
#       )
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
#   # Print the plot and the summary data frame with empirical p-values
#   print(p)
#   print(df_all_summary)
# }
# 
# check_se(
#   est_theta = 0.2,
#   thetaU = 0.3,
#   N = 200000,
#   prop_invalid = 0.5,
#   overlap = 1,
#   sctfc = FALSE
# )

####################################################################

## Example usage:
# analyze_all_scenarios(
#   est_theta = 0.2,
#   thetaU = 0.5,
#   N = 1e+05,
#   prop_invalid = 0.1,
#   overlap = 0.7,
#   sctfc = TRUE
# )




# check_se <- function(est_theta,
#                      thetaU,
#                      N,
#                      prop_invalid,
#                      overlap,
#                      sctfc = FALSE) {
#   # Define the specific scenario filenames
#   scenarios <- c("BI", "DN", "BN", "DI")
#   df_all <- data.frame()
#   df_all_summary <- data.frame()
# 
#   for (scenario in scenarios) {
#     # Construct the file path based on input parameters for each scenario
#     file_path <- sprintf(
#       "./results/%s_est_theta%g_thetaU%g_N%s_prop_invalid%g_overlap%g.rda",
#       scenario,
#       est_theta,
#       thetaU,
#       format(N, scientific = sctfc),
#       prop_invalid,
#       overlap
#     )
#     # Load the data
#     load(file_path)
#     est <- na.omit(est)
# 
#     # Calculate estimated SE using the standard deviation across iterations
#     estimated_se <- apply(est[, c(
#       "MR-Mix","IB-Mix"
#     )], 2, sd)
# 
#     # Extract built-in SE values across iterations to create an empirical distribution
#     builtin_se_values <- est[, c(
#       "MR-Mix_se","IB-Mix_se"
#     )]
# 
#     # Calculate the mean of built-in SEs
#     builtin_se_means <- colMeans(builtin_se_values, na.rm = TRUE)
#     #colMeans(builtin_se_values, na.rm = TRUE); apply(builtin_se_values, 2, median, na.rm = TRUE)
# 
#     # Calculate empirical p-values for estimated SEs
#     p_values <- sapply(1:ncol(builtin_se_values), function(i) {
#       mean(builtin_se_values[, i] >= estimated_se[i])
#     })
# 
#     methods <- c(
#       "MR-Mix", "IB-Mix"
#     )
# 
#     # Create a data frame for the current scenario with full SE distribution
#     df_scenario <- data.frame(
#       Method = rep(methods, each = nrow(builtin_se_values)),
#       Builtin_SE = as.vector(builtin_se_values),
#       Scenario = scenario
#     )
# 
#     # Create a summary data frame for mean built-in SEs, estimated SEs, and p-values
#     df_scenario_summary <- data.frame(
#       Method = methods,
#       Mean_Builtin_SE = builtin_se_means,
#       Estimated_SE = estimated_se,
#       Empirical_P = p_values,
#       Scenario = scenario
#     )
# 
#     # Combine the data into the main data frames
#     df_all <- rbind(df_all, df_scenario)
#     df_all_summary <- rbind(df_all_summary, df_scenario_summary)
#   }
# 
#   # Plot the comparison of SEs using box plots and overlaying mean and estimated SE points
#   p <- ggplot(df_all, aes(x = Method, y = Builtin_SE, fill = Scenario)) +
#     geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#     geom_point(
#       data = df_all_summary,
#       aes(x = Method, y = Mean_Builtin_SE, color = "Mean Builtin SE"),
#       size = 3,
#       shape = 16
#     ) +
#     geom_point(
#       data = df_all_summary,
#       aes(x = Method, y = Estimated_SE, color = "Estimated SE"),
#       size = 3,
#       shape = 18
#     ) +
#     labs(title = "Comparison of Built-in SE Distribution and Estimated SE", y = "Standard Error (SE)", x = "Method") +
#     scale_color_manual(
#       name = "SE Type",
#       values = c(
#         "Mean Builtin SE" = "blue",
#         "Estimated SE" = "red"
#       )
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     facet_wrap( ~ Scenario, ncol = 2)  # Creates a 2x2 grid of plots
# 
#   # Print the plot and the summary data frame with empirical p-values
#   print(p)
#   print(df_all_summary)
# }


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# Install required packages if necessary
packages <- c("ggplot2", "cowplot", "MASS")
# Uncomment next line if you need to install them
# sapply(packages, function(pkg) if (!require(pkg, character.only = TRUE)) install.packages(pkg))
library(ggplot2)
library(cowplot)
library(MASS)

set.seed(42)

# 1. Generate clustered data
n1 <- 400; n2 <- 200; n3 <- 200

# Large cluster around (0,0)
clust1 <- mvrnorm(n1, mu = c(0,0), Sigma = matrix(c(0.08, 0, 0, 0.08), 2))
# Small cluster (1,0)
clust2 <- mvrnorm(n2, mu = c(1,0), Sigma = matrix(c(0.06, 0, 0, 0.06), 2))
# Small cluster (1,1)
clust3 <- mvrnorm(n3, mu = c(1,1), Sigma = matrix(c(0.06, 0, 0, 0.06), 2))

dat <- as.data.frame(rbind(clust1, clust2, clust3))
names(dat) <- c("x","y")
dat$cluster <- factor(c(rep("Main (0,0)", n1), rep("(1,0)", n2), rep("(1,1)", n3)))

# 2. 2D scatter with density contours
scatter <- ggplot(dat, aes(x=x, y=y, color=cluster)) +
  geom_point(size=1.3, alpha=0.6, show.legend=FALSE) +
  stat_density2d(aes(fill = ..level..), geom = "polygon", color=NA, alpha=0.14, show.legend=FALSE) +
  scale_color_manual(values = c("#4257b2", "#ff6e54", "#38b96d")) +
  scale_fill_gradient(low = "#dddddd", high = "#7e57c2") +
  labs(title = "2D Mode Detection: True clusters") +
  xlim(-0.5, 1.6) + ylim(-0.5, 1.6) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(1,1,1,1)
  )

# 3. x-marginal: distribution inverted below
marginal <- ggplot(dat, aes(x=x)) +
  geom_histogram(aes(y=-..density..), bins=30, fill="#4257b2", alpha=0.5, color=NA) +
  geom_density(aes(y=-..density..), color="#d81b60", size=1) +
  labs(x="x", y="Marginal density") +
  xlim(-0.5, 1.6) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.y = element_text(vjust=0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face="bold"),
    plot.margin = margin(0,1,1,1)
  )

# 4. Compose: scatter above, x-marginal below
final_plot <- plot_grid(
  scatter, marginal, ncol=1, align="v", rel_heights = c(1.2, 0.75)
)

# # Add combined caption using cowplot::ggdraw
# library(grid)
# plot_with_caption <- ggdraw(final_plot) +
#   draw_label(
#     "Projected (inverted) x-axis marginal suggests bimodality,\nbut 2D clustering finds correct modes.",
#     x = 0.5, y = 0.08, hjust = 0.5, fontface = "italic", size=12
#   )

# Save as high-res PNG/PDf
ggsave("2D_mode_vs_marginal.png", final_plot, width=6, height=7, dpi=300)
ggsave("2D_mode_vs_marginal.pdf", final_plot, width=6, height=7)

# Print to R display
print(plot_with_caption)