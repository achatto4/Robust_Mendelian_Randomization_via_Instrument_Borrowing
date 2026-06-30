#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

if (!requireNamespace("xtable", quietly = TRUE)) {
  stop("Package 'xtable' is required. Install with: install.packages('xtable')")
}

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "../../../results/simulation_results/mode_comp"
out_dir <- if (length(args) >= 2) args[2] else file.path(results_dir, "NULL_reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading NULL summaries from: ", normalizePath(results_dir, mustWork = FALSE))

summary_files <- list.files(
  path = results_dir,
  pattern = "_NULL_.*_type1_summary\\.csv$",
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No files matched '*_NULL_*_type1_summary.csv' in: ", results_dir)
}

raw <- bind_rows(lapply(summary_files, readr::read_csv, show_col_types = FALSE))

required_cols <- c(
  "scenario", "theta", "theta_alt", "phi", "alpha", "n_rep",
  "type1_mrmode", "type1_ib_mode", "thetaU", "N", "prop_invalid", "overlap"
)
missing <- setdiff(required_cols, names(raw))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

n_before <- nrow(raw)

dat <- raw %>%
  mutate(
    scenario = factor(scenario, levels = c("BI", "BN", "DI", "DN")),
    theta_alt = as.numeric(theta_alt),
    alpha = as.numeric(alpha),
    N = as.numeric(N),
    prop_invalid = as.numeric(prop_invalid),
    overlap = as.numeric(overlap),
    type1_mrmode = as.numeric(type1_mrmode),
    type1_ib_mode = as.numeric(type1_ib_mode)
  ) %>%
  filter(
    !is.na(scenario), !is.na(alpha), !is.na(theta_alt), !is.na(N),
    !is.na(prop_invalid), !is.na(overlap),
    !is.na(type1_mrmode), !is.na(type1_ib_mode),
    is.finite(alpha), is.finite(theta_alt), is.finite(N),
    is.finite(prop_invalid), is.finite(overlap),
    is.finite(type1_mrmode), is.finite(type1_ib_mode)
  ) %>%
  mutate(
    alpha_pct = 100 * alpha,
    mrmode_pct = 100 * type1_mrmode,
    ibmode_pct = 100 * type1_ib_mode,
    ib_dev_pp = ibmode_pct - alpha_pct,
    mr_dev_pp = mrmode_pct - alpha_pct
  ) %>%
  arrange(scenario, theta_alt, alpha, N, prop_invalid, overlap)

n_after <- nrow(dat)
message("Rows loaded: ", n_before, " | rows used after filtering: ", n_after, " | rows dropped: ", n_before - n_after)

if (n_after == 0) {
  stop("All rows were filtered out; no usable data for tables/plots.")
}

long_dat <- dat %>%
  transmute(
    scenario, theta_alt, alpha, alpha_pct, N, prop_invalid, overlap,
    MRMode = mrmode_pct,
    `IB-Mode` = ibmode_pct
  ) %>%
  pivot_longer(cols = c("MRMode", "IB-Mode"), names_to = "method", values_to = "type1_pct")

# -----------------------------
# Concise LaTeX tables
# -----------------------------

# Table 1: headline calibration, pooled across all settings.
calib_tbl <- long_dat %>%
  group_by(method, alpha) %>%
  summarise(
    nominal_pct = 100 * first(alpha),
    observed_mean_pct = mean(type1_pct),
    deviation_pp = observed_mean_pct - nominal_pct,
    q10_pct = quantile(type1_pct, 0.10),
    q90_pct = quantile(type1_pct, 0.90),
    settings = n(),
    .groups = "drop"
  ) %>%
  arrange(alpha, method)

# Table 2: IB-Mode robustness to non-zero secondary effect, pooled across scenario/N/invalid/overlap.
ib_theta_tbl <- dat %>%
  group_by(alpha, theta_alt) %>%
  summarise(
    nominal_pct = 100 * first(alpha),
    observed_mean_pct = mean(ibmode_pct),
    deviation_pp = observed_mean_pct - nominal_pct,
    q10_pct = quantile(ibmode_pct, 0.10),
    q90_pct = quantile(ibmode_pct, 0.90),
    settings = n(),
    .groups = "drop"
  ) %>%
  arrange(alpha, theta_alt)

write_xtable <- function(df, file, caption, label, digits = 2) {
  tab <- xtable::xtable(df, caption = caption, label = label, digits = digits)
  print(
    tab,
    file = file,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    caption.placement = "top",
    table.placement = "htbp"
  )
}

write_xtable(
  calib_tbl,
  file.path(out_dir, "NULL_table_main_calibration.tex"),
  caption = "Observed type-I error (percent) by nominal significance level (percent), pooled across all settings.",
  label = "tab:null_main_calibration"
)

write_xtable(
  ib_theta_tbl,
  file.path(out_dir, "NULL_table_ibmode_theta_alt.tex"),
  caption = "IB-Mode type-I error (percent) versus secondary trait effect $\\theta_{alt}$.",
  label = "tab:null_ibmode_theta_alt"
)

# -----------------------------
# Intuitive, minimal figures
# -----------------------------

theme_pub <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Figure 1: main calibration (both methods) with spread bands.
fig1_df <- long_dat %>%
  group_by(method, alpha) %>%
  summarise(
    nominal_pct = 100 * first(alpha),
    mean_pct = mean(type1_pct),
    q10_pct = quantile(type1_pct, 0.10),
    q90_pct = quantile(type1_pct, 0.90),
    .groups = "drop"
  )

fig1 <- ggplot(fig1_df, aes(x = nominal_pct, y = mean_pct, color = method, group = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.6) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = q10_pct, ymax = q90_pct), width = 0.7, alpha = 0.7) +
  scale_x_continuous(breaks = c(20, 10, 5)) +
  labs(
    title = "Type-I Error Calibration at phi = 1, theta = 0",
    subtitle = "Error bars show 10th-90th percentile across all simulation settings",
    x = "Nominal significance level (%)",
    y = "Observed type-I error (%)",
    color = NULL
  ) +
  theme_pub

# Figure 2: scenario-wise IB-Mode calibration.
fig2_df <- dat %>%
  group_by(scenario, alpha) %>%
  summarise(
    nominal_pct = first(alpha_pct),
    mean_pct = mean(ibmode_pct),
    q10_pct = quantile(ibmode_pct, 0.10),
    q90_pct = quantile(ibmode_pct, 0.90),
    .groups = "drop"
  )

fig2 <- ggplot(fig2_df, aes(x = nominal_pct, y = mean_pct, group = scenario, color = scenario)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.6) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = q10_pct, ymax = q90_pct), width = 0.7, alpha = 0.7) +
  scale_x_continuous(breaks = c(20, 10, 5)) +
  labs(
    title = "IB-Mode Calibration Across Scenarios",
    x = "Nominal significance level (%)",
    y = "Observed type-I error (%)",
    color = "Scenario"
  ) +
  theme_pub

# Figure 3: sample size + invalid-IV effect via deviation heatmap.
fig3_df <- dat %>%
  group_by(alpha, N, prop_invalid) %>%
  summarise(
    mean_dev_pp = mean(ib_dev_pp),
    .groups = "drop"
  )

fig3 <- ggplot(fig3_df, aes(x = factor(N), y = factor(prop_invalid), fill = mean_dev_pp)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", mean_dev_pp)), size = 3.1) +
  facet_wrap(~ alpha, nrow = 1, labeller = label_both) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    name = "Deviation\n(pp)"
  ) +
  labs(
    title = "IB-Mode Deviation (Observed - Nominal)",
    subtitle = "Averaged over scenario, overlap, and secondary-effect settings",
    x = "Sample size (N)",
    y = "Invalid IV proportion"
  ) +
  theme_pub

# Figure 4: robustness to non-zero secondary-trait effect.
# For visualization, focus on calibration-relevant alpha levels and include 0.1 when available.
preferred_alpha <- c(0.005, 0.01, 0.05, 0.1)
plot_alpha <- preferred_alpha[preferred_alpha %in% sort(unique(dat$alpha))]
if (length(plot_alpha) == 0) {
  plot_alpha <- sort(unique(dat$alpha[dat$alpha < 1]))
}
if (length(plot_alpha) == 0) {
  plot_alpha <- sort(unique(dat$alpha))
}

fig4_df <- dat %>%
  filter(alpha %in% plot_alpha) %>%
  group_by(alpha, theta_alt) %>%
  summarise(
    nominal_pct = first(alpha_pct),
    mean_pct = mean(ibmode_pct),
    q10_pct = quantile(ibmode_pct, 0.10),
    q90_pct = quantile(ibmode_pct, 0.90),
    .groups = "drop"
  )

hline_df <- distinct(fig4_df, alpha, nominal_pct)

fig4 <- ggplot(fig4_df, aes(x = theta_alt, y = mean_pct, color = factor(alpha), group = alpha)) +
  geom_hline(
    data = hline_df,
    aes(yintercept = nominal_pct, color = factor(alpha)),
    linetype = "dotted",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = q10_pct, ymax = q90_pct), width = 0.02, alpha = 0.7) +
  labs(
    title = "IB-Mode Robustness to Non-zero Secondary-Trait Effect",
    subtitle = "Dotted lines mark nominal alpha levels",
    x = expression(theta[alt]),
    y = "Observed type-I error (%)",
    color = expression(alpha)
  ) +
  theme_pub

# Figure 5: scatter-grid style calibration plot (similar to requested style), no NA rows.
# --- ORIGINAL VERSION (all four scenarios BI/BN/DI/DN, all overlaps) -- COMMENTED OUT per reviewer (F2) ---
# fig5 <- dat %>%
#   filter(alpha %in% plot_alpha) %>%
#   mutate(type1_ib_mode_plot = pmax(type1_ib_mode, min(plot_alpha) / 10)) %>%
#   ggplot(aes(x = alpha, y = type1_ib_mode_plot)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
#   geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
#   geom_point(
#     aes(color = factor(N), shape = factor(prop_invalid)),
#     alpha = 0.65,
#     size = 2.0,
#     position = position_jitter(width = 0.0004, height = 0.0)
#   ) +
#   facet_grid(theta_alt ~ scenario, labeller = label_both) +
#   scale_x_continuous(
#     trans = "log10",
#     breaks = plot_alpha,
#     labels = scales::label_number(accuracy = 0.001)
#   ) +
#   scale_y_continuous(
#     trans = "log10",
#     breaks = plot_alpha,
#     labels = scales::label_number(accuracy = 0.001)
#   ) +
#   labs(
#     title = "IB-Mode Type-I Error Calibration (phi = 1, theta = 0)",
#     subtitle = "Panels by scenario and secondary-trait effect; both axes use log scale",
#     x = "Nominal significance level",
#     y = "Observed type-I error",
#     color = "Sample size (N)",
#     shape = "Invalid IV proportion"
#   ) +
#   theme_pub

# --- REVISED VERSION (F2): restrict to InSIDE-violated scenarios only (BN, DN) at D_ov = 0.75,
#     i.e. the same scenario family as main-text Figure 3. Scenario codes:
#     B=balanced / D=directional pleiotropy ; I=InSIDE satisfied / N=noInside (InSIDE violated).
fig5 <- dat %>%
  filter(alpha %in% plot_alpha) %>%
  filter(scenario %in% c("BN", "DN")) %>%   # keep only InSIDE-violated scenarios
  filter(abs(overlap - 0.75) < 1e-8) %>%    # degree of invalid-instrument overlap D_ov = 0.75
  mutate(type1_ib_mode_plot = pmax(type1_ib_mode, min(plot_alpha) / 10)) %>%
  ggplot(aes(x = alpha, y = type1_ib_mode_plot)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  geom_point(
    aes(color = factor(N), shape = factor(prop_invalid)),
    alpha = 0.65,
    size = 2.0,
    position = position_jitter(width = 0.0004, height = 0.0)
  ) +
  facet_grid(theta_alt ~ scenario, labeller = label_both) +
  scale_x_continuous(
    trans = "log10",
    breaks = plot_alpha,
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = plot_alpha,
    labels = scales::label_number(accuracy = 0.001)
  ) +
  labs(
    title = "IB-Mode Type-I Error Calibration (phi = 1, theta = 0)",
    subtitle = "InSIDE-violated scenarios (BN, DN), invalid-instrument overlap D_ov = 0.75; both axes log scale",
    x = "Nominal significance level",
    y = "Observed type-I error",
    color = "Sample size (N)",
    shape = "Invalid IV proportion"
  ) +
  theme_pub

ggsave(file.path(out_dir, "NULL_fig1_calibration_methods.png"), fig1, width = 8.5, height = 5.2, dpi = 320)
ggsave(file.path(out_dir, "NULL_fig2_calibration_scenarios.png"), fig2, width = 8.5, height = 5.2, dpi = 320)
ggsave(file.path(out_dir, "NULL_fig3_deviation_heatmap_N_invalid.png"), fig3, width = 10.0, height = 4.8, dpi = 320)
ggsave(file.path(out_dir, "NULL_fig4_robustness_theta_alt.png"), fig4, width = 8.5, height = 5.2, dpi = 320)
ggsave(file.path(out_dir, "NULL_fig5_scatter_grid_calibration.png"), fig5, width = 13.5, height = 9.5, dpi = 320)

# Figure 3-style robustness panel by scenario and alpha.
fig_theta_scenario_df <- dat %>%
  filter(alpha %in% plot_alpha) %>%
  group_by(scenario, theta_alt, alpha) %>%
  summarise(
    nominal_pct = first(alpha_pct),
    mean_pct = pmax(mean(ibmode_pct), 0.001),
    q10_pct = pmax(quantile(ibmode_pct, 0.10), 0.001),
    q90_pct = pmax(quantile(ibmode_pct, 0.90), 0.001),
    .groups = "drop"
  )

fig_theta_scenario <- ggplot(
  fig_theta_scenario_df,
  aes(x = theta_alt, y = mean_pct, color = factor(alpha), group = alpha)
) +
  geom_hline(
    aes(yintercept = nominal_pct, color = factor(alpha)),
    linetype = "dotted",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.0) +
  geom_errorbar(aes(ymin = q10_pct, ymax = q90_pct), width = 0.02, alpha = 0.75) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_x_continuous(breaks = sort(unique(fig_theta_scenario_df$theta_alt))) +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "IB-Mode Remains Calibrated with Non-zero Secondary-Trait Effect",
    subtitle = "Dotted lines are nominal alpha levels; y-axis on log scale",
    x = expression(theta[alt]),
    y = "Observed type-I error (%)",
    color = expression(alpha)
  ) +
  theme_pub

ggsave(file.path(out_dir, "NULL_fig3_ibmode_thetaAlt_robustness.png"), fig_theta_scenario, width = 10.8, height = 5.4, dpi = 320)

# Figure 6: overlap sensitivity curves (IB-Mode only).
fig6_df <- dat %>%
  group_by(scenario, alpha, overlap) %>%
  summarise(
    mean_pct = mean(ibmode_pct),
    q10_pct = quantile(ibmode_pct, 0.10),
    q90_pct = quantile(ibmode_pct, 0.90),
    nominal_pct = first(alpha_pct),
    .groups = "drop"
  )

fig6 <- ggplot(fig6_df, aes(x = overlap, y = mean_pct, color = factor(alpha), group = alpha)) +
  geom_hline(
    aes(yintercept = nominal_pct, color = factor(alpha)),
    linetype = "dotted",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(ymin = q10_pct, ymax = q90_pct, fill = factor(alpha)),
    alpha = 0.12,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.9) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_x_continuous(breaks = sort(unique(fig6_df$overlap))) +
  labs(
    title = "Overlap Sensitivity Curves (IB-Mode)",
    subtitle = "Mean observed type-I error vs overlap; shaded bands show 10th-90th percentile",
    x = "Sample overlap",
    y = "Observed type-I error (%)",
    color = expression(alpha)
  ) +
  theme_pub

ggsave(file.path(out_dir, "NULL_fig6_overlap_sensitivity_curves.png"), fig6, width = 10.5, height = 4.8, dpi = 320)

message("Done. Wrote concise LaTeX tables and high-signal figures to: ", normalizePath(out_dir, mustWork = FALSE))
