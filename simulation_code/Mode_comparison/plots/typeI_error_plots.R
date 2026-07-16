#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "../../../results/simulation_results/mode_comp"

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

# fig1 (type-I calibration) and fig4 (robustness to theta_alt) are the
# type-I figures used in the paper; built above as ggplot objects — export as needed.
