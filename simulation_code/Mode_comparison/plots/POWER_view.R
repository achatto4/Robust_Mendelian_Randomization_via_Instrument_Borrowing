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
results_dir <- if (length(args) >= 1) args[1] else "../../../results/simulation_results/mode_comp_power"
out_dir <- if (length(args) >= 2) args[2] else file.path(results_dir, "POWER_reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading POWER summaries from: ", normalizePath(results_dir, mustWork = FALSE))

summary_files <- list.files(
  path = results_dir,
  pattern = "_POWER_.*_power_summary\\.csv$",
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No files matched '*_POWER_*_power_summary.csv' in: ", results_dir)
}

raw <- bind_rows(lapply(summary_files, readr::read_csv, show_col_types = FALSE))

required_cols <- c(
  "scenario", "theta", "theta_alt", "phi", "alpha", "n_rep",
  "power_mrmode", "power_ib_mode", "thetaU", "N", "prop_invalid",
  "overlap", "NxNy_alt_ratio"
)
missing <- setdiff(required_cols, names(raw))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

n_before <- nrow(raw)

dat <- raw %>%
  mutate(
    scenario = factor(scenario, levels = c("BI", "BN", "DI", "DN")),
    theta = as.numeric(theta),
    theta_alt = as.numeric(theta_alt),
    alpha = as.numeric(alpha),
    N = as.numeric(N),
    prop_invalid = as.numeric(prop_invalid),
    overlap = as.numeric(overlap),
    NxNy_alt_ratio = as.numeric(NxNy_alt_ratio),
    power_mrmode = as.numeric(power_mrmode),
    power_ib_mode = as.numeric(power_ib_mode)
  ) %>%
  filter(
    !is.na(scenario), !is.na(theta), !is.na(theta_alt), !is.na(alpha),
    !is.na(N), !is.na(prop_invalid), !is.na(overlap), !is.na(NxNy_alt_ratio),
    !is.na(power_mrmode), !is.na(power_ib_mode),
    is.finite(theta), is.finite(theta_alt), is.finite(alpha),
    is.finite(N), is.finite(prop_invalid), is.finite(overlap),
    is.finite(NxNy_alt_ratio), is.finite(power_mrmode), is.finite(power_ib_mode)
  ) %>%
  mutate(
    power_mrmode_pct = 100 * power_mrmode,
    power_ibmode_pct = 100 * power_ib_mode
  ) %>%
  arrange(scenario, theta, theta_alt, alpha, N, prop_invalid, overlap, NxNy_alt_ratio)

n_after <- nrow(dat)
message("Rows loaded: ", n_before, " | rows used after filtering: ", n_after, " | rows dropped: ", n_before - n_after)
if (n_after == 0) stop("All rows were filtered out; no usable data for tables/plots.")

# Preferred alpha set for plots
preferred_alpha <- c(0.005, 0.01, 0.05, 0.1)
plot_alpha <- preferred_alpha[preferred_alpha %in% sort(unique(dat$alpha))]
if (length(plot_alpha) == 0) plot_alpha <- sort(unique(dat$alpha))

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

# -----------------------------
# Tables
# -----------------------------

# Table 1: Power vs auxiliary sample size (NxNy_alt_ratio) for IB-Mode
power_aux_tbl <- dat %>%
  group_by(theta, alpha, NxNy_alt_ratio) %>%
  summarise(
    mean_power_pct = mean(power_ibmode_pct),
    q10_pct = quantile(power_ibmode_pct, 0.10),
    q90_pct = quantile(power_ibmode_pct, 0.90),
    settings = n(),
    .groups = "drop"
  ) %>%
  arrange(theta, alpha, NxNy_alt_ratio)

write_xtable(
  power_aux_tbl,
  file.path(out_dir, "POWER_table_ibmode_auxN.tex"),
  caption = "IB-Mode power (percent) versus auxiliary outcome sample-size ratio $N_{alt}/N$ (summarized over scenarios; fixed $N=10^5$, invalid-IV proportion 0.5, overlap 0.75).",
  label = "tab:power_ibmode_auxN"
)

# -----------------------------
# Figures
# -----------------------------

theme_pub <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Figure 1: Power vs auxiliary sample size ratio (NxNy_alt_ratio)
fig1_df <- dat %>%
  filter(alpha %in% plot_alpha) %>%
  group_by(theta, scenario, alpha, NxNy_alt_ratio) %>%
  summarise(
    mean_power = mean(power_ibmode_pct),
    q10 = quantile(power_ibmode_pct, 0.10),
    q90 = quantile(power_ibmode_pct, 0.90),
    .groups = "drop"
  ) %>%
  mutate(Nalt_over_N = 1 / NxNy_alt_ratio)

fig1 <- ggplot(fig1_df, aes(x = Nalt_over_N, y = mean_power, color = factor(alpha), group = alpha)) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.9) +
  facet_grid(theta ~ scenario, labeller = label_both) +
  scale_x_continuous(
    trans = "log10",
    breaks = sort(unique(fig1_df$Nalt_over_N)),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    title = "IB-Mode Power vs Auxiliary Outcome Sample Size",
    subtitle = "Fixed N = 1e5, invalid-IV proportion = 0.5, overlap = 0.75; x-axis on log scale",
    x = expression(N[alt]/N),
    y = "Power (%)",
    color = expression(alpha)
  ) +
  theme_pub

# Figure 2: Power vs secondary-trait effect (theta_alt)
fig2_df <- dat %>%
  filter(alpha %in% plot_alpha) %>%
  group_by(theta, scenario, alpha, theta_alt) %>%
  summarise(
    mean_power = mean(power_ibmode_pct),
    q10 = quantile(power_ibmode_pct, 0.10),
    q90 = quantile(power_ibmode_pct, 0.90),
    .groups = "drop"
  )

fig2 <- ggplot(fig2_df, aes(x = theta_alt, y = mean_power, color = factor(alpha), group = alpha)) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.9) +
  facet_grid(theta ~ scenario, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(fig2_df$theta_alt))) +
  labs(
    title = "IB-Mode Power vs Secondary-Trait Effect",
    subtitle = "Fixed N = 1e5, invalid-IV proportion = 0.5, overlap = 0.75; averaged over auxiliary sample size",
    x = expression(theta[alt]),
    y = "Power (%)",
    color = expression(alpha)
  ) +
  theme_pub

# Figure 3: Heatmap of power vs primary effect and auxiliary sample size
fig3_df <- dat %>%
  filter(alpha %in% plot_alpha) %>%
  group_by(theta, alpha, NxNy_alt_ratio) %>%
  summarise(mean_power = mean(power_ibmode_pct), .groups = "drop")

fig3 <- ggplot(fig3_df, aes(x = factor(NxNy_alt_ratio), y = factor(theta), fill = mean_power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", mean_power)), size = 3.1) +
  facet_wrap(~ alpha, nrow = 1, labeller = label_both) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "Power (%)") +
  labs(
    title = "IB-Mode Power vs Auxiliary Sample Size",
    subtitle = "Averaged over scenario and secondary-trait effect; fixed N = 1e5, invalid-IV proportion = 0.5, overlap = 0.75",
    x = expression(N/N[alt]),
    y = expression(theta)
  ) +
  theme_pub

# Save figures

ggsave(file.path(out_dir, "POWER_fig1_power_vs_auxN.png"), fig1, width = 12.5, height = 6.5, dpi = 320)
ggsave(file.path(out_dir, "POWER_fig2_power_vs_thetaAlt.png"), fig2, width = 12.5, height = 6.5, dpi = 320)
ggsave(file.path(out_dir, "POWER_fig3_heatmap_power_auxN.png"), fig3, width = 10.5, height = 5.2, dpi = 320)

message("Done. Wrote POWER tables and figures to: ", normalizePath(out_dir, mustWork = FALSE))
