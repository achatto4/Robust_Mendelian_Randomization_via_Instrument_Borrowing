# ------------------ Setup ------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(viridis)
library(ggsci)

# Path to your CSV files
data_dir <- "RDA_results/"  # directory of final_results_*.csv produced by the exposure_*.R scripts

csv_files <- list.files(
  path = data_dir, pattern = "^final_results_.*\\.csv$",
  full.names = TRUE
)
if (length(csv_files) == 0) stop("No CSV files found in the specified directory.")

data_list <- setNames(
  lapply(csv_files, read.csv, stringsAsFactors = FALSE),
  str_remove(basename(csv_files), "\\.csv$")
)

# ------------------ Exposure Abbreviation Map ------------------
exposure_map <- c(
  "BFP" = "BF", "bmi" = "BMI", "bw" = "BW",
  "dbp" = "DBP", "drinkpw" = "Alcohol", "eversmok" = "Smoking",
  "FG" = "FG", "hdl" = "HDL", "height" = "Height",
  "ldl" = "LDL", "logTG" = "TG", "sbp" = "SBP"
)

# ------------------ Outcome Normalization Map ------------------
outcome_map <- c(
  "cad" = "CAD", "cad_mvp" = "CAD",
  "t2d" = "T2D", "t2d_mvp" = "T2D",
  "stroke" = "Stroke", "stk" = "Stroke", "stk_mvp" = "Stroke",
  "asthma" = "Asthma", "astm" = "Asthma", "astm_mvp" = "Asthma"
)

# ------------------ Causal Pair Categories ------------------
causal_pairs <- list(
  considered_causal = c(
    "SBP -> CAD", "DBP -> CAD", "Smoking -> CAD",
    "SBP -> Stroke", "DBP -> Stroke", "Smoking -> Stroke",
    "LDL -> CAD", "Smoking -> T2D", "LDL -> Stroke"
  ),
  supported_by_literature = c(
    "BMI -> CAD", "TG -> CAD", "BMI -> T2D",
    "BF -> CAD", "BF -> T2D", "FG -> T2D",
    "Height -> CAD", "BMI -> Stroke", "BF -> Stroke",
    "Smoking -> Asthma"
  ),
  unknown_conflict = c(
    "SBP -> T2D", "HDL -> T2D", "DBP -> T2D", "TG -> T2D",
    "BW -> CAD", "BW -> T2D", "TG -> Stroke", "BMI -> Asthma",
    "FG -> CAD", "FG -> Stroke", "BF -> Asthma", "BW -> Stroke",
    "Alcohol -> Stroke", "Height -> Stroke", "LDL -> T2D",
    "Alcohol -> T2D", "Alcohol -> CAD"
  ),
  implausible_unsupported = c(
    "Alcohol -> Asthma", "Height -> Asthma", "SBP -> Asthma",
    "DBP -> Asthma", "FG -> Asthma", "TG -> Asthma",
    "Height -> T2D", "LDL -> Asthma", "BW -> Asthma",
    "HDL -> Asthma"
  ),
  considered_noncausal = c(
    "HDL -> CAD", "HDL -> Stroke"
  )
)

pair_categories <- imap_dfr(causal_pairs, ~ {
  tibble(
    exposure = str_extract(.x, "^[^ ]+"),
    outcome  = str_extract(.x, "[^ ]+$"),
    category = .y
  )
}) %>%
  mutate(category = recode(category,
                           considered_causal         = "causal",
                           supported_by_literature   = "causal",
                           considered_noncausal      = "noncausal",
                           unknown_conflict          = "correlated",
                           implausible_unsupported   = "unsupported"
  ))

pair_lookup <- pair_categories %>% rename(dict_category = category)

# ------------------ Combine & Normalize All Data ------------------
# NOTE: We will create *two* significance flags later: p<0.001 and p<0.05
bonf_threshold <- 0.001
p05_threshold  <- 0.05

df_all <- map2_dfr(data_list, names(data_list), ~ {
  exposure_raw <- str_remove(.y, "^final_results_")
  if (tolower(exposure_raw) %in% c("crp", "vitd")) return(NULL)
  
  exposure_clean <- dplyr::recode(exposure_raw, !!!exposure_map, .default = exposure_raw)
  
  .x %>%
    mutate(
      exposure = exposure_clean,
      outcome = {
        tmp <- stringr::str_to_lower(trait1)
        tmp <- stringr::str_replace(tmp, "_mvp$", "")
        mapped <- dplyr::recode(tmp, !!!outcome_map, .default = tmp)
        ifelse(mapped %in% c("CAD", "T2D", "Stroke", "Asthma"), mapped, toupper(mapped))
      },
      file = .y,
      method = case_when(
        str_trim(method) == "new_method_weighted" ~ "IB-MODE",
        str_trim(method) == "IB-MR-PRESSO" ~ "IB-PRESSO",
        TRUE ~ str_trim(method)
      )
    )
}) %>%
  filter(!str_detect(method, "(?i)IB-Mix|Simple mode")) %>%
  mutate(
    pval = suppressWarnings(as.numeric(pval))
  )

# ------------------ Keep only dictionary pairs & attach category ------------------
df_all <- df_all %>%
  inner_join(pair_lookup, by = c("exposure", "outcome")) %>%
  mutate(category = dict_category) %>%
  dplyr::select(-dict_category)

# Required facet order: noncausal, unsupp, corr, causal (internally: noncausal, unsupported, correlated, causal)
cat_levels <- c("noncausal", "unsupported", "correlated", "causal")
cat_short_map <- c(
  noncausal   = "noncausal",
  unsupported = "unsupp",
  correlated  = "corr",
  causal      = "causal"
)

df_all <- df_all %>%
  mutate(category = factor(as.character(category), levels = cat_levels))

# ------------------ Keep ONLY first auxiliary trait (preserve file order) ------------------
df_all <- df_all %>%
  group_by(exposure, outcome, method) %>%
  slice(1) %>%
  ungroup()

# ------------------ Plot + table maker ------------------
make_plot_and_tables <- function(df, thr, out_eps, out_sig_csv, out_pval_csv, plot_title_suffix) {
  
  df2 <- df %>%
    mutate(sig = !is.na(pval) & pval < thr)
  
  method_category_totals <- df2 %>%
    distinct(exposure, outcome, method, category) %>%
    count(method, category, name = "total_n")
  
  summary_counts <- df2 %>%
    group_by(method, category) %>%
    summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop") %>%
    left_join(method_category_totals, by = c("method", "category")) %>%
    mutate(perc_sig = 100 * n_sig / total_n)
  
  category_labels <- method_category_totals %>%
    group_by(category) %>%
    summarise(total_n = max(total_n, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      category = factor(as.character(category), levels = cat_levels),
      facet_lab = paste0(recode(as.character(category), !!!cat_short_map), " (N=", total_n, ")")
    ) %>%
    arrange(category) %>%
    dplyr::select(category, facet_lab)
  
  summary_counts_plot <- summary_counts %>%
    left_join(category_labels, by = "category") %>%
    mutate(facet_lab = factor(facet_lab, levels = category_labels$facet_lab)) %>%
    filter(!is.na(perc_sig) & perc_sig >= 0 & perc_sig <= 100)
  
  # ---- method ordering (IB-MODE with Weighted mode; IB-MR-PRESSO with MR-PRESSO) ----
  all_methods <- unique(df2$method)
  mode_methods <- c("IB-MODE", "Weighted mode")
  presso_methods <- c("IB-MR-PRESSO", "MR-PRESSO")
  remaining_methods <- setdiff(sort(all_methods), c(mode_methods, presso_methods))
  desired_method_order <- c(mode_methods, presso_methods, remaining_methods)
  
  summary_counts_plot$method <- factor(summary_counts_plot$method, levels = desired_method_order)
  
  group_breaks <- c(
    length(mode_methods) + 0.5,
    length(mode_methods) + length(presso_methods) + 0.5
  )
  
  if (nrow(summary_counts_plot) > 0) {
    
    plot_sig <- ggplot(summary_counts_plot, aes(x = method, y = perc_sig, fill = method)) +
      geom_col(width = 0.7, color = "gray30") +
      geom_vline(xintercept = group_breaks, linetype = "dotted", linewidth = 1) +
      facet_wrap(~ facet_lab, scales = "free_x", nrow = 1) +
      scale_y_continuous(
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = c(0, 0.05)),
        limits = c(0, 100)
      ) +
      labs(
        x = NULL,
        y = "Percentage of Significant Pairs (%)",
        title = paste0("Percentage of Significant Pairs by Method & Category", plot_title_suffix),
        fill = "Method"
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.grid.major.x = element_blank()
      ) +
      scale_fill_d3(
        palette = "category20",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          nrow = 2,
          byrow = TRUE
        )
      )
    
    ggsave(out_eps, plot = plot_sig, device = cairo_ps, width = 10, height = 5, units = "in")
    print(plot_sig)
    
  } else {
    message("No dictionary pairs found after normalization — plot skipped for threshold: ", thr)
    plot_sig <- NULL
  }
  
  # ---- wide tables for this threshold ----
  sig_table <- df2 %>%
    group_by(exposure, outcome, method, category) %>%
    summarise(
      any_sig = any(sig, na.rm = TRUE),
      min_p   = suppressWarnings(min(pval, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      min_p = ifelse(is.infinite(min_p), NA_real_, min_p),
      category = factor(as.character(category), levels = cat_levels)
    )
  
  sig_wide <- sig_table %>%
    mutate(`Significant?` = if_else(any_sig, "✓", "")) %>%
    dplyr::select(exposure, outcome, category, method, `Significant?`) %>%
    pivot_wider(names_from = method, values_from = `Significant?`) %>%
    arrange(category, outcome, exposure)
  
  pval_wide <- sig_table %>%
    dplyr::select(exposure, outcome, category, method, min_p) %>%
    pivot_wider(names_from = method, values_from = min_p) %>%
    arrange(category, outcome, exposure)
  
  write.csv(sig_wide, out_sig_csv, row.names = FALSE)
  write.csv(pval_wide, out_pval_csv, row.names = FALSE)
  
  invisible(list(plot = plot_sig, sig_wide = sig_wide, pval_wide = pval_wide))
}

# ------------------ MAKE BOTH PLOTS + RESULTS ------------------
# Final requested EPS names:
#  - Bonf (0.001): significant_pairs_plot_percent.eps
#  - p<0.05:       significant_pairs_plot_percent0.05.eps

# (1) Bonferroni threshold 0.001
out_bonf <- make_plot_and_tables(
  df = df_all,
  thr = bonf_threshold,
  out_eps = file.path(data_dir, "significant_pairs_plot_percent.eps"),
  out_sig_csv = file.path(data_dir, "significance_matrix_bonf0.001.csv"),
  out_pval_csv = file.path(data_dir, "min_pvalues_matrix_bonf0.001.csv"),
  plot_title_suffix = " (Bonferroni p < 0.001)"
)

# (2) Nominal threshold 0.05
out_p05 <- make_plot_and_tables(
  df = df_all,
  thr = p05_threshold,
  out_eps = file.path(data_dir, "significant_pairs_plot_percent0.05.eps"),
  out_sig_csv = file.path(data_dir, "significance_matrix_p0.05.csv"),
  out_pval_csv = file.path(data_dir, "min_pvalues_matrix_p0.05.csv"),
  plot_title_suffix = " (p < 0.05)"
)

# Optional: quick peek at one result table
print(head(out_bonf$sig_wide, 10))
#########DONE#########################
################################################################
# ------------------ Setup ------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(viridis)
library(ggsci)

# Path to your CSV files
data_dir <- "RDA_results/"  # directory of final_results_*.csv produced by the exposure_*.R scripts

csv_files <- list.files(
  path = data_dir, pattern = "^final_results_.*\\.csv$",
  full.names = TRUE
)
if (length(csv_files) == 0) stop("No CSV files found in the specified directory.")

data_list <- setNames(
  lapply(csv_files, read.csv, stringsAsFactors = FALSE),
  str_remove(basename(csv_files), "\\.csv$")
)

# ------------------ Exposure Abbreviation Map ------------------
exposure_map <- c(
  "BFP" = "BF", "bmi" = "BMI", "bw" = "BW",
  "dbp" = "DBP", "drinkpw" = "Alcohol", "eversmok" = "Smoking",
  "FG" = "FG", "hdl" = "HDL", "height" = "Height",
  "ldl" = "LDL", "logTG" = "TG", "sbp" = "SBP"
)

# ------------------ Outcome Normalization Map ------------------
outcome_map <- c(
  "cad" = "CAD", "cad_mvp" = "CAD",
  "t2d" = "T2D", "t2d_mvp" = "T2D",
  "stroke" = "Stroke", "stk" = "Stroke", "stk_mvp" = "Stroke",
  "asthma" = "Asthma", "astm" = "Asthma", "astm_mvp" = "Asthma"
)

# ------------------ Causal Pair Categories ------------------
causal_pairs <- list(
  considered_causal = c(
    "SBP -> CAD", "DBP -> CAD", "Smoking -> CAD",
    "SBP -> Stroke", "DBP -> Stroke", "Smoking -> Stroke",
    "LDL -> CAD", "Smoking -> T2D", "LDL -> Stroke"
  ),
  supported_by_literature = c(
    "BMI -> CAD", "TG -> CAD", "BMI -> T2D",
    "BF -> CAD", "BF -> T2D", "FG -> T2D",
    "Height -> CAD", "BMI -> Stroke", "BF -> Stroke",
    "Smoking -> Asthma"
  ),
  unknown_conflict = c(
    "SBP -> T2D", "HDL -> T2D", "DBP -> T2D", "TG -> T2D",
    "BW -> CAD", "BW -> T2D", "TG -> Stroke", "BMI -> Asthma",
    "FG -> CAD", "FG -> Stroke", "BF -> Asthma", "BW -> Stroke",
    "Alcohol -> Stroke", "Height -> Stroke", "LDL -> T2D",
    "Alcohol -> T2D", "Alcohol -> CAD"
  ),
  implausible_unsupported = c(
    "Alcohol -> Asthma", "Height -> Asthma", "SBP -> Asthma",
    "DBP -> Asthma", "FG -> Asthma", "TG -> Asthma",
    "Height -> T2D", "LDL -> Asthma", "BW -> Asthma",
    "HDL -> Asthma"
  ),
  considered_noncausal = c(
    "HDL -> CAD", "HDL -> Stroke"
  )
)

pair_categories <- imap_dfr(causal_pairs, ~ {
  tibble(
    exposure = str_extract(.x, "^[^ ]+"),
    outcome  = str_extract(.x, "[^ ]+$"),
    category = .y
  )
}) %>%
  mutate(category = recode(category,
                           considered_causal         = "causal",
                           supported_by_literature   = "causal",
                           considered_noncausal      = "noncausal",
                           unknown_conflict          = "correlated",
                           implausible_unsupported   = "unsupported"
  ))

pair_lookup <- pair_categories %>% rename(dict_category = category)

# ------------------ Combine & Normalize All Data ------------------
# NOTE: We will create *two* significance flags later: p<0.001 and p<0.05
bonf_threshold <- 0.001
p05_threshold  <- 0.05

df_all <- map2_dfr(data_list, names(data_list), ~ {
  exposure_raw <- str_remove(.y, "^final_results_")
  if (tolower(exposure_raw) %in% c("crp", "vitd")) return(NULL)
  
  exposure_clean <- dplyr::recode(exposure_raw, !!!exposure_map, .default = exposure_raw)
  
  .x %>%
    mutate(
      exposure = exposure_clean,
      outcome = {
        tmp <- stringr::str_to_lower(trait1)
        tmp <- stringr::str_replace(tmp, "_mvp$", "")
        mapped <- dplyr::recode(tmp, !!!outcome_map, .default = tmp)
        ifelse(mapped %in% c("CAD", "T2D", "Stroke", "Asthma"), mapped, toupper(mapped))
      },
      file = .y,
      method = case_when(
        str_trim(method) == "new_method_weighted" ~ "IB-MODE",
        str_trim(method) == "IB-MR-PRESSO" ~ "IB-PRESSO",
        TRUE ~ str_trim(method)
      )
    )
}) %>%
  filter(!str_detect(method, "(?i)IB-Mix|Simple mode")) %>%
  mutate(
    pval = suppressWarnings(as.numeric(pval))
  )

# ------------------ Keep only dictionary pairs & attach category ------------------
df_all <- df_all %>%
  inner_join(pair_lookup, by = c("exposure", "outcome")) %>%
  mutate(category = dict_category) %>%
  dplyr::select(-dict_category)

# Required facet order: noncausal, unsupp, corr, causal (internally: noncausal, unsupported, correlated, causal)
cat_levels <- c("noncausal", "unsupported", "correlated", "causal")
cat_short_map <- c(
  noncausal   = "noncausal",
  unsupported = "unsupp",
  correlated  = "corr",
  causal      = "causal"
)

df_all <- df_all %>%
  mutate(category = factor(as.character(category), levels = cat_levels))

# ------------------ Keep ONLY first auxiliary trait (preserve file order) ------------------
df_all <- df_all %>%
  group_by(exposure, outcome, method) %>%
  group_modify(~ {
    if (grepl("^IB", .y$method)) {
      if (nrow(.x) >= 2) .x[2, ] else .x[1, ]   # fallback if only 1 row
    } else {
      .x[1, ]
    }
  }) %>%
  ungroup()

# ------------------ Plot + table maker ------------------
make_plot_and_tables <- function(df, thr, out_eps, out_sig_csv, out_pval_csv, plot_title_suffix) {
  
  df2 <- df %>%
    mutate(sig = !is.na(pval) & pval < thr)
  
  method_category_totals <- df2 %>%
    distinct(exposure, outcome, method, category) %>%
    count(method, category, name = "total_n")
  
  summary_counts <- df2 %>%
    group_by(method, category) %>%
    summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop") %>%
    left_join(method_category_totals, by = c("method", "category")) %>%
    mutate(perc_sig = 100 * n_sig / total_n)
  
  category_labels <- method_category_totals %>%
    group_by(category) %>%
    summarise(total_n = max(total_n, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      category = factor(as.character(category), levels = cat_levels),
      facet_lab = paste0(recode(as.character(category), !!!cat_short_map), " (N=", total_n, ")")
    ) %>%
    arrange(category) %>%
    dplyr::select(category, facet_lab)
  
  summary_counts_plot <- summary_counts %>%
    left_join(category_labels, by = "category") %>%
    mutate(facet_lab = factor(facet_lab, levels = category_labels$facet_lab)) %>%
    filter(!is.na(perc_sig) & perc_sig >= 0 & perc_sig <= 100)
  
  # ---- method ordering (IB-MODE with Weighted mode; IB-MR-PRESSO with MR-PRESSO) ----
  all_methods <- unique(df2$method)
  mode_methods <- c("IB-MODE", "Weighted mode")
  presso_methods <- c("IB-MR-PRESSO", "MR-PRESSO")
  remaining_methods <- setdiff(sort(all_methods), c(mode_methods, presso_methods))
  desired_method_order <- c(mode_methods, presso_methods, remaining_methods)
  
  summary_counts_plot$method <- factor(summary_counts_plot$method, levels = desired_method_order)
  
  group_breaks <- c(
    length(mode_methods) + 0.5,
    length(mode_methods) + length(presso_methods) + 0.5
  )
  
  if (nrow(summary_counts_plot) > 0) {
    
    plot_sig <- ggplot(summary_counts_plot, aes(x = method, y = perc_sig, fill = method)) +
      geom_col(width = 0.7, color = "gray30") +
      geom_vline(xintercept = group_breaks, linetype = "dotted", linewidth = 1) +
      facet_wrap(~ facet_lab, scales = "free_x", nrow = 1) +
      scale_y_continuous(
        breaks = scales::breaks_pretty(n = 5),
        expand = expansion(mult = c(0, 0.05)),
        limits = c(0, 100)
      ) +
      labs(
        x = NULL,
        y = "Percentage of Significant Pairs (%)",
        title = paste0("Percentage of Significant Pairs by Method & Category", plot_title_suffix),
        fill = "Method"
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        panel.grid.major.x = element_blank()
      ) +
      scale_fill_d3(
        palette = "category20",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          nrow = 2,
          byrow = TRUE
        )
      )
    
    ggsave(out_eps, plot = plot_sig, device = cairo_ps, width = 10, height = 5, units = "in")
    print(plot_sig)
    
  } else {
    message("No dictionary pairs found after normalization — plot skipped for threshold: ", thr)
    plot_sig <- NULL
  }
  
  # ---- wide tables for this threshold ----
  sig_table <- df2 %>%
    group_by(exposure, outcome, method, category) %>%
    summarise(
      any_sig = any(sig, na.rm = TRUE),
      min_p   = suppressWarnings(min(pval, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      min_p = ifelse(is.infinite(min_p), NA_real_, min_p),
      category = factor(as.character(category), levels = cat_levels)
    )
  
  sig_wide <- sig_table %>%
    mutate(`Significant?` = if_else(any_sig, "✓", "")) %>%
    dplyr::select(exposure, outcome, category, method, `Significant?`) %>%
    pivot_wider(names_from = method, values_from = `Significant?`) %>%
    arrange(category, outcome, exposure)
  
  pval_wide <- sig_table %>%
    dplyr::select(exposure, outcome, category, method, min_p) %>%
    pivot_wider(names_from = method, values_from = min_p) %>%
    arrange(category, outcome, exposure)
  
  write.csv(sig_wide, out_sig_csv, row.names = FALSE)
  write.csv(pval_wide, out_pval_csv, row.names = FALSE)
  
  invisible(list(plot = plot_sig, sig_wide = sig_wide, pval_wide = pval_wide))
}
# ------------------ MAKE BOTH PLOTS + RESULTS ------------------
# (1) Bonferroni threshold 0.001
out_bonf <- make_plot_and_tables(
  df = df_all,
  thr = bonf_threshold,
  out_eps = file.path(data_dir, "significant_pairs_plot_percent_IB2.eps"),
  out_sig_csv = file.path(data_dir, "significance_matrix_bonf0.001_IB2.csv"),
  out_pval_csv = file.path(data_dir, "min_pvalues_matrix_bonf0.001_IB2.csv"),
  plot_title_suffix = " (Bonferroni p < 0.001)"
)

# (2) Nominal threshold 0.05
out_p05 <- make_plot_and_tables(
  df = df_all,
  thr = p05_threshold,
  out_eps = file.path(data_dir, "significant_pairs_plot_percent0.05_IB2.eps"),
  out_sig_csv = file.path(data_dir, "significance_matrix_p0.05_IB2.csv"),
  out_pval_csv = file.path(data_dir, "min_pvalues_matrix_p0.05_IB2.csv"),
  plot_title_suffix = " (p < 0.05)"
)

# Optional: quick peek at one result table
print(head(out_bonf$sig_wide, 10))
################################################################
################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Calculate Z² as beta² / se²
df_all <- df_all %>%
  mutate(Z2 = (b^2) / (se^2))

# Prepare data for IB-MODE and Weighted mode, causal scenarios
plot_wide <- df_all %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    IBMODE_Z2m1 = `IB-MODE` - 1,
    WMode_Z2m1  = `Weighted mode` - 1,
    pair_label  = paste(exposure, outcome, sep = " → "),
    above_diag  = IBMODE_Z2m1 > WMode_Z2m1
  )

# Count points above diagonal
n_above <- sum(plot_wide$above_diag)
n_total <- nrow(plot_wide)
prop_above <- round(100 * n_above / n_total, 1)

p <- ggplot(plot_wide, aes(x = WMode_Z2m1, y = IBMODE_Z2m1)) +
  geom_point(aes(color = above_diag), size = 2.5) +
  scale_color_manual(
    values = c("TRUE" = "royalblue", "FALSE" = "firebrick"),
    labels = c("TRUE" = "Above y = x", "FALSE" = "Below y = x"),
    name = "IB-MODE > Weighted mode"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", size = 1) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 3,
    fontface = "plain",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.25,
    point.padding = 0.4,
    max.overlaps = Inf
  ) +
  theme_bw(base_size = 10) +
  labs(
    x = expression("Weighted mode  " ~ (Z^2 - 1)),
    y = expression("IB-MODE  " ~ (Z^2 - 1)),
    title = expression("Comparison of " ~ Z^2 - 1 ~ " statistics: IB-MODE vs Weighted mode (causal pairs)")
  ) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0.0)), # transparent legend box
    legend.key = element_rect(fill = alpha("white", 0.0)),        # transparent legend keys
    legend.text = element_text(size = 10, color = "grey40"),      # smaller, grey text
    legend.title = element_text(size = 11, color = "grey20"),     # smaller, light title
    legend.spacing.x = unit(0.5, "cm"),                           # tighten items
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0)           # reduce margin
  )

print(p)
ggsave(
  file.path(data_dir, "IBMODE_vs_Weightedmode_Z2minus1_plot_pubworthy_subtle.eps"),
  plot = p,
  width = 8, height = 6, device = cairo_ps
)


##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Epsilon to avoid log(0)
epsilon <- 1e-8

# Prepare data for IB-MODE and Weighted mode, causal scenarios
plot_var <- df_all %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(mean_var = mean(se^2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_var) %>%
  mutate(
    IBMODE_logvar = log(`IB-MODE` + epsilon),
    WMode_logvar  = log(`Weighted mode` + epsilon),
    pair_label    = paste(exposure, outcome, sep = " → "),
    above_diag    = IBMODE_logvar > WMode_logvar
  ) %>%
  filter(is.finite(IBMODE_logvar), is.finite(WMode_logvar))

# Count points where IB-MODE has higher log variance
n_above <- sum(plot_var$above_diag)
n_total <- nrow(plot_var)
prop_above <- round(100 * n_above / n_total, 1)

# Subtle publication-style log-variance plot
p_logvar <- ggplot(plot_var, aes(x = WMode_logvar, y = IBMODE_logvar)) +
  geom_point(aes(color = above_diag), size = 2.5) +
  scale_color_manual(
    values = c("TRUE" = "royalblue", "FALSE" = "firebrick"),
    labels = c("TRUE" = "IB-MODE > Weighted mode", "FALSE" = "IB-MODE ≤ Weighted mode"),
    name = "Log variance comparison"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", size = 1) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 3,
    fontface = "plain",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.25,
    point.padding = 0.4,
    max.overlaps = Inf
  ) +
  theme_bw(base_size = 10) +
  labs(
    x = expression("Weighted mode log-variance (" ~ log(se^2) ~ ")"),
    y = expression("IB-MODE log-variance (" ~ log(se^2) ~ ")"),
    title = expression("Log-variance comparison: IB-MODE vs Weighted mode (causal pairs)")
  )  +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0.0)),
    legend.key = element_rect(fill = alpha("white", 0.0)),
    legend.text = element_text(size = 10, color = "grey40"),
    legend.title = element_text(size = 11, color = "grey20"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0)
  )

print(p_logvar)

ggsave(
  file.path(data_dir, "IBMODE_vs_Weightedmode_logvariance_plot_pubworthy_subtle.eps"),
  plot = p_logvar,
  width = 8, height = 6, device = cairo_ps
)

##############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

epsilon <- 1e-8  # Avoid division by zero

plot_var <- df_all %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(mean_var = mean(se^2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_var) %>%
  mutate(
    WMode_var  = `Weighted mode` + epsilon,
    IBMODE_var = `IB-MODE` + epsilon,
    var_ratio  = WMode_var / IBMODE_var,                           # *** Here: WMODE / IB-MODE ***
    log_WMode_var = log10(WMode_var),
    pair_label = paste(exposure, outcome, sep = " → ")
  ) %>%
  filter(is.finite(log_WMode_var), is.finite(var_ratio))

p_ratio <- ggplot(plot_var, aes(x = log_WMode_var, y = var_ratio)) +
  geom_point(size = 2.5, color = "royalblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", size = 1) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 3,
    fontface = "plain",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.25,
    point.padding = 0.4,
    max.overlaps = Inf
  ) +
  theme_bw(base_size = 10) +
  labs(
    x = expression("Weighted mode log-variance  " ~ log[10](se^2)),
    y = "Relative efficiency: Var(Weighted mode) / Var(IB-MODE)",
    title = expression("Relative efficiency: IB-MODE vs Weighted mode (causal pairs)")
  ) +
  theme(
    legend.position = "none"
  )

print(p_ratio)

ggsave(
  file.path(data_dir, "IBMODE_vs_Weightedmode_variance_ratio_plot.eps"),
  plot = p_ratio,
  width = 8, height = 6, device = cairo_ps
)

#########################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Prepare data
plot_wide <- df_all %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_WeightedMode = log(`Weighted mode`),
    ratio_ibmode = `IB-MODE` / `Weighted mode`,
    ratio_wmode = `Weighted mode` / `IB-MODE`,
    pair_label  = paste(exposure, outcome, sep = " → ")
  ) %>%
  select(logZ2_WeightedMode, ratio_ibmode, ratio_wmode, pair_label)

# Gather ratios into long format for facets:
plot_long <- plot_wide %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "ratio_type",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_type = factor(ratio_type, levels = c("ratio_ibmode","ratio_wmode"),
                        labels = c("IB-MODE / Weighted mode", "Weighted mode / IB-MODE")),
    above_diag = ratio_value > 1
  )

# Panel plot (now with fixed y axis scale)
p_panel <- ggplot(plot_long, aes(x = logZ2_WeightedMode, y = ratio_value)) +
  geom_point(aes(color = above_diag), size = 2.5) +
  scale_color_manual(
    values = c("TRUE" = "royalblue", "FALSE" = "firebrick"),
    labels = c("TRUE" = "Ratio > 1", "FALSE" = "Ratio ≤ 1"),
    name = "Ratio > 1"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", size = 1) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 2.8,
    fontface = "plain",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.2,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  facet_wrap(~ratio_type, scales = "fixed") +
  theme_bw(base_size = 10) +
  labs(
    x = expression("log " ~ Z^2 ~ " [Weighted mode]"),
    y = expression("Ratio"),
    title = expression("Comparison of Z^2 Ratios vs. log " ~ Z^2 ~ " (Weighted mode)")
  ) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 0)),
    legend.text = element_text(size = 10, color = "grey40"),
    legend.title = element_text(size = 11, color = "grey20"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0)
  )

print(p_panel)

# # Optionally save:
# ggsave("panel_z2ratios_vs_logz2.eps", plot = p_panel, width = 10, height = 5, device = cairo_ps)
#################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Prepare data
plot_wide <- df_all %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category %in% c("unsupported","noncausal"), method %in% c("IB-MR-PRESSO", "MR-PRESSO")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_MRPRESSO = log(`MR-PRESSO`),
    ratio_ibmrpresso = `IB-MR-PRESSO` / `MR-PRESSO`,
    ratio_mrpresso = `MR-PRESSO` / `IB-MR-PRESSO`,
    pair_label  = paste(exposure, outcome, sep = " → ")
  ) %>%
  select(logZ2_MRPRESSO, ratio_ibmrpresso, ratio_mrpresso, pair_label)

# Gather ratios into long format for facets:
plot_long <- plot_wide %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "ratio_type",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_type = factor(ratio_type,
                        levels = c("ratio_ibmrpresso","ratio_mrpresso"),
                        labels = c("IB-MR-PRESSO / MR-PRESSO", "MR-PRESSO / IB-MR-PRESSO")
    ),
    above_diag = ratio_value > 1
  )

# Panel plot (now with fixed y axis scale)
p_panel <- ggplot(plot_long, aes(x = logZ2_MRPRESSO, y = ratio_value)) +
  geom_point(aes(color = above_diag), size = 2.5) +
  scale_color_manual(
    values = c("TRUE" = "royalblue", "FALSE" = "firebrick"),
    labels = c("TRUE" = "Ratio > 1", "FALSE" = "Ratio ≤ 1"),
    name = "Ratio > 1"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", size = 1) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 2.8,
    fontface = "plain",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.2,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  facet_wrap(~ratio_type, scales = "fixed") +
  theme_bw(base_size = 10) +
  labs(
    x = expression("log " ~ Z^2 ~ " [MR-PRESSO]"),
    y = expression("Ratio"),
    title = expression("Comparison of Z^2 Ratios vs. log " ~ Z^2 ~ " (MR-PRESSO)")
  ) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 0)),
    legend.text = element_text(size = 10, color = "grey40"),
    legend.title = element_text(size = 11, color = "grey20"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0)
  )

print(p_panel)

# # Optionally save:
# ggsave("panel_z2ratios_vs_logz2_mrpresso.eps", plot = p_panel, width = 10, height = 5, device = cairo_ps)

#################################################################################

# Load the packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ---- Filter to Significant Trait Pairs in Any Method ----
# Identify trait pairs that are significant for at least one method
signif_pairs <- df_all %>%
  group_by(exposure, outcome) %>%
  summarise(any_signif = any(pval < 0.005, na.rm = TRUE), .groups = "drop") %>%
  filter(any_signif) %>%
  select(exposure, outcome)

# Filter df_all to keep only those trait pairs
df_all_sig <- df_all %>%
  inner_join(signif_pairs, by = c("exposure", "outcome"))

# ---- Data Preparation: IB-MODE vs WEIGHTED MODE ----
plot_wide_mode <- df_all_sig %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_ref = log(`Weighted mode`),
    ratio_ib_vs_ref = `IB-MODE` / `Weighted mode`,
    ratio_ref_vs_ib = `Weighted mode` / `IB-MODE`,
    pair_label  = paste(exposure, outcome, sep = " → "),
    method_pair = "IB-MODE vs Weighted mode"
  ) %>%
  select(logZ2_ref, ratio_ib_vs_ref, ratio_ref_vs_ib, pair_label, method_pair)

# ---- Data Preparation: IB-PRESSO vs MR-PRESSO ----
plot_wide_presso <- df_all_sig %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category == "causal", method %in% c("IB-PRESSO", "MR-PRESSO")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_ref = log(`MR-PRESSO`),
    ratio_ib_vs_ref = `IB-PRESSO` / `MR-PRESSO`,
    ratio_ref_vs_ib = `MR-PRESSO` / `IB-PRESSO`,
    pair_label  = paste(exposure, outcome, sep = " → "),
    method_pair = "IB-PRESSO vs MR-PRESSO"
  ) %>%
  select(logZ2_ref, ratio_ib_vs_ref, ratio_ref_vs_ib, pair_label, method_pair)

# ---- Combine and Long Format for Facets ----
plot_wide_all <- bind_rows(plot_wide_mode, plot_wide_presso)
plot_long_all <- plot_wide_all %>%
  pivot_longer(
    cols = c(ratio_ib_vs_ref, ratio_ref_vs_ib),
    names_to = "ratio_type",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_type = factor(
      ratio_type,
      levels = c("ratio_ib_vs_ref", "ratio_ref_vs_ib"),
      labels = c("IB method / Reference", "Reference / IB method")
    ),
    method_pair = factor(
      method_pair,
      levels = c("IB-MODE vs Weighted mode", "IB-PRESSO vs MR-PRESSO")
    ),
    ratio_group = case_when(
      abs(ratio_value - 1) < 1e-8 ~ "Ratio = 1",
      ratio_value > 1 ~ "Ratio > 1",
      ratio_value < 1 ~ "Ratio < 1"
    ),
    ratio_group = factor(ratio_group, levels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"))
  )

# ---- Plot ----
p_unified <- ggplot(plot_long_all, aes(x = logZ2_ref, y = ratio_value)) +
  geom_point(aes(color = ratio_group), size = 2.8, alpha = 0.85, stroke=0.2) +
  scale_color_manual(
    values = c("Ratio > 1" = "#2166AC",       # blue
               "Ratio < 1" = "#B2182B",       # red
               "Ratio = 1" = "#F4A51C"),      # orange/gold
    breaks = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    labels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    name = expression(bold("Z² Ratio Class"))
  ) +
  geom_hline(yintercept = 1, linetype = "longdash", color = "grey50", size = 0.7) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 2.7,
    fontface = "italic",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.18,
    point.padding = 0.28,
    max.overlaps = 12,
    family = "Helvetica"
  ) +
  facet_grid(method_pair ~ ratio_type, scales = "free_y") +
  labs(
    x = expression(bold("log ") * Z^2 * " [Reference method]"),
    y = expression(bold("Z² Ratio")),
    title = "Comparative Z² Ratios between IB and Reference Methods",
    subtitle = "Each point: Exposure → Outcome pair. Rows: Method; Columns: Ratio direction."
  ) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(1.3, "lines"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 7)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 0)),
    legend.text = element_text(size = 13, color = "grey30"),
    legend.title = element_text(size = 14, color = "grey10", face="bold"),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
    strip.background = element_rect(fill = "grey96", color = "grey70", linewidth = 0.8),
    strip.text.x = element_text(face = "bold", size = 15, margin = margin(b=2)),
    strip.text.y = element_text(face = "bold", size = 15, angle = 270, margin = margin(r=2))
  )

print(p_unified)
ggsave("unified_2x2_ibmode_ibmrpresso_causal.eps", plot = p_unified, width = 11, height = 10, device = cairo_ps)

#################################################################################
# NEW CODE: 2x1 plot with only IB/Reference ratio
#################################################################################

# ---- Prepare data for 2x1 plot (IB/Reference only) ----
plot_long_ib_only <- plot_wide_all %>%
  pivot_longer(
    cols = ratio_ib_vs_ref,
    names_to = "ratio_type",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_type = factor(
      ratio_type,
      levels = "ratio_ib_vs_ref",
      labels = "IB method / Reference"
    ),
    method_pair = factor(
      method_pair,
      levels = c("IB-MODE vs Weighted mode", "IB-PRESSO vs MR-PRESSO")
    ),
    ratio_group = case_when(
      abs(ratio_value - 1) < 1e-8 ~ "Ratio = 1",
      ratio_value > 1 ~ "Ratio > 1",
      ratio_value < 1 ~ "Ratio < 1"
    ),
    ratio_group = factor(ratio_group, levels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"))
  )

# ---- 2x1 Plot: IB/Reference only ----
p_ib_only <- ggplot(plot_long_ib_only, aes(x = logZ2_ref, y = ratio_value)) +
  geom_point(aes(color = ratio_group), size = 3.2, alpha = 0.85, stroke=0.2) +
  scale_color_manual(
    values = c("Ratio > 1" = "#2166AC",       # blue
               "Ratio < 1" = "#B2182B",       # red
               "Ratio = 1" = "#F4A51C"),      # orange/gold
    breaks = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    labels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    name = expression(bold("Z² Ratio Class"))
  ) +
  geom_hline(yintercept = 1, linetype = "longdash", color = "grey50", size = 0.7) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 3.0,
    fontface = "italic",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.20,
    point.padding = 0.30,
    max.overlaps = 15,
    family = "Helvetica"
  ) +
  facet_wrap(~ method_pair, ncol = 1, scales = "free_y") +
  labs(
    x = expression(bold("log ") * Z^2 * " [Reference method]"),
    y = expression(bold("Z² Ratio (IB method / Reference)")),
    title = "IB Methods vs Reference Methods: Z² Ratios"
  ) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 0.7,
    panel.spacing = unit(1.5, "lines"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 7)),
    plot.subtitle = element_text(size = 12.5, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 0)),
    legend.text = element_text(size = 13, color = "grey30"),
    legend.title = element_text(size = 14, color = "grey10", face="bold"),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
    strip.background = element_rect(fill = "grey96", color = "grey70", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 15, margin = margin(t=2, b=2))
  )

print(p_ib_only)
ggsave("unified_2x1_IB_vs_reference_only.eps", plot = p_ib_only, width = 8, height = 10, device = cairo_ps)

#############################################
############################################

library(dplyr)
library(tidyr)

df_causal <- df_all_sig %>%
  filter(category == "causal")

# IB-MODE vs Weighted mode
mode_tab <- df_causal %>%
  filter(method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(stat = mean(b^2 / se^2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = stat) %>%
  mutate(IB_Mode = 100 * (`IB-MODE` - `Weighted mode`) / `Weighted mode`) %>%
  group_by(exposure) %>%
  summarise(IB_Mode = mean(IB_Mode, na.rm = TRUE))

# IB-MR-PRESSO vs MR-PRESSO
presso_tab <- df_causal %>%
  filter(method %in% c("IB-MR-PRESSO", "MR-PRESSO")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(stat = mean(b^2 / se^2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = stat) %>%
  mutate(IB_PRESSO = 100 * (`IB-MR-PRESSO` - `MR-PRESSO`) / `MR-PRESSO`) %>%
  group_by(exposure) %>%
  summarise(IB_PRESSO = mean(IB_PRESSO, na.rm = TRUE))

# Final table
final_table <- full_join(mode_tab, presso_tab, by = "exposure")

# Order exposures
final_table$exposure <- factor(final_table$exposure,
                               levels = c("FG","TG","BF","BMI","DBP","Height","LDL","SBP","Smoking"))

final_table <- arrange(final_table, exposure)

print(final_table)

# ------------------ Setup ------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(viridis)
library(ggsci)
library(ggrepel)

# Path to your CSV files
data_dir <- "RDA_results/"  # directory of final_results_*.csv produced by the exposure_*.R scripts

csv_files <- list.files(
  path = data_dir, pattern = "^final_results_.*\\.csv$",
  full.names = TRUE
)
if (length(csv_files) == 0) stop("No CSV files found in the specified directory.")

data_list <- setNames(
  lapply(csv_files, read.csv, stringsAsFactors = FALSE),
  str_remove(basename(csv_files), "\\.csv$")
)

# ------------------ Manual override for IB methods with [2] values ------------------
ib_mode_2_values <- tribble(
  ~exposure, ~outcome, ~estimate, ~ci_lower, ~ci_upper,
  # HDL
  "HDL", "Asthma", -0.07, -0.16, 0.03,
  "HDL", "CAD", -0.06, -0.13, 0.01,
  "HDL", "eGFR", 0.04, 0.01, 0.07,
  "HDL", "Stroke", -0.13, -0.22, -0.05,
  "HDL", "T2D", -0.08, -0.13, -0.04,
  # LDL
  "LDL", "HTN", 0.08, 0.04, 0.13,
  "LDL", "Asthma", -0.01, -0.12, 0.10,
  "LDL", "CAD", 0.22, 0.15, 0.30,
  "LDL", "eGFR", 0.06, 0.03, 0.09,
  "LDL", "Stroke", 0.08, -0.02, 0.18,
  "LDL", "T2D", -0.03, -0.10, 0.03,
  # Triglyceride
  "TG", "HTN", 0.17, 0.12, 0.21,
  "TG", "Asthma", -0.06, -0.13, 0.02,
  "TG", "CAD", 0.27, 0.20, 0.34,
  "TG", "eGFR", -0.04, -0.06, -0.01,
  "TG", "Stroke", 0.09, 0.02, 0.17,
  "TG", "T2D", 0.10, 0.02, 0.18,
  # SBP
  "SBP", "Asthma", 0.00, -0.02, 0.03,
  "SBP", "CAD", 0.03, 0.02, 0.05,
  "SBP", "eGFR", -0.01, -0.02, -0.01,
  "SBP", "Stroke", 0.05, 0.03, 0.07,
  "SBP", "T2D", 0.01, -0.00, 0.03,
  # DBP
  "DBP", "Asthma", -0.00, -0.04, 0.03,
  "DBP", "CAD", 0.04, 0.02, 0.06,
  "DBP", "eGFR", -0.02, -0.03, -0.01,
  "DBP", "Stroke", 0.07, 0.04, 0.10,
  "DBP", "T2D", 0.02, -0.01, 0.04,
  # BMI
  "BMI", "HDLC", -0.56, -0.68, -0.43,
  "BMI", "HTN", 1.01, 0.88, 1.13,
  "BMI", "LDLC", -0.28, -0.40, -0.16,
  "BMI", "Asthma", 0.22, -0.11, 0.55,
  "BMI", "CAD", 0.61, 0.44, 0.78,
  "BMI", "eGFR", -0.13, -0.19, -0.07,
  "BMI", "Stroke", 0.12, -0.10, 0.35,
  "BMI", "T2D", 0.93, 0.79, 1.08,
  # BFP
  "BF", "Asthma", 0.44, -0.02, 0.91,
  "BF", "CAD", 0.56, 0.29, 0.84,
  "BF", "Stroke", 0.42, 0.09, 0.75,
  "BF", "T2D", 1.15, 0.96, 1.35,
  # HEIGHT
  "Height", "Asthma", 0.07, -0.13, 0.28,
  "Height", "CAD", -0.19, -0.31, -0.07,
  "Height", "Stroke", -0.17, -0.41, 0.06,
  "Height", "T2D", 0.18, -0.02, 0.38,
  # EVERSMOK
  "Smoking", "Asthma", -0.24, -1.48, 1.01,
  "Smoking", "CAD", 0.57, -0.69, 1.83,
  "Smoking", "Stroke", 1.41, -0.33, 3.14,
  "Smoking", "T2D", 1.51, 0.18, 2.84,
  # DRINKPW
  "Alcohol", "Asthma", -0.01, -0.25, 0.22,
  "Alcohol", "CAD", 0.22, 0.04, 0.40,
  "Alcohol", "Stroke", 0.01, -0.22, 0.24,
  "Alcohol", "T2D", 0.59, 0.09, 1.09,
  # FG
  "FG", "Asthma", -0.09, -0.31, 0.12,
  "FG", "CAD", 0.22, 0.08, 0.36,
  "FG", "Stroke", 0.04, -0.18, 0.27,
  "FG", "T2D", 0.63, 0.49, 0.78,
  # BW
  "BW", "Asthma", 0.09, -0.35, 0.53,
  "BW", "CAD", 0.10, -0.32, 0.53,
  "BW", "Stroke", -0.08, -0.60, 0.44,
  "BW", "T2D", -0.56, -0.99, -0.13
)

ib_presso_2_values <- tribble(
  ~exposure, ~outcome, ~estimate, ~ci_lower, ~ci_upper,
  # HDL
  "HDL", "Asthma", 0.05, -0.00, 0.09,
  "HDL", "CAD", -0.30, -0.36, -0.25,
  "HDL", "eGFR", 0.04, 0.02, 0.06,
  "HDL", "Stroke", -0.08, -0.13, -0.02,
  "HDL", "T2D", -0.23, -0.27, -0.18,
  # LDL
  "LDL", "HTN", 0.10, 0.07, 0.14,
  "LDL", "Asthma", -0.04, -0.09, 0.01,
  "LDL", "CAD", 0.49, 0.44, 0.53,
  "LDL", "eGFR", 0.06, 0.04, 0.08,
  "LDL", "Stroke", 0.18, 0.13, 0.24,
  "LDL", "T2D", -0.03, -0.07, 0.01,
  # Triglyceride
  "TG", "HTN", 0.32, 0.27, 0.36,
  "TG", "Asthma", -0.06, -0.11, -0.01,
  "TG", "CAD", 0.41, 0.36, 0.45,
  "TG", "eGFR", -0.01, -0.03, 0.02,
  "TG", "Stroke", 0.13, 0.07, 0.18,
  "TG", "T2D", 0.12, 0.06, 0.17,
  # SBP
  "SBP", "Asthma", 0.00, -0.00, 0.01,
  "SBP", "CAD", 0.04, 0.04, 0.05,
  "SBP", "eGFR", -0.01, -0.01, -0.01,
  "SBP", "Stroke", 0.03, 0.03, 0.04,
  "SBP", "T2D", 0.02, 0.02, 0.03,
  # DBP
  "DBP", "Asthma", -0.00, -0.01, 0.01,
  "DBP", "CAD", 0.06, 0.05, 0.07,
  "DBP", "eGFR", -0.01, -0.02, -0.01,
  "DBP", "Stroke", 0.04, 0.03, 0.05,
  "DBP", "T2D", 0.02, 0.01, 0.02,
  # BMI
  "BMI", "HDLC", -0.53, -0.58, -0.47,
  "BMI", "HTN", 0.71, 0.64, 0.78,
  "BMI", "LDLC", -0.21, -0.24, -0.18,
  "BMI", "Asthma", 0.23, 0.15, 0.31,
  "BMI", "CAD", 0.39, 0.32, 0.47,
  "BMI", "eGFR", -0.09, -0.12, -0.05,
  "BMI", "Stroke", 0.21, 0.13, 0.30,
  "BMI", "T2D", 0.92, 0.86, 0.98,
  # BFP
  "BF", "Asthma", 0.28, 0.18, 0.38,
  "BF", "CAD", 0.43, 0.33, 0.53,
  "BF", "Stroke", 0.31, 0.21, 0.42,
  "BF", "T2D", 1.18, 1.10, 1.26,
  # HEIGHT
  "Height", "Asthma", -0.00, -0.05, 0.04,
  "Height", "CAD", -0.13, -0.17, -0.09,
  "Height", "Stroke", -0.05, -0.09, 0.00,
  "Height", "T2D", -0.02, -0.06, 0.01,
  # EVERSMOK
  "Smoking", "Asthma", 0.26, -0.08, 0.60,
  "Smoking", "CAD", 0.59, 0.23, 0.94,
  "Smoking", "Stroke", 0.68, 0.27, 1.09,
  "Smoking", "T2D", 0.20, -0.17, 0.57,
  # DRINKPW
  "Alcohol", "Asthma", -0.10, -0.26, 0.07,
  "Alcohol", "CAD", 0.09, -0.06, 0.24,
  "Alcohol", "Stroke", 0.07, -0.20, 0.34,
  "Alcohol", "T2D", -0.01, -0.56, 0.54,
  # FG
  "FG", "Asthma", -0.09, -0.19, 0.01,
  "FG", "CAD", 0.21, 0.12, 0.31,
  "FG", "Stroke", 0.09, -0.03, 0.20,
  "FG", "T2D", 1.09, 0.92, 1.26,
  # BW
  "BW", "Asthma", 0.05, -0.06, 0.15,
  "BW", "CAD", -0.16, -0.26, -0.05,
  "BW", "Stroke", -0.08, -0.20, 0.04,
  "BW", "T2D", -0.27, -0.37, -0.16
)

# Calculate p-values from confidence intervals (approximate)
calculate_pval <- function(est, lower, upper) {
  se <- (upper - lower) / (2 * 1.96)
  z <- abs(est / se)
  pval <- 2 * (1 - pnorm(z))
  return(pval)
}

ib_mode_2_values <- ib_mode_2_values %>%
  mutate(pval = calculate_pval(estimate, ci_lower, ci_upper),
         method = "IB-MODE")

ib_presso_2_values <- ib_presso_2_values %>%
  mutate(pval = calculate_pval(estimate, ci_lower, ci_upper),
         method = "IB-MR-PRESSO")

# Combine the manual values
manual_ib_values <- bind_rows(ib_mode_2_values, ib_presso_2_values)

# ------------------ Exposure Abbreviation Map ------------------
exposure_map <- c(
  "BFP" = "BF", "bmi" = "BMI", "bw" = "BW",
  "dbp" = "DBP", "drinkpw" = "Alcohol", "eversmok" = "Smoking",
  "FG" = "FG", "hdl" = "HDL", "height" = "Height",
  "ldl" = "LDL", "logTG" = "TG", "sbp" = "SBP"
)

# ------------------ Outcome Normalization Map ------------------
outcome_map <- c(
  "cad" = "CAD", "cad_mvp" = "CAD",
  "t2d" = "T2D", "t2d_mvp" = "T2D",
  "stroke" = "Stroke", "stk" = "Stroke", "stk_mvp" = "Stroke",
  "asthma" = "Asthma", "astm" = "Asthma", "astm_mvp" = "Asthma"
)

# ------------------ Causal Pair Categories ------------------
causal_pairs <- list(
  considered_causal = c(
    "SBP -> CAD", "DBP -> CAD", "Smoking -> CAD",
    "SBP -> Stroke", "DBP -> Stroke", "Smoking -> Stroke",
    "LDL -> CAD", "Smoking -> T2D", "LDL -> Stroke"
  ),
  supported_by_literature = c(
    "BMI -> CAD", "TG -> CAD", "BMI -> T2D",
    "BF -> CAD", "BF -> T2D", "FG -> T2D",
    "Height -> CAD", "BMI -> Stroke", "BF -> Stroke",
    "Smoking -> Asthma"
  ),
  unknown_conflict = c(
    "SBP -> T2D", "HDL -> T2D", "DBP -> T2D", "TG -> T2D",
    "BW -> CAD", "BW -> T2D", "TG -> Stroke", "BMI -> Asthma",
    "FG -> CAD", "FG -> Stroke", "BF -> Asthma", "BW -> Stroke",
    "Alcohol -> Stroke", "Height -> Stroke", "LDL -> T2D",
    "Alcohol -> T2D", "Alcohol -> CAD"
  ),
  implausible_unsupported = c(
    "Alcohol -> Asthma", "Height -> Asthma", "SBP -> Asthma",
    "DBP -> Asthma", "FG -> Asthma", "TG -> Asthma",
    "Height -> T2D", "LDL -> Asthma", "BW -> Asthma",
    "HDL -> Asthma"
  ),
  considered_noncausal = c(
    "HDL -> CAD", "HDL -> Stroke"
  )
)

pair_categories <- imap_dfr(causal_pairs, ~ {
  tibble(
    exposure = str_extract(.x, "^[^ ]+"),
    outcome  = str_extract(.x, "[^ ]+$"),
    category = .y
  )
}) %>%
  mutate(category = recode(category,
                           considered_causal         = "causal",
                           supported_by_literature   = "causal",
                           considered_noncausal      = "noncausal",
                           unknown_conflict          = "correlated",
                           implausible_unsupported   = "unsupported"
  ))

pair_lookup <- pair_categories %>% rename(dict_category = category)

# ------------------ Combine & Normalize All Data ------------------
bonf_threshold <- 0.001
p05_threshold  <- 0.05

df_all <- map2_dfr(data_list, names(data_list), ~ {
  exposure_raw <- str_remove(.y, "^final_results_")
  if (tolower(exposure_raw) %in% c("crp", "vitd")) return(NULL)
  
  exposure_clean <- dplyr::recode(exposure_raw, !!!exposure_map, .default = exposure_raw)
  
  .x %>%
    mutate(
      exposure = exposure_clean,
      outcome = {
        tmp <- stringr::str_to_lower(trait1)
        tmp <- stringr::str_replace(tmp, "_mvp$", "")
        mapped <- dplyr::recode(tmp, !!!outcome_map, .default = tmp)
        ifelse(mapped %in% c("CAD", "T2D", "Stroke", "Asthma"), mapped, toupper(mapped))
      },
      file = .y,
      method = if_else(str_trim(method) == "new_method_weighted", "IB-MODE", str_trim(method))
    )
}) %>%
  filter(!str_detect(method, "(?i)IB-Mix|Simple mode")) %>%
  mutate(
    pval = suppressWarnings(as.numeric(pval))
  )

# ------------------ Keep only dictionary pairs & attach category ------------------
df_all <- df_all %>%
  inner_join(pair_lookup, by = c("exposure", "outcome")) %>%
  mutate(category = dict_category) %>%
  dplyr::select(-dict_category)

# ------------------ Replace ONLY IB method values with manual [2] values ------------------
# Keep non-IB methods as is
df_non_ib <- df_all %>%
  filter(!method %in% c("IB-MODE", "IB-MR-PRESSO"))

# Add manual IB values
df_ib <- manual_ib_values %>%
  left_join(pair_lookup, by = c("exposure", "outcome")) %>%
  mutate(
    category = dict_category,
    trait1 = outcome,
    b = estimate,
    se = (ci_upper - ci_lower) / (2 * 1.96)
  ) %>%
  dplyr::select(exposure, outcome, method, b, se, pval, category, trait1)

# Combine
df_all <- bind_rows(df_non_ib, df_ib) %>%
  filter(!is.na(category))

cat_levels <- c("noncausal", "unsupported", "correlated", "causal")
cat_short_map <- c(
  noncausal   = "noncausal",
  unsupported = "unsupp",
  correlated  = "corr",
  causal      = "causal"
)

df_all <- df_all %>%
  mutate(category = factor(as.character(category), levels = cat_levels))

# ============================================================
# Z² RATIO PLOT CODE STARTS HERE
# ============================================================

# ---- Filter to Significant Trait Pairs in Any Method ----
signif_pairs <- df_all %>%
  group_by(exposure, outcome) %>%
  summarise(any_signif = any(pval < 0.005, na.rm = TRUE), .groups = "drop") %>%
  filter(any_signif) %>%
  dplyr::select(exposure, outcome)  # Use dplyr::select explicitly

# Filter df_all to keep only those trait pairs
df_all_sig <- df_all %>%
  inner_join(signif_pairs, by = c("exposure", "outcome"))

# ---- Data Preparation: IB-MODE vs WEIGHTED MODE ----
plot_wide_mode <- df_all_sig %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_ref = log(`Weighted mode`),
    ratio_ib_vs_ref = `IB-MODE` / `Weighted mode`,
    ratio_ref_vs_ib = `Weighted mode` / `IB-MODE`,
    pair_label  = paste(exposure, outcome, sep = " → "),
    method_pair = "IB-MODE vs Weighted mode"
  ) %>%
  dplyr::select(logZ2_ref, ratio_ib_vs_ref, ratio_ref_vs_ib, pair_label, method_pair)  # Use dplyr::select

# ---- Data Preparation: IB-MR-PRESSO vs MR-PRESSO ----
plot_wide_presso <- df_all_sig %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category == "causal", method %in% c("IB-MR-PRESSO", "MR-PRESSO")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_ref = log(`MR-PRESSO`),
    ratio_ib_vs_ref = `IB-MR-PRESSO` / `MR-PRESSO`,
    ratio_ref_vs_ib = `MR-PRESSO` / `IB-MR-PRESSO`,
    pair_label  = paste(exposure, outcome, sep = " → "),
    method_pair = "IB-MR-PRESSO vs MR-PRESSO"
  ) %>%
  dplyr::select(logZ2_ref, ratio_ib_vs_ref, ratio_ref_vs_ib, pair_label, method_pair)  # Use dplyr::select

# ---- Combine and Long Format for Facets ----
plot_wide_all <- bind_rows(plot_wide_mode, plot_wide_presso)
plot_long_all <- plot_wide_all %>%
  pivot_longer(
    cols = c(ratio_ib_vs_ref, ratio_ref_vs_ib),
    names_to = "ratio_type",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_type = factor(
      ratio_type,
      levels = c("ratio_ib_vs_ref", "ratio_ref_vs_ib"),
      labels = c("IB method / Reference", "Reference / IB method")
    ),
    method_pair = factor(
      method_pair,
      levels = c("IB-MODE vs Weighted mode", "IB-MR-PRESSO vs MR-PRESSO")
    ),
    ratio_group = case_when(
      abs(ratio_value - 1) < 1e-8 ~ "Ratio = 1",
      ratio_value > 1 ~ "Ratio > 1",
      ratio_value < 1 ~ "Ratio < 1"
    ),
    ratio_group = factor(ratio_group, levels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"))
  )

# ---- Plot ----
p_unified_IB2 <- ggplot(plot_long_all, aes(x = logZ2_ref, y = ratio_value)) +
  geom_point(aes(color = ratio_group), size = 2.8, alpha = 0.85, stroke=0.2) +
  scale_color_manual(
    values = c("Ratio > 1" = "#2166AC",       # blue
               "Ratio < 1" = "#B2182B",       # red
               "Ratio = 1" = "#F4A51C"),      # orange/gold
    breaks = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    labels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    name = expression(bold("Z² Ratio Class"))
  ) +
  geom_hline(yintercept = 1, linetype = "longdash", color = "grey50", size = 0.7) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 2.7,
    fontface = "italic",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.18,
    point.padding = 0.28,
    max.overlaps = 12,
    family = "Helvetica"
  ) +
  facet_grid(method_pair ~ ratio_type, scales = "free_y") +
  labs(
    x = expression(bold("log ") * Z^2 * " [Reference method]"),
    y = expression(bold("Z² Ratio")),
    title = "Comparative Z² Ratios between IB and Reference Methods",
    subtitle = "Each point: Exposure → Outcome pair. IB methods use second-priority auxiliary traits."
  ) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(1.3, "lines"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 7)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 0)),
    legend.text = element_text(size = 13, color = "grey30"),
    legend.title = element_text(size = 14, color = "grey10", face="bold"),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
    strip.background = element_rect(fill = "grey96", color = "grey70", linewidth = 0.8),
    strip.text.x = element_text(face = "bold", size = 15, margin = margin(b=2)),
    strip.text.y = element_text(face = "bold", size = 15, angle = 270, margin = margin(r=2))
  )

print(p_unified_IB2)
ggsave(file.path(data_dir, "unified_2x2_ibmode_ibmrpresso_causal_IB2.eps"), 
       plot = p_unified_IB2, width = 11, height = 10, device = cairo_ps)

#################################################################################
# OPTIONAL: 2x1 plot with only Reference/IB ratio
#################################################################################

# ---- Prepare data for 2x1 plot (Reference/IB only) ----
plot_long_ref_only <- plot_wide_all %>%
  pivot_longer(
    cols = ratio_ref_vs_ib,
    names_to = "ratio_type",
    values_to = "ratio_value"
  ) %>%
  mutate(
    ratio_type = factor(
      ratio_type,
      levels = "ratio_ref_vs_ib",
      labels = "Reference / IB method"
    ),
    method_pair = factor(
      method_pair,
      levels = c("IB-MODE vs Weighted mode", "IB-MR-PRESSO vs MR-PRESSO")
    ),
    ratio_group = case_when(
      abs(ratio_value - 1) < 1e-8 ~ "Ratio = 1",
      ratio_value > 1 ~ "Ratio > 1",
      ratio_value < 1 ~ "Ratio < 1"
    ),
    ratio_group = factor(ratio_group, levels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"))
  )

# ---- 2x1 Plot: IB/Reference only (rewritten version) ----
p_ref_only_IB2 <- ggplot(plot_long_ib_only, aes(x = logZ2_ref, y = ratio_value)) +
  geom_point(aes(color = ratio_group), size = 3.2, alpha = 0.85, stroke=0.2) +
  scale_color_manual(
    values = c("Ratio > 1" = "#2166AC",       # blue
               "Ratio < 1" = "#B2182B",       # red
               "Ratio = 1" = "#F4A51C"),      # orange/gold
    breaks = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    labels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"),
    name = expression(bold("Z² Ratio Class"))
  ) +
  geom_hline(yintercept = 1, linetype = "longdash", color = "grey50", size = 0.7) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 3.0,
    fontface = "italic",
    color = "grey20",
    segment.color = "grey80",
    box.padding = 0.20,
    point.padding = 0.30,
    max.overlaps = 15,
    family = "Helvetica"
  ) +
  facet_wrap(~ method_pair, ncol = 1, scales = "free_y") +
  labs(
    x = expression(bold("log ") * Z^2 * " [Reference method]"),
    y = expression(bold("Z² Ratio (IB method / Reference)")),
    title = "IB Methods vs Reference Methods: Z² Ratios"
  ) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 0.7,
    panel.spacing = unit(1.5, "lines"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 7)),
    plot.subtitle = element_text(size = 12.5, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 0)),
    legend.text = element_text(size = 13, color = "grey30"),
    legend.title = element_text(size = 14, color = "grey10", face="bold"),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = -8, r = 0, b = 0, l = 0),
    strip.background = element_rect(fill = "grey96", color = "grey70", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 15, margin = margin(t=2, b=2))
  )

print(p_ref_only_IB2)
ggsave(file.path(data_dir, "unified_2x1_IB_vs_reference_only_IB2.eps"), 
       plot = p_ref_only_IB2, width = 8, height = 10, device = cairo_ps)
#########################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(grid)

# ============================================================

# 1. EFFICACY TABLE (IB-MODE EFFICIENCY GAINS ONLY)

# ============================================================

tab_eff <- data.frame(
  Exposure = c("FG","TG","bf","bmi","dbp","height","ldl","sbp","smoking"),
  Gain     = c(-21.1,-14.9,169.0,122.0,58.0,-10.1,450.0,4.0,300.0)
)
tab_eff$Color <- ifelse(tab_eff$Gain >= 0, "green3", "red3")

# Reshape into 3×3 matrix

table_matrix <- matrix(
  paste0(tab_eff$Exposure, "\n", sprintf("%.1f", tab_eff$Gain)),
  nrow = 3, byrow = TRUE
)
color_matrix <- matrix(tab_eff$Color, nrow = 3, byrow = TRUE)

# Prominent table grob

gt <- tableGrob(
  table_matrix,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(
      fg_params = list(col = as.vector(color_matrix), fontsize = 20, fontface = "bold"),
      bg_params = list(fill = "#E8E8E8", col = "black", lwd = 2)
    ),
    base_size = 20
  )
)
gt$gp <- gpar(lwd = 3)   # thick outer border

# ============================================================

# 2. ORIGINAL ANALYSIS (SIGNIFICANT PAIRS → Z² RATIOS)

# ============================================================

signif_pairs <- df_all %>%
  group_by(exposure, outcome) %>%
  summarise(any_signif = any(pval < 0.005, na.rm = TRUE), .groups = "drop") %>%
  filter(any_signif) %>%
  select(exposure, outcome)

df_all_sig <- df_all %>%
  inner_join(signif_pairs, by = c("exposure", "outcome"))

plot_wide_mode <- df_all_sig %>%
  mutate(Z2 = (b^2) / (se^2)) %>%
  filter(category == "causal", method %in% c("IB-MODE", "Weighted mode")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(Z2 = mean(Z2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = Z2) %>%
  mutate(
    logZ2_ref = log(`Weighted mode`),
    ratio_ib_vs_ref = `IB-MODE` / `Weighted mode`,
    pair_label = paste(exposure, outcome, sep = " → ")
  ) %>%
  select(logZ2_ref, ratio_ib_vs_ref, pair_label) %>%
  filter(is.finite(logZ2_ref), is.finite(ratio_ib_vs_ref))

plot_single <- plot_wide_mode %>%
  mutate(
    ratio_group = case_when(
      abs(ratio_ib_vs_ref - 1) < 1e-8 ~ "Ratio = 1",
      ratio_ib_vs_ref > 1 ~ "Ratio > 1",
      ratio_ib_vs_ref < 1 ~ "Ratio < 1"
    ),
    ratio_group = factor(ratio_group, levels = c("Ratio > 1", "Ratio = 1", "Ratio < 1"))
  )

# ============================================================

# 3. MAIN SCATTER PLOT (BIGGER AND BOLDER)

# ============================================================

p_mode <- ggplot(plot_single, aes(x = logZ2_ref, y = ratio_ib_vs_ref)) +
  geom_point(aes(color = ratio_group), size = 7, alpha = 0.95, stroke = 1.8) +
  scale_color_manual(
    values = c("Ratio > 1" = "#2166AC",
               "Ratio < 1" = "#B2182B",
               "Ratio = 1" = "#F4A51C")
  ) +
  geom_hline(yintercept = 1, linetype = "longdash", color = "grey20", size = 2) +
  ggrepel::geom_text_repel(
    aes(label = pair_label),
    size = 6, fontface = "italic", color = "grey20",
    segment.size = 1.2,
    segment.color = "grey40",
    box.padding = 0.35, point.padding = 0.45,
    max.overlaps = 25
  ) +
  labs(
    x = expression(bold("log ") * Z^2 * " [Weighted mode reference]"),
    y = expression(bold("Z² Ratio (IB-MODE / Weighted mode)")),
    title = "Z² Ratio: IB-MODE vs Weighted mode"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text  = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1.5),
    panel.border = element_rect(size = 2),
    legend.position = "right",
    legend.text = element_text(size = 18, face = "bold"),
    aspect.ratio = 1
  )

# ============================================================

# 4. INSERT PROMINENT TABLE INTO PLOT

# ============================================================

p_final <- p_mode +
  annotation_custom(
    grob = gt,
    xmin = max(plot_single$logZ2_ref) - 4,
    xmax = max(plot_single$logZ2_ref) + 0.5,
    ymin = max(plot_single$ratio_ib_vs_ref) - 4.5,
    ymax = max(plot_single$ratio_ib_vs_ref) + 0.5
  )

p_final


# ============================================================
# 5. SAVE THE FINAL PLOT
# ============================================================
ggsave(
  filename = "IB_MODE_efficiency_plot.png",
  plot = p_final,
  width = 10,
  height = 10,
  dpi = 320
)


###################################causal gains

# Load required packages
library(dplyr)
library(tidyr)
# Read in your dataset (assuming it's saved as "data.txt" with tab separation)
df <- read.csv(file.choose(), header = TRUE)
df$Z.2.1=df$Z.2.1+1
# ---- IB-mode vs Weighted mode ----
eff_gain_mode <- df %>%
  filter(method %in% c("IB-mode", "Weighted mode")) %>%
  select(exposure, trait1, method, Z.2.1) %>%
  pivot_wider(names_from = method, values_from = Z.2.1) %>%
  mutate(
    perc_gain_IBmode = 
      ( `IB-mode` - `Weighted mode` ) / `Weighted mode` * 100
  ) %>%
  group_by(exposure) %>%
  summarise(perc_gain_IBmode = mean(perc_gain_IBmode, na.rm = TRUE), .groups = "drop")

# ---- IB-MR-PRESSO vs MR-PRESSO ----
eff_gain_presso <- df %>%
  filter(method %in% c("IB-MR-PRESSO", "MR-PRESSO")) %>%
  select(exposure, trait1, method, Z.2.1) %>%
  pivot_wider(names_from = method, values_from = Z.2.1) %>%
  mutate(
    perc_gain_IBpresso = 
      ( `IB-MR-PRESSO` - `MR-PRESSO` ) / `MR-PRESSO` * 100
  ) %>%
  group_by(exposure) %>%
  summarise(perc_gain_IBpresso = mean(perc_gain_IBpresso, na.rm = TRUE), .groups = "drop")

# ---- Merge into one table ----
eff_gain_table <- full_join(eff_gain_mode, eff_gain_presso, by = "exposure")

print(eff_gain_table)

