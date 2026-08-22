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
        str_trim(method) == "IB-Mode" ~ "IB-MODE",
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
  
  # ---- method ordering (IB-MODE with Weighted mode; IB-PRESSO with MR-PRESSO) ----
  all_methods <- unique(df2$method)
  mode_methods <- c("IB-MODE", "Weighted mode")
  presso_methods <- c("IB-PRESSO", "MR-PRESSO")
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
        str_trim(method) == "IB-Mode" ~ "IB-MODE",
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
  
  # ---- method ordering (IB-MODE with Weighted mode; IB-PRESSO with MR-PRESSO) ----
  all_methods <- unique(df2$method)
  mode_methods <- c("IB-MODE", "Weighted mode")
  presso_methods <- c("IB-PRESSO", "MR-PRESSO")
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

#################################################################################
# 2x1 plot: IB vs reference-only ratio
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

# IB-PRESSO vs MR-PRESSO
presso_tab <- df_causal %>%
  filter(method %in% c("IB-PRESSO", "MR-PRESSO")) %>%
  group_by(exposure, outcome, method) %>%
  summarise(stat = mean(b^2 / se^2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = stat) %>%
  mutate(IB_PRESSO = 100 * (`IB-PRESSO` - `MR-PRESSO`) / `MR-PRESSO`) %>%
  group_by(exposure) %>%
  summarise(IB_PRESSO = mean(IB_PRESSO, na.rm = TRUE))

# Final table
final_table <- full_join(mode_tab, presso_tab, by = "exposure")

# Order exposures
final_table$exposure <- factor(final_table$exposure,
                               levels = c("FG","TG","BF","BMI","DBP","Height","LDL","SBP","Smoking"))

final_table <- arrange(final_table, exposure)

print(final_table)
