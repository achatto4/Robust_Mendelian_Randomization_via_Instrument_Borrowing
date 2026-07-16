
# =========================================================
# Load libraries
# =========================================================
library(dplyr)
library(tidyr)
library(readr)
library(kableExtra)
library(knitr)

# =========================================================
# Define base directory
# =========================================================
base_dir <- "RDA_results/"  # IB-Mode RDA outputs (final_results_*.csv)

# =========================================================
# Helper: Process a single results CSV file
# =========================================================
process_file <- function(file_path) {
  cat("\nProcessing file:", basename(file_path), "\n")
  df <- readr::read_csv(file_path, show_col_types = FALSE)
  
  # Remove '_mvp' from trait1 for readability
  df <- df %>% mutate(trait1 = gsub("_mvp$", "", trait1))
  
  # IB methods
  df_ib <- df %>%
    dplyr::filter(method %in% c("IB-Mode","IB-MR-PRESSO")) %>%
    dplyr::mutate(
      method_label = ifelse(method == "IB-Mode", "IB-MODE", "IB-MR-PRESSO"),
      lower = b - 1.96*se,
      upper = b + 1.96*se,
      ib_cell = sprintf("\\makecell{%0.2f \\\\ (%0.2f, %0.2f)}", b, lower, upper)  # <<< Only one backslash!
    ) %>%
    dplyr::mutate(col_id = method_label) %>%
    dplyr::select(trait1, col_id, ib_cell) %>%
    tidyr::pivot_wider(names_from = col_id, values_from = ib_cell)
  
  # Main MR methods
  main_methods <- c("Inverse variance weighted","MR Egger","MR-ConMix","MR-Mix","MR-cML",
                    "Weighted median","Weighted mode","MR-PRESSO")
  
  df_main <- df %>%
    dplyr::filter(method %in% main_methods) %>%
    dplyr::mutate(
      lower = b - 1.96*se,
      upper = b + 1.96*se,
      cell = sprintf("\\makecell{%0.2f \\\\ (%0.2f, %0.2f)}", b, lower, upper)  # <<< Only one backslash!
    ) %>%
    dplyr::mutate(col_id = dplyr::case_when(
      method == "Inverse variance weighted" ~ "IVW",
      method == "MR Egger" ~ "Egger",
      method == "MR-ConMix" ~ "ConMix",
      method == "MR-Mix" ~ "MRMix",
      method == "MR-cML" ~ "MRcML",
      method == "Weighted median" ~ "Median",
      method == "Weighted mode" ~ "Mode",
      method == "MR-PRESSO" ~ "PRESSO",
      TRUE ~ method
    )) %>%
    dplyr::select(trait1, col_id, cell) %>%
    tidyr::pivot_wider(names_from = col_id, values_from = cell)
  
  # Merge main and IB methods
  final <- dplyr::full_join(df_main, df_ib, by="trait1")
  return(final)
}

# =========================================================
# Build group table
# =========================================================
make_group_table <- function(files, label){
  table_list <- list()
  header <- data.frame(trait1 = label, IVW="",Egger="",ConMix="",MRMix="",
                       MRcML="",Median="",Mode="",PRESSO="",
                       `IB-MODE`="", `IB-MR-PRESSO`="", check.names = FALSE, stringsAsFactors = FALSE)
  table_list[[1]] <- header
  
  # For each file (exposure):
  for (file in files) {
    # Get exposure from filename, upper case
    expo <- toupper(gsub("^final_results_|\\.csv$", "", basename(file)))
    # Make a "row" labelling the exposure (bold LaTeX for prettiness, blank in all other columns)
    exp_row <- header[1,]
    exp_row[1,] <- ""
    exp_row[1,1] <- paste0("\\textbf{", expo, "}")  # Bold upper case
    table_list[[length(table_list) + 1]] <- exp_row
    
    df <- process_file(file.path(base_dir, file))
    if(nrow(df) > 0)
      table_list[[length(table_list) + 1]] <- df
  }
  final_table <- bind_rows(table_list)
  return(final_table)
}

# =========================================================
# Define file groups
# =========================================================
lipid_files      <- c("final_results_hdl.csv","final_results_ldl.csv","final_results_logTG.csv")
bp_files         <- c("final_results_sbp.csv","final_results_dbp.csv")
anthro_files     <- c("final_results_bmi.csv","final_results_BFP.csv","final_results_height.csv")
lifestyle_files  <- c("final_results_eversmok.csv","final_results_drinkpw.csv")
glucose_other_files <- c("final_results_FG.csv","final_results_vitD.csv","final_results_bw.csv")

# =========================================================
# Build all tables
# =========================================================
lipid_table      <- make_group_table(lipid_files, "Lipids")
bp_table         <- make_group_table(bp_files, "Blood pressure")
anthro_table     <- make_group_table(anthro_files, "Anthropometrics")
lifestyle_table  <- make_group_table(lifestyle_files, "Lifestyle")
glucose_table    <- make_group_table(glucose_other_files, "Glucose / Other")

# =========================================================
# Print LaTeX tables
# =========================================================
print_table <- function(tbl, caption){
  kable(tbl,
        format = "latex",
        booktabs = TRUE,
        escape = FALSE,
        caption = caption) %>%
    kable_styling(latex_options = "scale_down", font_size = 9)
}

print_table(lipid_table, "MR estimates for lipid traits (with IB-MODE and IB-MR-PRESSO)")
print_table(bp_table, "MR estimates for blood pressure traits")
print_table(anthro_table, "MR estimates for anthropometric traits")
print_table(lifestyle_table, "MR estimates for lifestyle traits")
print_table(glucose_table, "MR estimates for glucose/other traits")
