
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
    dplyr::filter(method %in% c("new_method_weighted","IB-MR-PRESSO")) %>%
    dplyr::mutate(
      method_label = ifelse(method == "new_method_weighted", "IB-MODE", "IB-MR-PRESSO"),
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

###############################################################################

library(dplyr)
library(tidyr)
library(stringr)

# Read csv file
df <- read.csv("RDA_results_presso/final_results_dbp.csv", stringsAsFactors = FALSE)

mrpresso_method <- "MR-PRESSO"
ib_methods <- c("IB-MR-PRESSO", "new_method_weighted")
traits <- unique(df$trait1)
latex_rows <- c()

for (trait in traits) {
  sub <- df %>% filter(trait1 == trait)
  mrpresso_row <- sub %>% filter(method == mrpresso_method & (trait2 == "NA" | is.na(trait2)))
  ib_rows <- sub %>% filter(method %in% ib_methods & trait2 != "NA" & !is.na(trait2)) %>%
    arrange(factor(method, levels=ib_methods), trait2)
  n_ib <- nrow(ib_rows)
  if (n_ib < 3) {
    ib_rows <- bind_rows(ib_rows, as.data.frame(matrix(NA, ncol=ncol(ib_rows), nrow=3-n_ib, dimnames = list(NULL, colnames(ib_rows)))))
  }
  # MR-PRESSO: beta (beta-1.96*se, beta+1.96*se)
  if (nrow(mrpresso_row) > 0) {
    beta <- mrpresso_row$b[1]
    se <- mrpresso_row$se[1]
    lower <- beta - 1.96*se
    upper <- beta + 1.96*se
    mrpresso_str <- sprintf("%.2f (%.2f, %.2f)", beta, lower, upper)
    z_ref <- abs(beta/se)
  } else {
    mrpresso_str <- ""
    z_ref <- NA
  }
  ib_cells <- c()
  for (i in 1:3) {
    ir <- ib_rows[i, ]
    if (!is.na(ir$method)) {
      beta <- ir$b
      se <- ir$se
      lower <- beta - 1.96 * se
      upper <- beta + 1.96 * se
      effect_ci_line <- sprintf("%.2f (%.2f, %.2f)", beta, lower, upper)
      subtrait <- ir$trait2
      z <- abs(beta/se)
      eff_str <- if (!is.na(z_ref) && !is.na(z)) {
        sprintf("%.0f\\%%", 100 * (z^2 - 1)/(z_ref^2 - 1))
      } else {
        ""
      }
      # Use \makecell with two lines, and escape percent for LaTeX safety.
      cell <- sprintf("\\makecell{%s \\\\ (%s) (%s)}", effect_ci_line, subtrait, eff_str)
    } else {
      cell <- ""
    }
    ib_cells <- c(ib_cells, cell)
  }
  latex_rows <- c(latex_rows, sprintf("%s & %s & %s & %s & %s \\\\", trait, mrpresso_str, ib_cells[1], ib_cells[2], ib_cells[3]))
}

table_string <- paste0(
  "\\begin{tabular}{l",
  strrep(">{\\centering\\arraybackslash}m{3.5cm}", 4), # wider columns, vertical alignment
  "}\n\\hline\nTrait1 & MR-PRESSO & IB\\_PRESSO(1) & IB\\_PRESSO(2) & IB\\_PRESSO(3)\\\\\n\\hline\n",
  paste(latex_rows, collapse="\n"),
  "\n\\hline\n\\end{tabular}\n"
)
cat(table_string)


###########################################################################################

library(dplyr)
library(tidyr)
library(stringr)

# --------- List category files ------------
base_dir <- "RDA_results/"  # IB-Mode RDA outputs (final_results_*.csv)
lipid_files      <- c("final_results_hdl.csv","final_results_ldl.csv","final_results_logTG.csv")
bp_files         <- c("final_results_sbp.csv","final_results_dbp.csv")
anthro_files     <- c("final_results_bmi.csv","final_results_BFP.csv","final_results_height.csv")
lifestyle_files  <- c("final_results_eversmok.csv","final_results_drinkpw.csv")
glucose_other_files <- c("final_results_FG.csv","final_results_vitD.csv","final_results_bw.csv")

categories <- list(
  "Lipids"        = lipid_files,
  "Blood Pressure"= bp_files,
  "Anthropometrics"= anthro_files,
  "Lifestyle"     = lifestyle_files,
  "Glucose / Other"= glucose_other_files
)

# --- Elegant Cell Formatter ---
format_cell <- function(b, se, trait2) {
  if (is.na(b) | is.na(se)) return("")
  lower <- b - 1.96 * se
  upper <- b + 1.96 * se
  trait2_line <- ifelse(trait2 == "NA" | is.na(trait2), "", trait2)
  sprintf("\\makecell{%.2f \\\\ (%.2f, %.2f) \\\\ %s}", b, lower, upper, trait2_line)
}

# --- Elegant Table function ---
make_category_table <- function(files, category_label) {
  dfs <- lapply(file.path(base_dir, files), function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$Exposure <- toupper(gsub("^final_results_|\\.csv$", "", basename(f)))
    d
  })
  df <- bind_rows(dfs)
  traits <- unique(df$trait1)
  latex_rows <- c()
  for (trait in traits) {
    sub <- df %>% filter(trait1 == trait)
    expo <- sub$Exposure[1]
    # MR-PRESSO row (main trait2 only)
    mrpresso_row <- sub %>% filter(method=="MR-PRESSO", trait2 == "NA" | is.na(trait2))
    mrpresso_cell <- if (nrow(mrpresso_row) > 0)
      format_cell(mrpresso_row$b[1], mrpresso_row$se[1], trait)
    else ""
    # IB-PRESSO top 3 (rank by abs(Z), show trait2)
    ib_presso_rows <- sub %>% filter(method=="IB-MR-PRESSO", trait2 != "NA" & !is.na(trait2)) %>%
      arrange(desc(abs(b/se)))
    ib_presso_cells <- rep("",3)
    for (i in 1:3) if(i <= nrow(ib_presso_rows)) {
      ir <- ib_presso_rows[i, ]
      ib_presso_cells[i] <- format_cell(ir$b, ir$se, ir$trait2)
    }
    # Weighted Mode row
    weighted_row <- sub %>% filter(method=="Weighted mode", trait2 == "NA" | is.na(trait2))
    weighted_cell <- if (nrow(weighted_row)>0)
      format_cell(weighted_row$b[1], weighted_row$se[1], trait)
    else ""
    # IB-MODE top 3 (rank by abs(Z), show trait2)
    ibmode_rows <- sub %>% filter(method %in% c("IB-MODE","new_method_weighted"), trait2 != "NA" & !is.na(trait2)) %>%
      arrange(desc(abs(b/se)))
    ibmode_cells <- rep("",3)
    for (i in 1:3) if(i <= nrow(ibmode_rows)) {
      ir <- ibmode_rows[i, ]
      ibmode_cells[i] <- format_cell(ir$b, ir$se, ir$trait2)
    }
    # Compose full row for LaTeX
    row <- c(expo, trait, mrpresso_cell, ib_presso_cells, weighted_cell, ibmode_cells)
    latex_rows <- c(latex_rows, paste(row, collapse = " & "), "\\\\")
  }
  # Elegant header
  header <- paste(
    "\\textbf{Exposure}", "\\textbf{Trait}",
    "\\textbf{MR-PRESSO}", "\\textbf{IB-PRESSO 1}", "\\textbf{IB-PRESSO 2}", "\\textbf{IB-PRESSO 3}",
    "\\textbf{Weighted Mode}", "\\textbf{IB-MODE 1}", "\\textbf{IB-MODE 2}", "\\textbf{IB-MODE 3}",
    sep = " & "
  )
  table_string <- paste0(
    "\n% ============================\n",
    "% ", category_label, " MR Estimates (Publication Table)\n",
    "% ============================\n",
    "\\begin{table}[ht]\n\\centering\n",
    "\\caption{Mendelian randomization estimates for ", category_label, " (beta, confidence interval, and component/intermediate trait)}\n",
    "\\scriptsize\n",
    "\\begin{tabular}{l l",
    strrep(">{\\centering\\arraybackslash}m{3.3cm}", 8),
    "}\n\\hline\n", header, " \\\\\n\\hline\n",
    paste(latex_rows, collapse = "\n"),
    "\n\\hline\n\\end{tabular}\n\\end{table}\n"
  )
  return(table_string)
}

# ----------- Loop across categories and print elegant LaTeX code -----------
for (catname in names(categories)) {
  tbl <- make_category_table(categories[[catname]], catname)
  cat(tbl)
}


############################################################################

library(dplyr)
library(tidyr)
library(stringr)

# --------- List category files ------------
base_dir <- "RDA_results_presso/"  # IB-PRESSO RDA outputs (final_results_*.csv)
lipid_files      <- c("final_results_hdl.csv","final_results_ldl.csv","final_results_logTG.csv")
bp_files         <- c("final_results_sbp.csv","final_results_dbp.csv")
anthro_files     <- c("final_results_bmi.csv","final_results_BFP.csv","final_results_height.csv")
lifestyle_files  <- c("final_results_eversmok.csv","final_results_drinkpw.csv")
glucose_other_files <- c("final_results_FG.csv","final_results_vitD.csv","final_results_bw.csv")

categories <- list(
  "Lipids"        = lipid_files,
  "Blood Pressure"= bp_files,
  "Anthropometrics"= anthro_files,
  "Lifestyle"     = lifestyle_files,
  "Glucose / Other"= glucose_other_files
)

mrpresso_method <- "MR-PRESSO"
ib_methods <- c("IB-MR-PRESSO", "new_method_weighted")

format_cell <- function(beta, se, subtrait=NULL, eff_str=NULL) {
  if (is.na(beta) | is.na(se)) return("")
  lower <- beta - 1.96*se
  upper <- beta + 1.96*se
  ci <- sprintf("%.2f (%.2f, %.2f)", beta, lower, upper)
  # Multi-line cell using makecell
  if (!is.null(subtrait) && !is.null(eff_str)) {
    return(sprintf("\\makecell{%s \\\\ (%s) (%s)}", ci, subtrait, eff_str))
  } else if (!is.null(subtrait)) {
    return(sprintf("\\makecell{%s \\\\ (%s)}", ci, subtrait))
  } else {
    return(ci)
  }
}

make_category_table <- function(files, category_label) {
  dfs <- lapply(file.path(base_dir, files), function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    # Exposure name from filename
    d$Exposure <- toupper(gsub("^final_results_|\\.csv$", "", basename(f)))
    d
  })
  df <- bind_rows(dfs)
  
  traits <- unique(df$trait1)
  latex_rows <- c()
  
  for (trait in traits) {
    sub <- df %>% filter(trait1 == trait)
    expo <- sub$Exposure[1]
    
    # MR-PRESSO: main effect row (trait2 NA)
    mrpresso_row <- sub %>% filter(method == mrpresso_method & (trait2 == "NA" | is.na(trait2)))
    if (nrow(mrpresso_row) > 0) {
      beta <- mrpresso_row$b[1]
      se <- mrpresso_row$se[1]
      mrpresso_str <- format_cell(beta, se)
      z_ref <- abs(beta/se)
    } else {
      mrpresso_str <- ""
      z_ref <- NA
    }
    
    # IB methods, trait2 non-NA: show up to 3, with both method and trait2 visible
    ib_rows <- sub %>% filter(method %in% ib_methods & trait2 != "NA" & !is.na(trait2)) %>%
      arrange(factor(method, levels=ib_methods), trait2)
    n_ib <- nrow(ib_rows)
    if (n_ib < 3) {
      ib_rows <- bind_rows(ib_rows, as.data.frame(matrix(NA, nrow=3-n_ib, ncol=ncol(ib_rows), 
                                                         dimnames=list(NULL, colnames(ib_rows)))))
    }
    ib_cells <- c()
    for (i in 1:3) {
      ir <- ib_rows[i, ]
      if (!is.na(ir$method)) {
        beta <- ir$b
        se <- ir$se
        subtrait <- ir$trait2
        z <- abs(beta/se)
        eff_str <- if (!is.na(z_ref) && !is.na(z)) {
          sprintf("%.0f\\%%", 100 * (z^2 - 1)/(z_ref^2 - 1))
        } else {
          ""
        }
        cell <- format_cell(beta, se, subtrait, eff_str)
      } else {
        cell <- ""
      }
      ib_cells <- c(ib_cells, cell)
    }
    
    row <- sprintf("%s & %s & %s & %s & %s & %s \\\\", expo, trait, mrpresso_str, ib_cells[1], ib_cells[2], ib_cells[3])
    latex_rows <- c(latex_rows, row)
  }
  
  table_string <- paste0(
    "\n% ============================\n",
    "% ", category_label, " MR Estimates (Publication Table)\n",
    "% ============================\n",
    "\\begin{table}[ht]\n\\centering\n",
    "\\caption{Mendelian randomization estimates for ", category_label, " (beta, confidence interval, and component/intermediate trait)}\n",
    "\\scriptsize\n",
    "\\begin{tabular}{l l",
    strrep(">{\\centering\\arraybackslash}m{3.5cm}", 4),
    "}\n\\hline\n",
    "Exposure & Trait & MR-PRESSO & IB\\_PRESSO(1) & IB\\_PRESSO(2) & IB\\_PRESSO(3) \\\\\n\\hline\n",
    paste(latex_rows, collapse = "\n"),
    "\n\\hline\n\\end{tabular}\n\\end{table}\n"
  )
  return(table_string)
}

# ----------- Loop across categories and print LaTeX code -----------
for (catname in names(categories)) {
  tbl <- make_category_table(categories[[catname]], catname)
  cat(tbl)
}


#############################################

library(dplyr)
library(stringr)

# ---------- CONFIG ----------
base_dir <- "RDA_results_presso/"  # IB-PRESSO RDA outputs (final_results_*.csv)
lipid_files      <- c("final_results_hdl.csv","final_results_ldl.csv","final_results_logTG.csv")
bp_files         <- c("final_results_sbp.csv","final_results_dbp.csv")
anthro_files     <- c("final_results_bmi.csv","final_results_BFP.csv","final_results_height.csv")
lifestyle_files  <- c("final_results_eversmok.csv","final_results_drinkpw.csv")
glucose_other_files <- c("final_results_FG.csv","final_results_vitD.csv","final_results_bw.csv")

categories <- list(
  "Lipids"          = lipid_files,
  "Blood Pressure"  = bp_files,
  "Anthropometrics" = anthro_files,
  "Lifestyle"       = lifestyle_files,
  "Glucose / Other" = glucose_other_files
)

# ---------- HELPERS ----------
clean_name <- function(x) {
  x <- gsub("_mvp$","",x)
  x <- gsub("_"," ",x)
  tools::toTitleCase(x)
}

format_main_cell <- function(beta,se) {
  if (is.na(beta) | is.na(se)) return("")
  lower <- beta - 1.96*se
  upper <- beta + 1.96*se
  sprintf("\\makecell{%.2f \\\\ (%.2f, %.2f)}", beta, lower, upper)
}

format_ib_cell <- function(beta,se,subtrait) {
  if (is.na(beta) | is.na(se)) return("")
  lower <- beta - 1.96*se
  upper <- beta + 1.96*se
  sprintf("\\makecell{%.2f \\\\ (%.2f, %.2f) \\\\ (%s)}",
          beta, lower, upper, subtrait)
}

make_comparison_table <- function(files, category_label, main_method, ib_method, table_title) {
  dfs <- lapply(file.path(base_dir, files), function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$Exposure <- toupper(gsub("^final_results_|\\.csv$", "", basename(f)))
    d
  })
  df <- bind_rows(dfs)
  df$trait1 <- clean_name(df$trait1)
  df$trait2 <- clean_name(df$trait2)
  df$Exposure <- clean_name(df$Exposure)
  
  # group order (file blocks)
  group_order <- unique(df$Exposure)
  latex_rows <- c()
  
  for (g in group_order) {
    subg <- df %>% filter(Exposure == g)
    traits <- unique(subg$trait1)
    
    latex_rows <- c(latex_rows, sprintf("\\multicolumn{5}{l}{\\textbf{%s}} \\\\", g))
    
    for (trait in traits) {
      sub <- subg %>% filter(trait1 == trait)
      
      # main method row
      main_row <- sub %>% filter(method == main_method & (is.na(trait2) | trait2=="NA"))
      if (nrow(main_row)>0) {
        main_cell <- format_main_cell(main_row$b[1], main_row$se[1])
      } else main_cell <- ""
      
      # only pull IB from the specified ib_method
      ib_rows <- sub %>% 
        filter(method == ib_method & !is.na(trait2) & trait2!="NA") %>%
        arrange(trait2)
      if (nrow(ib_rows)<3) {
        ib_rows <- bind_rows(
          ib_rows,
          as.data.frame(matrix(NA, nrow=3-nrow(ib_rows), ncol=ncol(ib_rows), 
                               dimnames=list(NULL, colnames(ib_rows))))
        )
      }
      ib_cells <- c()
      for (i in 1:3) {
        if (!is.na(ib_rows$method[i])) {
          ib_cells <- c(ib_cells, format_ib_cell(ib_rows$b[i], ib_rows$se[i], ib_rows$trait2[i]))
        } else {
          ib_cells <- c(ib_cells, "")
        }
      }
      row <- sprintf("%s & %s & %s & %s & %s \\\\",
                     trait, main_cell, ib_cells[1], ib_cells[2], ib_cells[3])
      latex_rows <- c(latex_rows, row)
    }
  }
  
  ib_header <- if (main_method=="Weighted mode") {
    c("IB-MODE1","IB-MODE2","IB-MODE3")
  } else {
    c("IB-PRESSO1","IB-PRESSO2","IB-PRESSO3")
  }
  
  paste0(
    "\\begin{table}[ht]\n\\centering\n",
    "\\caption{",table_title," for ",category_label,"}\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{1.3}\n",
    "\\begin{tabular}{l",
    strrep(">{\\centering\\arraybackslash}m{3.8cm}", 4),
    "}\n\\hline\n",
    "Outcome & ",main_method," & ",paste(ib_header, collapse=" & ")," \\\\\n\\hline\n",
    paste(latex_rows, collapse="\n"),
    "\n\\hline\n\\end{tabular}\n\\end{table}\n"
  )
}

# ---------- RUN ----------
for (catname in names(categories)) {
  cat(make_comparison_table(categories[[catname]], catname, 
                            "Weighted mode", "new_method_weighted",
                            "Weighted mode vs IB-MODE"))
  
  cat(make_comparison_table(categories[[catname]], catname, 
                            "MR-PRESSO", "IB-MR-PRESSO",
                            "MR-PRESSO vs IB-PRESSO"))
}
