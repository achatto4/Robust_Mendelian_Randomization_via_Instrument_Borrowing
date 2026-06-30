#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "../../../results/simulation_results/mode_comp"
alpha_arg <- if (length(args) >= 2) args[2] else "0.001,0.005,0.01,0.05,0.1"
update_rda <- if (length(args) >= 3) tolower(args[3]) %in% c("1", "true", "yes", "y") else FALSE

alpha_levels <- unique(suppressWarnings(as.numeric(strsplit(alpha_arg, ",")[[1]])))
alpha_levels <- sort(alpha_levels[is.finite(alpha_levels) & alpha_levels > 0 & alpha_levels <= 1])
if (length(alpha_levels) == 0) {
  stop("No valid alpha levels provided. Example: 0.001,0.01,0.05")
}

rda_files <- list.files(
  path = results_dir,
  pattern = "NULL.*\\.rda$",
  full.names = TRUE
)

if (length(rda_files) == 0) {
  stop("No NULL .rda files found in: ", results_dir)
}

message("Rewriting type-I summary CSVs from .rda files in: ", normalizePath(results_dir, mustWork = FALSE))
message("Alpha levels: ", paste(alpha_levels, collapse = ", "))

for (f in rda_files) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)

  if (!exists("est", envir = env, inherits = FALSE)) {
    warning("Skipping (missing 'est'): ", basename(f))
    next
  }

  est <- get("est", envir = env, inherits = FALSE)
  if (!(is.matrix(est) || is.data.frame(est))) {
    warning("Skipping (est is not matrix/data.frame): ", basename(f))
    next
  }
  est <- as.matrix(est)

  cn <- colnames(est)
  idx_mode <- if (!is.null(cn) && "MRMode" %in% cn) match("MRMode", cn) else 4
  idx_ib <- if (!is.null(cn) && "mode_new" %in% cn) match("mode_new", cn) else 5
  idx_mode_se <- if (!is.null(cn) && "MRMode_se" %in% cn) match("MRMode_se", cn) else 6
  idx_ib_se <- if (!is.null(cn) && "mode_new_se" %in% cn) match("mode_new_se", cn) else 7

  mode_p <- rep(NA_real_, nrow(est))
  ib_mode_p <- rep(NA_real_, nrow(est))

  valid_mode <- !is.na(est[, idx_mode]) & !is.na(est[, idx_mode_se]) & est[, idx_mode_se] > 0
  valid_ib <- !is.na(est[, idx_ib]) & !is.na(est[, idx_ib_se]) & est[, idx_ib_se] > 0

  mode_p[valid_mode] <- 2 * pnorm(-abs(est[valid_mode, idx_mode] / est[valid_mode, idx_mode_se]))
  ib_mode_p[valid_ib] <- 2 * pnorm(-abs(est[valid_ib, idx_ib] / est[valid_ib, idx_ib_se]))

  meta_cols <- c("scenario", "theta", "theta_alt", "phi", "thetaU", "N", "prop_invalid", "overlap")
  meta <- as.list(setNames(rep(NA, length(meta_cols)), meta_cols))
  n_rep_val <- nrow(est)

  if (exists("type1_summary", envir = env, inherits = FALSE)) {
    ts0 <- get("type1_summary", envir = env, inherits = FALSE)
    if (is.data.frame(ts0) && nrow(ts0) > 0) {
      for (nm in meta_cols) {
        if (nm %in% names(ts0)) meta[[nm]] <- ts0[[nm]][1]
      }
      if ("n_rep" %in% names(ts0) && is.finite(as.numeric(ts0$n_rep[1]))) {
        n_rep_val <- as.numeric(ts0$n_rep[1])
      }
    }
  }

  type1_summary_new <- do.call(rbind, lapply(alpha_levels, function(a) {
    data.frame(
      scenario = meta$scenario,
      theta = as.numeric(meta$theta),
      theta_alt = as.numeric(meta$theta_alt),
      phi = as.numeric(meta$phi),
      alpha = a,
      n_rep = n_rep_val,
      n_valid_mode = sum(!is.na(mode_p)),
      n_valid_ib_mode = sum(!is.na(ib_mode_p)),
      type1_mrmode = mean(mode_p < a, na.rm = TRUE),
      type1_ib_mode = mean(ib_mode_p < a, na.rm = TRUE),
      thetaU = as.numeric(meta$thetaU),
      N = as.numeric(meta$N),
      prop_invalid = as.numeric(meta$prop_invalid),
      overlap = as.numeric(meta$overlap)
    )
  }))

  out_csv <- sub("\\.rda$", "_type1_summary.csv", f)
  write.csv(type1_summary_new, out_csv, row.names = FALSE)

  if (update_rda) {
    assign("type1_summary", type1_summary_new, envir = env)
    obj_names <- ls(envir = env, all.names = TRUE)
    save(list = obj_names, file = f, envir = env)
  }

  message("Rewrote: ", basename(out_csv))
}

message("Done.")
