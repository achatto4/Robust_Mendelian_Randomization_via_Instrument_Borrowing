# Read MR-PRESSO's outlier-corrected estimate.
#
# MRPRESSO 1.0 fits the corrected model only inside
#     if (DISTORTIONtest & OUTLIERtest) { ... }
# so callers must pass DISTORTIONtest = TRUE for the "Outlier-corrected" row to
# be populated. Falls back to the "Raw" row when no correction was made.

.presso_corrected <- function(presso) {
  mm <- presso[["Main MR results"]]
  oc <- mm[mm[["MR Analysis"]] == "Outlier-corrected", , drop = FALSE]
  rw <- mm[mm[["MR Analysis"]] == "Raw", , drop = FALSE]
  b_oc <- suppressWarnings(as.numeric(oc[["Causal Estimate"]][1]))
  corrected <- length(b_oc) == 1L && !is.na(b_oc)
  row <- if (corrected) oc else rw
  oi <- presso[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]
  list(b          = suppressWarnings(as.numeric(row[["Causal Estimate"]][1])),
       se         = suppressWarnings(as.numeric(row[["Sd"]][1])),
       n_outliers = if (corrected && is.numeric(oi)) length(oi) else 0L,
       corrected  = corrected)
}
