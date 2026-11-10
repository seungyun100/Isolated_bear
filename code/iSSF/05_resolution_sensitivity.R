# Optional raster-resolution sensitivity: 10, 50, 100, 200 m.
# Requires preprocessed, identical choice sets at each resolution.
# Data preparation from native rasters is not reimplemented here.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})
source(file.path("code", "01_config.R"))
source(file.path("code", "iSSF", "01_issf_common.R"))
create_issf_directories()

reference_fits_file <- file.path(issf_model_dir, "population_fits.rds")
reference_index_file <- file.path(issf_model_dir, "track_index.rds")
if (!all(file.exists(c(reference_fits_file, reference_index_file, issf_data_file)))) {
  stop("Run the final 100-m iSSF before the resolution sensitivity analysis.")
}
reference <- readRDS(reference_fits_file)
ref <- readRDS(issf_data_file)
check_issf_data(ref, strict=TRUE)
track_index <- readRDS(reference_index_file)

choice_key <- function(x) {
  paste(x$id_season, x$step_id_id_year, as.integer(x$case_), sep="::")
}
ref_counts <- table(choice_key(ref))

resolution_results <- list()
counter <- 1L
for (resolution in c(10L, 50L, 100L, 200L)) {
  path <- if (resolution == 100L) issf_data_file else {
    file.path(issf_private_dir, paste0("FINAL_iSSF_data_", resolution, "m_filtered.RDS"))
  }
  if (!file.exists(path)) {
    if (resolution == 100L) stop("Missing reference resolution data.")
    message("Skipping resolution ", resolution, " m; input unavailable: ", path)
    next
  }
  dat <- readRDS(path)
  if (resolution != 100L) {
    # Interaction names in issf_terms are model terms, not literal data columns.
    required_fields <- c("id_season", "step_id_id_year", "case_", "sl_", "log_sl_", "cos_ta_",
                         "z_end_ELE", "z_end_ELE2", "z_end_VRM", "z_end_VRM2", "z_end_PBUD",
                         "z_end_PDBF", "z_end_PDBF2", "z_log_end_DistMainR",
                         "z_log_end_DistMinorR", "z_log_end_DistEdgeAbs", "z_start_VRM",
                         "z_start_PDBF", "z_log_start_DistMinorR")
    if (!all(required_fields %in% names(dat))) stop("Missing predictor columns at ", resolution, " m.")
    if (nrow(dat) != nrow(ref) || !identical(table(choice_key(dat)), ref_counts)) {
      stop("Choice sets differ at ", resolution, " m. Resolution comparisons require identical steps.")
    }
  }
  if (resolution == 100L) {
    pooled <- reference
  } else {
    # First fit one track as sequential smoke test.
    first_track <- track_index$id_season[1]
    smoke <- fit_issf_track(dat[dat$id_season == first_track, , drop=FALSE])
    fits <- setNames(vector("list", nrow(track_index)), track_index$id_season)
    fits[[first_track]] <- smoke
    for (i in seq_len(nrow(track_index))) {
      id_i <- track_index$id_season[i]
      if (is.null(fits[[id_i]])) fits[[id_i]] <- fit_issf_track(dat[dat$id_season == id_i, , drop=FALSE])
    }
    pooled <- pool_to_population(fits, track_index)$population
  }
  for (key in issf_groups$key) {
    r <- reference[[key]]$beta[issf_terms]
    b <- pooled[[key]]$beta[issf_terms]
    resolution_results[[counter]] <- data.frame(
      resolution_m=resolution, group=key, term=issf_terms,
      beta=as.numeric(b), beta_reference_100m=as.numeric(r),
      abs_difference=abs(b-r), direction_agrees=sign(b) == sign(r),
      relative_change=abs(b-r) / pmax(abs(r), 1e-8), row.names=NULL
    )
    counter <- counter + 1L
  }
  cat("Completed", resolution, "m resolution.\n")
}
if (!length(resolution_results)) stop("No sensitivity data available.")
comparison <- dplyr::bind_rows(resolution_results)
write.csv(comparison, file.path(issf_result_dir, "resolution_sensitivity_comparison.csv"), row.names=FALSE)
cat("Sensitivity comparison saved. Missing optional resolutions were not fabricated.\n")
