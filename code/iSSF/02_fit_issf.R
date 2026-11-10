# Final 100-m iSSF: 17-term model fitted to 41 seasonal tracks.
# Pool tracks within individuals and individuals within sex-season groups.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})
source(file.path("code", "01_config.R"))
source(file.path("code", "iSSF", "01_issf_common.R"))
create_issf_directories()
if (!file.exists(issf_data_file)) stop("Missing restricted input: ", issf_data_file)

set.seed(20261102L)
df <- readRDS(issf_data_file)
qc <- check_issf_data(df, strict=TRUE)
write.csv(as.data.frame(qc), file.path(issf_result_dir, "sample_qc.csv"), row.names=FALSE)

track_index <- df |>
  distinct(Sex, Season, real_id, id_season) |>
  arrange(.data$Sex, .data$Season, .data$real_id, .data$id_season)
if (nrow(track_index) != 41L || anyDuplicated(track_index$id_season)) {
  stop("Expected 41 uniquely identified seasonal tracks.")
}

# Sequential smoke test before fitting all tracks.
first_track <- track_index$id_season[1]
smoke <- fit_issf_track(df[df$id_season == first_track, , drop=FALSE])
if (length(smoke$beta) != 17L) stop("Sequential model smoke test failed.")
cat("Sequential smoke test passed:", first_track, "\n")

track_fits <- setNames(vector("list", nrow(track_index)), track_index$id_season)
track_fits[[first_track]] <- smoke
for (i in seq_len(nrow(track_index))) {
  id_i <- track_index$id_season[i]
  if (is.null(track_fits[[id_i]])) {
    track_fits[[id_i]] <- fit_issf_track(df[df$id_season == id_i, , drop=FALSE])
  }
  cat(sprintf("Fitted %02d/%02d %s\n", i, nrow(track_index), id_i))
}

pooled <- pool_to_population(track_fits, track_index)
stopifnot(setequal(names(pooled$population), issf_groups$key))
for (key in issf_groups$key) {
  if (length(pooled$population[[key]]$beta) != 17L) stop("Incomplete group: ", key)
}

saveRDS(track_fits, file.path(issf_model_dir, "track_fits.rds"))
saveRDS(track_index, file.path(issf_model_dir, "track_index.rds"))
saveRDS(pooled$individual, file.path(issf_model_dir, "individual_fits.rds"))
saveRDS(pooled$population, file.path(issf_model_dir, "population_fits.rds"))
beta_list <- lapply(pooled$population, function(x) x$beta)
saveRDS(beta_list, file.path(issf_result_dir, "population_beta_for_simulation.rds"))
write.csv(make_coef_table(pooled$population),
          file.path(issf_result_dir, "population_coefficients_100m.csv"), row.names=FALSE)
write.csv(make_coef_table(pooled$individual),
          file.path(issf_result_dir, "individual_coefficients_100m.csv"), row.names=FALSE)

# Optional exact comparison with a previously saved reference analysis.
if (file.exists(issf_ref_beta_file)) {
  reference <- readRDS(issf_ref_beta_file)
  if (!setequal(names(reference), names(beta_list))) stop("Reference group mismatch.")
  differences <- vapply(names(beta_list), function(key) {
    if (!setequal(names(reference[[key]]), names(beta_list[[key]]))) stop("Reference term mismatch.")
    max(abs(as.numeric(reference[[key]][issf_terms]) - beta_list[[key]][issf_terms]))
  }, numeric(1))
  write.csv(data.frame(group=names(differences), max_abs_diff=as.numeric(differences)),
            file.path(issf_result_dir, "reference_coefficient_qc.csv"), row.names=FALSE)
  if (any(differences > 1e-6)) stop("Re-fitted beta differs from saved reference: see QC file.")
}
cat("iSSF fitting completed: 41 tracks and four population-level groups.\n")
