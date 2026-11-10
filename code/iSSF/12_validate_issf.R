# Leave-one-individual-out (LOIO) simulation validation of the final 100-m iSSF.
# Uses training-individual coefficients and empirical movement steps.
# No manuscript or supplementary figures are generated here.

suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
})
source(file.path("R", "00_config.R"))
source(file.path("R", "10_issf_common.R"))
source(file.path("R", "20_simulation_functions.R"))
create_issf_directories()

input_files <- c(
  issf_data_file, issf_scaling_file,
  file.path(issf_model_dir, "individual_fits.rds"),
  file.path(issf_model_dir, "track_index.rds")
)
if (!all(file.exists(input_files))) stop("Run R/11_fit_issf.R and provide private inputs first.")
df <- readRDS(issf_data_file)
check_issf_data(df, strict=TRUE)
sc <- read_issf_scaling()
cache <- read_sim_environment()
individual <- readRDS(file.path(issf_model_dir, "individual_fits.rds"))
track_index <- readRDS(file.path(issf_model_dir, "track_index.rds"))

n_trajectories <- 200L
n_null <- 999L
n_candidates <- 11L
max_attempts <- 10L
seed <- 20261102L
used <- df[df$case_ == 1L, , drop=FALSE]

# Precompute positive visit-intensity rank bins once per held-out track.
# Zero-visit cells are a separate class, excluded from the positive-bin rho.
make_visit_rank <- function(visits) {
  positive <- which(visits > 0L)
  if (!length(positive)) return(NULL)
  bins <- integer(length(visits))
  bins[positive] <- pmax(1L, pmin(10L, ceiling(
    rank(visits[positive], ties.method="average") / length(positive) * 10
  )))
  list(bins=bins, area=tabulate(bins[positive], nbins=10L))
}

area_adjusted_rho <- function(rank_object, test_cells) {
  if (is.null(rank_object)) return(NA_real_)
  bins <- rank_object$bins
  hits <- test_cells[!is.na(test_cells) & test_cells >= 1L & test_cells <= length(bins)]
  if (!length(hits)) return(NA_real_)
  use <- tabulate(bins[hits], nbins=10L)
  ok <- rank_object$area > 0L
  if (sum(ok) < 3L) return(NA_real_)
  ratio <- use[ok] / rank_object$area[ok]
  if (stats::sd(ratio) == 0) return(NA_real_)
  suppressWarnings(stats::cor(which(ok), ratio, method="spearman"))
}

make_visit_counts <- function(simulations, template) {
  all_cells <- unlist(lapply(simulations, function(traj) {
    if (nrow(traj) < 2L) return(integer(0))
    cells <- terra::cellFromXY(template, as.matrix(traj[-1L, c("x", "y")]))
    cells[!is.na(cells)]
  }), use.names=FALSE)
  tabulate(all_cells, nbins=terra::ncell(template))
}

# Sequential one-track/one-trajectory smoke test before the full validation.
first <- track_index[1, , drop=FALSE]
first_individuals <- unique(track_index$real_id[
  track_index$Sex == first$Sex & track_index$Season == first$Season &
    track_index$real_id != first$real_id
])
first_fit <- pool_ivw(individual[paste(first$Sex, first$Season, first_individuals, sep="__")])
first_obs <- used[used$id_season == first$id_season, , drop=FALSE]
first_obs <- first_obs[order(first_obs$t1_), , drop=FALSE]
first_pool <- movement_pool_for_group(df, first$Sex, first$Season, first$real_id)
smoke <- simulate_one_track(
  c(first_obs$x1_[1], first_obs$y1_[1]), first_fit$beta, first_pool,
  cache, sc, n_steps=3L, seed=seed
)
if (nrow(smoke) < 2L) stop("LOIO one-track simulation smoke test failed.")
cat("LOIO sequential smoke test passed.\n")

qc_rows <- vector("list", nrow(track_index))
for (i in seq_len(nrow(track_index))) {
  row <- track_index[i, , drop=FALSE]
  sex <- row$Sex; season <- row$Season; held <- row$real_id
  group_ids <- unique(track_index$real_id[
    track_index$Sex == sex & track_index$Season == season & track_index$real_id != held
  ])
  if (length(group_ids) < 1L) stop("LOIO fold has no training individuals.")
  keys <- paste(sex, season, group_ids, sep="__")
  beta <- pool_ivw(individual[keys])$beta
  movement <- movement_pool_for_group(df, sex, season, held)
  track_used <- used[used$id_season == row$id_season, , drop=FALSE]
  track_used <- track_used[order(track_used$t1_), , drop=FALSE]
  n_steps <- nrow(track_used)
  start <- c(track_used$x1_[1], track_used$y1_[1])
  simulations <- vector("list", n_trajectories)
  for (j in seq_len(n_trajectories)) {
    simulations[[j]] <- simulate_one_track(
      start, beta, movement, cache, sc, n_steps=n_steps,
      n_candidates=n_candidates, max_attempts=max_attempts,
      seed=seed + i * 100000L + j
    )
  }
  completion <- vapply(simulations, function(x) max(x$step) / n_steps, numeric(1))
  visits <- make_visit_counts(simulations, cache$template)
  rank_object <- make_visit_rank(visits)
  obs_xy <- as.matrix(track_used[, c("x2_", "y2_")])
  obs_cells <- terra::cellFromXY(cache$template, obs_xy)
  observed_rho <- area_adjusted_rho(rank_object, obs_cells)

  max_distance <- max(sqrt((obs_xy[, 1] - start[1])^2 +
                           (obs_xy[, 2] - start[2])^2), na.rm=TRUE)
  null_rho <- rep(NA_real_, n_null)
  if (is.finite(max_distance) && max_distance > 0) {
    set.seed(seed + i * 100000L + 99000L)
    for (b in seq_len(n_null)) {
      radius <- max_distance * sqrt(runif(n_steps))
      theta <- runif(n_steps, -pi, pi)
      null_xy <- cbind(start[1] + radius * cos(theta),
                       start[2] + radius * sin(theta))
      null_cells <- terra::cellFromXY(cache$template, null_xy)
      null_rho[b] <- area_adjusted_rho(rank_object, null_cells)
    }
  }
  valid_null <- null_rho[is.finite(null_rho)]
  empirical_p <- if (is.finite(observed_rho) && length(valid_null)) {
    (1 + sum(valid_null >= observed_rho)) / (1 + length(valid_null))
  } else NA_real_
  qc_rows[[i]] <- data.frame(
    Sex=sex, Season=season, heldout_id=held, id_season=row$id_season,
    n_training_individuals=length(group_ids), n_test_steps=n_steps,
    n_sim=n_trajectories, n_full=sum(completion == 1),
    mean_completion=mean(completion), min_completion=min(completion),
    spearman_rho=observed_rho, n_valid_null=length(valid_null), empirical_p=empirical_p
  )
  cat(sprintf("LOIO %02d/%02d: %s, rho=%s\n", i, nrow(track_index),
              row$id_season, format(observed_rho, digits=3)))
  write.csv(dplyr::bind_rows(qc_rows[seq_len(i)]),
            file.path(issf_loio_dir, "loio_validation_incremental.csv"), row.names=FALSE)
}
qc <- dplyr::bind_rows(qc_rows)
write.csv(qc, file.path(issf_loio_dir, "loio_validation_by_track.csv"), row.names=FALSE)
summary <- qc |>
  dplyr::group_by(.data$Sex, .data$Season) |>
  dplyr::summarise(mean_rho=mean(.data$spearman_rho, na.rm=TRUE),
                   n_tracks=dplyr::n(), n_valid=sum(is.finite(.data$spearman_rho)),
                   .groups="drop")
write.csv(summary, file.path(issf_loio_dir, "loio_validation_summary.csv"), row.names=FALSE)
cat("LOIO simulation validation completed.\n")
