# Final iSSF connectivity simulations based on RSF-defined core habitat.
# Four sex-season groups; 1000 trajectories each; 11 candidates per step.
# Each trajectory is simulated to 4x its base horizon; 1x and 2x are truncations.
# This script produces numerical/spatial analysis outputs, not manuscript figures.

suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(sf)
})
source(file.path("R", "00_config.R"))
source(file.path("R", "10_issf_common.R"))
source(file.path("R", "20_simulation_functions.R"))
create_sim_directories()

n_sim <- 1000L
n_candidates <- 11L
max_attempts <- 10L
seed <- 42L
core_min_area_km2 <- 93.1
expected_core_count <- 4L
expected_core_area_km2 <- 2638.42
expected_occupied_overlap_km2 <- 764.84

# Use --smoke-only for a rapid, sequential preflight on private inputs.
args <- commandArgs(trailingOnly=TRUE)
smoke_only <- "--smoke-only" %in% args
if (length(setdiff(args, "--smoke-only"))) stop("Unsupported command-line argument.")

required <- c(
  issf_data_file, issf_scaling_file,
  file.path(issf_result_dir, "population_beta_for_simulation.rds"),
  file.path(model_raster_dir, "core_suitable_binary_bins8_10_min93.1km2_100m.tif"),
  file.path(model_raster_dir, "core_suitable_patch_ID_bins8_10_min93.1km2_100m.tif"),
  file.path(study_area_dir, "mcp_camera_gps.shp")
)
if (!all(file.exists(required))) {
  stop("Missing simulation input(s):\n", paste(required[!file.exists(required)], collapse="\n"))
}
df <- readRDS(issf_data_file)
check_issf_data(df, strict=TRUE)
sc <- read_issf_scaling()
cache <- read_sim_environment()
beta_list <- readRDS(file.path(issf_result_dir, "population_beta_for_simulation.rds"))
if (!setequal(names(beta_list), issf_groups$key)) stop("Four sex-season beta vectors required.")
for (key in issf_groups$key) {
  if (!setequal(names(beta_list[[key]]), issf_terms) || any(!is.finite(beta_list[[key]]))) {
    stop("Invalid 17-term population beta for ", key)
  }
}

# New (final) RSF habitat definition. No legacy core habitat or legacy trajectories.
core_file <- required[4]
patch_file <- required[5]
core_r <- terra::rast(core_file)
patch_r <- terra::rast(patch_file)
if (!terra::compareGeom(core_r, patch_r, stopOnError=FALSE)) {
  stop("RSF core binary/patch rasters have different geometry.")
}
if (any(abs(terra::res(core_r) - 100) > 1e-6)) stop("Core raster resolution must be 100 m.")
if (!terra::same.crs(core_r, cache$template)) stop("Core/environment CRS mismatch.")
# Resample categorical patch IDs only, never bilinearly interpolate labels.
if (!terra::compareGeom(core_r, cache$template, stopOnError=FALSE)) {
  core_r <- terra::resample(core_r, cache$template, method="near")
  patch_r <- terra::resample(patch_r, cache$template, method="near")
}
core_values <- as.numeric(terra::values(core_r, mat=FALSE))
patch_cells <- as.integer(terra::values(patch_r, mat=FALSE))
core_idx <- which(is.finite(core_values) & core_values > 0)
patch_ids <- sort(unique(patch_cells[!is.na(patch_cells) & patch_cells > 0]))
if (length(patch_ids) != expected_core_count) stop("Expected four final RSF core patches.")
core_area <- length(core_idx) * prod(terra::res(cache$template)) / 1e6
if (abs(core_area - expected_core_area_km2) > 2) {
  stop("Core total area is inconsistent with final RSF definition: ", core_area)
}

mcp <- terra::vect(required[6])
if (!terra::same.crs(mcp, cache$template)) mcp <- terra::project(mcp, terra::crs(cache$template))
mcp_core <- terra::mask(core_r, mcp, touches=FALSE)
mcp_cells <- which(is.finite(terra::values(mcp_core, mat=FALSE)) &
                     as.numeric(terra::values(mcp_core, mat=FALSE)) > 0)
occupied_area <- length(mcp_cells) * prod(terra::res(cache$template)) / 1e6
if (abs(occupied_area - expected_occupied_overlap_km2) > 2) {
  stop("Occupied core area failed QC: ", occupied_area)
}
occupied_ids <- sort(unique(patch_cells[mcp_cells]))
occupied_ids <- occupied_ids[!is.na(occupied_ids) & occupied_ids > 0L]
target_ids <- setdiff(patch_ids, occupied_ids)
if (length(occupied_ids) != 1L || length(target_ids) != 3L) {
  stop("Expected one occupied core and three target cores.")
}

complete_env <- Reduce(`&`, lapply(cache$values, is.finite))
start_cells <- mcp_cells[complete_env[mcp_cells] & patch_cells[mcp_cells] %in% occupied_ids]
if (!length(start_cells)) stop("No valid simulation starting cells.")
start_xy <- terra::xyFromCell(cache$template, start_cells)
patch_qc <- data.frame(
  patch_id=patch_ids,
  n_cells=vapply(patch_ids, function(id) sum(patch_cells == id, na.rm=TRUE), integer(1)),
  occupied=patch_ids %in% occupied_ids,
  target=patch_ids %in% target_ids
)
patch_qc$area_km2 <- patch_qc$n_cells * prod(terra::res(cache$template)) / 1e6
write.csv(patch_qc, file.path(sim_summary_dir, "core_patch_qc.csv"), row.names=FALSE)
cat("Final RSF core QC passed; occupied IDs:", occupied_ids,
    "; target IDs:", paste(target_ids, collapse=", "), "\n")

# Sequential smoke tests for each sex-season; reject invalid starts/proposals.
for (g in seq_len(nrow(issf_groups))) {
  sex <- issf_groups$Sex[g]; season <- issf_groups$Season[g]
  movement <- movement_pool_for_group(df, sex, season)
  smoke <- simulate_one_track(
    start_xy[1, ], beta_list[[issf_groups$key[g]]], movement, cache, sc,
    n_steps=5L, n_candidates=n_candidates, max_attempts=max_attempts,
    seed=seed + g, patch_cells=patch_cells, target_ids=target_ids
  )
  if (nrow(smoke) < 2L) stop("Sequential simulation smoke test failed: ", issf_groups$key[g])
}
cat("All four sequential trajectory smoke tests passed.\n")
if (smoke_only) {
  cat("Smoke-only check complete; no full simulations run.\n")
} else {
  manifest <- list(
    n_sim=n_sim, n_candidates=n_candidates, max_attempts=max_attempts,
    seed=seed, horizons=simulation_horizons,
    core_min_area_km2=core_min_area_km2, occupied_ids=occupied_ids,
    target_ids=target_ids, start_n=length(start_cells),
    beta=beta_list, scaling=sc,
    core_file=basename(core_file), patch_file=basename(patch_file),
    input_md5=tools::md5sum(c(issf_data_file, core_file, patch_file))
  )
  manifest_file <- file.path(sim_out_dir, "simulation_manifest.rds")
  if (file.exists(manifest_file)) {
    old <- readRDS(manifest_file)
    if (!isTRUE(all.equal(old, manifest, check.attributes=FALSE))) {
      stop("Existing output was created with different simulation parameters. Use a fresh results directory.")
    }
  } else {
    existing_traj <- list.files(sim_trajectory_dir, pattern="^trajectory_.*\\.rds$", recursive=TRUE)
    if (length(existing_traj)) stop("Untracked old trajectory results found; do not mix simulation runs.")
    saveRDS(manifest, manifest_file)
  }

  summaries <- list()
  counter <- 1L
  for (g in seq_len(nrow(issf_groups))) {
    sex <- issf_groups$Sex[g]; season <- issf_groups$Season[g]
    key <- issf_groups$key[g]
    beta <- beta_list[[key]]
    movement <- movement_pool_for_group(df, sex, season)
    horizons <- simulation_horizons[[season]]
    max_steps <- max(horizons)
    out_dir <- file.path(sim_trajectory_dir, key)
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    set.seed(seed + g * 1000L)
    chosen_start <- sample(seq_len(nrow(start_xy)), size=n_sim, replace=TRUE)
    visit_counts <- integer(terra::ncell(cache$template))
    group_rows <- vector("list", n_sim * length(horizons))
    for (i in seq_len(n_sim)) {
      path <- file.path(out_dir, sprintf("trajectory_%04d.rds", i))
      if (file.exists(path)) {
        traj <- readRDS(path)
        if (!is.data.frame(traj) || !all(c("step", "x", "y", "patch_id") %in% names(traj))) {
          stop("Corrupted resumed trajectory: ", path)
        }
        if (abs(traj$x[1] - start_xy[chosen_start[i], 1]) > 1e-5 ||
            abs(traj$y[1] - start_xy[chosen_start[i], 2]) > 1e-5) {
          stop("Resumed trajectory start does not match the fixed simulation manifest.")
        }
      } else {
        traj <- simulate_one_track(
          start_xy[chosen_start[i], ], beta, movement, cache, sc,
          n_steps=max_steps, n_candidates=n_candidates,
          max_attempts=max_attempts, seed=seed + g * 10000000L + i,
          patch_cells=patch_cells, target_ids=target_ids
        )
        saveRDS(traj, path)
      }
      for (h in seq_along(horizons)) {
        h_max <- horizons[h]
        section <- traj[traj$step > 0 & traj$step <= h_max, , drop=FALSE]
        hit <- which(section$patch_id %in% target_ids)
        group_rows[[i * length(horizons) - length(horizons) + h]] <- data.frame(
          Sex=sex, Season=season, simulation_id=i, multiplier=c(1L, 2L, 4L)[h],
          requested_steps=h_max, completed_steps=max(traj$step),
          completed_within_horizon=min(max(traj$step), h_max),
          reached_target=length(hit) > 0L,
          first_arrival_step=if (length(hit)) section$step[hit[1]] else NA_integer_,
          first_target_id=if (length(hit)) section$patch_id[hit[1]] else NA_integer_
        )
      }
      primary <- traj[traj$step > 0 & traj$step <= horizons[1], c("x", "y"), drop=FALSE]
      if (nrow(primary)) {
        cells <- terra::cellFromXY(cache$template, as.matrix(primary))
        cells <- cells[!is.na(cells)]
        visit_counts <- visit_counts + tabulate(cells, nbins=length(visit_counts))
      }
      if (i %% 50L == 0L) cat(key, "completed", i, "of", n_sim, "trajectories\n")
    }
    group_table <- dplyr::bind_rows(group_rows)
    write.csv(group_table, file.path(sim_summary_dir, paste0(key, "_connectivity_by_trajectory.csv")), row.names=FALSE)
    summaries[[counter]] <- group_table
    counter <- counter + 1L
    visits_r <- cache$template
    terra::values(visits_r) <- visit_counts
    names(visits_r) <- "visit_count_primary_horizon"
    terra::writeRaster(visits_r, file.path(sim_raster_dir, paste0(key, "_visit_counts_100m.tif")),
                       overwrite=TRUE, datatype="INT4S")
  }
  all_tracks <- dplyr::bind_rows(summaries)
  summary <- all_tracks |>
    dplyr::group_by(.data$Sex, .data$Season, .data$multiplier, .data$requested_steps) |>
    dplyr::summarise(n_trajectories=dplyr::n(), n_arrived=sum(.data$reached_target),
                     proportion_arrived=mean(.data$reached_target),
                     mean_steps_completed=mean(.data$completed_within_horizon),
                     .groups="drop")
  if (any(summary$n_trajectories != n_sim)) stop("Simulation trajectory-count QC failed.")
  write.csv(all_tracks, file.path(sim_summary_dir, "all_trajectory_connectivity.csv"), row.names=FALSE)
  write.csv(summary, file.path(sim_summary_dir, "connectivity_summary.csv"), row.names=FALSE)
  cat("Connectivity simulations completed. Inspect arrival results and completion QC.\n")
}
