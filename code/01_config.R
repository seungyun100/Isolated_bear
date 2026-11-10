# Centralized paths for Shikoku RSF, iSSF and connectivity analysis.
# Run from the repository root, or set SHIKOKU_PROJECT_ROOT.

project_dir <- normalizePath(
  Sys.getenv("SHIKOKU_PROJECT_ROOT", unset=Sys.getenv("RSF_PROJECT_ROOT", unset=".")),
  winslash="/", mustWork=TRUE
)
private_dir <- file.path(project_dir, "data_private")
native_var_dir <- file.path(private_dir, "variables")
study_area_dir <- file.path(private_dir, "study_area")
gps_source_file <- file.path(private_dir, "gps", "df_step_all.RDS")

rsf_output_dir <- file.path(project_dir, "results", "rsf")
base_dir <- file.path(rsf_output_dir, "01_preparation")
base_raster_dir <- file.path(base_dir, "rasters")
selected_raster_dir <- file.path(base_raster_dir, "selected_predictors")
temp_focal_dir <- file.path(base_raster_dir, "temp_focal_candidates")
base_vector_dir <- file.path(base_dir, "vectors")
base_table_dir <- file.path(base_dir, "tables")
base_rds_dir <- file.path(base_dir, "rds")

model_dir <- file.path(rsf_output_dir, "02_model")
model_raster_dir <- file.path(model_dir, "rasters")
model_vector_dir <- file.path(model_dir, "vectors")
model_table_dir <- file.path(model_dir, "tables")
model_rds_dir <- file.path(model_dir, "rds")
temp_full_dir <- file.path(model_dir, "temp_full_predictors")
figure_dir <- file.path(rsf_output_dir, "figures")

create_rsf_directories <- function() {
  directories <- c(base_dir, base_raster_dir, selected_raster_dir,
                   temp_focal_dir, base_vector_dir, base_table_dir, base_rds_dir,
                   model_dir, model_raster_dir, model_vector_dir,
                   model_table_dir, model_rds_dir, temp_full_dir, figure_dir)
  invisible(lapply(directories, dir.create, recursive=TRUE, showWarnings=FALSE))
}

# iSSF and connectivity simulation paths.
issf_private_dir <- file.path(private_dir, "issf")
issf_data_file <- file.path(issf_private_dir, "FINAL_iSSF_data_100m_filtered.RDS")
issf_scaling_file <- file.path(issf_private_dir, "FINAL_scaling_parameters_reference_10m.RDS")
issf_ref_coef_file <- file.path(issf_private_dir, "FINAL_population_coefficients_100m.RDS")
issf_ref_beta_file <- file.path(issf_private_dir, "FINAL_iSSF_simulation_beta_100m.RDS")
issf_out_dir <- file.path(project_dir, "results", "issf")
issf_model_dir <- file.path(issf_out_dir, "models")
issf_result_dir <- file.path(issf_out_dir, "results")
issf_figure_dir <- file.path(issf_out_dir, "figures")
issf_loio_dir <- file.path(issf_out_dir, "validation")
sim_env_dir <- file.path(private_dir, "simulation", "environment_100m")
sim_out_dir <- file.path(project_dir, "results", "simulation")
sim_trajectory_dir <- file.path(sim_out_dir, "trajectories")
sim_summary_dir <- file.path(sim_out_dir, "summaries")
sim_raster_dir <- file.path(sim_out_dir, "rasters")

create_issf_directories <- function() {
  directories <- c(issf_model_dir, issf_result_dir, issf_figure_dir, issf_loio_dir)
  invisible(lapply(directories, dir.create, recursive=TRUE, showWarnings=FALSE))
}
create_sim_directories <- function() {
  directories <- c(sim_out_dir, sim_trajectory_dir, sim_summary_dir, sim_raster_dir)
  invisible(lapply(directories, dir.create, recursive=TRUE, showWarnings=FALSE))
}
