# Centralized paths and analysis settings for Shikoku landscape RSF.
# Run scripts from the repository root, or set RSF_PROJECT_ROOT.

project_dir <- normalizePath(Sys.getenv("RSF_PROJECT_ROOT", unset="."),
                             winslash="/", mustWork=TRUE)
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
