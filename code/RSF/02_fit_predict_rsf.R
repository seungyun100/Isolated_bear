# Shikoku landscape RSF | 02: model selection and habitat mapping
#
# PURPOSE
#   1. Load Script-01 outputs
#   2. Extract predictors for fixed 20,000-cell sample
#   3. Fit and evaluate landscape RSF
#   4. Compare linear vs quadratic DBF
#   5. Re-evaluate VRM and Edge with quadratic DBF fixed
#   6. Select final model
#   7. Predict relative landscape suitability across Shikoku
#   8. Create ABC-calibrated equal-area suitability classes
#   9. Evaluate GPS distribution and bootstrap threshold
#  10. Define suitable landscape and >=93.1 km2 core patches
#
# IMPORTANT
#   - DBF focal scale was selected linearly in Script 01.
#   - DBF functional form is evaluated here.
#   - DistEdge structure:
#       absolute distance
#       outside-forest indicator
#       absolute distance x outside interaction
#   - All outputs are saved under:
#       results/rsf/02_model
# 1. Packages
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(car)

gc()

# Project-relative paths (configure in code/01_config.R)
source(file.path("code", "01_config.R"))
create_rsf_directories()

# 3. Load Script-01 Stage-1 objects
#
# IMPORTANT:
# SpatRaster objects saved inside RDS are NOT reused.
# template_100m is reconstructed from TIFF.
stage1_file <- file.path(base_rds_dir, "RSF_stage1_objects.rds")

stopifnot(
  file.exists(stage1_file)
)

stage1 <- readRDS(stage1_file)

required_stage1_names <- c(
  "A",
  "B",
  "C",
  "ABC",
  "sample_base_100m",
  "focal_file_table",
  "best_scales_100m",
  "selected_focal_table_100m",
  "DistMainR_file",
  "DistMinorR_file",
  "DistEdge_file"
)

if (!all(required_stage1_names %in% names(stage1))) {
  stop("Stage-1 RDS is missing required objects.")
}

A <- stage1$A
B <- stage1$B
C <- stage1$C
ABC <- stage1$ABC

sample_base_100m <- stage1$sample_base_100m
focal_file_table <- stage1$focal_file_table
best_scales_100m <- stage1$best_scales_100m
selected_focal_table_100m <- stage1$selected_focal_table_100m

stopifnot(
  nrow(sample_base_100m) == 20000
)

template_100m_file <- file.path(base_raster_dir, "template_100m.tif")

stopifnot(
  file.exists(template_100m_file)
)

template_100m <- rast(template_100m_file)

cat("\nTemplate reconstructed from disk:\n")
print(template_100m)

cat(
  "Fixed sample rows:",
  nrow(sample_base_100m),
  "\n"
)

# 4. Predictor file paths
get_selected_file <- function(variable) {
  f <- selected_focal_table_100m$file[
    selected_focal_table_100m$variable == variable
  ]

  if (length(f) != 1) {
    stop(
      paste("Selected predictor file not unique:", variable)
    )
  }

  if (!file.exists(f)) {
    stop(
      paste("Selected predictor file missing:", f)
    )
  }

  f
}

BUD_file_100m <- get_selected_file("BUD")
DBF_file_100m <- get_selected_file("DBF")
ENF_file_100m <- get_selected_file("ENF")
VRM_file_100m <- get_selected_file("VRM")
ELE_file_100m <- get_selected_file("ELE")

DistMainR_file_100m <- stage1$DistMainR_file
DistMinorR_file_100m <- stage1$DistMinorR_file
DistEdge_file_100m <- stage1$DistEdge_file

predictor_files <- c(
  BUD=BUD_file_100m,
  DBF=DBF_file_100m,
  ENF=ENF_file_100m,
  VRM=VRM_file_100m,
  ELE=ELE_file_100m,
  DistMainR=DistMainR_file_100m,
  DistMinorR=DistMinorR_file_100m,
  DistEdge_signed=DistEdge_file_100m
)

if (!all(file.exists(predictor_files))) {
  print(
    predictor_files[
      !file.exists(predictor_files)
    ]
  )

  stop("Some predictor files are missing.")
}

cat("\nPredictor files:\n")
print(predictor_files)

# 5. Convert fixed sample cells to raster row / column
sample_rc <- rowColFromCell(template_100m, sample_base_100m$cell)

stopifnot(
  nrow(sample_rc) == 20000,
  ncol(sample_rc) == 2
)

sample_rows <- sample_rc[, 1]
sample_cols <- sample_rc[, 2]

cat(
  "\nSample row/column indices prepared:",
  nrow(sample_rc),
  "cells\n"
)

# 6. Memory-safe raster extraction function
#
# Reads only raster rows containing sampled cells.
# Does NOT load the full raster into RAM.
extract_cells_by_rows <- function(file, variable, sample_rows, sample_cols, template) {
  cat("\n========================================\n")
  cat("Extracting:", variable, "\n")
  cat("File:", file, "\n")

  if (!file.exists(file)) {
    stop(
      paste("File not found:", file)
    )
  }

  rr <- rast(file)

  geom_ok <- compareGeom(rr, template, stopOnError=FALSE)

  cat("Geometry match:", geom_ok, "\n")

  if (!geom_ok) {
    stop(
      paste("Geometry mismatch:", variable)
    )
  }

  unique_rows <- sort(
    unique(sample_rows)
  )

  out <- rep(
    NA_real_,
    length(sample_rows)
  )

  cat(
    "Raster rows to read:",
    length(unique_rows),
    "of",
    nrow(rr),
    "\n"
  )

  terra::readStart(rr)

  for (i in seq_along(unique_rows)) {
    row_i <- unique_rows[i]

    row_values <- terra::readValues(rr, row=row_i, nrows=1, mat=FALSE)

    idx <- which(sample_rows == row_i)

    out[idx] <- row_values[
      sample_cols[idx]
    ]

    if (i %% 100 == 0 || i == length(unique_rows)) {
      cat(
        "  completed:",
        i,
        "/",
        length(unique_rows),
        "rows\n"
      )
    }
  }

  terra::readStop(rr)

  rm(rr)
  gc()

  if (length(out) != length(sample_rows)) {
    stop(
      paste("Unexpected extracted length:", variable)
    )
  }

  cat(
    "NA count:",
    sum(is.na(out)),
    "\n"
  )

  if (anyNA(out)) {
    stop(
      paste("NA detected:", variable)
    )
  }

  cat("Completed:", variable, "\n")

  out
}

# 7. Extract predictors for fixed 20,000-cell sample
sample_model_100m <- sample_base_100m

sample_model_100m$BUD <- extract_cells_by_rows(
  predictor_files[["BUD"]],
  "BUD",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$DBF <- extract_cells_by_rows(
  predictor_files[["DBF"]],
  "DBF",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$ENF <- extract_cells_by_rows(
  predictor_files[["ENF"]],
  "ENF",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$VRM <- extract_cells_by_rows(
  predictor_files[["VRM"]],
  "VRM",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$ELE <- extract_cells_by_rows(
  predictor_files[["ELE"]],
  "ELE",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$DistMainR <- extract_cells_by_rows(
  predictor_files[["DistMainR"]],
  "DistMainR",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$DistMinorR <- extract_cells_by_rows(
  predictor_files[["DistMinorR"]],
  "DistMinorR",
  sample_rows,
  sample_cols,
  template_100m
)

sample_model_100m$DistEdge_signed <- extract_cells_by_rows(
  predictor_files[["DistEdge_signed"]],
  "DistEdge_signed",
  sample_rows,
  sample_cols,
  template_100m
)

stopifnot(
  nrow(sample_model_100m) == 20000
)

na_summary <- colSums(
  is.na(sample_model_100m)
)

print(na_summary)

if (any(na_summary > 0)) {
  stop("NA detected in sample_model_100m.")
}

cat("\nAll predictors extracted successfully.\n")

# 8. Save extracted sample IMMEDIATELY
#
# If RStudio is restarted later, load this file instead of
# repeating raster extraction.
sample_predictor_file <- file.path(model_rds_dir, "sample_predictors_100m.rds")

saveRDS(sample_model_100m, sample_predictor_file)

write.csv(
  sample_model_100m,
  file.path(model_table_dir, "sample_predictors_100m.csv"),
  row.names=FALSE
)

cat("Saved model sample:\n", sample_predictor_file, "\n")

# 9. Transform road and edge distances
sample_model_100m <- sample_model_100m %>%
  mutate(
    log_DistMainR=log1p(.data$DistMainR),
    log_DistMinorR=log1p(.data$DistMinorR),
    abs_DistEdge=abs(.data$DistEdge_signed)
  )

# 10. Determine DistEdge sign using LULC
#
# Forest-center status is reconstructed over the whole
# 100-m template, not only the forest-dominated A/C sample.
lulc_file <- file.path(native_var_dir, "lulc_1820_shikoku.tif")

stopifnot(
  file.exists(lulc_file)
)

lulc_native <- rast(lulc_file)

if (!same.crs(lulc_native, template_100m)) {
  stop("LULC CRS does not match template_100m.")
}

lulc_crop <- crop(lulc_native, template_100m, snap="out")

lulc_center_100m <- resample(lulc_crop, template_100m, method="near", threads=TRUE)

forest_classes <- c(6, 7, 8, 9, 11)

forest_center_100m <- lulc_center_100m %in% forest_classes

DistEdge_sign_raster <- rast(DistEdge_file_100m)

positive_cells <- ifel(DistEdge_sign_raster > 0, 1, NA)

negative_cells <- ifel(DistEdge_sign_raster < 0, 1, NA)

positive_forest_cells <- ifel(DistEdge_sign_raster > 0 & forest_center_100m, 1, NA)

negative_forest_cells <- ifel(DistEdge_sign_raster < 0 & forest_center_100m, 1, NA)

n_positive <- global(positive_cells, "sum", na.rm=TRUE)[1, 1]

n_negative <- global(negative_cells, "sum", na.rm=TRUE)[1, 1]

n_positive_forest <- global(positive_forest_cells, "sum", na.rm=TRUE)[1, 1]

n_negative_forest <- global(negative_forest_cells, "sum", na.rm=TRUE)[1, 1]

positive_forest_prop <- n_positive_forest / n_positive
negative_forest_prop <- n_negative_forest / n_negative

edge_sign_QC <- data.frame(
  edge_sign=c("positive", "negative"),
  n_cells=c(n_positive, n_negative),
  forest_cells=c(n_positive_forest, n_negative_forest),
  forest_proportion=c(positive_forest_prop, negative_forest_prop)
)

print(edge_sign_QC)

write.csv(
  edge_sign_QC,
  file.path(model_table_dir, "DistEdge_sign_LULC_QC.csv"),
  row.names=FALSE
)

if (!all(is.finite(edge_sign_QC$forest_proportion))) {
  stop("Could not determine DistEdge sign.")
}

if (positive_forest_prop > negative_forest_prop) {
  inside_sign <- "positive"
  outside_sign <- "negative"
} else {
  inside_sign <- "negative"
  outside_sign <- "positive"
}

cat("\nDistEdge inside sign:", inside_sign, "\n")

cat("DistEdge outside sign:", outside_sign, "\n")

sample_model_100m <- sample_model_100m %>%
  mutate(
    outside_forest=case_when(
      .data$DistEdge_signed == 0 ~ 0L,
      .env$outside_sign == "positive" & .data$DistEdge_signed > 0 ~ 1L,
      .env$outside_sign == "negative" & .data$DistEdge_signed < 0 ~ 1L,
      TRUE ~ 0L
    )
  )

print(
  table(sample_model_100m$outside_forest)
)

write.csv(
  data.frame(inside_sign=inside_sign, outside_sign=outside_sign),
  file.path(model_table_dir, "DistEdge_sign_definition.csv"),
  row.names=FALSE
)

# 11. Correlation among candidate environmental variables
cor_vars_100m <- c("BUD", "DBF", "ENF", "VRM", "ELE", "log_DistMainR", "log_DistMinorR", "abs_DistEdge")

cor_all_100m <- cor(sample_model_100m[, cor_vars_100m], use="complete.obs", method="pearson")

print(
  round(cor_all_100m, 2)
)

write.csv(
  round(cor_all_100m, 4),
  file.path(model_table_dir, "correlation_all_predictors_100m.csv"),
  row.names=TRUE
)

cor_idx <- which(
  abs(cor_all_100m) > 0.6 &
    abs(cor_all_100m) < 1,
  arr.ind=TRUE
)

high_cor_pairs_100m <- data.frame(
  var1=rownames(cor_all_100m)[cor_idx[, 1]],
  var2=colnames(cor_all_100m)[cor_idx[, 2]],
  r=cor_all_100m[cor_idx]
) %>%
  filter(.data$var1 < .data$var2) %>%
  arrange(
    desc(
      abs(.data$r)
    )
  )

print(high_cor_pairs_100m)

write.csv(
  high_cor_pairs_100m,
  file.path(model_table_dir, "high_correlation_pairs_100m.csv"),
  row.names=FALSE
)

# 12. Scaling parameters
#
# Scaling is based only on the fixed 20,000-cell model sample.
scale_vars_100m <- c("BUD", "DBF", "ENF", "VRM", "ELE", "log_DistMainR", "log_DistMinorR", "abs_DistEdge")

scaling_params_100m <- data.frame(
  variable=scale_vars_100m,
  mean=sapply(sample_model_100m[, scale_vars_100m], mean),
  sd=sapply(sample_model_100m[, scale_vars_100m], sd),
  row.names=NULL
)

if (
  any(!is.finite(scaling_params_100m$sd)) ||
  any(scaling_params_100m$sd <= 0)
) {
  stop("Invalid scaling SD detected.")
}

get_scale_mean <- function(variable) {
  scaling_params_100m$mean[
    scaling_params_100m$variable == variable
  ]
}

get_scale_sd <- function(variable) {
  scaling_params_100m$sd[
    scaling_params_100m$variable == variable
  ]
}

sample_model_100m <- sample_model_100m %>%
  mutate(
    BUD_z=(.data$BUD - get_scale_mean("BUD")) /
      get_scale_sd("BUD"),
    DBF_z=(.data$DBF - get_scale_mean("DBF")) /
      get_scale_sd("DBF"),
    ENF_z=(.data$ENF - get_scale_mean("ENF")) /
      get_scale_sd("ENF"),
    VRM_z=(.data$VRM - get_scale_mean("VRM")) /
      get_scale_sd("VRM"),
    ELE_z=(.data$ELE - get_scale_mean("ELE")) /
      get_scale_sd("ELE"),
    log_DistMainR_z=(.data$log_DistMainR - get_scale_mean("log_DistMainR")) /
      get_scale_sd("log_DistMainR"),
    log_DistMinorR_z=(.data$log_DistMinorR - get_scale_mean("log_DistMinorR")) /
      get_scale_sd("log_DistMinorR"),
    abs_DistEdge_z=(.data$abs_DistEdge - get_scale_mean("abs_DistEdge")) /
      get_scale_sd("abs_DistEdge")
  )

write.csv(
  scaling_params_100m,
  file.path(model_table_dir, "scaling_params_100m.csv"),
  row.names=FALSE
)

# 13. DBF / ENF / ELE confounding sensitivity
#
# DBF quadratic structure is used here.
m_DBF <- glm(
  presence ~ DBF_z + I(DBF_z^2),
  data=sample_model_100m,
  family=binomial()
)

m_ENF <- glm(
  presence ~ ENF_z,
  data=sample_model_100m,
  family=binomial()
)

m_ELE <- glm(
  presence ~ ELE_z,
  data=sample_model_100m,
  family=binomial()
)

m_DBF_ENF <- glm(
  presence ~ DBF_z + I(DBF_z^2) + ENF_z,
  data=sample_model_100m,
  family=binomial()
)

m_DBF_ELE <- glm(
  presence ~ DBF_z + I(DBF_z^2) + ELE_z,
  data=sample_model_100m,
  family=binomial()
)

m_ENF_ELE <- glm(
  presence ~ ENF_z + ELE_z,
  data=sample_model_100m,
  family=binomial()
)

m_DBF_ENF_ELE <- glm(
  presence ~ DBF_z + I(DBF_z^2) + ENF_z + ELE_z,
  data=sample_model_100m,
  family=binomial()
)

extract_coef <- function(model, model_name) {
  cc <- summary(model)$coefficients

  data.frame(
    model=model_name,
    variable=rownames(cc),
    beta=cc[, "Estimate"],
    SE=cc[, "Std. Error"],
    z=cc[, "z value"],
    p=cc[, "Pr(>|z|)"],
    row.names=NULL
  ) %>%
    filter(
      .data$variable != "(Intercept)"
    )
}

confounding_results <- bind_rows(
  extract_coef(m_DBF, "DBF quadratic"),
  extract_coef(m_ENF, "ENF"),
  extract_coef(m_ELE, "ELE"),
  extract_coef(m_DBF_ENF, "DBF quadratic + ENF"),
  extract_coef(m_DBF_ELE, "DBF quadratic + ELE"),
  extract_coef(m_ENF_ELE, "ENF + ELE"),
  extract_coef(m_DBF_ENF_ELE, "DBF quadratic + ENF + ELE")
)

print(
  tibble::as_tibble(confounding_results),
  n=Inf,
  width=Inf
)

write.csv(
  confounding_results,
  file.path(model_table_dir, "DBF_ENF_ELE_sensitivity_100m.csv"),
  row.names=FALSE
)

# 14. Create 4-fold spatial k-means CV folds
set.seed(20261025)

xy_scaled_100m <- scale(
  sample_model_100m[, c("x", "y")]
)

km_100m <- kmeans(xy_scaled_100m, centers=4, nstart=100)

sample_model_100m$fold <- km_100m$cluster

stopifnot(
  length(unique(sample_model_100m$fold)) == 4
)

fold_summary <- sample_model_100m %>%
  count(fold, zone) %>%
  tidyr::pivot_wider(names_from=zone, values_from=n, values_fill=0) %>%
  mutate(total=.data$A + .data$C)

print(fold_summary)

write.csv(
  fold_summary,
  file.path(model_table_dir, "spatial_CV_fold_summary_100m.csv"),
  row.names=FALSE
)

# 15. CV functions
auc_rank <- function(obs, pred) {
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)

  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  r <- rank(pred, ties.method="average")

  (
    sum(r[obs == 1]) -
      n1 * (n1 + 1) / 2
  ) /
    (n1 * n0)
}

evaluate_model_cv <- function(formula, model_name, data) {
  metric_list <- vector("list", 4)

  coef_list <- vector("list", 4)

  for (fold_i in 1:4) {
    train <- data %>%
      filter(.data$fold != .env$fold_i)

    test <- data %>%
      filter(.data$fold == .env$fold_i)

    mod <- glm(
      formula,
      data=train,
      family=binomial()
    )

    pred <- predict(mod, newdata=test, type="response")

    if (any(!is.finite(pred))) {
      stop(
        paste("Non-finite prediction:", model_name, "fold", fold_i)
      )
    }

    pred_clip <- pmin(
      pmax(pred, 1e-15),
      1 - 1e-15
    )

    metric_list[[fold_i]] <- data.frame(
      model=model_name,
      fold=fold_i,
      n_test=nrow(test),
      AUC=auc_rank(test$presence, pred),
      logloss=-mean(
        test$presence * log(pred_clip) +
          (1 - test$presence) * log(1 - pred_clip)
      ),
      Brier=mean(
        (test$presence - pred)^2
      )
    )

    cc <- coef(mod)

    coef_list[[fold_i]] <- data.frame(
      model=model_name,
      fold=fold_i,
      variable=names(cc),
      beta=as.numeric(cc),
      row.names=NULL
    )
  }

  list(
    metrics=bind_rows(metric_list),
    coefficients=bind_rows(coef_list)
  )
}

summarise_cv <- function(cv_results) {
  cv_results %>%
    group_by(.data$model) %>%
    summarise(
      mean_AUC=mean(.data$AUC),
      sd_AUC=sd(.data$AUC),
      mean_logloss=mean(.data$logloss),
      sd_logloss=sd(.data$logloss),
      mean_Brier=mean(.data$Brier),
      sd_Brier=sd(.data$Brier),
      .groups="drop"
    )
}

# 16. Linear vs quadratic DBF comparison
form_linear_DBF <- presence ~
  BUD_z +
  DBF_z +
  log_DistMainR_z +
  log_DistMinorR_z

form_quadratic_DBF <- presence ~
  BUD_z +
  DBF_z +
  I(DBF_z^2) +
  log_DistMainR_z +
  log_DistMinorR_z

m_linear_DBF <- glm(
  form_linear_DBF,
  data=sample_model_100m,
  family=binomial()
)

m_quadratic_DBF <- glm(
  form_quadratic_DBF,
  data=sample_model_100m,
  family=binomial()
)

DBF_AIC <- data.frame(
  model=c("linear_DBF", "quadratic_DBF"),
  AIC=c(
    AIC(m_linear_DBF),
    AIC(m_quadratic_DBF)
  )
)

DBF_AIC$delta_AIC <- DBF_AIC$AIC -
  min(DBF_AIC$AIC)

DBF_LRT <- anova(m_linear_DBF, m_quadratic_DBF, test="Chisq")

cv_linear_DBF <- evaluate_model_cv(form_linear_DBF, "linear_DBF", sample_model_100m)

cv_quadratic_DBF <- evaluate_model_cv(form_quadratic_DBF, "quadratic_DBF", sample_model_100m)

DBF_cv_metrics <- bind_rows(cv_linear_DBF$metrics, cv_quadratic_DBF$metrics)

DBF_cv_summary <- summarise_cv(DBF_cv_metrics)

cat("\nDBF AIC:\n")

print(DBF_AIC)

cat("\nDBF LRT:\n")

print(DBF_LRT)

cat("\nDBF CV:\n")

print(DBF_cv_summary)

write.csv(
  DBF_AIC,
  file.path(model_table_dir, "DBF_linear_vs_quadratic_AIC.csv"),
  row.names=FALSE
)

write.csv(
  DBF_cv_metrics,
  file.path(model_table_dir, "DBF_linear_vs_quadratic_CV_folds.csv"),
  row.names=FALSE
)

write.csv(
  DBF_cv_summary,
  file.path(model_table_dir, "DBF_linear_vs_quadratic_CV_summary.csv"),
  row.names=FALSE
)

# 17. Candidate models with quadratic DBF fixed
form_base <- presence ~
  BUD_z +
  DBF_z +
  I(DBF_z^2) +
  log_DistMainR_z +
  log_DistMinorR_z

form_VRM <- presence ~
  BUD_z +
  DBF_z +
  I(DBF_z^2) +
  VRM_z +
  log_DistMainR_z +
  log_DistMinorR_z

form_Edge <- presence ~
  BUD_z +
  DBF_z +
  I(DBF_z^2) +
  log_DistMainR_z +
  log_DistMinorR_z +
  abs_DistEdge_z * outside_forest

form_VRM_Edge <- presence ~
  BUD_z +
  DBF_z +
  I(DBF_z^2) +
  VRM_z +
  log_DistMainR_z +
  log_DistMinorR_z +
  abs_DistEdge_z * outside_forest

# 18. Sequential smoke test
smoke_train <- sample_model_100m %>%
  filter(.data$fold != 1)

smoke_test <- sample_model_100m %>%
  filter(.data$fold == 1)

smoke_forms <- list(base=form_base, VRM=form_VRM, Edge=form_Edge, VRM_Edge=form_VRM_Edge)

for (model_i in names(smoke_forms)) {
  smoke_mod <- glm(
    smoke_forms[[model_i]],
    data=smoke_train,
    family=binomial()
  )

  smoke_pred <- predict(smoke_mod, newdata=smoke_test, type="response")

  if (
    length(smoke_pred) != nrow(smoke_test) ||
    any(!is.finite(smoke_pred))
  ) {
    stop(
      paste("Smoke test failed:", model_i)
    )
  }
}

cat("\nCandidate-model smoke test passed.\n")

# 19. Candidate-model spatial CV
cv_base <- evaluate_model_cv(form_base, "base", sample_model_100m)

cv_VRM <- evaluate_model_cv(form_VRM, "VRM", sample_model_100m)

cv_Edge <- evaluate_model_cv(form_Edge, "Edge", sample_model_100m)

cv_VRM_Edge <- evaluate_model_cv(form_VRM_Edge, "VRM_Edge", sample_model_100m)

candidate_cv_metrics <- bind_rows(
  cv_base$metrics,
  cv_VRM$metrics,
  cv_Edge$metrics,
  cv_VRM_Edge$metrics
)

candidate_cv_summary <- summarise_cv(candidate_cv_metrics) %>%
  arrange(
    .data$mean_logloss,
    .data$mean_Brier,
    desc(.data$mean_AUC)
  )

print(
  tibble::as_tibble(candidate_cv_metrics),
  n=Inf,
  width=Inf
)

print(
  tibble::as_tibble(candidate_cv_summary),
  n=Inf,
  width=Inf
)

write.csv(
  candidate_cv_metrics,
  file.path(model_table_dir, "candidate_model_spatial_CV_folds.csv"),
  row.names=FALSE
)

write.csv(
  candidate_cv_summary,
  file.path(model_table_dir, "candidate_model_spatial_CV_summary.csv"),
  row.names=FALSE
)

# 20. Candidate coefficient stability
candidate_cv_coef <- bind_rows(
  cv_base$coefficients,
  cv_VRM$coefficients,
  cv_Edge$coefficients,
  cv_VRM_Edge$coefficients
)

candidate_coef_stability <- candidate_cv_coef %>%
  filter(
    .data$variable != "(Intercept)"
  ) %>%
  group_by(.data$model, .data$variable) %>%
  summarise(
    mean_beta=mean(.data$beta),
    min_beta=min(.data$beta),
    max_beta=max(.data$beta),
    n_positive=sum(.data$beta > 0),
    n_negative=sum(.data$beta < 0),
    .groups="drop"
  )

print(
  tibble::as_tibble(candidate_coef_stability),
  n=Inf,
  width=Inf
)

write.csv(
  candidate_cv_coef,
  file.path(model_table_dir, "candidate_model_spatial_CV_coefficients.csv"),
  row.names=FALSE
)

write.csv(
  candidate_coef_stability,
  file.path(model_table_dir, "candidate_model_coefficient_stability.csv"),
  row.names=FALSE
)

# 21. Full-sample candidate models and AIC
candidate_formulas <- list(
  base=form_base,
  VRM=form_VRM,
  Edge=form_Edge,
  VRM_Edge=form_VRM_Edge
)

candidate_models <- lapply(
  candidate_formulas,
  function(ff) {
    glm(
      ff,
      data=sample_model_100m,
      family=binomial()
    )
  }
)

candidate_AIC <- data.frame(
  model=names(candidate_models),
  AIC=sapply(candidate_models, AIC),
  row.names=NULL
)

candidate_AIC$delta_AIC <- candidate_AIC$AIC -
  min(candidate_AIC$AIC)

candidate_AIC <- candidate_AIC %>%
  arrange(.data$AIC)

print(candidate_AIC)

write.csv(
  candidate_AIC,
  file.path(model_table_dir, "candidate_model_AIC.csv"),
  row.names=FALSE
)

# 22. Final-model selection
#
# Default:
# lowest mean spatial-CV log-loss.
#
# Change manually if interpretation / stability warrants it:
# "base", "VRM", "Edge", "VRM_Edge"
final_model_override <- "base"  # Fixed a priori to reproduce the reported final RSF

valid_model_names <- c("base", "VRM", "Edge", "VRM_Edge")

if (
  !is.na(final_model_override) &&
  !final_model_override %in% valid_model_names
) {
  stop("Invalid final_model_override.")
}

if (is.na(final_model_override)) {
  final_model_name <- candidate_cv_summary$model[
    which.min(candidate_cv_summary$mean_logloss)
  ]
} else {
  final_model_name <- final_model_override
}

final_formula <- candidate_formulas[[final_model_name]]
m_final_100m <- candidate_models[[final_model_name]]

cat("\n========================================\n")

cat("FINAL MODEL:", final_model_name, "\n")

cat("========================================\n")

print(final_formula)

print(
  summary(m_final_100m)
)

# 23. Final coefficient table
cc_final <- summary(m_final_100m)$coefficients

final_coef <- data.frame(
  variable=rownames(cc_final),
  beta=cc_final[, "Estimate"],
  SE=cc_final[, "Std. Error"],
  z=cc_final[, "z value"],
  p=cc_final[, "Pr(>|z|)"],
  row.names=NULL
)

print(final_coef)

write.csv(
  final_coef,
  file.path(model_table_dir, "final_model_coefficients_100m.csv"),
  row.names=FALSE
)

# 24. Final CV outputs
final_cv_result <- switch(
  final_model_name,
  base=cv_base,
  VRM=cv_VRM,
  Edge=cv_Edge,
  VRM_Edge=cv_VRM_Edge
)

final_cv_metrics <- final_cv_result$metrics

final_cv_summary <- summarise_cv(final_cv_metrics)

final_cv_coef <- final_cv_result$coefficients

final_coef_stability <- final_cv_coef %>%
  filter(
    .data$variable != "(Intercept)"
  ) %>%
  group_by(.data$variable) %>%
  summarise(
    mean_beta=mean(.data$beta),
    min_beta=min(.data$beta),
    max_beta=max(.data$beta),
    n_positive=sum(.data$beta > 0),
    n_negative=sum(.data$beta < 0),
    .groups="drop"
  )

write.csv(
  final_cv_metrics,
  file.path(model_table_dir, "final_model_spatial_CV_folds.csv"),
  row.names=FALSE
)

write.csv(
  final_cv_summary,
  file.path(model_table_dir, "final_model_spatial_CV_summary.csv"),
  row.names=FALSE
)

write.csv(
  final_coef_stability,
  file.path(model_table_dir, "final_model_coefficient_stability.csv"),
  row.names=FALSE
)

# 25. Final VIF / GVIF
final_vif <- car::vif(m_final_100m)

print(final_vif)

if (is.matrix(final_vif)) {
  final_vif_table <- data.frame(
    variable=rownames(final_vif),
    final_vif,
    row.names=NULL,
    check.names=FALSE
  )
} else {
  final_vif_table <- data.frame(
    variable=names(final_vif),
    VIF=as.numeric(final_vif)
  )
}

write.csv(
  final_vif_table,
  file.path(model_table_dir, "final_model_VIF_GVIF.csv"),
  row.names=FALSE
)

# 26. DBF quadratic vertex
final_beta <- coef(m_final_100m)

if (
  all(
    c("DBF_z", "I(DBF_z^2)") %in% names(final_beta)
  )
) {
  DBF_vertex_z <- -final_beta["DBF_z"] /
    (
      2 *
        final_beta["I(DBF_z^2)"]
    )

  DBF_vertex_percentile <- mean(sample_model_100m$DBF_z <= DBF_vertex_z, na.rm=TRUE) * 100

  DBF_vertex_raw <- get_scale_mean("DBF") +
    DBF_vertex_z *
    get_scale_sd("DBF")

  DBF_vertex_table <- data.frame(
    vertex_z=as.numeric(DBF_vertex_z),
    vertex_raw=as.numeric(DBF_vertex_raw),
    percentile=DBF_vertex_percentile
  )
} else {
  DBF_vertex_table <- data.frame(vertex_z=NA_real_, vertex_raw=NA_real_, percentile=NA_real_)
}

print(DBF_vertex_table)

write.csv(
  DBF_vertex_table,
  file.path(model_table_dir, "DBF_quadratic_vertex.csv"),
  row.names=FALSE
)

# 27. DBF conditional response curve
DBF_seq <- seq(
  quantile(sample_model_100m$DBF_z, 0.01),
  quantile(sample_model_100m$DBF_z, 0.99),
  length.out=300
)

DBF_curve <- data.frame(
  BUD_z=0,
  DBF_z=DBF_seq,
  VRM_z=0,
  log_DistMainR_z=0,
  log_DistMinorR_z=0,
  abs_DistEdge_z=0,
  outside_forest=0
)

DBF_curve$fitted_score <- predict(m_final_100m, newdata=DBF_curve, type="response")

DBF_curve$DBF_raw <- get_scale_mean("DBF") +
  DBF_curve$DBF_z *
  get_scale_sd("DBF")

write.csv(
  DBF_curve,
  file.path(model_table_dir, "DBF_conditional_response_curve.csv"),
  row.names=FALSE
)

# 28. Main-road sign-reversal diagnostic
m_main_only <- glm(
  presence ~ log_DistMainR_z,
  data=sample_model_100m,
  family=binomial()
)

m_main_BUD <- glm(
  presence ~ log_DistMainR_z + BUD_z,
  data=sample_model_100m,
  family=binomial()
)

m_main_DBF <- glm(
  presence ~
    log_DistMainR_z +
    DBF_z +
    I(DBF_z^2),
  data=sample_model_100m,
  family=binomial()
)

m_main_BUD_DBF <- glm(
  presence ~
    log_DistMainR_z +
    BUD_z +
    DBF_z +
    I(DBF_z^2),
  data=sample_model_100m,
  family=binomial()
)

mainroad_sign_diagnostic <- data.frame(
  model=c(
    "MainR only",
    "MainR + BUD",
    "MainR + DBF quadratic",
    "MainR + BUD + DBF quadratic",
    "Final model"
  ),
  beta_MainR=c(
    coef(m_main_only)["log_DistMainR_z"],
    coef(m_main_BUD)["log_DistMainR_z"],
    coef(m_main_DBF)["log_DistMainR_z"],
    coef(m_main_BUD_DBF)["log_DistMainR_z"],
    coef(m_final_100m)["log_DistMainR_z"]
  )
)

print(mainroad_sign_diagnostic)

write.csv(
  mainroad_sign_diagnostic,
  file.path(model_table_dir, "MainRoad_sign_reversal_diagnostic.csv"),
  row.names=FALSE
)

# 29. Save model objects before mapping
saveRDS(
  sample_model_100m,
  file.path(model_rds_dir, "sample_model_100m_FINAL.rds")
)

saveRDS(
  scaling_params_100m,
  file.path(model_rds_dir, "scaling_params_100m_FINAL.rds")
)

saveRDS(
  m_final_100m,
  file.path(model_rds_dir, "landscape_RSF_final_100m.rds")
)

saveRDS(
  list(
    final_model_name=final_model_name,
    final_formula=final_formula,
    inside_sign=inside_sign,
    outside_sign=outside_sign,
    candidate_cv_summary=candidate_cv_summary,
    candidate_AIC=candidate_AIC
  ),
  file.path(model_rds_dir, "final_model_selection_objects.rds")
)

# 30. Full-Shikoku 100-m template
xmin_full <- floor(
  xmin(lulc_native) / 100
) * 100

xmax_full <- ceiling(
  xmax(lulc_native) / 100
) * 100

ymin_full <- floor(
  ymin(lulc_native) / 100
) * 100

ymax_full <- ceiling(
  ymax(lulc_native) / 100
) * 100

template_full_100m <- rast(
  xmin=xmin_full,
  xmax=xmax_full,
  ymin=ymin_full,
  ymax=ymax_full,
  resolution=100,
  crs=crs(lulc_native)
)

names(template_full_100m) <- "template_full_100m"

values(template_full_100m) <- NA_real_

writeRaster(
  template_full_100m,
  file.path(model_raster_dir, "template_full_Shikoku_100m.tif"),
  overwrite=TRUE
)

# 31. Full-Shikoku land mask
land_fine <- ifel(
  is.na(lulc_native),
  NA,
  ifel(lulc_native == 0, 0, 1)
)

land_prop_100m <- resample(land_fine, template_full_100m, method="average", threads=TRUE)

land_mask_100m <- ifel(land_prop_100m >= 0.5, 1, NA)

names(land_mask_100m) <- "land"

writeRaster(
  land_mask_100m,
  file.path(model_raster_dir, "Shikoku_land_mask_100m.tif"),
  overwrite=TRUE,
  datatype="INT1U"
)

# 32. Find selected native focal rasters
best_scale_lookup <- setNames(best_scales_100m$scale, best_scales_100m$variable)

get_selected_native_file <- function(variable) {
  sc <- best_scale_lookup[[variable]]

  f <- focal_file_table$native_file[
    focal_file_table$variable == variable &
      focal_file_table$scale == sc
  ]

  if (length(f) != 1) {
    stop(
      paste("Selected native focal file not unique:", variable)
    )
  }

  if (!file.exists(f)) {
    stop(
      paste("Selected native focal file missing:", f)
    )
  }

  f
}

selected_native_files <- c(
  BUD=get_selected_native_file("BUD"),
  DBF=get_selected_native_file("DBF"),
  ENF=get_selected_native_file("ENF"),
  VRM=get_selected_native_file("VRM"),
  ELE=get_selected_native_file("ELE")
)

# 33. Helper for full-Shikoku raster generation
build_full_raster <- function(source_file, output_file, method, output_name) {
  rr <- rast(source_file)

  if (!same.crs(rr, template_full_100m)) {
    stop(
      paste("CRS mismatch:", source_file)
    )
  }

  rr_crop <- crop(rr, template_full_100m, snap="out")

  rr_full <- resample(rr_crop, template_full_100m, method=method, threads=TRUE)

  names(rr_full) <- output_name

  writeRaster(rr_full, output_file, overwrite=TRUE)

  rm(rr, rr_crop, rr_full)

  gc()

  output_file
}

# 34. Generate full-Shikoku predictor rasters
full_BUD_file <- build_full_raster(
  selected_native_files[["BUD"]],
  file.path(temp_full_dir, "BUD_100m.tif"),
  "average",
  "BUD"
)

full_DBF_file <- build_full_raster(
  selected_native_files[["DBF"]],
  file.path(temp_full_dir, "DBF_100m.tif"),
  "average",
  "DBF"
)

full_VRM_file <- build_full_raster(
  selected_native_files[["VRM"]],
  file.path(temp_full_dir, "VRM_100m.tif"),
  "average",
  "VRM"
)

DistMainR_native_file <- file.path(
  native_var_dir,
  "distance",
  "DistMainR_notunnel_whole_mainland.tif"
)

DistMinorR_native_file <- file.path(
  native_var_dir,
  "distance",
  "DistMinorR_notunnel_whole_mainland.tif"
)

DistEdge_native_file <- file.path(native_var_dir, "distance", "Dist_edge1820.tif")

full_Main_file <- build_full_raster(
  DistMainR_native_file,
  file.path(temp_full_dir, "DistMainR_100m.tif"),
  "bilinear",
  "DistMainR"
)

full_Minor_file <- build_full_raster(
  DistMinorR_native_file,
  file.path(temp_full_dir, "DistMinorR_100m.tif"),
  "bilinear",
  "DistMinorR"
)

full_Edge_file <- build_full_raster(
  DistEdge_native_file,
  file.path(temp_full_dir, "DistEdge_signed_100m.tif"),
  "bilinear",
  "DistEdge_signed"
)

# 35. Load full-Shikoku predictor rasters
BUD_full <- rast(full_BUD_file)

DBF_full <- rast(full_DBF_file)

VRM_full <- rast(full_VRM_file)

Main_full <- rast(full_Main_file)

Minor_full <- rast(full_Minor_file)

Edge_full <- rast(full_Edge_file)

# 36. Standardize full-Shikoku predictors
BUD_z_full <- (
  BUD_full -
    get_scale_mean("BUD")
) /
  get_scale_sd("BUD")

DBF_z_full <- (
  DBF_full -
    get_scale_mean("DBF")
) /
  get_scale_sd("DBF")

VRM_z_full <- (
  VRM_full -
    get_scale_mean("VRM")
) /
  get_scale_sd("VRM")

log_Main_full <- log1p(Main_full)

log_Minor_full <- log1p(Minor_full)

log_Main_z_full <- (
  log_Main_full -
    get_scale_mean("log_DistMainR")
) /
  get_scale_sd("log_DistMainR")

log_Minor_z_full <- (
  log_Minor_full -
    get_scale_mean("log_DistMinorR")
) /
  get_scale_sd("log_DistMinorR")

abs_Edge_full <- abs(Edge_full)

abs_Edge_z_full <- (
  abs_Edge_full -
    get_scale_mean("abs_DistEdge")
) /
  get_scale_sd("abs_DistEdge")

if (outside_sign == "positive") {
  outside_full <- ifel(Edge_full > 0, 1, 0)
} else {
  outside_full <- ifel(Edge_full < 0, 1, 0)
}

names(outside_full) <- "outside_forest"

# 37. Manual full-Shikoku prediction
final_beta <- coef(m_final_100m)

lp_full <- BUD_z_full * 0 +
  final_beta["(Intercept)"]

add_term <- function(lp, beta_name, raster_term) {
  if (beta_name %in% names(final_beta)) {
    lp <- lp +
      raster_term *
      final_beta[[beta_name]]
  }

  lp
}

lp_full <- add_term(lp_full, "BUD_z", BUD_z_full)

lp_full <- add_term(lp_full, "DBF_z", DBF_z_full)

lp_full <- add_term(
  lp_full,
  "I(DBF_z^2)",
  DBF_z_full^2
)

lp_full <- add_term(lp_full, "VRM_z", VRM_z_full)

lp_full <- add_term(lp_full, "log_DistMainR_z", log_Main_z_full)

lp_full <- add_term(lp_full, "log_DistMinorR_z", log_Minor_z_full)

lp_full <- add_term(lp_full, "abs_DistEdge_z", abs_Edge_z_full)

lp_full <- add_term(lp_full, "outside_forest", outside_full)

if ("abs_DistEdge_z:outside_forest" %in% names(final_beta)) {
  lp_full <- lp_full +
    abs_Edge_z_full *
    outside_full *
    final_beta[["abs_DistEdge_z:outside_forest"]]
}

if ("outside_forest:abs_DistEdge_z" %in% names(final_beta)) {
  lp_full <- lp_full +
    abs_Edge_z_full *
    outside_full *
    final_beta[["outside_forest:abs_DistEdge_z"]]
}

suitability_full <- 1 /
  (
    1 +
      exp(-lp_full)
  )

suitability_full <- mask(suitability_full, land_mask_100m)

names(suitability_full) <- "relative_landscape_suitability"

full_suitability_file <- file.path(
  model_raster_dir,
  "relative_landscape_suitability_full_Shikoku_100m.tif"
)

writeRaster(suitability_full, full_suitability_file, overwrite=TRUE)

print(
  global(
    suitability_full,
    c("min", "mean", "max"),
    na.rm=TRUE
  )
)

# 38. A+B+C suitability raster
ABC_v <- vect(ABC)

suitability_ABC <- crop(suitability_full, ABC_v)

suitability_ABC <- mask(suitability_ABC, ABC_v, touches=FALSE)

names(suitability_ABC) <- "relative_landscape_suitability"

ABC_suitability_file <- file.path(
  model_raster_dir,
  "relative_landscape_suitability_ABC_100m.tif"
)

writeRaster(suitability_ABC, ABC_suitability_file, overwrite=TRUE)

# 39. A+B+C equal-area cut points
ABC_values <- values(suitability_ABC, mat=FALSE, na.rm=TRUE)

if (length(ABC_values) == 0) {
  stop("No valid A+B+C suitability values.")
}

q_probs <- seq(0, 1, by=0.1)

equal_area_cutpoints <- as.numeric(
  quantile(ABC_values, probs=q_probs, na.rm=TRUE, names=FALSE, type=7)
)

if (
  length(unique(equal_area_cutpoints)) !=
  length(equal_area_cutpoints)
) {
  stop("Duplicated equal-area cut points detected.")
}

cutpoint_table <- data.frame(quantile=q_probs, cutpoint=equal_area_cutpoints)

print(cutpoint_table)

write.csv(
  cutpoint_table,
  file.path(model_table_dir, "ABC_equal_area_bin_cutpoints.csv"),
  row.names=FALSE
)

# 40. A+B+C equal-area bins
ABC_rcl <- cbind(equal_area_cutpoints[-11], equal_area_cutpoints[-1], 1:10)

ABC_bin <- classify(suitability_ABC, rcl=ABC_rcl, include.lowest=TRUE, right=TRUE)

names(ABC_bin) <- "suitability_bin"

writeRaster(
  ABC_bin,
  file.path(model_raster_dir, "relative_suitability_equal_area_bins_ABC_100m.tif"),
  overwrite=TRUE,
  datatype="INT1U"
)

ABC_bin_freq <- as.data.frame(
  freq(ABC_bin)
)

ABC_bin_summary <- ABC_bin_freq %>%
  transmute(
    bin=.data$value,
    n_cells=.data$count,
    area_km2=.data$count * 0.01,
    landscape_percent=.data$count /
      sum(.data$count) *
      100
  )

print(ABC_bin_summary)

write.csv(
  ABC_bin_summary,
  file.path(model_table_dir, "ABC_equal_area_bin_area_summary.csv"),
  row.names=FALSE
)

# 41. Apply ABC cut points to full Shikoku
full_bin_breaks <- c(-1e30, equal_area_cutpoints[2:10], 1e30)

full_bin_rcl <- cbind(full_bin_breaks[-11], full_bin_breaks[-1], 1:10)

Shikoku_bin <- classify(suitability_full, rcl=full_bin_rcl, include.lowest=TRUE, right=TRUE)

names(Shikoku_bin) <- "suitability_bin"

Shikoku_bin_file <- file.path(
  model_raster_dir,
  "relative_suitability_ABC_calibrated_bins_full_Shikoku_100m.tif"
)

writeRaster(Shikoku_bin, Shikoku_bin_file, overwrite=TRUE, datatype="INT1U")

Shikoku_bin_freq <- as.data.frame(
  freq(Shikoku_bin)
)

Shikoku_bin_summary <- Shikoku_bin_freq %>%
  transmute(
    bin=.data$value,
    n_cells=.data$count,
    area_km2=.data$count * 0.01,
    Shikoku_percent=.data$count /
      sum(.data$count) *
      100
  )

print(Shikoku_bin_summary)

write.csv(
  Shikoku_bin_summary,
  file.path(model_table_dir, "Shikoku_ABC_calibrated_bin_area_summary.csv"),
  row.names=FALSE
)

# 42. Load restricted GPS-step data (never publish the underlying records)
if (!file.exists(gps_source_file)) {
  stop("Private GPS-step file missing: ", gps_source_file)
}
df_step_all <- readRDS(gps_source_file)
required_gps_cols <- c("id", "case_", "t2_", "x2_", "y2_")
if (!all(required_gps_cols %in% names(df_step_all))) {
  stop("Private GPS-step file has unexpected columns.")
}

# 43. Actual GPS endpoints
gps_used <- df_step_all %>%
  filter(.data$case_ == TRUE) %>%
  distinct(id, t2_, x2_, y2_, .keep_all=TRUE) %>%
  select(id, t2_, x2_, y2_) %>%
  filter(
    !is.na(.data$x2_),
    !is.na(.data$y2_)
  ) %>%
  mutate(
    real_id=sub(
      "_[0-9]{4}$",
      "",
      .data$id
    )
  )

cat(
  "GPS used fixes:",
  nrow(gps_used),
  "\n"
)

cat(
  "Real individuals:",
  n_distinct(gps_used$real_id),
  "\n"
)

# 44. Extract suitability and bins to GPS
gps_points_v <- vect(
  gps_used[, c("x2_", "y2_")],
  geom=c("x2_", "y2_"),
  crs=crs(suitability_full)
)

gps_score_extract <- terra::extract(suitability_full, gps_points_v)

gps_bin_extract <- terra::extract(Shikoku_bin, gps_points_v)

score_cols <- setdiff(
  names(gps_score_extract),
  "ID"
)

bin_cols <- setdiff(
  names(gps_bin_extract),
  "ID"
)

if (
  length(score_cols) != 1 ||
  length(bin_cols) != 1
) {
  stop("Unexpected GPS extraction structure.")
}

gps_used$score <- gps_score_extract[[score_cols]]
gps_used$bin <- as.integer(gps_bin_extract[[bin_cols]])

gps_used <- gps_used %>%
  filter(
    !is.na(.data$score),
    !is.na(.data$bin)
  )

gps_bin_summary <- gps_used %>%
  count(bin, name="GPS_n") %>%
  mutate(
    GPS_percent=.data$GPS_n /
      sum(.data$GPS_n) *
      100
  )

print(gps_bin_summary)

write.csv(
  gps_bin_summary,
  file.path(model_table_dir, "GPS_distribution_ABC_equal_area_bins.csv"),
  row.names=FALSE
)

saveRDS(
  gps_used,
  file.path(model_rds_dir, "GPS_with_final_suitability_and_bins.rds")
)

# 45. 90% GPS coverage threshold
get_pooled_threshold <- function(bin_values, target=0.90) {
  coverage <- sapply(
    1:10,
    function(bin_min) {
      mean(bin_values >= bin_min)
    }
  )

  eligible <- which(coverage >= target)

  if (length(eligible) == 0) {
    return(NA_integer_)
  }

  max(eligible)
}

observed_GPS_threshold <- get_pooled_threshold(gps_used$bin, target=0.90)

cat("Observed pooled GPS threshold:", observed_GPS_threshold, "\n")

# 46. Individual bootstrap
set.seed(20261102)

n_boot <- 1000
target_GPS <- 0.90

real_ids <- unique(gps_used$real_id)

n_real_ids <- length(real_ids)

bootstrap_threshold <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  sampled_ids <- sample(real_ids, size=n_real_ids, replace=TRUE)

  sampled_bins <- unlist(
    lapply(
      sampled_ids,
      function(id_i) {
        gps_used$bin[
          gps_used$real_id == id_i
        ]
      }
    ),
    use.names=FALSE
  )

  threshold_bin <- get_pooled_threshold(sampled_bins, target=target_GPS)

  bootstrap_threshold[[b]] <- data.frame(
    bootstrap=b,
    threshold_bin=threshold_bin,
    n_locations=length(sampled_bins)
  )
}

bootstrap_threshold <- bind_rows(bootstrap_threshold)

threshold_distribution <- bootstrap_threshold %>%
  count(threshold_bin, name="n") %>%
  mutate(
    proportion=.data$n /
      sum(.data$n),
    percent=.data$proportion *
      100
  ) %>%
  arrange(.data$threshold_bin)

threshold_summary <- data.frame(
  observed_threshold=observed_GPS_threshold,
  median_bootstrap_threshold=median(bootstrap_threshold$threshold_bin, na.rm=TRUE),
  q025=as.numeric(
    quantile(bootstrap_threshold$threshold_bin, 0.025, type=1, na.rm=TRUE)
  ),
  q975=as.numeric(
    quantile(bootstrap_threshold$threshold_bin, 0.975, type=1, na.rm=TRUE)
  )
)

print(threshold_distribution)

print(threshold_summary)

write.csv(
  threshold_distribution,
  file.path(model_table_dir, "GPS_individual_bootstrap_90pct_threshold_distribution.csv"),
  row.names=FALSE
)

write.csv(
  threshold_summary,
  file.path(model_table_dir, "GPS_individual_bootstrap_90pct_threshold_summary.csv"),
  row.names=FALSE
)

saveRDS(
  bootstrap_threshold,
  file.path(model_rds_dir, "GPS_individual_bootstrap_90pct_threshold.rds")
)

# 47. Final suitable-bin threshold
#
# Conservative threshold = lower 2.5% bootstrap threshold.
suitable_min_bin <- as.integer(threshold_summary$q025)
if (suitable_min_bin != 8L) {
  stop("Bootstrap cutoff differs from the manuscript's bins 8-10. Investigate inputs.")
}

if (!suitable_min_bin %in% 1:10) {
  stop("Invalid suitable_min_bin.")
}

cat(
  "\nFinal suitable bins:",
  paste(suitable_min_bin:10, collapse=", "),
  "\n"
)

# 48. Suitable landscape raster
suitable_100m <- ifel(Shikoku_bin >= suitable_min_bin, 1, NA)

names(suitable_100m) <- "suitable"

suitable_file <- file.path(
  model_raster_dir,
  paste0("suitable_landscape_bins", suitable_min_bin, "_10_100m.tif")
)

writeRaster(suitable_100m, suitable_file, overwrite=TRUE, datatype="INT1U")

suitable_cells <- global(suitable_100m, "notNA")[1, 1]

suitable_area_km2 <- suitable_cells * 0.01

cat("Suitable area:", suitable_area_km2, "km2\n")

# 49. Contiguous suitable patches
patch_100m <- patches(suitable_100m, directions=8, zeroAsNA=TRUE)

names(patch_100m) <- "patch_id"

patch_freq <- as.data.frame(
  freq(patch_100m)
)

patch_summary <- data.frame(patch_id=patch_freq$value, n_cells=patch_freq$count)

patch_summary$area_km2 <- patch_summary$n_cells * 0.01

patch_summary <- patch_summary[
  order(patch_summary$area_km2, decreasing=TRUE),
]

row.names(patch_summary) <- NULL

write.csv(
  patch_summary,
  file.path(model_table_dir, "suitable_patch_summary_100m.csv"),
  row.names=FALSE
)

# 50. Core suitable habitat >=93.1 km2
core_threshold_km2 <- 93.1

core_patch_table <- patch_summary %>%
  filter(.data$area_km2 >= core_threshold_km2) %>%
  arrange(
    desc(.data$area_km2)
  )

if (nrow(core_patch_table) == 0) {
  stop("No suitable patch >= 93.1 km2.")
}

core_patch_table$core_rank <- seq_len(
  nrow(core_patch_table)
)

print(core_patch_table)

write.csv(
  core_patch_table,
  file.path(model_table_dir, "core_suitable_patch_summary_min93.1km2.csv"),
  row.names=FALSE
)

core_patch_ids <- core_patch_table$patch_id

core_patch_id_100m <- subst(patch_100m, from=core_patch_ids, to=core_patch_ids, others=NA)

names(core_patch_id_100m) <- "patch_id"

core_patch_file <- file.path(
  model_raster_dir,
  paste0("core_suitable_patch_ID_bins", suitable_min_bin, "_10_min93.1km2_100m.tif")
)

writeRaster(core_patch_id_100m, core_patch_file, overwrite=TRUE, datatype="INT4S")

core_binary_100m <- ifel(
  !is.na(core_patch_id_100m),
  1,
  NA
)

names(core_binary_100m) <- "core_suitable"

core_binary_file <- file.path(
  model_raster_dir,
  paste0("core_suitable_binary_bins", suitable_min_bin, "_10_min93.1km2_100m.tif")
)

writeRaster(core_binary_100m, core_binary_file, overwrite=TRUE, datatype="INT1U")

# 51. Core polygons
core_polygon <- as.polygons(core_patch_id_100m, dissolve=TRUE, na.rm=TRUE)

names(core_polygon) <- "patch_id"

core_polygon_file <- file.path(
  model_vector_dir,
  paste0("core_suitable_patches_bins", suitable_min_bin, "_10_min93.1km2_100m.gpkg")
)

writeVector(core_polygon, core_polygon_file, filetype="GPKG", overwrite=TRUE)

# 52. Occupied-range overlap summary
A_v <- vect(A)

suitable_within_A <- mask(suitable_100m, A_v, touches=FALSE)

core_within_A <- mask(core_binary_100m, A_v, touches=FALSE)

suitable_A_cells <- global(suitable_within_A, "notNA")[1, 1]

core_A_cells <- global(core_within_A, "notNA")[1, 1]

core_total_cells <- global(core_binary_100m, "notNA")[1, 1]

suitable_A_area_km2 <- suitable_A_cells * 0.01
core_A_area_km2 <- core_A_cells * 0.01
core_total_area_km2 <- core_total_cells * 0.01

core_outside_A_area_km2 <- core_total_area_km2 -
  core_A_area_km2

# 53. Final habitat summary
final_habitat_summary <- data.frame(
  final_model=final_model_name,
  suitable_min_bin=suitable_min_bin,
  suitable_max_bin=10,
  suitable_area_km2=suitable_area_km2,
  suitable_within_A_km2=suitable_A_area_km2,
  core_threshold_km2=core_threshold_km2,
  n_core_patches=nrow(core_patch_table),
  core_total_area_km2=core_total_area_km2,
  core_within_A_km2=core_A_area_km2,
  core_outside_A_km2=core_outside_A_area_km2
)

print(final_habitat_summary)

write.csv(
  final_habitat_summary,
  file.path(model_table_dir, "final_suitable_core_summary.csv"),
  row.names=FALSE
)

# Cross-check against the published landscape classification.
# These checks should flag changed inputs or analysis environments.
if (nrow(core_patch_table) != 4L) {
  stop("Core-patch count differs from the four patches reported in the manuscript.")
}
if (abs(suitable_area_km2 - 2883.04) > 1.0) {
  stop("Suitable-landscape area differs from the reported approximately 2,883 km2.")
}

# 54. Save final analysis objects
saveRDS(
  list(
    final_model_name=final_model_name,
    final_formula=final_formula,
    final_coefficients=final_coef,
    scaling_params=scaling_params_100m,
    inside_sign=inside_sign,
    outside_sign=outside_sign,
    equal_area_cutpoints=equal_area_cutpoints,
    suitable_min_bin=suitable_min_bin,
    core_threshold_km2=core_threshold_km2,
    core_patch_table=core_patch_table,
    final_habitat_summary=final_habitat_summary
  ),
  file.path(model_rds_dir, "RSF_final_analysis_objects.rds")
)

# 55. Final QC
required_final_outputs <- c(
  full_suitability_file,
  ABC_suitability_file,
  Shikoku_bin_file,
  suitable_file,
  core_patch_file,
  core_binary_file,
  core_polygon_file,
  file.path(model_rds_dir, "landscape_RSF_final_100m.rds"),
  file.path(model_rds_dir, "RSF_final_analysis_objects.rds")
)

missing_final_outputs <- required_final_outputs[
  !file.exists(required_final_outputs)
]

if (length(missing_final_outputs) > 0) {
  print(missing_final_outputs)

  stop("Script 02 finished with missing outputs.")
}

# 56. Final console summary
cat("\n========================================\n")

cat("RSF Script 02 completed successfully.\n")

cat("========================================\n")

cat("\nFinal model:", final_model_name, "\n")

cat("\nCandidate CV summary:\n")

print(candidate_cv_summary)

cat("\nCandidate coefficient stability:\n")

print(candidate_coef_stability)

cat("\nCandidate AIC:\n")

print(candidate_AIC)

cat("\nFinal coefficients:\n")

print(final_coef)

cat("\nFinal CV:\n")

print(final_cv_summary)

cat("\nDBF vertex:\n")

print(DBF_vertex_table)

cat("\nABC cut points:\n")

print(cutpoint_table)

cat("\nGPS bootstrap threshold:\n")

print(threshold_summary)

cat("\nFinal suitable/core summary:\n")

print(final_habitat_summary)

cat("========================================\n")

gc()
