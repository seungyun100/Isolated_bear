# Shikoku landscape RSF | 01: study-area and predictor preparation
#
# PURPOSE
#   1. Define A / B / C study zones
#   2. Build 100-m analysis grid
#   3. Build forest-dominated fitting domain
#   4. Draw fixed 10,000 A + 10,000 C sample
#   5. Prepare focal environmental rasters
#   6. Select focal scale for BUD / DBF / ENF / VRM / ELE
#   7. Prepare Main-road / Minor-road / signed Edge rasters
#   8. Retain all inputs required by Script 02
#
# IMPORTANT
#   - DBF focal scale is selected using a linear univariate model.
#   - Quadratic DBF structure is evaluated later in Script 02.
#   - DistEdge remains SIGNED here.
#   - Edge parameterization is evaluated later in Script 02.
#   - Outputs are saved under results/rsf/.
# 1. Packages
library(terra)
library(sf)
library(dplyr)

gc()

# Project-relative paths (configure in R/00_config.R)
source(file.path("R", "00_config.R"))
create_rsf_directories()

# 3. Original input files
area_file <- file.path(study_area_dir, "mcp_camera_gps.shp")

lulc_file <- file.path(native_var_dir, "lulc_1820_shikoku.tif")

DistMainR_file <- file.path(
  native_var_dir,
  "distance",
  "DistMainR_notunnel_whole_mainland.tif"
)

DistMinorR_file <- file.path(
  native_var_dir,
  "distance",
  "DistMinorR_notunnel_whole_mainland.tif"
)

DistEdge_file <- file.path(native_var_dir, "distance", "Dist_edge1820.tif")

stopifnot(file.exists(area_file))
stopifnot(file.exists(lulc_file))
stopifnot(file.exists(DistMainR_file))
stopifnot(file.exists(DistMinorR_file))
stopifnot(file.exists(DistEdge_file))

# 4. Focal-scale candidates
focal_scales <- c(25, 50, 100, 150, 250, 500, 1000, 2000, 4000)

focal_variables <- c("BUD", "DBF", "ENF", "VRM", "ELE")

# 5. Native focal-file helper functions
find_native_file <- function(filename) {
  hits <- list.files(native_var_dir, recursive=TRUE, full.names=TRUE)

  hits <- hits[
    basename(hits) == filename
  ]

  if (length(hits) == 0) {
    stop(
      paste("Cannot find file:", filename)
    )
  }

  if (length(hits) > 1) {
    print(hits)

    stop(
      paste("Multiple files with identical filename found:", filename)
    )
  }

  hits
}

get_native_focal_file <- function(variable, scale) {
  filename <- switch(
    variable,
    BUD=paste0("BUD1820_focal_", scale, ".tif"),
    DBF=paste0("DBF1820_focal_", scale, ".tif"),
    ENF=paste0("ENF1820_focal_", scale, ".tif"),
    VRM=paste0("VRM_focal_", scale, ".tif"),
    ELE=paste0("ELE_focal_", scale, ".tif")
  )

  find_native_file(filename)
}

# 6. Define A / B / C
A <- st_read(area_file, quiet=TRUE) %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(32653)

buf_8.9 <- st_buffer(A, dist=8900)

buf_17.8 <- st_buffer(A, dist=17800)

B <- st_make_valid(
  st_difference(buf_8.9, A)
)

C <- st_make_valid(
  st_difference(buf_17.8, buf_8.9)
)

area_summary_100m <- data.frame(
  zone=c("A", "B", "C"),
  area_km2=c(
    as.numeric(st_area(A)),
    as.numeric(st_area(B)),
    as.numeric(st_area(C))
  ) / 1e6
)

print(area_summary_100m)

write.csv(
  area_summary_100m,
  file.path(base_table_dir, "area_summary_100m.csv"),
  row.names=FALSE
)

st_write(
  st_sf(zone="A", geometry=A),
  file.path(base_vector_dir, "zone_A_occupied.gpkg"),
  delete_dsn=TRUE,
  quiet=TRUE
)

st_write(
  st_sf(zone="B", geometry=B),
  file.path(base_vector_dir, "zone_B_uncertainty.gpkg"),
  delete_dsn=TRUE,
  quiet=TRUE
)

st_write(
  st_sf(zone="C", geometry=C),
  file.path(base_vector_dir, "zone_C_comparison.gpkg"),
  delete_dsn=TRUE,
  quiet=TRUE
)

ABC <- st_buffer(A, dist=17800) %>%
  st_make_valid()

st_write(
  st_sf(zone="ABC", geometry=ABC),
  file.path(base_vector_dir, "study_extent_ABC.gpkg"),
  delete_dsn=TRUE,
  quiet=TRUE
)

# 7. Create 100-m analysis template
bbox_C <- st_bbox(C)

xmin_100 <- floor(
  as.numeric(bbox_C["xmin"]) / 100
) * 100

xmax_100 <- ceiling(
  as.numeric(bbox_C["xmax"]) / 100
) * 100

ymin_100 <- floor(
  as.numeric(bbox_C["ymin"]) / 100
) * 100

ymax_100 <- ceiling(
  as.numeric(bbox_C["ymax"]) / 100
) * 100

template_100m <- rast(
  xmin=xmin_100,
  xmax=xmax_100,
  ymin=ymin_100,
  ymax=ymax_100,
  resolution=100,
  crs="EPSG:32653"
)

names(template_100m) <- "template_100m"

values(template_100m) <- NA_real_

print(template_100m)
print(res(template_100m))
print(origin(template_100m))
print(ext(template_100m))

writeRaster(
  template_100m,
  file.path(base_raster_dir, "template_100m.tif"),
  overwrite=TRUE
)

# 8. Forest proportion at 100 m
lulc <- rast(lulc_file)

if (!same.crs(lulc, template_100m)) {
  stop("LULC CRS does not match template_100m.")
}

lulc_crop <- crop(lulc, template_100m, snap="out")

forest_classes <- c(6, 7, 8, 9, 11)

forest_fine <- ifel(
  lulc_crop == 0,
  NA,
  ifel(lulc_crop %in% forest_classes, 1, 0)
)

names(forest_fine) <- "forest"

forest_prop_100m <- resample(forest_fine, template_100m, method="average", threads=TRUE)

names(forest_prop_100m) <- "forest_prop"

stopifnot(
  compareGeom(forest_prop_100m, template_100m, stopOnError=FALSE)
)

writeRaster(
  forest_prop_100m,
  file.path(base_raster_dir, "forest_prop_100m.tif"),
  overwrite=TRUE
)

print(
  global(
    forest_prop_100m,
    c("min", "mean", "max"),
    na.rm=TRUE
  )
)

# 9. A / C forest-dominated fitting domain
A_v <- vect(A)
B_v <- vect(B)
C_v <- vect(C)

forest_A_100m <- mask(forest_prop_100m, A_v, touches=FALSE)

forest_C_100m <- mask(forest_prop_100m, C_v, touches=FALSE)

forest_threshold <- 0.5

eligible_A_100m <- forest_A_100m >= forest_threshold
eligible_C_100m <- forest_C_100m >= forest_threshold

domain_A_100m <- ifel(eligible_A_100m, 1, NA)

domain_C_100m <- ifel(eligible_C_100m, 1, NA)

eligible_domain_100m <- cover(domain_A_100m, domain_C_100m)

eligible_summary_100m <- data.frame(
  zone=c("A", "C"),
  eligible_cells=c(
    global(domain_A_100m, "sum", na.rm=TRUE)[1, 1],
    global(domain_C_100m, "sum", na.rm=TRUE)[1, 1]
  )
)

eligible_summary_100m$area_km2 <-
  eligible_summary_100m$eligible_cells * 0.01

print(eligible_summary_100m)

write.csv(
  eligible_summary_100m,
  file.path(base_table_dir, "eligible_cells_100m.csv"),
  row.names=FALSE
)

writeRaster(
  domain_A_100m,
  file.path(base_raster_dir, "domain_A_forest50_100m.tif"),
  overwrite=TRUE,
  datatype="INT1U"
)

writeRaster(
  domain_C_100m,
  file.path(base_raster_dir, "domain_C_forest50_100m.tif"),
  overwrite=TRUE,
  datatype="INT1U"
)

# 10. Extract eligible A / C cell numbers
A_cells_100m <- which(
  values(eligible_A_100m, mat=FALSE) == 1
)

C_cells_100m <- which(
  values(eligible_C_100m, mat=FALSE) == 1
)

cat(
  "A eligible cells:",
  length(A_cells_100m),
  "\n"
)

cat(
  "C eligible cells:",
  length(C_cells_100m),
  "\n"
)

# 11. Fixed balanced sample
# A = 10,000
# C = 10,000
n_sample <- 10000

stopifnot(
  length(A_cells_100m) >= n_sample
)

stopifnot(
  length(C_cells_100m) >= n_sample
)

set.seed(20261025)

A_sample_cells_100m <- sample(A_cells_100m, size=n_sample, replace=FALSE)

C_sample_cells_100m <- sample(C_cells_100m, size=n_sample, replace=FALSE)

A_xy <- xyFromCell(template_100m, A_sample_cells_100m)

C_xy <- xyFromCell(template_100m, C_sample_cells_100m)

A_sample_100m <- data.frame(
  cell=A_sample_cells_100m,
  x=A_xy[, 1],
  y=A_xy[, 2],
  zone="A",
  presence=1L
)

C_sample_100m <- data.frame(
  cell=C_sample_cells_100m,
  x=C_xy[, 1],
  y=C_xy[, 2],
  zone="C",
  presence=0L
)

sample_base_100m <- bind_rows(A_sample_100m, C_sample_100m)

sample_base_100m$forest_prop <-
  values(forest_prop_100m, mat=FALSE)[sample_base_100m$cell]

stopifnot(
  nrow(sample_base_100m) == 20000
)

stopifnot(
  sum(sample_base_100m$presence == 1) == 10000
)

stopifnot(
  sum(sample_base_100m$presence == 0) == 10000
)

if (anyNA(sample_base_100m)) {
  stop("NA detected in sample_base_100m.")
}

print(
  table(sample_base_100m$zone)
)

print(
  summary(sample_base_100m$forest_prop)
)

saveRDS(
  sample_base_100m,
  file.path(base_rds_dir, "sample_base_100m.rds")
)

write.csv(
  sample_base_100m,
  file.path(base_table_dir, "sample_base_100m.csv"),
  row.names=FALSE
)

# 12. Build native focal-file table
focal_file_table <- bind_rows(
  lapply(
    focal_variables,
    function(v) {
      bind_rows(
        lapply(
          focal_scales,
          function(sc) {
            f <- get_native_focal_file(v, sc)

            data.frame(variable=v, scale=sc, native_file=f, stringsAsFactors=FALSE)
          }
        )
      )
    }
  )
)

print(focal_file_table, row.names=FALSE)

write.csv(
  focal_file_table,
  file.path(base_table_dir, "native_focal_files_100m.csv"),
  row.names=FALSE
)

# 13. Native focal raster geometry QC
focal_native_info <- bind_rows(
  lapply(
    seq_len(
      nrow(focal_file_table)
    ),
    function(i) {
      rr <- rast(focal_file_table$native_file[i])

      data.frame(
        variable=focal_file_table$variable[i],
        scale=focal_file_table$scale[i],
        crs_ok=same.crs(rr, template_100m),
        res_x=res(rr)[1],
        res_y=res(rr)[2],
        xmin=xmin(rr),
        xmax=xmax(rr),
        ymin=ymin(rr),
        ymax=ymax(rr)
      )
    }
  )
)

print(focal_native_info)

if (!all(focal_native_info$crs_ok)) {
  stop("Some focal rasters have a CRS different from EPSG:32653.")
}

write.csv(
  focal_native_info,
  file.path(base_table_dir, "native_focal_geometry_100m.csv"),
  row.names=FALSE
)

# 14. Convert all focal rasters to common 100-m grid
# Temporary outputs only
focal_file_table$output_file <- NA_character_

for (i in seq_len(
  nrow(focal_file_table)
)) {
  v <- focal_file_table$variable[i]
  sc <- focal_file_table$scale[i]

  native_file <-
    focal_file_table$native_file[i]

  output_file <- file.path(
    temp_focal_dir,
    paste0(v, "_focal_", sc, "_100m.tif")
  )

  cat("Creating:", v, "| focal =", sc, "m | resolution = 100 m\n")

  rr <- rast(native_file)

  if (!same.crs(rr, template_100m)) {
    stop(
      paste("CRS mismatch:", v, sc)
    )
  }

  rr_crop <- crop(rr, template_100m, snap="out")

  rr_100m <- resample(rr_crop, template_100m, method="average", threads=TRUE)

  names(rr_100m) <- v

  if (!compareGeom(rr_100m, template_100m, stopOnError=FALSE)) {
    stop(
      paste("Geometry mismatch:", v, sc)
    )
  }

  writeRaster(rr_100m, output_file, overwrite=TRUE)

  focal_file_table$output_file[i] <-
    output_file

  rm(rr, rr_crop, rr_100m)

  gc()
}

write.csv(
  focal_file_table,
  file.path(base_table_dir, "focal_candidate_files_100m.csv"),
  row.names=FALSE
)

# 15. Prepare Main-road distance raster at 100 m
DistMainR_native <- rast(DistMainR_file)

if (!same.crs(DistMainR_native, template_100m)) {
  stop("DistMainR CRS does not match template_100m.")
}

DistMainR_crop <- crop(DistMainR_native, template_100m, snap="out")

DistMainR_100m <- resample(DistMainR_crop, template_100m, method="bilinear", threads=TRUE)

names(DistMainR_100m) <- "DistMainR"

stopifnot(
  compareGeom(DistMainR_100m, template_100m, stopOnError=FALSE)
)

writeRaster(
  DistMainR_100m,
  file.path(selected_raster_dir, "DistMainR_100m.tif"),
  overwrite=TRUE
)

# 16. Prepare Minor-road distance raster at 100 m
DistMinorR_native <- rast(DistMinorR_file)

if (!same.crs(DistMinorR_native, template_100m)) {
  stop("DistMinorR CRS does not match template_100m.")
}

DistMinorR_crop <- crop(DistMinorR_native, template_100m, snap="out")

DistMinorR_100m <- resample(DistMinorR_crop, template_100m, method="bilinear", threads=TRUE)

names(DistMinorR_100m) <- "DistMinorR"

stopifnot(
  compareGeom(DistMinorR_100m, template_100m, stopOnError=FALSE)
)

writeRaster(
  DistMinorR_100m,
  file.path(selected_raster_dir, "DistMinorR_100m.tif"),
  overwrite=TRUE
)

# 17. Prepare SIGNED forest-edge distance raster at 100 m
#
# IMPORTANT:
# Script 01 does NOT decide whether positive or negative means
# inside/outside forest.
#
# Script 02 will:
#   - determine the sign using land-cover information
#   - calculate absolute distance
#   - create inside/outside indicator
#   - evaluate abs(distance) x outside interaction
DistEdge_native <- rast(DistEdge_file)

if (!same.crs(DistEdge_native, template_100m)) {
  stop("DistEdge CRS does not match template_100m.")
}

DistEdge_crop <- crop(DistEdge_native, template_100m, snap="out")

DistEdge_100m <- resample(DistEdge_crop, template_100m, method="bilinear", threads=TRUE)

names(DistEdge_100m) <- "DistEdge_signed"

stopifnot(
  compareGeom(DistEdge_100m, template_100m, stopOnError=FALSE)
)

edge_range_100m <- global(
  DistEdge_100m,
  c("min", "max"),
  na.rm=TRUE
)

print(edge_range_100m)

if (edge_range_100m[1, 1] >= 0 || edge_range_100m[1, 2] <= 0) {
  stop("DistEdge does not appear to be a signed-distance raster.")
}

DistEdge_output_file <- file.path(selected_raster_dir, "DistEdge_signed_100m.tif")

writeRaster(DistEdge_100m, DistEdge_output_file, overwrite=TRUE)

sample_DistEdge_100m <-
  values(DistEdge_100m, mat=FALSE)[sample_base_100m$cell]

if (anyNA(sample_DistEdge_100m)) {
  stop("NA detected in DistEdge at sampled A/C cells.")
}

edge_sample_summary_100m <- data.frame(
  n=length(sample_DistEdge_100m),
  min=min(sample_DistEdge_100m),
  q25=as.numeric(
    quantile(sample_DistEdge_100m, 0.25)
  ),
  median=median(sample_DistEdge_100m),
  mean=mean(sample_DistEdge_100m),
  q75=as.numeric(
    quantile(sample_DistEdge_100m, 0.75)
  ),
  max=max(sample_DistEdge_100m),
  n_negative=sum(sample_DistEdge_100m < 0),
  n_zero=sum(sample_DistEdge_100m == 0),
  n_positive=sum(sample_DistEdge_100m > 0)
)

print(edge_sample_summary_100m)

write.csv(
  edge_sample_summary_100m,
  file.path(base_table_dir, "DistEdge_signed_sample_QC_100m.csv"),
  row.names=FALSE
)

# 18. Focal-scale selection function
#
# IMPORTANT:
# Each scale is evaluated using a simple LINEAR univariate model.
#
# DBF focal-scale selection also remains linear.
# DBF quadratic structure is evaluated later in Script 02.
test_focal_scales_100m <- function(variable, file_table, sample_base) {
  temp <- file_table[
    file_table$variable == variable,
  ]

  temp <- temp[
    order(temp$scale),
  ]

  results <- vector(
    "list",
    nrow(temp)
  )

  for (i in seq_len(
    nrow(temp)
  )) {
    sc <- temp$scale[i]

    rr <- rast(temp$output_file[i])

    if (!compareGeom(rr, template_100m, stopOnError=FALSE)) {
      stop(
        paste("Geometry mismatch:", variable, sc)
      )
    }

    x <- values(rr, mat=FALSE)[sample_base$cell]

    if (anyNA(x)) {
      stop(
        paste("NA detected:", variable, sc)
      )
    }

    dat <- data.frame(presence=sample_base$presence, x=x)

    mod <- glm(
      presence ~ x,
      data=dat,
      family=binomial()
    )

    cc <- summary(mod)$coefficients

    results[[i]] <- data.frame(
      variable=variable,
      scale=sc,
      n=nrow(dat),
      beta=unname(
        coef(mod)["x"]
      ),
      SE=cc[
        "x",
        "Std. Error"
      ],
      z=cc[
        "x",
        "z value"
      ],
      p=cc[
        "x",
        "Pr(>|z|)"
      ],
      AIC=AIC(mod),
      stringsAsFactors=FALSE
    )

    rm(rr, x, dat, mod)

    gc()
  }

  bind_rows(results)
}

# 19. Run focal-scale selection
scale_results_100m <- bind_rows(
  lapply(
    focal_variables,
    function(v) {
      cat("\nTesting focal scales:", v, "\n")

      test_focal_scales_100m(v, focal_file_table, sample_base_100m)
    }
  )
)

scale_results_100m <- scale_results_100m %>%
  group_by(.data$variable) %>%
  mutate(
    min_AIC=min(.data$AIC),
    delta_AIC=.data$AIC - .data$min_AIC
  ) %>%
  ungroup() %>%
  arrange(.data$variable, .data$scale)

print(
  tibble::as_tibble(scale_results_100m),
  n=Inf,
  width=Inf
)

write.csv(
  scale_results_100m,
  file.path(base_table_dir, "focal_scale_results_100m.csv"),
  row.names=FALSE
)

# 20. Best and competitive focal scales
best_scales_100m <- scale_results_100m %>%
  group_by(.data$variable) %>%
  slice_min(order_by=.data$AIC, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  arrange(.data$variable)

competitive_scales_100m <- scale_results_100m %>%
  filter(.data$delta_AIC <= 2) %>%
  arrange(.data$variable, .data$delta_AIC, .data$scale)

cat("\nBest focal scales:\n")

print(
  tibble::as_tibble(best_scales_100m),
  n=Inf,
  width=Inf
)

cat("\nCompetitive focal scales, delta AIC <= 2:\n")

print(
  tibble::as_tibble(competitive_scales_100m),
  n=Inf,
  width=Inf
)

write.csv(
  best_scales_100m,
  file.path(base_table_dir, "best_focal_scales_100m.csv"),
  row.names=FALSE
)

write.csv(
  competitive_scales_100m,
  file.path(base_table_dir, "competitive_focal_scales_100m.csv"),
  row.names=FALSE
)

# 21. Copy SELECTED focal rasters to permanent output folder
best_scale_lookup_100m <- setNames(best_scales_100m$scale, best_scales_100m$variable)

selected_focal_table_100m <- bind_rows(
  lapply(
    focal_variables,
    function(v) {
      sc <- best_scale_lookup_100m[
        v
      ]

      source_file <- focal_file_table$output_file[
        focal_file_table$variable == v &
          focal_file_table$scale == sc
      ]

      if (length(source_file) != 1) {
        stop(
          paste("Selected focal file not uniquely found:", v, sc)
        )
      }

      destination_file <- file.path(
        selected_raster_dir,
        paste0(v, "_focal_", sc, "_100m.tif")
      )

      ok <- file.copy(from=source_file, to=destination_file, overwrite=TRUE)

      if (!ok) {
        stop(
          paste("Failed to copy selected raster:", v)
        )
      }

      data.frame(variable=v, selected_scale=sc, file=destination_file, stringsAsFactors=FALSE)
    }
  )
)

print(selected_focal_table_100m)

write.csv(
  selected_focal_table_100m,
  file.path(base_table_dir, "selected_predictor_rasters_100m.csv"),
  row.names=FALSE
)

# 22. Geometry QC for all permanent predictors
permanent_predictor_files <- c(
  selected_focal_table_100m$file,
  file.path(selected_raster_dir, "DistMainR_100m.tif"),
  file.path(selected_raster_dir, "DistMinorR_100m.tif"),
  file.path(selected_raster_dir, "DistEdge_signed_100m.tif")
)

stopifnot(
  all(
    file.exists(permanent_predictor_files)
  )
)

predictor_geometry_QC <- bind_rows(
  lapply(
    permanent_predictor_files,
    function(f) {
      rr <- rast(f)

      data.frame(
        file=basename(f),
        geometry_match=compareGeom(rr, template_100m, stopOnError=FALSE),
        crs_match=same.crs(rr, template_100m),
        resolution_x=res(rr)[1],
        resolution_y=res(rr)[2],
        stringsAsFactors=FALSE
      )
    }
  )
)

print(predictor_geometry_QC)

if (!all(predictor_geometry_QC$geometry_match)) {
  stop("Some permanent predictor rasters do not match template_100m.")
}

if (!all(predictor_geometry_QC$crs_match)) {
  stop("Some permanent predictor rasters have a CRS mismatch.")
}

write.csv(
  predictor_geometry_QC,
  file.path(base_table_dir, "permanent_predictor_geometry_QC.csv"),
  row.names=FALSE
)

# 23. Save Stage-1 objects required by Script 02
stage1_objects_100m <- list(
  A=A,
  B=B,
  C=C,
  ABC=ABC,
  template_100m=template_100m,
  forest_threshold=forest_threshold,
  sample_base_100m=sample_base_100m,
  focal_scales=focal_scales,
  focal_variables=focal_variables,
  focal_file_table=focal_file_table,
  scale_results_100m=scale_results_100m,
  best_scales_100m=best_scales_100m,
  competitive_scales_100m=competitive_scales_100m,
  selected_focal_table_100m=selected_focal_table_100m,
  DistMainR_file=file.path(selected_raster_dir, "DistMainR_100m.tif"),
  DistMinorR_file=file.path(selected_raster_dir, "DistMinorR_100m.tif"),
  DistEdge_file=file.path(selected_raster_dir, "DistEdge_signed_100m.tif"),
  edge_range_100m=edge_range_100m
)

saveRDS(
  stage1_objects_100m,
  file.path(base_rds_dir, "RSF_stage1_objects.rds")
)

# 24. Save small supporting objects separately
saveRDS(
  sample_base_100m,
  file.path(base_rds_dir, "sample_base_100m.rds")
)

saveRDS(
  scale_results_100m,
  file.path(base_rds_dir, "focal_scale_results_100m.rds")
)

saveRDS(
  best_scales_100m,
  file.path(base_rds_dir, "best_focal_scales_100m.rds")
)

saveRDS(
  selected_focal_table_100m,
  file.path(base_rds_dir, "selected_predictor_rasters_100m.rds")
)

# 25. Remove TEMPORARY focal candidate rasters
#
# All candidate-scale rasters were generated only for focal-scale
# selection. Selected rasters have already been copied to
# selected_predictors.
temp_files <- list.files(temp_focal_dir, full.names=TRUE)

if (length(temp_files) > 0) {
  delete_ok <- file.remove(temp_files)

  if (!all(delete_ok)) {
    warning("Some temporary focal candidate files could not be deleted.")
  }
}

remaining_temp_files <- list.files(temp_focal_dir, full.names=TRUE)

cat(
  "Remaining temporary focal files:",
  length(remaining_temp_files),
  "\n"
)

# 26. Final QC
required_outputs <- c(
  file.path(base_raster_dir, "template_100m.tif"),
  file.path(base_raster_dir, "forest_prop_100m.tif"),
  file.path(selected_raster_dir, "DistMainR_100m.tif"),
  file.path(selected_raster_dir, "DistMinorR_100m.tif"),
  file.path(selected_raster_dir, "DistEdge_signed_100m.tif"),
  selected_focal_table_100m$file,
  file.path(base_rds_dir, "RSF_stage1_objects.rds")
)

missing_outputs <- required_outputs[
  !file.exists(required_outputs)
]

if (length(missing_outputs) > 0) {
  print(missing_outputs)

  stop("Script 01 finished with missing required outputs.")
}

cat("\n========================================\n")

cat("RSF Script 01 completed successfully.\n")

cat("========================================\n")

cat("Output root:\n", base_dir, "\n\n")

cat(
  "Fixed sample size:",
  nrow(sample_base_100m),
  "\n"
)

cat(
  "A sample:",
  sum(sample_base_100m$presence == 1),
  "\n"
)

cat(
  "C sample:",
  sum(sample_base_100m$presence == 0),
  "\n"
)

cat("\nSelected focal scales:\n")

print(
  best_scales_100m[
    ,
    c("variable", "scale", "AIC", "delta_AIC")
  ]
)

cat("\nPermanent predictor rasters:\n")

print(permanent_predictor_files)

cat("\nSigned edge-distance range:\n")

print(edge_range_100m)

cat(
  "\nStage-1 object:\n",
  file.path(base_rds_dir, "RSF_stage1_objects.rds"),
  "\n"
)

cat("========================================\n")

gc()
