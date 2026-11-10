# Shikoku landscape RSF
# Manuscript Fig. 5: habitat covariates across ten RSF suitability bins
#
# Generates Fig. 5 only; no supplementary or exploratory figures.
# 1. Packages
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)

gc()

# Project-relative paths (configure in code/01_config.R)
source(file.path("code", "01_config.R"))
create_rsf_directories()
predictor_dir <- selected_raster_dir

# 3. Input files
bin_raster_file <- file.path(
  model_raster_dir,
  "relative_suitability_equal_area_bins_ABC_100m.tif"
)

BUD_file <- file.path(predictor_dir, "BUD_focal_4000_100m.tif")

DBF_file <- file.path(predictor_dir, "DBF_focal_4000_100m.tif")

ENF_file <- file.path(predictor_dir, "ENF_focal_2000_100m.tif")

VRM_file <- file.path(predictor_dir, "VRM_focal_4000_100m.tif")

ELE_file <- file.path(predictor_dir, "ELE_focal_4000_100m.tif")

DistMainR_file <- file.path(predictor_dir, "DistMainR_100m.tif")

DistMinorR_file <- file.path(predictor_dir, "DistMinorR_100m.tif")

DistEdge_file <- file.path(predictor_dir, "DistEdge_signed_100m.tif")

required_files <- c(
  bin_raster_file,
  BUD_file,
  DBF_file,
  ENF_file,
  VRM_file,
  ELE_file,
  DistMainR_file,
  DistMinorR_file,
  DistEdge_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  cat("\nMissing files:\n")
  print(missing_files)
  stop("Required files are missing.")
}

cat("\nAll required files found.\n")

# 4. Read rasters
bin_raster <- rast(bin_raster_file)

BUD_r <- rast(BUD_file)
DBF_r <- rast(DBF_file)
ENF_r <- rast(ENF_file)
VRM_r <- rast(VRM_file)
ELE_r <- rast(ELE_file)

DistMainR_r <- rast(DistMainR_file)
DistMinorR_r <- rast(DistMinorR_file)
DistEdge_r <- rast(DistEdge_file)

cat("\nRasters loaded.\n")

# 5. Geometry QC
#
# Predictor rasters may cover a larger extent.
# CRS, resolution, and origin must match.
predictor_list <- list(
  BUD=BUD_r,
  DBF=DBF_r,
  ENF=ENF_r,
  VRM=VRM_r,
  ELE=ELE_r,
  DistMainR=DistMainR_r,
  DistMinorR=DistMinorR_r,
  DistEdge=DistEdge_r
)

geometry_qc <- data.frame(
  variable=names(predictor_list),
  same_crs=sapply(
    predictor_list,
    function(rr) same.crs(rr, bin_raster)
  ),
  same_resolution=sapply(
    predictor_list,
    function(rr) isTRUE(all.equal(res(rr), res(bin_raster)))
  ),
  same_origin=sapply(
    predictor_list,
    function(rr) isTRUE(all.equal(origin(rr), origin(bin_raster)))
  )
)

cat("\nGeometry QC:\n")
print(geometry_qc)

if (
  !all(geometry_qc$same_crs) ||
  !all(geometry_qc$same_resolution) ||
  !all(geometry_qc$same_origin)
) {
  stop(
    paste(
      "Some predictor rasters do not match",
      "CRS, resolution, or origin of the bin raster."
    )
  )
}

cat("\nCRS, resolution, and origin all match.\n")

# 6. Crop and mask helper
crop_mask_to_bin <- function(rr, rr_name, bin_raster) {
  cat("Processing:", rr_name, "\n")

  rr_crop <- crop(rr, bin_raster, snap="near")

  rr_mask <- mask(rr_crop, bin_raster)

  if (!compareGeom(rr_mask, bin_raster, stopOnError=FALSE)) {
    stop(
      paste("Geometry mismatch after crop/mask:", rr_name)
    )
  }

  rr_mask
}

# 7. Crop and mask predictors
BUD_m <- crop_mask_to_bin(BUD_r, "BUD", bin_raster)

DBF_m <- crop_mask_to_bin(DBF_r, "DBF", bin_raster)

ENF_m <- crop_mask_to_bin(ENF_r, "ENF", bin_raster)

VRM_m <- crop_mask_to_bin(VRM_r, "VRM", bin_raster)

ELE_m <- crop_mask_to_bin(ELE_r, "ELE", bin_raster)

DistMainR_m <- crop_mask_to_bin(DistMainR_r, "DistMainR", bin_raster)

DistMinorR_m <- crop_mask_to_bin(DistMinorR_r, "DistMinorR", bin_raster)

DistEdge_m <- crop_mask_to_bin(DistEdge_r, "DistEdge", bin_raster)

cat("\nAll predictors processed successfully.\n")

# 8. Build raster stack
cov_stack <- c(
  bin_raster,
  ELE_m,
  VRM_m,
  DBF_m,
  ENF_m,
  DistEdge_m,
  BUD_m,
  DistMainR_m,
  DistMinorR_m
)

names(cov_stack) <- c("bin", "ELE", "VRM", "DBF", "ENF", "DistEdge", "BUD", "DistMainR", "DistMinorR")

# 9. Convert to dataframe
cov_df <- as.data.frame(cov_stack, na.rm=FALSE)

cov_df <- cov_df %>%
  filter(
    !is.na(.data$bin)
  ) %>%
  mutate(
    bin=as.integer(.data$bin)
  )

if (!all(cov_df$bin %in% 1:10)) {
  stop("Unexpected RSF-bin values detected.")
}

cat(
  "\nNumber of cells in study extent:",
  nrow(cov_df),
  "\n"
)

cat("\nCells per bin:\n")
print(table(cov_df$bin))

# 10. Long format
cov_long <- cov_df %>%
  select(bin, ELE, VRM, DBF, ENF, DistEdge, BUD, DistMainR, DistMinorR) %>%
  pivot_longer(cols=-bin, names_to="variable", values_to="value") %>%
  filter(
    !is.na(.data$value)
  )

# 11. Variable order and labels
variable_levels <- c("ELE", "VRM", "DBF", "ENF", "DistEdge", "BUD", "DistMainR", "DistMinorR")

variable_labels <- c(
  ELE="Elevation (m)",
  VRM="VRM",
  DBF="Deciduous broadleaf\nforest (proportion)",
  ENF="Evergreen needleleaf\nforest (proportion)",
  DistEdge="Distance to\nforest edge (m)",
  BUD="Built-up area\n(proportion)",
  DistMainR="Distance to\nmain road (m)",
  DistMinorR="Distance to\nminor road (m)"
)

cov_long <- cov_long %>%
  mutate(
    variable=factor(.data$variable, levels=variable_levels),
    variable_label=factor(
      variable_labels[
        as.character(.data$variable)
      ],
      levels=variable_labels[
        variable_levels
      ]
    )
  )

# 12. Study-extent mean
mean_table <- cov_long %>%
  group_by(.data$variable, .data$variable_label) %>%
  summarise(
    mean_value=mean(.data$value, na.rm=TRUE),
    .groups="drop"
  )

cat("\nStudy-extent means:\n")
print(tibble::as_tibble(mean_table), n=Inf, width=Inf)

# 13. Bin summaries
#
# whiskers = 5th and 95th percentiles
# box      = 25th and 75th percentiles
# middle   = median
sum_table <- cov_long %>%
  group_by(.data$variable, .data$variable_label, .data$bin) %>%
  summarise(
    n=dplyr::n(),
    q05=as.numeric(
      quantile(.data$value, 0.05, na.rm=TRUE)
    ),
    q25=as.numeric(
      quantile(.data$value, 0.25, na.rm=TRUE)
    ),
    q50=as.numeric(
      quantile(.data$value, 0.50, na.rm=TRUE)
    ),
    q75=as.numeric(
      quantile(.data$value, 0.75, na.rm=TRUE)
    ),
    q95=as.numeric(
      quantile(.data$value, 0.95, na.rm=TRUE)
    ),
    .groups="drop"
  )

# 14. Suitable-bin shading
#
# Suitable landscape = bins 8-10
#
# Bin centers:
# 1, 2, ..., 10
#
# Therefore shading begins halfway between bins 7 and 8,
# and extends exactly to the right plotting boundary.
suitable_min_bin <- 8

shade_df <- data.frame(xmin=suitable_min_bin - 0.5, xmax=10.5)

# 15. Figure
p_cov_bins <- ggplot() +

  # ----------------------------------------------------------
# Suitable landscape shading
# ----------------------------------------------------------
geom_rect(
  data=shade_df,
  aes(xmin=.data$xmin, xmax=.data$xmax, ymin=-Inf, ymax=Inf),
  inherit.aes=FALSE,
  fill="grey75",
  alpha=0.22
) +

  # ----------------------------------------------------------
# Study-extent mean
# ----------------------------------------------------------
geom_hline(
  data=mean_table,
  aes(yintercept=.data$mean_value),
  inherit.aes=FALSE,
  linetype="dashed",
  linewidth=0.55,
  color="grey40"
) +

  # ----------------------------------------------------------
# 5th-95th percentile whisker
# ----------------------------------------------------------
geom_segment(
  data=sum_table,
  aes(x=.data$bin, xend=.data$bin, y=.data$q05, yend=.data$q95),
  linewidth=0.40,
  color="grey35"
) +

  # ----------------------------------------------------------
# Lower whisker cap
# ----------------------------------------------------------
geom_segment(
  data=sum_table,
  aes(x=.data$bin - 0.13, xend=.data$bin + 0.13, y=.data$q05, yend=.data$q05),
  linewidth=0.40,
  color="grey35"
) +

  # ----------------------------------------------------------
# Upper whisker cap
# ----------------------------------------------------------
geom_segment(
  data=sum_table,
  aes(x=.data$bin - 0.13, xend=.data$bin + 0.13, y=.data$q95, yend=.data$q95),
  linewidth=0.40,
  color="grey35"
) +

  # ----------------------------------------------------------
# IQR box
# ----------------------------------------------------------
geom_rect(
  data=sum_table,
  aes(xmin=.data$bin - 0.27, xmax=.data$bin + 0.27, ymin=.data$q25, ymax=.data$q75),
  inherit.aes=FALSE,
  fill="white",
  color="grey30",
  linewidth=0.55
) +

  # ----------------------------------------------------------
# Median
# ----------------------------------------------------------
geom_segment(
  data=sum_table,
  aes(x=.data$bin - 0.27, xend=.data$bin + 0.27, y=.data$q50, yend=.data$q50),
  linewidth=0.75,
  color="grey10"
) +

  # ----------------------------------------------------------
# Panels
# ----------------------------------------------------------
facet_grid(
  rows=vars(variable_label),
  scales="free_y",
  switch="y"
) +

  # ----------------------------------------------------------
# X axis
#
# 0.5 and 10.5 are the exact outer edges of bins 1 and 10.
# expand = c(0,0) removes additional ggplot padding.
# ----------------------------------------------------------
scale_x_continuous(
  breaks=1:10,
  limits=c(0.5, 10.5),
  expand=c(0, 0)
) +

  labs(
    x="Landscape suitability bin (low to high)",
    y=NULL
  ) +

  theme_bw(base_size=14) +

  theme(
    # --------------------------------------------------------
    # Panel
    # --------------------------------------------------------
    panel.background=element_rect(fill="white", color=NA),

    panel.border=element_rect(fill=NA, color="grey55", linewidth=0.55),

    # --------------------------------------------------------
    # Grid
    # --------------------------------------------------------
    panel.grid.major.x=element_blank(),

    panel.grid.minor=element_blank(),

    panel.grid.major.y=element_line(color="grey87", linewidth=0.45),

    # --------------------------------------------------------
    # Variable labels
    # --------------------------------------------------------
    strip.background=element_blank(),

    strip.text.y.left=element_text(
      angle=0,
      size=14,
      color="grey15",
      lineheight=0.95,
      margin=margin(r=18)
    ),

    strip.placement="outside",

    # --------------------------------------------------------
    # Axis
    # --------------------------------------------------------
    axis.text.x=element_text(size=12, color="grey20"),

    axis.text.y=element_text(size=11, color="grey20"),

    axis.title.x=element_text(
      size=14,
      color="grey15",
      margin=margin(t=8)
    ),

    # --------------------------------------------------------
    # Panel spacing
    # --------------------------------------------------------
    panel.spacing.y=grid::unit(0.12, "lines"),

    plot.margin=margin(t=8, r=12, b=8, l=8)
  )

# 16. Display
print(p_cov_bins)

# 17. Save summary tables
summary_file <- file.path(figure_dir, "Fig5_covariate_bins_summary.csv")

mean_file <- file.path(figure_dir, "Fig5_study_extent_means.csv")

write.csv(sum_table, summary_file, row.names=FALSE)

write.csv(mean_table, mean_file, row.names=FALSE)

# 18. Save PNG and PDF
figure_file_png <- file.path(figure_dir, "Fig5_habitat_covariates_by_bin.png")

figure_file_pdf <- file.path(figure_dir, "Fig5_habitat_covariates_by_bin.pdf")

ggsave(
  filename=figure_file_png,
  plot=p_cov_bins,
  width=9.4,
  height=12.8,
  units="in",
  dpi=400,
  bg="white"
)

ggsave(
  filename=figure_file_pdf,
  plot=p_cov_bins,
  width=9.4,
  height=12.8,
  units="in",
  bg="white"
)

# 19. Final QC
stopifnot(
  file.exists(figure_file_png),
  file.exists(figure_file_pdf),
  file.exists(summary_file),
  file.exists(mean_file)
)

cat("\n========================================\n")

cat("Manuscript Figure 5 completed successfully.\n")

cat("========================================\n")

cat("\nPNG:\n", figure_file_png, "\n")

cat("\nPDF:\n", figure_file_pdf, "\n")

cat("\nBin summary:\n", summary_file, "\n")

cat("\nStudy-extent means:\n", mean_file, "\n")

cat("========================================\n")

gc()
