# Manuscript Figure 6: population-level iSSF log-RSS by sex and season.
# Two-stage IVW coefficient covariance is propagated to confidence intervals.
# No other manuscript or supplementary figures are generated.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
source(file.path("code", "01_config.R"))
source(file.path("code", "iSSF", "01_issf_common.R"))
create_issf_directories()

population_file <- file.path(issf_model_dir, "population_fits.rds")
if (!all(file.exists(c(issf_data_file, issf_scaling_file, population_file)))) {
  stop("Run code/iSSF/02_fit_issf.R and provide private iSSF input files first.")
}
df <- readRDS(issf_data_file)
check_issf_data(df, strict=TRUE)
sc <- readRDS(issf_scaling_file)
population <- readRDS(population_file)

spec <- list(
  ELE=list(raw="end_ELE", scale="end_ELE", term="z_end_ELE", quadratic="z_end_ELE2", log=FALSE, label="Elevation (m)"),
  VRM=list(raw="end_VRM", scale="end_VRM", term="z_end_VRM", quadratic="z_end_VRM2", log=FALSE, label="VRM"),
  PDBF=list(raw="end_PDBF", scale="end_PDBF", term="z_end_PDBF", quadratic="z_end_PDBF2", log=FALSE, label="Deciduous broadleaf forest"),
  PBUD=list(raw="end_PBUD", scale="end_PBUD", term="z_end_PBUD", quadratic=NULL, log=FALSE, label="Built-up area"),
  DistMainR=list(raw="end_DistMainR", scale="log_end_DistMainR", term="z_log_end_DistMainR", quadratic=NULL, log=TRUE, label="Distance to main road (m)"),
  DistMinorR=list(raw="end_DistMinorR", scale="log_end_DistMinorR", term="z_log_end_DistMinorR", quadratic=NULL, log=TRUE, label="Distance to minor road (m)"),
  DistEdgeAbs=list(raw="end_DistEdgeAbs", scale="log_end_DistEdgeAbs", term="z_log_end_DistEdgeAbs", quadratic=NULL, log=TRUE, label="Distance to forest edge (m)")
)

curve_for_group <- function(sex, season, variable) {
  key <- paste(sex, season, sep="__")
  fit <- population[[key]]
  sp <- spec[[variable]]
  if (is.null(fit)) stop("Missing sex-season group: ", key)
  if (!sp$raw %in% names(df)) stop("Missing raw iSSF endpoint variable: ", sp$raw)
  use_values <- df[[sp$raw]][df$case_ == 1L & df$Sex == sex & df$Season == season]
  group_values <- df[[sp$raw]][df$Sex == sex & df$Season == season]
  group_values <- group_values[is.finite(group_values)]
  use_values <- use_values[is.finite(use_values)]
  if (!length(use_values) || !length(group_values)) stop("Missing observed/available predictor range.")
  x <- seq(as.numeric(stats::quantile(group_values, 0.01)),
           as.numeric(stats::quantile(group_values, 0.99)), length.out=200L)
  ref <- mean(use_values)
  transform_z <- function(v) {
    z <- if (sp$log) log(pmax(v, 0.0001)) else v
    (z - sc[[sp$scale]][["mean"]]) / sc[[sp$scale]][["sd"]]
  }
  z <- transform_z(x)
  z_ref <- transform_z(ref)
  terms <- sp$term
  contrasts <- matrix(z - z_ref, ncol=1L)
  if (!is.null(sp$quadratic)) {
    terms <- c(terms, sp$quadratic)
    contrasts <- cbind(z - z_ref, z^2 - z_ref^2)
  }
  beta <- fit$beta[terms]
  vc <- fit$vcov[terms, terms, drop=FALSE]
  if (any(!is.finite(beta)) || any(!is.finite(vc))) stop("Invalid coefficient covariance.")
  effect <- as.numeric(contrasts %*% beta)
  variance <- rowSums((contrasts %*% vc) * contrasts)
  se <- sqrt(pmax(variance, 0))
  data.frame(Sex=sex, Season=season, variable=variable,
             x=x, reference=ref, log_RSS=effect,
             lower=effect - 1.96 * se, upper=effect + 1.96 * se)
}

curves <- dplyr::bind_rows(lapply(issf_groups$key, function(group_key) {
  parts <- strsplit(group_key, "__", fixed=TRUE)[[1]]
  dplyr::bind_rows(lapply(names(spec), function(v) curve_for_group(parts[1], parts[2], v)))
}))
write.csv(curves, file.path(issf_result_dir, "figure6_logRSS_curves.csv"), row.names=FALSE)

plot_order <- c("ELE", "VRM", "PDBF", "DistEdgeAbs", "PBUD", "DistMainR", "DistMinorR")
for (season_i in c("Mating", "Hyperphagia")) {
  season_data <- curves[curves$Season == season_i, , drop=FALSE]
  panels <- lapply(plot_order, function(variable_i) {
    dat <- season_data[season_data$variable == variable_i, , drop=FALSE]
    ggplot(dat, aes(x=x, y=log_RSS, colour=Sex, fill=Sex)) +
      geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.15, linewidth=0, colour=NA) +
      geom_line(linewidth=0.8) +
      labs(x=spec[[variable_i]]$label, y=NULL) +
      theme_bw(base_size=10) +
      theme(panel.grid.minor=element_blank(), legend.position="bottom")
  })
  # Seven habitat panels plus a blank eighth panel, as used in the manuscript layout.
  p <- patchwork::wrap_plots(c(panels, list(patchwork::plot_spacer())), ncol=4L, guides="collect") +
    patchwork::plot_annotation(title=paste("iSSF log-RSS:", season_i)) &
    theme(legend.position="bottom")
  prefix <- file.path(issf_figure_dir, paste0("Fig6_", tolower(season_i), "_logRSS"))
  ggsave(paste0(prefix, ".png"), p, width=14, height=7.2, dpi=400, bg="white")
  ggsave(paste0(prefix, ".pdf"), p, width=14, height=7.2, bg="white")
}
cat("Figure 6 completed for both seasons.\n")
