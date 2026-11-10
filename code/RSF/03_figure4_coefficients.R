# Manuscript Figure 4: coefficients of the final landscape RSF.
# This script does not construct maps or supplementary figures.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
source(file.path("R", "00_config.R"))
create_rsf_directories()

model_file <- file.path(model_rds_dir, "landscape_RSF_final_100m.rds")
selection_file <- file.path(model_rds_dir, "RSF_final_analysis_objects.rds")
stopifnot(file.exists(model_file), file.exists(selection_file))
model <- readRDS(model_file)
selection <- readRDS(selection_file)
stopifnot(identical(as.character(selection$final_model_name), "base"))

coef_mat <- summary(model)$coefficients
labels <- c(
  "BUD_z"="Built-up area",
  "DBF_z"="Deciduous broadleaf forest",
  "I(DBF_z^2)"="Deciduous broadleaf forest²",
  "log_DistMainR_z"="Distance to main road",
  "log_DistMinorR_z"="Distance to minor road"
)
plot_data <- data.frame(
  term=rownames(coef_mat), estimate=coef_mat[, "Estimate"],
  se=coef_mat[, "Std. Error"], row.names=NULL
) %>%
  filter(.data$term != "(Intercept)") %>%
  mutate(
    lower=.data$estimate - 1.96 * .data$se,
    upper=.data$estimate + 1.96 * .data$se,
    label=unname(labels[.data$term]),
    label=factor(.data$label, levels=rev(unname(labels)))
  )
if (anyNA(plot_data$label) || nrow(plot_data) != 5L) {
  stop("The final RSF model does not contain the five expected coefficients.")
}

p <- ggplot(plot_data, aes(x=estimate, y=label)) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  geom_errorbar(aes(xmin=lower, xmax=upper), width=0.15, orientation="y", linewidth=0.6) +
  geom_point(size=2.5) +
  labs(x="Standardized coefficient estimate", y=NULL) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank())

file_prefix <- file.path(figure_dir, "Fig4_landscape_RSF_coefficients")
ggsave(paste0(file_prefix, ".png"), p, width=6.5, height=4.8, dpi=400, bg="white")
ggsave(paste0(file_prefix, ".pdf"), p, width=6.5, height=4.8, bg="white")
stopifnot(file.exists(paste0(file_prefix, ".png")), file.exists(paste0(file_prefix, ".pdf")))
cat("Figure 4 completed.\n")
