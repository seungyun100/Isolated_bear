# Common model specification and fixed-effect, two-stage IVW tools.
# Source after code/01_config.R. Does not load private data on its own.

issf_terms <- c(
  "sl_", "log_sl_", "cos_ta_", "log_sl_:cos_ta_",
  "z_end_ELE", "z_end_ELE2", "z_end_VRM", "z_end_VRM2",
  "z_end_PBUD", "z_end_PDBF", "z_end_PDBF2",
  "z_log_end_DistMainR", "z_log_end_DistMinorR", "z_log_end_DistEdgeAbs",
  "log_sl_:z_start_VRM", "log_sl_:z_start_PDBF", "log_sl_:z_log_start_DistMinorR"
)

issf_formula <- case_ ~ sl_ + log_sl_ + cos_ta_ + log_sl_:cos_ta_ +
  z_end_ELE + z_end_ELE2 + z_end_VRM + z_end_VRM2 + z_end_PBUD +
  z_end_PDBF + z_end_PDBF2 + z_log_end_DistMainR +
  z_log_end_DistMinorR + z_log_end_DistEdgeAbs +
  log_sl_:z_start_VRM + log_sl_:z_start_PDBF +
  log_sl_:z_log_start_DistMinorR + strata(step_id_id_year)

issf_groups <- data.frame(
  Sex=c("Female", "Female", "Male", "Male"),
  Season=c("Mating", "Hyperphagia", "Mating", "Hyperphagia")
)
issf_groups$key <- paste(issf_groups$Sex, issf_groups$Season, sep="__")

issf_expected_qc <- data.frame(
  Sex=c("Female", "Female", "Male", "Male"),
  Season=c("Hyperphagia", "Mating", "Hyperphagia", "Mating"),
  n_individuals=c(7L, 7L, 5L, 2L),
  n_tracks=c(17L, 10L, 10L, 4L),
  n_used=c(8781L, 5261L, 3398L, 1384L)
)

check_issf_data <- function(df, strict=TRUE) {
  required <- c(
    "Sex", "Season", "real_id", "id_season", "case_", "step_id_id_year",
    "sl_", "log_sl_", "cos_ta_", "ta_", "x1_", "y1_", "x2_", "y2_",
    "t1_", "t2_", "z_end_ELE", "z_end_ELE2", "z_end_VRM", "z_end_VRM2",
    "z_end_PBUD", "z_end_PDBF", "z_end_PDBF2", "z_log_end_DistMainR",
    "z_log_end_DistMinorR", "z_log_end_DistEdgeAbs", "z_start_VRM",
    "z_start_PDBF", "z_log_start_DistMinorR"
  )
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("Missing iSSF columns: ", paste(missing, collapse=", "))
  if (anyNA(df[, required, drop=FALSE])) stop("Missing values in mandatory iSSF fields.")
  if (any(!df$case_ %in% c(0, 1, FALSE, TRUE))) stop("case_ must be binary.")
  if (!all(df$Sex %in% c("Female", "Male"))) stop("Invalid Sex values.")
  if (!all(df$Season %in% c("Mating", "Hyperphagia"))) stop("Invalid Season values.")

  used <- df[df$case_ == 1L, , drop=FALSE]
  if (!nrow(used) || any(used$sl_ <= 20)) {
    stop("Used steps must be present and strictly >20 m.")
  }
  if (any(!is.finite(df$sl_)) || any(!is.finite(df$log_sl_)) ||
      any(!is.finite(df$cos_ta_))) stop("Non-finite movement predictors.")

  # Every choice set must include one used step and ten available steps.
  group_key <- paste(df$id_season, df$step_id_id_year, sep="::")
  used_count <- rowsum(as.integer(df$case_ == 1L), group_key, reorder=FALSE)
  available_count <- tabulate(match(group_key, unique(group_key)), nbins=length(unique(group_key)))
  if (any(used_count[, 1] != 1L)) stop("Each step stratum must have exactly one used step.")
  if (any(available_count != 11L)) stop("Each step stratum must have 11 choices.")

  qc <- used |>
    dplyr::group_by(.data$Sex, .data$Season) |>
    dplyr::summarise(
      n_individuals=dplyr::n_distinct(.data$real_id),
      n_tracks=dplyr::n_distinct(.data$id_season), n_used=dplyr::n(), .groups="drop"
    ) |>
    dplyr::arrange(.data$Sex, .data$Season)
  expected <- issf_expected_qc |>
    dplyr::arrange(.data$Sex, .data$Season)
  if (strict && !isTRUE(all.equal(as.data.frame(qc), expected, check.attributes=FALSE))) {
    print(as.data.frame(qc)); print(expected)
    stop("Final iSSF sample QC differs from the verified 41-track dataset.")
  }
  if (strict && dplyr::n_distinct(used$id_season) != 41L) stop("Expected 41 tracks.")
  if (any((used |> dplyr::count(.data$id_season))$n < 80L)) stop("A retained track has <80 used steps.")
  invisible(qc)
}

# Term-wise inverse-variance weighting, with propagated cross-term covariance.
pool_ivw <- function(fits, terms=issf_terms) {
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("Cannot pool an empty fit list.")
  p <- length(terms)
  b <- matrix(NA_real_, nrow=length(fits), ncol=p, dimnames=list(NULL, terms))
  v <- b
  for (i in seq_along(fits)) {
    fi <- fits[[i]]
    if (!setequal(names(fi$beta), terms)) stop("Wrong coefficient names in fit ", i)
    if (!is.matrix(fi$vcov) || !all(terms %in% rownames(fi$vcov))) {
      stop("Invalid covariance matrix in fit ", i)
    }
    b[i, ] <- fi$beta[terms]
    v[i, ] <- diag(fi$vcov[terms, terms, drop=FALSE])
  }
  if (any(!is.finite(b)) || any(!is.finite(v) | v <= 0)) stop("Invalid IVW estimates/variances.")
  w <- 1 / v
  w <- sweep(w, 2L, colSums(w), "/")
  beta <- setNames(colSums(w * b), terms)
  vc <- matrix(0, p, p, dimnames=list(terms, terms))
  for (i in seq_along(fits)) {
    vc <- vc + (w[i, ] %o% w[i, ]) * fits[[i]]$vcov[terms, terms, drop=FALSE]
  }
  if (any(!is.finite(vc)) || any(diag(vc) <= 0)) stop("Invalid IVW covariance.")
  list(beta=beta, vcov=vc, se=sqrt(diag(vc)), n_contributors=length(fits))
}

fit_issf_track <- function(df) {
  if (sum(df$case_ == 1L) < 80L) stop("Track <80 used steps.")
  # Uses the same conditional logistic likelihood as amt::fit_issf().
  fit <- survival::clogit(issf_formula, data=df, method="efron")
  beta <- stats::coef(fit)
  cov <- stats::vcov(fit)
  if (!setequal(names(beta), issf_terms)) stop("Model did not estimate all 17 terms.")
  beta <- beta[issf_terms]
  cov <- cov[issf_terms, issf_terms, drop=FALSE]
  if (any(!is.finite(beta)) || any(!is.finite(cov)) || any(diag(cov) <= 0)) {
    stop("Track iSSF has singular/invalid estimates.")
  }
  list(beta=beta, vcov=cov, n_used=sum(df$case_ == 1L))
}

pool_to_population <- function(track_fits, track_index, terms=issf_terms) {
  individual <- list()
  individuals <- unique(track_index[, c("Sex", "Season", "real_id")])
  for (i in seq_len(nrow(individuals))) {
    sub <- individuals[i, , drop=FALSE]
    keep <- track_index$Sex == sub$Sex & track_index$Season == sub$Season &
      track_index$real_id == sub$real_id
    track_ids <- track_index$id_season[keep]
    k <- paste(sub$Sex, sub$Season, sub$real_id, sep="__")
    individual[[k]] <- pool_ivw(track_fits[track_ids], terms)
  }
  population <- list()
  for (i in seq_len(nrow(issf_groups))) {
    sex <- issf_groups$Sex[i]; season <- issf_groups$Season[i]
    sub <- individuals[individuals$Sex == sex & individuals$Season == season, , drop=FALSE]
    keys <- paste(sub$Sex, sub$Season, sub$real_id, sep="__")
    population[[issf_groups$key[i]]] <- pool_ivw(individual[keys], terms)
  }
  list(individual=individual, population=population)
}

make_coef_table <- function(fit_list) {
  dplyr::bind_rows(lapply(names(fit_list), function(key) {
    fit <- fit_list[[key]]
    data.frame(key=key, term=names(fit$beta), estimate=as.numeric(fit$beta),
               std_error=as.numeric(fit$se %||% sqrt(diag(fit$vcov))),
               variance=as.numeric(diag(fit$vcov)), row.names=NULL)
  }))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
