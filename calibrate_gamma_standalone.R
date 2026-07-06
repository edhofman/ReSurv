## Standalone calibration script for a revised Gamma simulation scenario.
##
## Purpose:
##   Try a small number of Gamma DGP parameter combinations and compare
##   CL, COX, and a fixed shallow-tree XGB. This script does not modify
##   ReSurv and does not use targets or ReSurvCV.
##
## Run:
##   Rscript calibrate_gamma_standalone.R
##
## Optional environment overrides:
##   GAMMA_CALIB_REPS=5
##   GAMMA_CALIB_PERIOD_EXPOSURE=200
##   GAMMA_CALIB_PERIOD_FREQUENCY=0.05
##   GAMMA_CALIB_XGB_NROUNDS=500
##   GAMMA_CALIB_DT_THREADS=auto        # default: all logical cores minus one
##   GAMMA_CALIB_XGB_NTHREADS=auto      # default: same as GAMMA_CALIB_DT_THREADS
##   GAMMA_CALIB_OUT_DIR=results/gamma_calibration_simple

setwd("~/random_scripts/resurv-final")

suppressPackageStartupMessages({
  library(data.table)
  library(ReSurv)
})

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) as.integer(default) else out
}

get_env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) as.numeric(default) else out
}

get_env_threads <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value) || identical(tolower(value), "auto")) {
    if (!is.null(default)) return(as.integer(default))
    cores <- parallel::detectCores(logical = TRUE)
    if (is.na(cores) || cores < 1L) cores <- 1L
    return(as.integer(max(1L, cores - 1L)))
  }
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) {
    cores <- parallel::detectCores(logical = TRUE)
    if (is.na(cores) || cores < 1L) cores <- 1L
    return(as.integer(max(1L, cores - 1L)))
  }
  as.integer(out)
}

configure_threads <- function(dt_threads, xgb_threads) {
  data.table::setDTthreads(as.integer(dt_threads))
  message(sprintf(
    "Using data.table threads=%d; XGB nthread=%d",
    data.table::getDTthreads(), as.integer(xgb_threads)
  ))
  invisible(TRUE)
}

## -------------------------------------------------------------------------
## Candidate DGP parameter combinations.
## -------------------------------------------------------------------------

make_candidates <- function() {
  data.table(
    candidate = c("G1_mild", "G2_base", "G3_stronger_break", "G4_compact_age", "G5_stronger_all"),
    beta_break = c(-0.8, -1.0, -1.2, -1.2, -1.4),
    age_mean = c(45, 45, 45, 45, 45),
    age_sd = c(5, 5, 5, 4, 5),
    age_min = c(30, 30, 30, 32, 30),
    age_max = c(60, 60, 60, 58, 60),
    age_cut1 = c(42, 42, 42, 42, 42),
    age_cut2 = c(49, 49, 49, 49, 49),
    age_low = c(-0.10, -0.15, -0.15, -0.10, -0.20),
    age_mid = c(0.00, 0.00, 0.00, 0.00, 0.00),
    age_high = c(0.20, 0.25, 0.30, 0.25, 0.35)
  )
}

## -------------------------------------------------------------------------
## Standalone revised Gamma simulator.
## -------------------------------------------------------------------------

rtfwd_inverse <- function(n, alpha, beta, lambda, k, b) {
  u <- stats::runif(n)
  (-log(u) / (beta^alpha * lambda^(alpha * k)) + b^(-alpha * k))^(1 / (-alpha * k))
}

simulate_gamma_candidate <- function(candidate_row,
                                     random_seed,
                                     years = 4L,
                                     time_unit = 1 / 360,
                                     period_exposure = 80L,
                                     period_frequency = 0.05,
                                     beta0 = 1.15129,
                                     beta1 = 1.95601,
                                     break_day = NULL) {
  set.seed(as.integer(random_seed))

  max_day <- as.integer(years / time_unit)
  if (is.null(break_day)) break_day <- floor(max_day / 2)

  ## Expected claims per accident day and claim type.
  lambda_n <- as.numeric(period_exposure) * as.numeric(period_frequency)

  ## Draw claim counts for each accident day and claim type.
  grid <- data.table(
    AP = rep(seq_len(max_day), times = 2L),
    claim_type = rep(c(0L, 1L), each = max_day)
  )
  grid[, n_claims := stats::rpois(.N, lambda = lambda_n)]

  claims <- grid[n_claims > 0L, .(
    AP = rep(AP, n_claims),
    claim_type = rep(claim_type, n_claims)
  )]

  if (nrow(claims) == 0L) {
    stop("No claims simulated. Increase period_exposure or period_frequency.")
  }

  ## Continuous accident time inside the accident day.
  claims[, AT := AP - 1 + stats::runif(.N)]

  ## Rounded numeric age. The models see this as a numeric continuous feature.
  claims[, age_cont := stats::rnorm(.N, mean = candidate_row$age_mean, sd = candidate_row$age_sd)]
  claims[, age := round(pmin(pmax(age_cont, candidate_row$age_min), candidate_row$age_max))]

  ## Threshold-based age effect. We do not pass the bands to the models.
  claims[, age_effect := data.table::fcase(
    age < candidate_row$age_cut1, candidate_row$age_low,
    age < candidate_row$age_cut2, candidate_row$age_mid,
    default = candidate_row$age_high
  )]

  ## Revised Gamma log-risk component:
  ##   - claim_type main effect;
  ##   - claim_type-specific structural break after break_day;
  ##   - threshold age effect, observed only through numeric rounded age.
  claims[, phi :=
    beta0 * as.numeric(claim_type == 0L) +
    beta1 * as.numeric(claim_type == 1L) +
    candidate_row$beta_break * as.numeric(claim_type == 1L & AP > break_day) +
    age_effect
  ]

  alpha <- 0.5
  beta <- 2 * 30
  k <- 1
  b <- max_day
  lambda_delay <- 0.1 * exp(claims$phi)^(1 / alpha)

  claims[, delay := rtfwd_inverse(.N, alpha = alpha, beta = beta,
                                  lambda = lambda_delay, k = k, b = b)]
  claims[, RT := AT + delay]
  claims[, RP := ceiling(RT)]
  claims[, claim_number := seq_len(.N)]

  out <- claims[, .(
    claim_number = claim_number,
    claim_type = claim_type,
    age = age,
    AP = as.integer(AP),
    RP = as.integer(RP)
  )]

  out[]
}

prepare_for_resurv <- function(dat, years = 4L) {
  ReSurv::IndividualDataPP(
    data = dat,
    id = NULL,
    categorical_features = "claim_type",
    continuous_features = c("AP", "age"),
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "days",
    output_time_granularity = "quarters",
    years = as.integer(years),
    continuous_features_spline = NULL,
    calendar_period_extrapolation = FALSE
  )
}

fixed_xgb_hparameters <- function(nrounds = 500L, nthread = max(1L, parallel::detectCores(logical = TRUE) - 1L)) {
  list(
    params = list(
      booster = "gbtree",
      eta = 0.05,
      max_depth = 2L,
      min_child_weight = 10,
      subsample = 0.85,
      lambda = 1,
      alpha = 0,
      nthread = as.integer(nthread)
    ),
    nrounds = as.integer(nrounds),
    early_stopping_rounds = 50L,
    print_every_n = 50L,
    verbose = FALSE
  )
}

actual_ibnr_count <- function(dat, years = 4L, time_unit = 1 / 360) {
  max_day <- as.integer(years / time_unit)
  z <- copy(as.data.table(dat))
  z[, DP := RP - AP + 1L]
  nrow(z[AP >= 1L & AP <= max_day & DP >= 1L & DP <= max_day & RP > max_day])
}

fit_and_score_one <- function(candidate_row,
                              rep_id,
                              years,
                              time_unit,
                              period_exposure,
                              period_frequency,
                              xgb_nrounds,
                              xgb_nthreads) {
  sim_seed <- as.integer(rep_id)
  fit_seed <- as.integer(100000L + rep_id)

  dat <- simulate_gamma_candidate(
    candidate_row = candidate_row,
    random_seed = sim_seed,
    years = years,
    time_unit = time_unit,
    period_exposure = period_exposure,
    period_frequency = period_frequency
  )

  ibnr_n <- actual_ibnr_count(dat, years = years, time_unit = time_unit)
  observed_n <- nrow(dat[RP <= as.integer(years / time_unit)])

  idata <- prepare_for_resurv(dat, years = years)

  message(sprintf(
    "candidate=%s rep=%d: n=%d observed_upper=%d actual_ibnr=%d",
    candidate_row$candidate, rep_id, nrow(dat), observed_n, ibnr_n
  ))

  cat("started fitting")
  system.time(fit_cox <- tryCatch(
    ReSurv::ReSurv(
      idata,
      hazard_model = "COX",
      random_seed = fit_seed,
      grouping_method = "probability",
      simplifier = TRUE
    ),
    error = function(e) e
  ))

cat("done with cox")
  system.time(fit_xgb <- tryCatch(
    ReSurv::ReSurv(
      idata,
      hazard_model = "XGB",
      random_seed = fit_seed,
      hparameters = fixed_xgb_hparameters(nrounds = xgb_nrounds, nthread = xgb_nthreads),
      grouping_method = "probability",
      simplifier = TRUE
    ),
    error = function(e) e
  ))
cat("done with xgb")

  fit_status <- data.table(
    candidate = candidate_row$candidate,
    rep_id = as.integer(rep_id),
    model = c("COX", "XGB"),
    ok = c(!inherits(fit_cox, "error"), !inherits(fit_xgb, "error")),
    error = c(
      if (inherits(fit_cox, "error")) conditionMessage(fit_cox) else NA_character_,
      if (inherits(fit_xgb, "error")) conditionMessage(fit_xgb) else NA_character_
    )
  )

  models <- list()
  if (!inherits(fit_cox, "error")) models$COX <- fit_cox
  if (!inherits(fit_xgb, "error")) models$XGB <- fit_xgb

  if (length(models) == 0L) {
    return(list(scores = data.table(), fit_status = fit_status))
  }

  score_obj <- ReSurv::Score_Reserving(
    models = models,
    newdata = dat,
    scoring_metrics = c("EI", "R-tot", "R-cell-wise", "R-cal-wise", "CRPS"),
    granularity = "output",
    chain_ladder = TRUE
  )

  scores <- data.table::rbindlist(unclass(score_obj), use.names = TRUE, fill = TRUE)
  scores[, `:=`(
    candidate = candidate_row$candidate,
    rep_id = as.integer(rep_id),
    sim_seed = sim_seed,
    fit_seed = fit_seed,
    n_claims = nrow(dat),
    observed_upper = observed_n,
    actual_ibnr = ibnr_n,
    beta_break = candidate_row$beta_break,
    age_mean = candidate_row$age_mean,
    age_sd = candidate_row$age_sd,
    age_min = candidate_row$age_min,
    age_max = candidate_row$age_max,
    age_cut1 = candidate_row$age_cut1,
    age_cut2 = candidate_row$age_cut2,
    age_low = candidate_row$age_low,
    age_mid = candidate_row$age_mid,
    age_high = candidate_row$age_high
  )]

  data.table::setcolorder(scores, c(
    "candidate", "rep_id", "sim_seed", "fit_seed", "model", "metric", "score",
    "n_claims", "observed_upper", "actual_ibnr",
    "beta_break", "age_mean", "age_sd", "age_min", "age_max",
    "age_cut1", "age_cut2", "age_low", "age_mid", "age_high"
  ))

  list(scores = scores[], fit_status = fit_status[])
}

summarize_scores <- function(scores) {
  scores[, .(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = stats::sd(score, na.rm = TRUE),
    median_score = stats::median(score, na.rm = TRUE),
    q025 = stats::quantile(score, 0.025, na.rm = TRUE, names = FALSE),
    q975 = stats::quantile(score, 0.975, na.rm = TRUE, names = FALSE),
    n = sum(!is.na(score))
  ), by = .(candidate, model, metric)]
}

make_metric_wide <- function(scores) {
  wide <- data.table::dcast(
    scores,
    candidate + rep_id + model + actual_ibnr + n_claims + observed_upper ~ metric,
    value.var = "score"
  )
  if ("EI" %in% names(wide)) wide[, EIR := EI - 1]
  wide[]
}

make_candidate_diagnostics <- function(wide_scores) {
  means <- wide_scores[, .(
    actual_ibnr_mean = mean(actual_ibnr, na.rm = TRUE),
    actual_ibnr_min = min(actual_ibnr, na.rm = TRUE),
    n_claims_mean = mean(n_claims, na.rm = TRUE),
    observed_upper_mean = mean(observed_upper, na.rm = TRUE),
    EI_mean = mean(EI, na.rm = TRUE),
    EIR_mean = mean(EIR, na.rm = TRUE),
    Rtot_mean = mean(`R-tot`, na.rm = TRUE),
    Rcell_mean = mean(`R-cell-wise`, na.rm = TRUE),
    Rcal_mean = mean(`R-cal-wise`, na.rm = TRUE),
    CRPS_mean = mean(CRPS, na.rm = TRUE)
  ), by = .(candidate, model)]

  d <- data.table::dcast(
    means,
    candidate + actual_ibnr_mean + actual_ibnr_min + n_claims_mean + observed_upper_mean ~ model,
    value.var = c("EIR_mean", "Rtot_mean", "Rcell_mean", "Rcal_mean", "CRPS_mean")
  )

  ## Selection flags. Lower is better for R-cell-wise, R-cal-wise, and CRPS.
  if (all(c("Rcell_mean_CL", "Rcell_mean_COX", "Rcell_mean_XGB") %in% names(d))) {
    d[, `:=`(
      cox_beats_cl_cell = Rcell_mean_COX < Rcell_mean_CL,
      xgb_beats_cl_cell = Rcell_mean_XGB < Rcell_mean_CL,
      xgb_beats_cox_cell = Rcell_mean_XGB < Rcell_mean_COX
    )]
  }

  if (all(c("Rcal_mean_CL", "Rcal_mean_COX", "Rcal_mean_XGB") %in% names(d))) {
    d[, `:=`(
      cox_beats_cl_cal = Rcal_mean_COX < Rcal_mean_CL,
      xgb_beats_cl_cal = Rcal_mean_XGB < Rcal_mean_CL,
      xgb_beats_cox_cal = Rcal_mean_XGB < Rcal_mean_COX
    )]
  }

  if (all(c("CRPS_mean_COX", "CRPS_mean_XGB") %in% names(d))) {
    d[, xgb_beats_cox_crps := CRPS_mean_XGB < CRPS_mean_COX]
  }

  d[]
}

## -------------------------------------------------------------------------
## Main loop.
## -------------------------------------------------------------------------

dt_threads <- get_env_threads("GAMMA_CALIB_DT_THREADS")
xgb_nthreads <- get_env_threads("GAMMA_CALIB_XGB_NTHREADS", default = dt_threads)
configure_threads(dt_threads = dt_threads, xgb_threads = xgb_nthreads)

reps <- seq_len(get_env_int("GAMMA_CALIB_REPS", 3L))
years <- get_env_int("GAMMA_CALIB_YEARS", 4L)
time_unit <- 1 / 360
period_exposure <- get_env_int("GAMMA_CALIB_PERIOD_EXPOSURE", 150L)
period_frequency <- get_env_num("GAMMA_CALIB_PERIOD_FREQUENCY", 0.05)
xgb_nrounds <- get_env_int("GAMMA_CALIB_XGB_NROUNDS", 500L)
out_dir <- Sys.getenv("GAMMA_CALIB_OUT_DIR", unset = "results/gamma_calibration_simple")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidates <- make_candidates()

candidates<-candidates[candidate =="G4_compact_age",]

all_scores <- list()
all_fit_status <- list()
kk <- 0L

for (ii in seq_len(nrow(candidates))) {
  for (rr in reps) {

    kk <- kk + 1L
    res <- tryCatch(
      fit_and_score_one(
        candidate_row = candidates[ii],
        rep_id = rr,
        years = years,
        time_unit = time_unit,
        period_exposure = period_exposure,
        period_frequency = period_frequency,
        xgb_nrounds = xgb_nrounds,
        xgb_nthreads = xgb_nthreads
      ),
      error = function(e) {
        list(
          scores = data.table(),
          fit_status = data.table(
            candidate = candidates$candidate[ii],
            rep_id = as.integer(rr),
            model = "whole_replication",
            ok = FALSE,
            error = conditionMessage(e)
          )
        )
      }
    )
    all_scores[[kk]] <- res$scores
    all_fit_status[[kk]] <- res$fit_status
    gc(verbose = FALSE)


    }
}

scores <- data.table::rbindlist(all_scores, use.names = TRUE, fill = TRUE)
fit_status <- data.table::rbindlist(all_fit_status, use.names = TRUE, fill = TRUE)

if (nrow(scores) == 0L) {
  data.table::fwrite(fit_status, file.path(out_dir, "gamma_calibration_fit_status.csv"))
  stop("No scores were produced. See gamma_calibration_fit_status.csv.")
}

score_summary <- summarize_scores(scores)
metric_wide <- make_metric_wide(scores)
diagnostics <- make_candidate_diagnostics(metric_wide)

## Save outputs.
data.table::fwrite(candidates, file.path(out_dir, "gamma_calibration_candidates.csv"))
data.table::fwrite(scores, file.path(out_dir, "gamma_calibration_scores.csv"))
data.table::fwrite(score_summary, file.path(out_dir, "gamma_calibration_score_summary.csv"))
data.table::fwrite(metric_wide, file.path(out_dir, "gamma_calibration_metric_wide.csv"))
data.table::fwrite(diagnostics, file.path(out_dir, "gamma_calibration_candidate_diagnostics.csv"))
data.table::fwrite(fit_status, file.path(out_dir, "gamma_calibration_fit_status.csv"))

cat("\nDone. Wrote files to:", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n\n")
cat("Candidate diagnostics:\n")
print(diagnostics)

