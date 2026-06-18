# XGBoost helper functions
#
# XGBoost data preprocessing and model fitting.
#
# @import xgboost
## xgboost ----

pkg.env$xgboost_pp <- function(X,
                               Y,
                               samples_TF = NULL,
                               training_test_split = .1) {

  if (!is.numeric(training_test_split) ||
      length(training_test_split) != 1L ||
      !is.finite(training_test_split) ||
      training_test_split <= 0 ||
      training_test_split > 1) {
    stop("`training_test_split` must be a number in (0, 1].", call. = FALSE)
  }

  if (!is.null(samples_TF)) {
    xy <- data.frame(X, Y, samples_TF = samples_TF, check.names = FALSE)
  } else {
    xy <- data.frame(X, Y, check.names = FALSE)
  }

  tmp <- xy[order(xy$DP_rev_i), , drop = FALSE]
  tmp$id <- seq_len(nrow(tmp))

  if (is.null(samples_TF)) {
    if (training_test_split == 1) {
      sampled_id <- tmp$id
    } else {
      n_sample <- ceiling(nrow(tmp) * training_test_split)
      n_sample <- max(1L, min(nrow(tmp), n_sample))
      sampled_id <- sample(tmp$id, size = n_sample, replace = FALSE)
    }

    samples_cn <- data.frame(id = sampled_id)

  } else {
    cond <- as.logical(tmp$samples_TF)

    if (anyNA(cond)) {
      stop("`samples_TF` must be coercible to TRUE/FALSE without NA values.",
           call. = FALSE)
    }

    samples_cn <- data.frame(id = tmp$id[cond])
    tmp$samples_TF <- NULL
  }

  make_efron_c <- function(z) {
    ave(
      seq_along(z),
      z,
      FUN = function(ind) (seq_along(ind) - 1) / length(ind)
    )
  }

  make_dmatrix <- function(df) {
    xgboost::xgb.DMatrix(
      data = as.matrix(df[, colnames(X), drop = FALSE]),
      label = df$I
    )
  }

  add_xgb_attrs <- function(dmat, df) {
    event_times <- unique(df$DP_rev_i)

    attr(dmat, "truncation") <- df$TR_i
    attr(dmat, "claim_arrival") <- df$DP_rev_i

    attr(dmat, "risk_sets") <- risks_in_the_tie(
      starts_i = df$TR_i,
      stops_i  = df$DP_rev_i,
      stops    = event_times
    )

    attr(dmat, "event_sets") <- events_in_the_tie(
      starts_i = df$TR_i,
      stops_i  = df$DP_rev_i,
      stops    = event_times
    )

    attr(dmat, "efron_c") <- df$efron_c

    tie_table <- table(df$DP_rev_i)

    attr(dmat, "tieid") <- unname(tie_table)

    attr(dmat, "groups") <- rep(
      as.integer(names(tie_table)),
      unname(tie_table)
    )

    dmat
  }

  tmp_train <- tmp[tmp$id %in% samples_cn$id, , drop = FALSE]
  tmp_train <- tmp_train[order(tmp_train$DP_rev_i), , drop = FALSE]
  tmp_train$efron_c <- make_efron_c(tmp_train$DP_rev_i)

  ds_train_m <- make_dmatrix(tmp_train)
  ds_train_m <- add_xgb_attrs(ds_train_m, tmp_train)

  if (training_test_split < 1) {
    tmp_test <- tmp[!(tmp$id %in% samples_cn$id), , drop = FALSE]
    tmp_test <- tmp_test[order(tmp_test$DP_rev_i), , drop = FALSE]
    tmp_test$efron_c <- make_efron_c(tmp_test$DP_rev_i)

    ds_test_m <- make_dmatrix(tmp_test)
    ds_test_m <- add_xgb_attrs(ds_test_m, tmp_test)

    return(list(
      ds_train_m = ds_train_m,
      ds_test_m  = ds_test_m,
      samples_cn = samples_cn
    ))
  }

  list(
    ds_train_m = ds_train_m,
    ds_test_m  = NULL,
    samples_cn = samples_cn
  )
}
pkg.env$fit_xgboost <- function(datads_pp,
                                hparameters = list()) {

  if (length(hparameters) == 0L) {
    hparameters <- list(
      params = list(
        booster = "gbtree",
        eta = .01,
        subsample = .5,
        alpha = 1,
        lambda = 1,
        min_child_weight = .2
      ),
      print_every_n = NULL,
      nrounds = 10,
      verbose = FALSE,
      early_stopping_rounds = 500
    )
  }

  evals <- list(train = datads_pp$ds_train_m)

  if (!is.null(datads_pp$ds_test_m)) {
    evals$eval <- datads_pp$ds_test_m
  }

  early_stopping_rounds <- hparameters$early_stopping_rounds

  if (is.null(datads_pp$ds_test_m)) {
    early_stopping_rounds <- NULL
  }

  out <- xgboost::xgb.train(
    params = hparameters$params,
    data = datads_pp$ds_train_m,
    obj = cox_loss_objective,
    nrounds = hparameters$nrounds,
    custom_metric = cox_evaluation_metrics,
    evals = evals,
    verbose = hparameters$verbose,
    print_every_n = hparameters$print_every_n,
    early_stopping_rounds = early_stopping_rounds,
    maximize = FALSE
  )

  out
}



# Cross-validation

