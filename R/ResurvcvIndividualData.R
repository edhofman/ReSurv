#' K-fold cross-validation of a ReSurv model
#'
#' @param IndividualDataPP An object of class \code{IndividualDataPP}.
#' @param model Character. Either \code{"NN"} or \code{"XGB"}.
#' @param hparameters_grid Named list defining the hyperparameter grid.
#' @param folds Integer. Number of folds.
#' @param random_seed Integer. Random seed.
#' @param continuous_features_scaling_method Character. Scaling method for continuous features.
#' @param print_every_n Integer. Passed to XGBoost.
#' @param nrounds Integer. Number of XGBoost boosting rounds.
#' @param early_stopping_rounds Integer. XGBoost early stopping.
#' @param epochs Integer. Number of NN epochs.
#' @param parallel Logical. Currently passed to the CV helper.
#' @param ncores Integer. Number of cores if parallel execution is used.
#' @param num_workers Deprecated for the native torch backend. Ignored.
#' @param verbose Logical. Print model fitting output.
#' @param verbose.cv Logical. Print CV progress.
#'
#' @return An object of class \code{ReSurvCV}.
#'
#' @export
ReSurvCV <- function(IndividualDataPP,
                     model,
                     hparameters_grid,
                     folds,
                     random_seed,
                     continuous_features_scaling_method = "minmax",
                     print_every_n = 1L,
                     nrounds = NULL,
                     early_stopping_rounds = NULL,
                     epochs = 1,
                     parallel = FALSE,
                     ncores = 1,
                     num_workers = 0,
                     verbose = FALSE,
                     verbose.cv = FALSE) {

  UseMethod("ReSurvCV")
}

#' @export
ReSurvCV.default <- function(IndividualDataPP,
                             model,
                             hparameters_grid,
                             folds,
                             random_seed,
                             continuous_features_scaling_method = "minmax",
                             print_every_n = 1L,
                             nrounds = NULL,
                             early_stopping_rounds = NULL,
                             epochs = 1,
                             parallel = FALSE,
                             ncores = 1,
                             num_workers = 0,
                             verbose = FALSE,
                             verbose.cv = FALSE) {

  stop("`IndividualDataPP` must be an object of class `IndividualDataPP`.",
       call. = FALSE)
}

#' @export
ReSurvCV.IndividualDataPP <- function(IndividualDataPP,
                                      model,
                                      hparameters_grid,
                                      folds,
                                      random_seed,
                                      continuous_features_scaling_method = "minmax",
                                      print_every_n = 1L,
                                      nrounds = NULL,
                                      early_stopping_rounds = NULL,
                                      epochs = 1,
                                      parallel = FALSE,
                                      ncores = 1,
                                      num_workers = 0,
                                      verbose = FALSE,
                                      verbose.cv = FALSE) {

  if (!is.character(model) || length(model) != 1L || !(model %in% c("NN", "XGB"))) {
    stop("`model` must be either \"NN\" or \"XGB\".", call. = FALSE)
  }

  if (!is.numeric(folds) || length(folds) != 1L ||
      !is.finite(folds) || folds < 2) {
    stop("`folds` must be a single integer greater than or equal to 2.",
         call. = FALSE)
  }

  folds <- as.integer(folds)

  n <- nrow(IndividualDataPP$training.data)

  if (folds > n) {
    stop("`folds` cannot exceed the number of training observations.",
         call. = FALSE)
  }

  if (!is.list(hparameters_grid) || length(hparameters_grid) == 0L) {
    stop("`hparameters_grid` must be a non-empty named list.",
         call. = FALSE)
  }

  if (is.null(names(hparameters_grid)) ||
      any(names(hparameters_grid) == "")) {
    stop("`hparameters_grid` must be a named list.",
         call. = FALSE)
  }

  set.seed(random_seed)

  ## Balanced random fold assignment.
  kfolds <- sample(rep(seq_len(folds), length.out = n))

  hparameters.f <- expand.grid(
    hparameters_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (model == "NN") {
    hparameters.f <- pkg.env$nn_hparameter_nodes_grid(
      hparameters.f,
      cv = TRUE
    )
  }

  out <- cbind(
    hparameters.f,
    train.lkh = NA_real_,
    test.lkh  = NA_real_,
    time      = NA_real_
  )

  if (model == "XGB") {
    out.cv <- pkg.env$xgboost_cv(
      IndividualDataPP = IndividualDataPP,
      folds = folds,
      kfolds = kfolds,
      random_seed = random_seed,
      print_every_n = print_every_n,
      nrounds = nrounds,
      early_stopping_rounds = early_stopping_rounds,
      hparameters.f = hparameters.f,
      out = out,
      parallel = parallel,
      ncores = ncores,
      verbose = verbose,
      verbose.cv = verbose.cv,
      continuous_features_scaling_method = continuous_features_scaling_method
    )
  }

  if (model == "NN") {
    out.cv <- pkg.env$deep_surv_cv(
      IndividualDataPP = IndividualDataPP,
      continuous_features_scaling_method = continuous_features_scaling_method,
      folds = folds,
      kfolds = kfolds,
      random_seed = random_seed,
      hparameters.f = hparameters.f,
      epochs = epochs,
      out = out,
      parallel = parallel,
      ncores = ncores,
      verbose = verbose,
      verbose.cv = verbose.cv
    )
  }

  best_id <- which(out.cv$test.lkh == min(out.cv$test.lkh, na.rm = TRUE))

  out <- list(
    out.cv = as.data.frame(out.cv),
    out.cv.best.oos = as.data.frame(out.cv[best_id, , drop = FALSE])
  )

  class(out) <- "ReSurvCV"

  out
}
