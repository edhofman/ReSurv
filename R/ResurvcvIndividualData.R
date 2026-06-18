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

  if (!is.character(model) || length(model) != 1L ||
      !(model %in% c("NN", "XGB"))) {
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

  grid <- hparameters_grid

  ## ------------------------------------------------------------------
  ## Hyperparameter grid normalization
  ## ------------------------------------------------------------------

  if (model == "NN") {

    ## Backward-compatible alias.
    if ("epsilon" %in% names(grid) && !("eps" %in% names(grid))) {
      grid$eps <- grid$epsilon
      grid$epsilon <- NULL
    }

    if ("epsilon" %in% names(grid) && "eps" %in% names(grid)) {
      stop(
        "Use only one of `eps` or `epsilon` in `hparameters_grid`. ",
        "The current torch backend uses `eps`.",
        call. = FALSE
      )
    }

    forbidden_nn <- intersect(
      names(grid),
      c("epochs", "verbose", "num_workers")
    )

    if (length(forbidden_nn) > 0L) {
      stop(
        "The following NN entries must be supplied as ReSurvCV arguments, ",
        "not inside `hparameters_grid`: ",
        paste(forbidden_nn, collapse = ", "),
        call. = FALSE
      )
    }

    required_nn <- c(
      "num_layers",
      "num_nodes",
      "activation",
      "optim",
      "lr",
      "xi",
      "eps"
    )

    missing_nn <- setdiff(required_nn, names(grid))

    if (length(missing_nn) > 0L) {
      stop(
        "Missing NN hyperparameters in `hparameters_grid`: ",
        paste(missing_nn, collapse = ", "),
        call. = FALSE
      )
    }

    hparameters.f <- expand.grid(
      grid,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    hparameters.f$num_layers <- as.integer(hparameters.f$num_layers)
    hparameters.f$num_nodes <- as.integer(hparameters.f$num_nodes)

    if (any(!is.finite(hparameters.f$num_layers)) ||
        any(hparameters.f$num_layers < 1L)) {
      stop("`num_layers` must contain positive integers.", call. = FALSE)
    }

    if (any(!is.finite(hparameters.f$num_nodes)) ||
        any(hparameters.f$num_nodes < 1L)) {
      stop("`num_nodes` must contain positive integers.", call. = FALSE)
    }

    ## Old encoding compatibility:
    ## activation = 1 -> LeakyReLU
    ## activation = 2 -> SELU
    activation_raw <- tolower(trimws(as.character(hparameters.f$activation)))

    activation_raw[activation_raw == "1"] <- "leakyrelu"
    activation_raw[activation_raw == "2"] <- "selu"
    activation_raw[activation_raw == "leaky-relu"] <- "leakyrelu"
    activation_raw[activation_raw == "leaky_relu"] <- "leakyrelu"
    activation_raw[activation_raw == "leakyrelu"] <- "leakyrelu"
    activation_raw[activation_raw == "relu"] <- "relu"
    activation_raw[activation_raw == "selu"] <- "selu"
    activation_raw[activation_raw == "tanh"] <- "tanh"

    allowed_activation <- c(
      "relu",
      "selu",
      "tanh",
      "leakyrelu"
    )

    bad_activation <- unique(
      activation_raw[!(activation_raw %in% allowed_activation)]
    )

    if (length(bad_activation) > 0L) {
      stop(
        "Unknown NN activation value(s): ",
        paste(bad_activation, collapse = ", "),
        ". Allowed values are relu, selu, tanh, leakyrelu, or old codes 1/2.",
        call. = FALSE
      )
    }

    hparameters.f$activation <- activation_raw

    ## Old encoding compatibility:
    ## optim = 1 -> Adam
    ## optim = 2 -> SGD
    optim_raw <- tolower(trimws(as.character(hparameters.f$optim)))

    optim_out <- optim_raw
    optim_out[optim_raw == "1"] <- "Adam"
    optim_out[optim_raw == "2"] <- "SGD"
    optim_out[optim_raw == "adam"] <- "Adam"
    optim_out[optim_raw == "sgd"] <- "SGD"
    optim_out[optim_raw == "adamw"] <- "AdamW"

    bad_optim <- unique(
      optim_raw[!(optim_raw %in% c("1", "2", "adam", "sgd", "adamw"))]
    )

    if (length(bad_optim) > 0L) {
      stop(
        "Unknown NN optimizer value(s): ",
        paste(bad_optim, collapse = ", "),
        ". Allowed values are Adam, SGD, AdamW, or old codes 1/2.",
        call. = FALSE
      )
    }

    hparameters.f$optim <- optim_out

    hparameters.f$lr <- as.numeric(hparameters.f$lr)
    hparameters.f$xi <- as.numeric(hparameters.f$xi)
    hparameters.f$eps <- as.numeric(hparameters.f$eps)

    if (any(!is.finite(hparameters.f$lr)) ||
        any(hparameters.f$lr <= 0)) {
      stop("`lr` must contain positive finite values.", call. = FALSE)
    }

    if (any(!is.finite(hparameters.f$xi)) ||
        any(hparameters.f$xi < 0) ||
        any(hparameters.f$xi > 1)) {
      stop("`xi` must lie in [0, 1].", call. = FALSE)
    }

    if (any(!is.finite(hparameters.f$eps)) ||
        any(hparameters.f$eps < 0)) {
      stop("`eps` must contain non-negative finite values.", call. = FALSE)
    }

    if (!("early_stopping" %in% names(hparameters.f))) {
      hparameters.f$early_stopping <- FALSE
    } else {
      hparameters.f$early_stopping <- as.logical(hparameters.f$early_stopping)
      hparameters.f$early_stopping[is.na(hparameters.f$early_stopping)] <- FALSE
    }

    if (!("patience" %in% names(hparameters.f))) {
      hparameters.f$patience <- max(1L, as.integer(epochs))
    } else {
      hparameters.f$patience <- as.integer(hparameters.f$patience)

      if (any(!is.finite(hparameters.f$patience)) ||
          any(hparameters.f$patience < 1L)) {
        stop("`patience` must contain positive integers.", call. = FALSE)
      }
    }

    ## Explicit node columns required by build_deepsurv_net().
    max_layers <- max(hparameters.f$num_layers, na.rm = TRUE)

    for (jj in seq_len(max_layers)) {
      hparameters.f[[paste0("node_", jj)]] <- NA_integer_
    }

    for (rr in seq_len(nrow(hparameters.f))) {
      nl <- hparameters.f$num_layers[rr]
      nn <- hparameters.f$num_nodes[rr]

      for (jj in seq_len(nl)) {
        hparameters.f[[paste0("node_", jj)]][rr] <- nn
      }
    }
  }

  if (model == "XGB") {

    forbidden_xgb <- intersect(
      names(grid),
      c("nrounds", "early_stopping_rounds", "print_every_n", "verbose")
    )

    if (length(forbidden_xgb) > 0L) {
      stop(
        "The following XGB entries must be supplied as ReSurvCV arguments, ",
        "not inside `hparameters_grid`: ",
        paste(forbidden_xgb, collapse = ", "),
        call. = FALSE
      )
    }

    hparameters.f <- expand.grid(
      grid,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    if (!("booster" %in% names(hparameters.f))) {
      hparameters.f$booster <- "gbtree"
    }

    if ("max_depth" %in% names(hparameters.f)) {
      hparameters.f$max_depth <- as.integer(hparameters.f$max_depth)
    }

    if ("min_child_weight" %in% names(hparameters.f)) {
      hparameters.f$min_child_weight <- as.numeric(hparameters.f$min_child_weight)
    }

    if ("eta" %in% names(hparameters.f)) {
      hparameters.f$eta <- as.numeric(hparameters.f$eta)
    }

    if ("subsample" %in% names(hparameters.f)) {
      hparameters.f$subsample <- as.numeric(hparameters.f$subsample)
    }

    if ("alpha" %in% names(hparameters.f)) {
      hparameters.f$alpha <- as.numeric(hparameters.f$alpha)
    }

    if ("lambda" %in% names(hparameters.f)) {
      hparameters.f$lambda <- as.numeric(hparameters.f$lambda)
    }
  }

  out <- cbind(
    hparameters.f,
    train.lkh = NA_real_,
    test.lkh  = NA_real_,
    time      = NA_real_
  )

  ## ------------------------------------------------------------------
  ## Run CV
  ## ------------------------------------------------------------------

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

  if (!("test.lkh" %in% names(out.cv)) ||
      !any(is.finite(out.cv$test.lkh))) {
    stop(
      "No finite validation likelihood was produced by cross-validation.",
      call. = FALSE
    )
  }

  best_id <- which.min(out.cv$test.lkh)

  hp_cols <- names(hparameters.f)

  best_row <- out.cv[
    best_id,
    ,
    drop = FALSE
  ]

  best_hp_table <- best_row[
    ,
    hp_cols,
    drop = FALSE
  ]

  if (model == "XGB") {

    best_params <- as.list(best_hp_table[1, , drop = TRUE])

    best_hparameters <- list(
      params = best_params,
      print_every_n = print_every_n,
      nrounds = nrounds,
      verbose = verbose,
      early_stopping_rounds = early_stopping_rounds
    )

  } else {

    best_hparameters <- as.data.frame(
      best_hp_table,
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    out.cv = as.data.frame(out.cv),
    out.cv.best.oos = as.data.frame(best_row),
    hparameters.best = best_hparameters
  )

  class(out) <- "ReSurvCV"

  out
}
