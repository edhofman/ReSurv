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

  stop(
    "`IndividualDataPP` must be an object of class `IndividualDataPP`.",
    call. = FALSE
  )
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

  ## ------------------------------------------------------------------
  ## Validation
  ## ------------------------------------------------------------------

  if (!is.character(model) ||
      length(model) != 1L ||
      !(model %in% c("NN", "XGB"))) {
    stop(
      "`model` must be either \"NN\" or \"XGB\".",
      call. = FALSE
    )
  }

  if (!is.numeric(folds) ||
      length(folds) != 1L ||
      !is.finite(folds) ||
      folds < 2) {
    stop(
      "`folds` must be a single integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  folds <- as.integer(folds)

  if (!is.numeric(epochs) ||
      length(epochs) != 1L ||
      !is.finite(epochs) ||
      epochs < 1) {
    stop(
      "`epochs` must be a positive integer.",
      call. = FALSE
    )
  }

  epochs <- as.integer(epochs)

  n <- nrow(IndividualDataPP$training.data)

  if (folds > n) {
    stop(
      "`folds` cannot exceed the number of training observations.",
      call. = FALSE
    )
  }

  if (!is.list(hparameters_grid) || length(hparameters_grid) == 0L) {
    stop(
      "`hparameters_grid` must be a non-empty named list.",
      call. = FALSE
    )
  }

  if (is.null(names(hparameters_grid)) ||
      any(names(hparameters_grid) == "")) {
    stop(
      "`hparameters_grid` must be a named list.",
      call. = FALSE
    )
  }

  if (isTRUE(parallel)) {
    warning(
      "`parallel = TRUE` is ignored in this single-process CV implementation.",
      call. = FALSE
    )
  }

  set.seed(random_seed)

  kfolds <- sample(
    rep(seq_len(folds), length.out = n)
  )

  grid <- hparameters_grid

  ## ------------------------------------------------------------------
  ## Hyperparameter grid: NN
  ## ------------------------------------------------------------------

  if (model == "NN") {

    if (!requireNamespace("torch", quietly = TRUE)) {
      stop(
        "Package `torch` is required for `model = \"NN\"`.",
        call. = FALSE
      )
    }

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

    if ("nodes" %in% names(grid)) {

      if ("num_layers" %in% names(grid) || "num_nodes" %in% names(grid)) {
        stop(
          "When using `nodes`, do not also provide `num_layers` or `num_nodes`. ",
          "`nodes` already defines the full NN architecture.",
          call. = FALSE
        )
      }

      nodes_list <- grid$nodes

      if (!is.list(nodes_list) || length(nodes_list) == 0L) {
        stop(
          "`nodes` must be a non-empty list, for example ",
          "list(c(4L), c(8L, 4L), c(16L, 8L, 4L)).",
          call. = FALSE
        )
      }

      nodes_list <- lapply(nodes_list, function(z) {
        z <- as.integer(z)

        if (length(z) == 0L ||
            any(is.na(z)) ||
            any(!is.finite(z)) ||
            any(z < 1L)) {
          stop(
            "Each element of `nodes` must be a positive integer vector.",
            call. = FALSE
          )
        }

        z
      })

      grid$nodes <- NULL

      required_nn <- c(
        "activation",
        "optim",
        "lr",
        "xi",
        "eps"
      )

      missing_nn <- setdiff(
        required_nn,
        names(grid)
      )

      if (length(missing_nn) > 0L) {
        stop(
          "Missing NN hyperparameters in `hparameters_grid`: ",
          paste(missing_nn, collapse = ", "),
          call. = FALSE
        )
      }

      base_grid <- expand.grid(
        grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      max_layers <- max(lengths(nodes_list))

      rows <- vector(
        "list",
        nrow(base_grid) * length(nodes_list)
      )

      rr <- 1L

      for (ii in seq_len(nrow(base_grid))) {

        for (aa in seq_along(nodes_list)) {

          arch <- nodes_list[[aa]]

          tmp_row <- base_grid[
            ii,
            ,
            drop = FALSE
          ]

          tmp_row$num_layers <- as.integer(length(arch))

          tmp_row$num_nodes <- paste(
            arch,
            collapse = "-"
          )

          for (jj in seq_len(max_layers)) {
            tmp_row[[paste0("node_", jj)]] <- if (jj <= length(arch)) {
              as.integer(arch[jj])
            } else {
              NA_integer_
            }
          }

          rows[[rr]] <- tmp_row
          rr <- rr + 1L
        }
      }

      hparameters.f <- do.call(
        rbind,
        rows
      )

      rownames(hparameters.f) <- NULL

    } else {

      required_nn <- c(
        "num_layers",
        "num_nodes",
        "activation",
        "optim",
        "lr",
        "xi",
        "eps"
      )

      missing_nn <- setdiff(
        required_nn,
        names(grid)
      )

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
        stop(
          "`num_layers` must contain positive integers.",
          call. = FALSE
        )
      }

      if (any(!is.finite(hparameters.f$num_nodes)) ||
          any(hparameters.f$num_nodes < 1L)) {
        stop(
          "`num_nodes` must contain positive integers.",
          call. = FALSE
        )
      }

      max_layers <- max(
        hparameters.f$num_layers,
        na.rm = TRUE
      )

      for (jj in seq_len(max_layers)) {
        hparameters.f[[paste0("node_", jj)]] <- NA_integer_
      }

      for (rr in seq_len(nrow(hparameters.f))) {

        nl <- hparameters.f$num_layers[rr]
        nn <- hparameters.f$num_nodes[rr]

        for (jj in seq_len(nl)) {
          hparameters.f[[paste0("node_", jj)]][rr] <- as.integer(nn)
        }
      }
    }

    hparameters.f$num_layers <- as.integer(hparameters.f$num_layers)

    node_cols <- grep(
      "^node_[0-9]+$",
      names(hparameters.f),
      value = TRUE
    )

    if (length(node_cols) == 0L) {
      stop(
        "NN hyperparameter grid must contain at least one `node_*` column.",
        call. = FALSE
      )
    }

    node_cols <- node_cols[
      order(as.integer(sub("^node_", "", node_cols)))
    ]

    for (cc in node_cols) {
      hparameters.f[[cc]] <- as.integer(hparameters.f[[cc]])
    }

    for (rr in seq_len(nrow(hparameters.f))) {

      nl <- hparameters.f$num_layers[rr]

      node_values <- unlist(
        hparameters.f[rr, node_cols, drop = FALSE],
        use.names = FALSE
      )

      if (any(is.na(node_values[seq_len(nl)])) ||
          any(node_values[seq_len(nl)] < 1L)) {
        stop(
          "Each active hidden layer must have a positive number of nodes.",
          call. = FALSE
        )
      }
    }

    activation_raw <- tolower(
      trimws(as.character(hparameters.f$activation))
    )

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

    optim_raw <- tolower(
      trimws(as.character(hparameters.f$optim))
    )

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
      stop(
        "`lr` must contain positive finite values.",
        call. = FALSE
      )
    }

    if (any(!is.finite(hparameters.f$xi)) ||
        any(hparameters.f$xi < 0) ||
        any(hparameters.f$xi > 1)) {
      stop(
        "`xi` must lie in [0, 1].",
        call. = FALSE
      )
    }

    if (any(!is.finite(hparameters.f$eps)) ||
        any(hparameters.f$eps < 0)) {
      stop(
        "`eps` must contain non-negative finite values.",
        call. = FALSE
      )
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
        stop(
          "`patience` must contain positive integers.",
          call. = FALSE
        )
      }
    }
  }

  ## ------------------------------------------------------------------
  ## Hyperparameter grid: XGB
  ## ------------------------------------------------------------------

  if (model == "XGB") {

    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop(
        "Package `xgboost` is required for `model = \"XGB\"`.",
        call. = FALSE
      )
    }

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

    if (is.null(nrounds)) {
      nrounds <- 10L
    }

    nrounds <- as.integer(nrounds)

    if (!is.finite(nrounds) || nrounds < 1L) {
      stop(
        "`nrounds` must be a positive integer.",
        call. = FALSE
      )
    }
  }

  rownames(hparameters.f) <- NULL

  ## ------------------------------------------------------------------
  ## Build design matrix once
  ## ------------------------------------------------------------------

  xy <- cv_method_design_matrix(
    IndividualDataPP = IndividualDataPP,
    model = model,
    continuous_features_scaling_method = continuous_features_scaling_method
  )

  X <- xy$X
  Y <- xy$Y
  ## ------------------------------------------------------------------
  ## Allocate output
  ## ------------------------------------------------------------------

  out.cv <- cbind(
    hparameters.f,
    train.lkh = NA_real_,
    test.lkh = NA_real_,
    time = NA_real_
  )

  rownames(out.cv) <- NULL

  ## ------------------------------------------------------------------
  ## XGB: fast cumulative-array Cox objective/evaluation
  ## ------------------------------------------------------------------

  if (model == "XGB") {

    xgb_folds <- vector("list", folds)

    for (ff in seq_len(folds)) {
      xgb_folds[[ff]] <- cv_make_xgb_fold(
        X = X,
        Y = Y,
        samples_TF = kfolds != ff
      )
    }

    for (hp in seq_len(nrow(hparameters.f))) {

      if (verbose.cv) {
        cat(
          as.character(Sys.time()),
          "Testing XGB hyperparameter combination",
          hp,
          "out of",
          nrow(hparameters.f),
          "\n"
        )
      }

      start_hp <- Sys.time()

      params <- cv_scalar_list(
        hparameters.f[hp, , drop = FALSE]
      )

      tmp.train.lkh <- numeric(folds)
      tmp.test.lkh  <- numeric(folds)

      for (ff in seq_len(folds)) {

        set.seed(as.integer(random_seed + hp * 10000L + ff))

        datads_pp <- xgb_folds[[ff]]

        use_early_stopping <- !is.null(early_stopping_rounds)

        evals <- list()

        early_stopping_rounds_eff <- NULL

        if (use_early_stopping && !is.null(datads_pp$ds_test_m)) {
          evals$eval <- datads_pp$ds_test_m
          early_stopping_rounds_eff <- early_stopping_rounds
        }

        model.out.k <- xgboost::xgb.train(
          params = params,
          data = datads_pp$ds_train_m,
          obj = cv_xgb_fast_loss_objective,
          nrounds = nrounds,
          custom_metric = if (length(evals) > 0L) {
            cv_xgb_fast_evaluation_metrics
          } else {
            NULL
          },
          evals = evals,
          verbose = verbose,
          print_every_n = print_every_n,
          early_stopping_rounds = early_stopping_rounds_eff,
          maximize = FALSE
        )

        best.it <- model.out.k$best_iteration

        if (is.null(best.it) ||
            length(best.it) == 0L ||
            is.na(best.it)) {
          best.it <- nrounds
        }

        pred_train <- cv_predict_xgb_at(
          model.out = model.out.k,
          dmat = datads_pp$ds_train_m,
          iteration = best.it
        )

        tmp.train.lkh[ff] <- cv_extract_metric_value(
          cv_xgb_fast_evaluation_metrics(
            pred_train,
            datads_pp$ds_train_m
          )
        )

        test_from_log <- NA_real_

        if (!is.null(model.out.k$evaluation_log) &&
            nrow(model.out.k$evaluation_log) > 0L) {
          test_from_log <- tryCatch(
            cv_get_xgb_eval_value(
              model.out = model.out.k,
              dataset = "eval",
              iteration = best.it
            ),
            error = function(e) {
              NA_real_
            }
          )
        }

        if (is.finite(test_from_log)) {

          tmp.test.lkh[ff] <- test_from_log

        } else {

          pred_test <- cv_predict_xgb_at(
            model.out = model.out.k,
            dmat = datads_pp$ds_test_m,
            iteration = best.it
          )

          tmp.test.lkh[ff] <- cv_extract_metric_value(
            cv_xgb_fast_evaluation_metrics(
              pred_test,
              datads_pp$ds_test_m
            )
          )
        }
      }

      out.cv[hp, "train.lkh"] <- mean(tmp.train.lkh, na.rm = TRUE)
      out.cv[hp, "test.lkh"]  <- mean(tmp.test.lkh, na.rm = TRUE)
      out.cv[hp, "time"] <- as.numeric(
        difftime(Sys.time(), start_hp, units = "mins")
      )
    }
  }

  ## ------------------------------------------------------------------
  ## NN: precomputed fold tensors and in-script torch fitting
  ## ------------------------------------------------------------------

  if (model == "NN") {

    nn_folds <- vector("list", folds)

    for (ff in seq_len(folds)) {
      nn_folds[[ff]] <- cv_make_nn_fold(
        X = X,
        Y = Y,
        samples_TF = kfolds != ff
      )
    }

    for (hp in seq_len(nrow(hparameters.f))) {

      if (verbose.cv) {
        cat(
          as.character(Sys.time()),
          "Testing NN hyperparameter combination",
          hp,
          "out of",
          nrow(hparameters.f),
          "\n"
        )
      }

      start_hp <- Sys.time()

      params <- cv_scalar_list(
        hparameters.f[hp, , drop = FALSE]
      )

      tmp.train.lkh <- numeric(folds)
      tmp.test.lkh  <- numeric(folds)

      for (ff in seq_len(folds)) {

        log_k <- cv_train_deepsurv_fold(
          fold_data = nn_folds[[ff]],
          params = params,
          epochs = epochs,
          verbose = verbose,
          seed = as.integer(random_seed + hp * 10000L + ff)
        )

        best.it <- which.min(log_k$val_loss)

        if (length(best.it) == 0L || is.na(best.it)) {
          best.it <- which.min(log_k$train_loss)
        }

        tmp.train.lkh[ff] <- log_k$train_loss[best.it]
        tmp.test.lkh[ff]  <- log_k$val_loss[best.it]
      }

      out.cv[hp, "train.lkh"] <- mean(tmp.train.lkh, na.rm = TRUE)
      out.cv[hp, "test.lkh"]  <- mean(tmp.test.lkh, na.rm = TRUE)
      out.cv[hp, "time"] <- as.numeric(
        difftime(Sys.time(), start_hp, units = "mins")
      )
    }
  }

  rownames(out.cv) <- NULL

  ## ------------------------------------------------------------------
  ## Select best row
  ## ------------------------------------------------------------------

  if (!("test.lkh" %in% names(out.cv)) ||
      !any(is.finite(out.cv$test.lkh))) {
    stop(
      "No finite validation likelihood was produced by cross-validation.",
      call. = FALSE
    )
  }

  best_id <- which.min(out.cv$test.lkh)

  out.cv.best <- out.cv[
    best_id,
    ,
    drop = FALSE
  ]

  rownames(out.cv.best) <- NULL

  ## ------------------------------------------------------------------
  ## Final hparameters object
  ## ------------------------------------------------------------------

  if (model == "XGB") {

    best_params <- cv_scalar_list(
      out.cv.best[, names(hparameters.f), drop = FALSE]
    )

    best_hparameters <- list(
      params = best_params,
      print_every_n = print_every_n,
      nrounds = nrounds,
      verbose = verbose,
      early_stopping_rounds = early_stopping_rounds
    )

  } else {

    node_cols <- grep(
      "^node_[0-9]+$",
      names(out.cv.best),
      value = TRUE
    )

    node_cols <- node_cols[
      order(as.integer(sub("^node_", "", node_cols)))
    ]

    node_values <- unlist(
      out.cv.best[
        ,
        node_cols,
        drop = FALSE
      ],
      use.names = FALSE
    )

    node_values <- as.integer(
      node_values[!is.na(node_values)]
    )

    if (length(node_values) == 0L) {
      stop(
        "Could not recover the selected NN architecture from `node_*` columns.",
        call. = FALSE
      )
    }

    best_hparameters <- list(
      num_layers = as.integer(length(node_values)),
      num_nodes = node_values,
      activation = as.character(out.cv.best$activation),
      optim = as.character(out.cv.best$optim),
      lr = as.numeric(out.cv.best$lr),
      xi = as.numeric(out.cv.best$xi),
      eps = as.numeric(out.cv.best$eps),
      early_stopping = if ("early_stopping" %in% names(out.cv.best)) {
        as.logical(out.cv.best$early_stopping)
      } else {
        FALSE
      },
      patience = if ("patience" %in% names(out.cv.best)) {
        as.integer(out.cv.best$patience)
      } else {
        max(1L, as.integer(epochs))
      },
      epochs = as.integer(epochs),
      verbose = isTRUE(verbose),
      num_workers = as.integer(num_workers)
    )
  }

  ## ------------------------------------------------------------------
  ## Return
  ## ------------------------------------------------------------------

  out <- list(
    out.cv = as.data.frame(out.cv),
    out.cv.best.oos = as.data.frame(out.cv.best),
    hparameters.best = best_hparameters
  )

  class(out) <- "ReSurvCV"

  out
}
