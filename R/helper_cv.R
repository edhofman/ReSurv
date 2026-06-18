pkg.env$cv_design_matrix <- function(IndividualDataPP,
                                     continuous_features_scaling_method = "minmax",
                                     remove_first_dummy = FALSE) {

  training_data <- IndividualDataPP$training.data

  categorical_features <- IndividualDataPP$data_information$categorical_features
  continuous_features  <- IndividualDataPP$data_information$continuous_features

  n <- nrow(training_data)

  is_baseline_model <- is.null(c(categorical_features, continuous_features))

  if (is_baseline_model) {
    X <- data.frame(intercept_1 = rep(1, n))
  } else {
    X_parts <- list()

    if (!is.null(categorical_features)) {
      X_parts[["categorical"]] <- pkg.env$model.matrix.creator(
        data = training_data,
        select_columns = categorical_features,
        remove_first_dummy = remove_first_dummy
      )
    }

    if (!is.null(continuous_features)) {
      scaler <- pkg.env$scaler(
        continuous_features_scaling_method = continuous_features_scaling_method
      )

      X_parts[["continuous"]] <- training_data |>
        dplyr::reframe(
          dplyr::across(
            dplyr::all_of(continuous_features),
            scaler
          )
        )
    }

    X <- do.call(cbind, X_parts)
    X <- as.data.frame(X)
  }

  Y <- training_data[, c("DP_rev_i", "I", "TR_i")]

  list(X = X, Y = Y)
}


pkg.env$xgboost_cv <- function(IndividualDataPP,
                               folds,
                               kfolds,
                               print_every_n = 1L,
                               nrounds = NULL,
                               verbose = FALSE,
                               early_stopping_rounds = NULL,
                               hparameters.f,
                               out,
                               verbose.cv = FALSE,
                               parallel = FALSE,
                               ncores = 1,
                               random_seed,
                               continuous_features_scaling_method = "minmax") {

  if (isTRUE(parallel)) {
    warning(
      "`parallel = TRUE` is not currently supported in the refactored CV path. ",
      "Running sequentially.",
      call. = FALSE
    )
  }

  xy <- pkg.env$cv_design_matrix(
    IndividualDataPP = IndividualDataPP,
    continuous_features_scaling_method = continuous_features_scaling_method,
    remove_first_dummy = TRUE
  )

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

    out[hp, c("train.lkh", "test.lkh", "time")] <- pkg.env$cv_xgboost(
      hp = hp,
      X = xy$X,
      Y = xy$Y,
      folds = folds,
      kfolds = kfolds,
      print_every_n = print_every_n,
      nrounds = nrounds,
      verbose = verbose,
      early_stopping_rounds = early_stopping_rounds,
      hparameters.f = hparameters.f,
      random_seed = random_seed
    )
  }

  out
}

pkg.env$get_xgb_eval_value <- function(model.out,
                                       dataset = c("train", "eval"),
                                       iteration = NULL) {
  dataset <- match.arg(dataset)

  evaluation_log <- model.out$evaluation_log

  if (is.null(evaluation_log) || nrow(evaluation_log) == 0L) {
    stop("XGBoost evaluation log is empty.", call. = FALSE)
  }

  if (is.null(iteration) || length(iteration) != 1L || is.na(iteration)) {
    iteration <- nrow(evaluation_log)
  }

  iteration <- min(iteration, nrow(evaluation_log))

  nm <- names(evaluation_log)

  candidates <- nm[
    grepl(paste0("^", dataset, "_"), nm) &
      grepl("log|partial|likelihood|lkh", nm, ignore.case = TRUE)
  ]

  if (length(candidates) == 0L) {
    candidates <- nm[
      grepl(paste0("^", dataset, "_"), nm) &
        nm != "iter"
    ]
  }

  if (length(candidates) == 0L) {
    stop(
      "Could not find an XGBoost evaluation metric for dataset `",
      dataset,
      "`. Available columns are: ",
      paste(nm, collapse = ", "),
      call. = FALSE
    )
  }

  evaluation_log[[candidates[1L]]][iteration]
}
pkg.env$cv_xgboost <- function(hp,
                               X,
                               Y,
                               folds,
                               kfolds,
                               print_every_n,
                               nrounds,
                               verbose,
                               early_stopping_rounds,
                               hparameters.f,
                               random_seed) {

  start <- Sys.time()

  hparameters <- list(
    params = as.list(hparameters.f[hp, , drop = FALSE]),
    print_every_n = print_every_n,
    nrounds = nrounds,
    verbose = verbose,
    early_stopping_rounds = early_stopping_rounds
  )

  tmp.train.lkh <- numeric(folds)
  tmp.test.lkh  <- numeric(folds)

  for (i in seq_len(folds)) {

    datads_pp <- pkg.env$xgboost_pp(
      X = X,
      Y = Y,
      samples_TF = kfolds != i
    )

    model.out.k <- pkg.env$fit_xgboost(
      datads_pp = datads_pp,
      hparameters = hparameters
    )

    best.it <- model.out.k$best_iteration

    if (is.null(best.it) || length(best.it) == 0L || is.na(best.it)) {
      best.it <- NULL
    }

    ## Score training fold directly from predictions.
    if (!is.null(best.it)) {
      pred_train <- tryCatch(
        predict(
          model.out.k,
          datads_pp$ds_train_m,
          iterationrange = c(1L, as.integer(best.it))
        ),
        error = function(e) {
          tryCatch(
            predict(
              model.out.k,
              datads_pp$ds_train_m,
              ntreelimit = as.integer(best.it)
            ),
            error = function(e2) {
              predict(model.out.k, datads_pp$ds_train_m)
            }
          )
        }
      )
    } else {
      pred_train <- predict(model.out.k, datads_pp$ds_train_m)
    }

    metric_train <- cox_evaluation_metrics(
      pred_train,
      datads_pp$ds_train_m
    )

    if (is.list(metric_train) && "value" %in% names(metric_train)) {
      tmp.train.lkh[i] <- as.numeric(metric_train$value)
    } else if (is.list(metric_train) && length(metric_train) >= 2L) {
      tmp.train.lkh[i] <- as.numeric(metric_train[[2L]])
    } else if (is.numeric(metric_train) && length(metric_train) == 1L) {
      tmp.train.lkh[i] <- as.numeric(metric_train)
    } else {
      stop(
        "Could not extract training CV metric from `cox_evaluation_metrics()`.",
        call. = FALSE
      )
    }

    ## Score validation fold directly from predictions.
    if (is.null(datads_pp$ds_test_m)) {
      tmp.test.lkh[i] <- NA_real_
    } else {

      if (!is.null(best.it)) {
        pred_test <- tryCatch(
          predict(
            model.out.k,
            datads_pp$ds_test_m,
            iterationrange = c(1L, as.integer(best.it))
          ),
          error = function(e) {
            tryCatch(
              predict(
                model.out.k,
                datads_pp$ds_test_m,
                ntreelimit = as.integer(best.it)
              ),
              error = function(e2) {
                predict(model.out.k, datads_pp$ds_test_m)
              }
            )
          }
        )
      } else {
        pred_test <- predict(model.out.k, datads_pp$ds_test_m)
      }

      metric_test <- cox_evaluation_metrics(
        pred_test,
        datads_pp$ds_test_m
      )

      if (is.list(metric_test) && "value" %in% names(metric_test)) {
        tmp.test.lkh[i] <- as.numeric(metric_test$value)
      } else if (is.list(metric_test) && length(metric_test) >= 2L) {
        tmp.test.lkh[i] <- as.numeric(metric_test[[2L]])
      } else if (is.numeric(metric_test) && length(metric_test) == 1L) {
        tmp.test.lkh[i] <- as.numeric(metric_test)
      } else {
        stop(
          "Could not extract validation CV metric from `cox_evaluation_metrics()`.",
          call. = FALSE
        )
      }
    }
  }

  time <- as.numeric(difftime(Sys.time(), start, units = "mins"))

  c(
    mean(tmp.train.lkh, na.rm = TRUE),
    mean(tmp.test.lkh, na.rm = TRUE),
    time
  )
}
pkg.env$nn_hparameter_nodes_grid <- function(hparameters, cv = FALSE) {

  if (!("num_layers" %in% names(hparameters))) {
    return(hparameters)
  }

  if (isTRUE(cv)) {
    max_layers <- max(hparameters$num_layers)

    node_names <- paste0("node_", seq_len(max_layers))

    hparameters <- hparameters |>
      dplyr::rowwise() |>
      dplyr::mutate(
        new = paste(
          paste(rep(num_nodes, num_layers), collapse = ","),
          paste(rep(NA_integer_, max_layers - num_layers), collapse = ","),
          sep = ","
        )
      ) |>
      tidyr::separate(
        new,
        into = node_names,
        sep = ",",
        convert = TRUE
      ) |>
      dplyr::ungroup()

    hparameters$num_nodes <- NULL

    return(as.data.frame(hparameters))
  }

  if (hparameters$num_layers == length(hparameters$num_nodes)) {
    for (i in seq_len(hparameters$num_layers)) {
      hparameters[[paste0("node_", i)]] <- hparameters$num_nodes[i]
    }
  } else if (length(hparameters$num_nodes) == 1L) {
    for (i in seq_len(hparameters$num_layers)) {
      hparameters[[paste0("node_", i)]] <- hparameters$num_nodes
    }
  } else {
    warning(
      "`num_nodes` was not supplied correctly. ",
      "Using the first value for all layers.",
      call. = FALSE
    )

    for (i in seq_len(hparameters$num_layers)) {
      hparameters[[paste0("node_", i)]] <- hparameters$num_nodes[1]
    }
  }

  hparameters[["num_nodes"]] <- NULL

  hparameters
}


pkg.env$deep_surv_cv <- function(IndividualDataPP,
                                 continuous_features_scaling_method,
                                 folds,
                                 kfolds,
                                 random_seed,
                                 verbose = FALSE,
                                 epochs,
                                 hparameters.f,
                                 out,
                                 parallel = FALSE,
                                 ncores = 1,
                                 verbose.cv = FALSE) {

  if (isTRUE(parallel)) {
    warning(
      "`parallel = TRUE` is not currently supported in the refactored CV path. ",
      "Running sequentially.",
      call. = FALSE
    )
  }

  xy <- pkg.env$cv_design_matrix(
    IndividualDataPP = IndividualDataPP,
    continuous_features_scaling_method = continuous_features_scaling_method,
    remove_first_dummy = FALSE
  )

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

    out[hp, c("train.lkh", "test.lkh", "time")] <- pkg.env$cv_deep_surv(
      hp = hp,
      X = xy$X,
      Y = xy$Y,
      folds = folds,
      kfolds = kfolds,
      random_seed = random_seed,
      verbose = verbose,
      epochs = epochs,
      hparameters.f = hparameters.f
    )
  }

  out
}


pkg.env$cv_deep_surv <- function(hp,
                                 X,
                                 Y,
                                 folds,
                                 kfolds,
                                 random_seed,
                                 verbose,
                                 epochs,
                                 hparameters.f) {

  start <- Sys.time()

  hparameters <- list(
    params = as.list(hparameters.f[hp, , drop = FALSE]),
    verbose = verbose,
    epochs = epochs
  )

  tmp.train.lkh <- numeric(folds)
  tmp.test.lkh  <- numeric(folds)

  for (i in seq_len(folds)) {
    datads_pp <- pkg.env$deep_surv_pp(
      X = X,
      Y = Y,
      samples_TF = kfolds != i
    )

    model.out.k <- pkg.env$fit_deep_surv(
      data = datads_pp,
      params = hparameters$params,
      verbose = hparameters$verbose,
      epochs = hparameters$epochs,
      num_workers = 0,
      seed = random_seed
    )

    best.it <- which.min(model.out.k$log$val_loss)

    if (length(best.it) == 0L || is.na(best.it)) {
      best.it <- which.min(model.out.k$log$train_loss)
    }

    tmp.train.lkh[i] <- model.out.k$log$train_loss[best.it]
    tmp.test.lkh[i]  <- model.out.k$log$val_loss[best.it]
  }

  time <- as.numeric(difftime(Sys.time(), start, units = "mins"))

  c(
    mean(tmp.train.lkh, na.rm = TRUE),
    mean(tmp.test.lkh, na.rm = TRUE),
    time
  )
}
