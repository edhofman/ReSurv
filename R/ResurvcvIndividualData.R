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
  ## Small in-script utilities
  ## ------------------------------------------------------------------

  scalar_list <- function(x) {
    x <- as.list(x[1L, , drop = FALSE])
    lapply(x, function(z) {
      if (length(z) == 1L) {
        z[[1L]]
      } else {
        z
      }
    })
  }

  extract_metric_value <- function(metric_object) {
    if (is.list(metric_object) && "value" %in% names(metric_object)) {
      return(as.numeric(metric_object$value))
    }

    if (is.list(metric_object) && length(metric_object) >= 2L) {
      return(as.numeric(metric_object[[2L]]))
    }

    if (is.numeric(metric_object) && length(metric_object) == 1L) {
      return(as.numeric(metric_object))
    }

    stop("Could not extract metric value.", call. = FALSE)
  }

  make_efron_c <- function(z) {
    ave(
      seq_along(z),
      z,
      FUN = function(ind) {
        (seq_along(ind) - 1) / length(ind)
      }
    )
  }

  get_xgb_eval_value <- function(model.out,
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

    iteration <- min(as.integer(iteration), nrow(evaluation_log))

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

    as.numeric(evaluation_log[[candidates[1L]]][iteration])
  }

  predict_xgb_at <- function(model.out,
                             dmat,
                             iteration = NULL) {

    if (is.null(iteration) || length(iteration) != 1L || is.na(iteration)) {
      return(predict(model.out, dmat))
    }

    iteration <- as.integer(iteration)

    tryCatch(
      predict(
        model.out,
        dmat,
        iterationrange = c(1L, iteration)
      ),
      error = function(e) {
        tryCatch(
          predict(
            model.out,
            dmat,
            ntreelimit = iteration
          ),
          error = function(e2) {
            predict(model.out, dmat)
          }
        )
      }
    )
  }

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

  training_data <- as.data.frame(IndividualDataPP$training.data)

  categorical_features <- IndividualDataPP$data_information$categorical_features
  continuous_features  <- IndividualDataPP$data_information$continuous_features

  if (length(categorical_features) == 0L) {
    categorical_features <- NULL
  }

  if (length(continuous_features) == 0L) {
    continuous_features <- NULL
  }

  remove_first_dummy <- identical(model, "XGB")

  is_baseline_model <- length(c(categorical_features, continuous_features)) == 0L

  if (is_baseline_model) {

    X <- data.frame(
      intercept_1 = rep(1, nrow(training_data))
    )

  } else {

    X_parts <- list()

    if (!is.null(categorical_features)) {

      if (!requireNamespace("fastDummies", quietly = TRUE)) {
        stop(
          "Package `fastDummies` is required for categorical features.",
          call. = FALSE
        )
      }

      X_tmp <- fastDummies::dummy_cols(
        training_data,
        select_columns = categorical_features,
        remove_selected_columns = TRUE,
        remove_first_dummy = remove_first_dummy
      )

      dummy_keep <- Reduce(
        `|`,
        lapply(
          categorical_features,
          function(z) {
            grepl(z, colnames(X_tmp))
          }
        )
      )

      X_cat <- X_tmp[
        ,
        colnames(X_tmp)[dummy_keep],
        drop = FALSE
      ]

      X_parts[["categorical"]] <- as.data.frame(X_cat)
    }

    if (!is.null(continuous_features)) {

      scale_one <- function(x) {
        x <- as.numeric(x)

        if (continuous_features_scaling_method == "minmax") {
          x_min <- min(x, na.rm = TRUE)
          x_max <- max(x, na.rm = TRUE)

          if (!is.finite(x_min) ||
              !is.finite(x_max) ||
              x_max == x_min) {
            return(rep(0, length(x)))
          }

          return(2 * (x - x_min) / (x_max - x_min) - 1)
        }

        if (continuous_features_scaling_method == "standard") {
          x_mean <- mean(x, na.rm = TRUE)
          x_sd   <- stats::sd(x, na.rm = TRUE)

          if (!is.finite(x_mean) ||
              !is.finite(x_sd) ||
              x_sd == 0) {
            return(rep(0, length(x)))
          }

          return((x - x_mean) / x_sd)
        }

        stop(
          "`continuous_features_scaling_method` must be either ",
          "\"minmax\" or \"standard\".",
          call. = FALSE
        )
      }

      X_cont <- as.data.frame(
        lapply(
          training_data[, continuous_features, drop = FALSE],
          scale_one
        )
      )

      names(X_cont) <- continuous_features

      X_parts[["continuous"]] <- X_cont
    }

    X <- as.data.frame(
      do.call(cbind, X_parts)
    )
  }

  Y <- training_data[
    ,
    c("DP_rev_i", "I", "TR_i"),
    drop = FALSE
  ]

  Y$DP_rev_i <- as.integer(Y$DP_rev_i)
  Y$I        <- as.integer(Y$I)
  Y$TR_i     <- as.integer(Y$TR_i)

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

    weighted_tabulate <- function(idx, weights, nbins) {

      idx <- as.integer(idx)
      weights <- as.numeric(weights)
      nbins <- as.integer(nbins)

      out <- numeric(nbins)

      keep <- !is.na(idx) &
        is.finite(idx) &
        idx >= 1L &
        idx <= nbins &
        !is.na(weights) &
        is.finite(weights)

      if (!any(keep)) {
        return(out)
      }

      tmp <- rowsum(
        x = matrix(weights[keep], ncol = 1L),
        group = idx[keep],
        reorder = FALSE
      )

      out[as.integer(rownames(tmp))] <- as.numeric(tmp[, 1L])

      out
    }

    xgb_fast_quantities <- function(preds, dtrain) {

      Ti <- as.integer(attr(dtrain, "truncation"))
      Ei <- as.integer(attr(dtrain, "claim_arrival"))
      efron_c <- as.numeric(attr(dtrain, "efron_c"))

      exp_p <- exp(preds)

      max_time <- max(c(Ei, Ti + 1L), na.rm = TRUE)

      start_idx <- Ti + 1L
      end_idx   <- Ei + 1L


      add <- weighted_tabulate(
        idx = start_idx,
        weights = exp_p,
        nbins = max_time + 1L
      )

      remove <- weighted_tabulate(
        idx = end_idx,
        weights = exp_p,
        nbins = max_time + 1L
      )
      risk_sum <- cumsum(add - remove)[seq_len(max_time)]

      event_sum <- weighted_tabulate(
        idx = Ei,
        weights = exp_p,
        nbins = max_time
      )

      denom <- risk_sum[Ei] - efron_c * event_sum[Ei]

      denom[!is.finite(denom) | denom <= 0] <- .Machine$double.eps

      list(
        Ti = Ti,
        Ei = Ei,
        efron_c = efron_c,
        exp_p = exp_p,
        denom = denom,
        max_time = max_time
      )
    }

    xgb_fast_evaluation_metrics <- function(preds, dtrain) {

      q <- xgb_fast_quantities(preds, dtrain)

      value <- (
        sum(log(q$denom)) -
          sum(preds)
      ) / length(preds)

      list(
        metric = "log_partial_likelihood",
        value = value
      )
    }

    xgb_fast_loss_objective <- function(preds, dtrain) {

      q <- xgb_fast_quantities(preds, dtrain)

      inv_denom <- 1 / q$denom

      alpha_i <- weighted_tabulate(
        idx = q$Ei,
        weights = inv_denom,
        nbins = q$max_time
      )

      beta_i <- weighted_tabulate(
        idx = q$Ei,
        weights = q$efron_c * inv_denom,
        nbins = q$max_time
      )

      gamma_i <- weighted_tabulate(
        idx = q$Ei,
        weights = inv_denom^2,
        nbins = q$max_time
      )

      omega_i <- weighted_tabulate(
        idx = q$Ei,
        weights = (1 - (1 - q$efron_c)^2) * inv_denom^2,
        nbins = q$max_time
      )

      cumsum_alpha <- c(0, cumsum(alpha_i))
      cumsum_gamma <- c(0, cumsum(gamma_i))

      alpha_i_lt <- cumsum_alpha[q$Ei + 1L] - cumsum_alpha[q$Ti + 1L]
      gamma_i_lt <- cumsum_gamma[q$Ei + 1L] - cumsum_gamma[q$Ti + 1L]

      beta_i_lt  <- beta_i[q$Ei]
      omega_i_lt <- omega_i[q$Ei]

      grad <- q$exp_p * (alpha_i_lt - beta_i_lt) - 1

      hess <- grad -
        (q$exp_p^2) * (gamma_i_lt - omega_i_lt) +
        1

      grad[!is.finite(grad)] <- 0
      hess[!is.finite(hess)] <- 1e-6

      list(
        grad = grad,
        hess = hess
      )
    }

    make_xgb_fold <- function(X, Y, samples_TF) {

      xy <- data.frame(
        X,
        Y,
        samples_TF = samples_TF,
        check.names = FALSE
      )

      tmp <- xy[
        order(xy$DP_rev_i),
        ,
        drop = FALSE
      ]

      tmp$id <- seq_len(nrow(tmp))

      cond <- as.logical(tmp$samples_TF)

      if (anyNA(cond)) {
        stop(
          "`samples_TF` must be coercible to TRUE/FALSE without NA values.",
          call. = FALSE
        )
      }

      tmp$samples_TF <- NULL

      x_cols <- colnames(X)

      make_dmatrix <- function(df) {
        xgboost::xgb.DMatrix(
          data = as.matrix(df[, x_cols, drop = FALSE]),
          label = df$I
        )
      }

      add_attrs <- function(dmat, df) {

        attr(dmat, "truncation") <- as.integer(df$TR_i)
        attr(dmat, "claim_arrival") <- as.integer(df$DP_rev_i)
        attr(dmat, "efron_c") <- as.numeric(make_efron_c(df$DP_rev_i))

        dmat
      }

      tmp_train <- tmp[cond, , drop = FALSE]
      tmp_test  <- tmp[!cond, , drop = FALSE]

      ds_train_m <- make_dmatrix(tmp_train)
      ds_train_m <- add_attrs(ds_train_m, tmp_train)

      ds_test_m <- make_dmatrix(tmp_test)
      ds_test_m <- add_attrs(ds_test_m, tmp_test)

      list(
        ds_train_m = ds_train_m,
        ds_test_m = ds_test_m
      )
    }

    xgb_folds <- vector("list", folds)

    for (ff in seq_len(folds)) {
      xgb_folds[[ff]] <- make_xgb_fold(
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

      params <- scalar_list(
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
          obj = xgb_fast_loss_objective,
          nrounds = nrounds,
          custom_metric = if (length(evals) > 0L) {
            xgb_fast_evaluation_metrics
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

        pred_train <- predict_xgb_at(
          model.out = model.out.k,
          dmat = datads_pp$ds_train_m,
          iteration = best.it
        )

        tmp.train.lkh[ff] <- extract_metric_value(
          xgb_fast_evaluation_metrics(
            pred_train,
            datads_pp$ds_train_m
          )
        )

        test_from_log <- NA_real_

        if (!is.null(model.out.k$evaluation_log) &&
            nrow(model.out.k$evaluation_log) > 0L) {
          test_from_log <- tryCatch(
            get_xgb_eval_value(
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

          pred_test <- predict_xgb_at(
            model.out = model.out.k,
            dmat = datads_pp$ds_test_m,
            iteration = best.it
          )

          tmp.test.lkh[ff] <- extract_metric_value(
            xgb_fast_evaluation_metrics(
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

    make_event_info <- function(y) {

      durations <- as.numeric(y[, 1])
      events    <- as.numeric(y[, 2])

      event_times <- sort(
        unique(durations[events == 1])
      )

      if (length(event_times) == 0L) {
        stop(
          "No events found in NN fold.",
          call. = FALSE
        )
      }

      event_counts <- as.integer(
        tabulate(
          match(durations[events == 1], event_times),
          nbins = length(event_times)
        )
      )

      list(
        event_times = event_times,
        event_counts = event_counts,
        n_events = sum(events == 1)
      )
    }

    make_nn_fold <- function(X, Y, samples_TF) {

      tmp_order <- order(Y$DP_rev_i)

      X_ordered <- as.data.frame(
        X[tmp_order, , drop = FALSE]
      )

      Y_ordered <- as.data.frame(
        Y[tmp_order, , drop = FALSE]
      )

      id_train <- as.logical(samples_TF[tmp_order])

      if (anyNA(id_train)) {
        stop(
          "`samples_TF` must be coercible to TRUE/FALSE without NA values.",
          call. = FALSE
        )
      }

      x_train <- as.matrix(X_ordered[id_train, , drop = FALSE])
      x_val   <- as.matrix(X_ordered[!id_train, , drop = FALSE])

      y_train <- as.matrix(Y_ordered[id_train, , drop = FALSE])
      y_val   <- as.matrix(Y_ordered[!id_train, , drop = FALSE])

      train_info <- make_event_info(y_train)
      val_info   <- make_event_info(y_val)

      list(
        x_train_t = torch::torch_tensor(
          x_train,
          dtype = torch::torch_float32()
        ),
        x_val_t = torch::torch_tensor(
          x_val,
          dtype = torch::torch_float32()
        ),
        dur_train_t = torch::torch_tensor(
          y_train[, 1],
          dtype = torch::torch_float32()
        ),
        event_train_t = torch::torch_tensor(
          y_train[, 2],
          dtype = torch::torch_float32()
        ),
        trunc_train_t = torch::torch_tensor(
          y_train[, 3],
          dtype = torch::torch_float32()
        ),
        dur_val_t = torch::torch_tensor(
          y_val[, 1],
          dtype = torch::torch_float32()
        ),
        event_val_t = torch::torch_tensor(
          y_val[, 2],
          dtype = torch::torch_float32()
        ),
        trunc_val_t = torch::torch_tensor(
          y_val[, 3],
          dtype = torch::torch_float32()
        ),
        train_event_times = train_info$event_times,
        train_event_counts = train_info$event_counts,
        train_n_events = train_info$n_events,
        val_event_times = val_info$event_times,
        val_event_counts = val_info$event_counts,
        val_n_events = val_info$n_events
      )
    }

    build_deepsurv_net <- function(input_dim, params) {

      activation_map <- list(
        "relu"       = torch::nn_relu,
        "selu"       = torch::nn_selu,
        "tanh"       = torch::nn_tanh,
        "leakyrelu"  = torch::nn_leaky_relu,
        "leaky_relu" = torch::nn_leaky_relu
      )

      activation <- tolower(
        as.character(params$activation)
      )

      act_fn <- activation_map[[activation]]

      if (is.null(act_fn)) {
        stop(
          "Unknown activation: ",
          params$activation,
          call. = FALSE
        )
      }

      layers <- list()
      in_features <- input_dim

      for (ll in seq_len(as.integer(params$num_layers))) {

        node_l <- as.integer(
          params[[paste0("node_", ll)]]
        )

        layers <- c(
          layers,
          list(
            torch::nn_linear(in_features, node_l),
            act_fn()
          )
        )

        in_features <- node_l
      }

      layers <- c(
        layers,
        list(
          torch::nn_linear(
            in_features,
            1L,
            bias = FALSE
          )
        )
      )

      do.call(
        torch::nn_sequential,
        layers
      )
    }

    cox_ph_loss_torch <- function(log_h,
                                  durations,
                                  events,
                                  truncation,
                                  event_times,
                                  event_counts,
                                  n_events) {

      log_h <- log_h$view(c(-1))

      exp_log_h <- torch::torch_exp(log_h)

      loss <- log_h$sum() * 0

      for (kk in seq_along(event_times)) {

        t_j <- event_times[kk]
        d_j <- event_counts[kk]

        if (d_j <= 0L) {
          next
        }

        event_ind <- ((durations == t_j) * (events == 1))$
          to(dtype = torch::torch_float32())

        risk_ind <- ((durations >= t_j) * (truncation < t_j))$
          to(dtype = torch::torch_float32())

        risk_sum <- torch::torch_sum(
          risk_ind * exp_log_h
        )

        event_sum_exp <- torch::torch_sum(
          event_ind * exp_log_h
        )

        event_sum_log <- torch::torch_sum(
          event_ind * log_h
        )

        efron_fraction <- torch::torch_tensor(
          seq(0, d_j - 1) / d_j,
          dtype = torch::torch_float32()
        )

        denominators <- risk_sum -
          efron_fraction * event_sum_exp

        loss <- loss +
          torch::torch_sum(torch::torch_log(denominators)) -
          event_sum_log
      }

      loss / n_events
    }

    nn_elastic_net_penalty_torch <- function(net,
                                             xi,
                                             eps) {

      l2 <- torch::torch_tensor(
        0,
        dtype = torch::torch_float32()
      )

      l1 <- torch::torch_tensor(
        0,
        dtype = torch::torch_float32()
      )

      for (p in net$parameters) {
        if (length(p$shape) >= 2L) {
          l2 <- l2 + torch::torch_sum(p^2)
          l1 <- l1 + torch::torch_sum(torch::torch_abs(p))
        }
      }

      eps * (xi * l2 + (1 - xi) * l1)
    }

    train_deepsurv_fold <- function(fold_data,
                                    params,
                                    epochs,
                                    verbose,
                                    seed) {

      torch::torch_manual_seed(
        as.integer(seed)
      )

      input_dim <- fold_data$x_train_t$shape[[2]]

      net <- build_deepsurv_net(
        input_dim = input_dim,
        params = params
      )

      optimizer <- switch(
        as.character(params$optim),
        "Adam" = torch::optim_adam(
          net$parameters,
          lr = as.numeric(params$lr)
        ),
        "SGD" = torch::optim_sgd(
          net$parameters,
          lr = as.numeric(params$lr)
        ),
        "AdamW" = if (exists(
          "optim_adamw",
          envir = asNamespace("torch")
        )) {
          torch::optim_adamw(
            net$parameters,
            lr = as.numeric(params$lr)
          )
        } else {
          torch::optim_adam(
            net$parameters,
            lr = as.numeric(params$lr)
          )
        },
        stop(
          "Unknown optimizer: ",
          params$optim,
          call. = FALSE
        )
      )

      train_losses <- numeric(0)
      val_losses   <- numeric(0)

      best_val <- Inf
      wait <- 0L

      for (epoch in seq_len(epochs)) {

        net$train()

        optimizer$zero_grad()

        train_log_h <- net(fold_data$x_train_t)$squeeze()

        cox_loss <- cox_ph_loss_torch(
          log_h = train_log_h,
          durations = fold_data$dur_train_t,
          events = fold_data$event_train_t,
          truncation = fold_data$trunc_train_t,
          event_times = fold_data$train_event_times,
          event_counts = fold_data$train_event_counts,
          n_events = fold_data$train_n_events
        )

        reg <- nn_elastic_net_penalty_torch(
          net = net,
          xi = as.numeric(params$xi),
          eps = as.numeric(params$eps)
        )

        total_loss <- cox_loss + reg

        total_loss$backward()

        optimizer$step()

        net$eval()

        torch::with_no_grad({

          train_log_h_eval <- net(fold_data$x_train_t)$squeeze()

          t_loss <- cox_ph_loss_torch(
            log_h = train_log_h_eval,
            durations = fold_data$dur_train_t,
            events = fold_data$event_train_t,
            truncation = fold_data$trunc_train_t,
            event_times = fold_data$train_event_times,
            event_counts = fold_data$train_event_counts,
            n_events = fold_data$train_n_events
          )

          val_log_h <- net(fold_data$x_val_t)$squeeze()

          v_loss <- cox_ph_loss_torch(
            log_h = val_log_h,
            durations = fold_data$dur_val_t,
            events = fold_data$event_val_t,
            truncation = fold_data$trunc_val_t,
            event_times = fold_data$val_event_times,
            event_counts = fold_data$val_event_counts,
            n_events = fold_data$val_n_events
          )
        })

        t_loss_val <- t_loss$item()
        v_loss_val <- v_loss$item()

        train_losses <- c(train_losses, t_loss_val)
        val_losses   <- c(val_losses, v_loss_val)

        if (isTRUE(params$early_stopping)) {

          if (v_loss_val < best_val) {

            best_val <- v_loss_val
            wait <- 0L

          } else {

            wait <- wait + 1L

            if (wait >= as.integer(params$patience)) {
              if (verbose) {
                message("Early stopping at epoch ", epoch)
              }
              break
            }
          }
        }

        if (verbose) {
          message(
            sprintf(
              "Epoch %d | train_loss: %.6f | val_loss: %.6f",
              epoch,
              t_loss_val,
              v_loss_val
            )
          )
        }
      }

      data.frame(
        train_loss = train_losses,
        val_loss = val_losses
      )
    }

    nn_folds <- vector("list", folds)

    for (ff in seq_len(folds)) {
      nn_folds[[ff]] <- make_nn_fold(
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

      params <- scalar_list(
        hparameters.f[hp, , drop = FALSE]
      )

      tmp.train.lkh <- numeric(folds)
      tmp.test.lkh  <- numeric(folds)

      for (ff in seq_len(folds)) {

        log_k <- train_deepsurv_fold(
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

    best_params <- scalar_list(
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
