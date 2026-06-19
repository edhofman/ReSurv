# Internal ReSurvCV method helpers
#
# These helpers support the current ReSurvCV.IndividualDataPP() implementation:
# design matrix construction, fast XGBoost CV folds, and torch fold training.
# They are internal and not exported.

cv_scalar_list <- function(x) {
  x <- as.list(x[1L, , drop = FALSE])
  lapply(x, function(z) {
    if (length(z) == 1L) {
      z[[1L]]
    } else {
      z
    }
  })
}

cv_extract_metric_value <- function(metric_object) {
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

cv_make_efron_c <- function(z) {
  ave(
    seq_along(z),
    z,
    FUN = function(ind) {
      (seq_along(ind) - 1) / length(ind)
    }
  )
}

cv_get_xgb_eval_value <- function(model.out,
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

cv_predict_xgb_at <- function(model.out,
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

cv_scale_one <- function(x, continuous_features_scaling_method) {
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

cv_method_design_matrix <- function(IndividualDataPP,
                                    model,
                                    continuous_features_scaling_method) {

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

      X_cont <- as.data.frame(
        lapply(
          training_data[, continuous_features, drop = FALSE],
          cv_scale_one,
          continuous_features_scaling_method = continuous_features_scaling_method
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

  list(X = X, Y = Y)
}

cv_weighted_tabulate <- function(idx, weights, nbins) {

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

cv_xgb_fast_quantities <- function(preds, dtrain) {

  Ti <- as.integer(attr(dtrain, "truncation"))
  Ei <- as.integer(attr(dtrain, "claim_arrival"))
  efron_c <- as.numeric(attr(dtrain, "efron_c"))

  exp_p <- exp(preds)

  max_time <- max(c(Ei, Ti + 1L), na.rm = TRUE)

  start_idx <- Ti + 1L
  end_idx   <- Ei + 1L


  add <- cv_weighted_tabulate(
    idx = start_idx,
    weights = exp_p,
    nbins = max_time + 1L
  )

  remove <- cv_weighted_tabulate(
    idx = end_idx,
    weights = exp_p,
    nbins = max_time + 1L
  )
  risk_sum <- cumsum(add - remove)[seq_len(max_time)]

  event_sum <- cv_weighted_tabulate(
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

cv_xgb_fast_evaluation_metrics <- function(preds, dtrain) {

  q <- cv_xgb_fast_quantities(preds, dtrain)

  value <- (
    sum(log(q$denom)) -
      sum(preds)
  ) / length(preds)

  list(
    metric = "log_partial_likelihood",
    value = value
  )
}

cv_xgb_fast_loss_objective <- function(preds, dtrain) {

  q <- cv_xgb_fast_quantities(preds, dtrain)

  inv_denom <- 1 / q$denom

  alpha_i <- cv_weighted_tabulate(
    idx = q$Ei,
    weights = inv_denom,
    nbins = q$max_time
  )

  beta_i <- cv_weighted_tabulate(
    idx = q$Ei,
    weights = q$efron_c * inv_denom,
    nbins = q$max_time
  )

  gamma_i <- cv_weighted_tabulate(
    idx = q$Ei,
    weights = inv_denom^2,
    nbins = q$max_time
  )

  omega_i <- cv_weighted_tabulate(
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

cv_make_xgb_fold <- function(X, Y, samples_TF) {

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
    attr(dmat, "efron_c") <- as.numeric(cv_make_efron_c(df$DP_rev_i))

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

cv_make_event_info <- function(y) {

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

cv_make_nn_fold <- function(X, Y, samples_TF) {

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

  train_info <- cv_make_event_info(y_train)
  val_info   <- cv_make_event_info(y_val)

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

cv_build_deepsurv_net <- function(input_dim, params) {

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

cv_cox_ph_loss_torch <- function(log_h,
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

cv_nn_elastic_net_penalty_torch <- function(net,
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

cv_train_deepsurv_fold <- function(fold_data,
                                   params,
                                   epochs,
                                   verbose,
                                   seed) {

  torch::torch_manual_seed(
    as.integer(seed)
  )

  input_dim <- fold_data$x_train_t$shape[[2]]

  net <- cv_build_deepsurv_net(
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

    cox_loss <- cv_cox_ph_loss_torch(
      log_h = train_log_h,
      durations = fold_data$dur_train_t,
      events = fold_data$event_train_t,
      truncation = fold_data$trunc_train_t,
      event_times = fold_data$train_event_times,
      event_counts = fold_data$train_event_counts,
      n_events = fold_data$train_n_events
    )

    reg <- cv_nn_elastic_net_penalty_torch(
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

      t_loss <- cv_cox_ph_loss_torch(
        log_h = train_log_h_eval,
        durations = fold_data$dur_train_t,
        events = fold_data$event_train_t,
        truncation = fold_data$trunc_train_t,
        event_times = fold_data$train_event_times,
        event_counts = fold_data$train_event_counts,
        n_events = fold_data$train_n_events
      )

      val_log_h <- net(fold_data$x_val_t)$squeeze()

      v_loss <- cv_cox_ph_loss_torch(
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
