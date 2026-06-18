# R torch implementation of DeepSurv for ReSurv
#
# Native R torch implementation of the NN backend.
# Contains: network builder, Efron loss, training loop, prediction, baseline hazards.

## Network builder ----

pkg.env$build_deepsurv_net <- function(input_dim, params) {

  activation_map <- list(
    "relu"       = torch::nn_relu,
    "selu"       = torch::nn_selu,
    "tanh"       = torch::nn_tanh,
    "leakyrelu"  = torch::nn_leaky_relu,
    "leaky_relu" = torch::nn_leaky_relu
  )

  activation <- tolower(as.character(params$activation))
  act_fn <- activation_map[[activation]]
  if (is.null(act_fn)) {
    stop("Unknown activation: ", params$activation)
  }

  layers <- list()
  in_features <- input_dim

  for (i in seq_len(params$num_layers)) {
    node_i <- as.integer(params[[paste0("node_", i)]])
    layers <- c(layers, list(
      torch::nn_linear(in_features, node_i),
      act_fn()
    ))
    in_features <- node_i
  }

  # Final output layer: single node, no bias, no activation

  layers <- c(layers, list(
    torch::nn_linear(in_features, 1L, bias = FALSE)
  ))

  do.call(torch::nn_sequential, layers)
}

## Efron partial likelihood loss with left truncation ----

pkg.env$cox_ph_loss_torch <- function(log_h,
                                      durations,
                                      events,
                                      truncation) {
  # log_h:      NN output phi(U_i, X_i)
  # durations:  DP_rev_i
  # events:     event indicator, usually I = 1 in ReSurv
  # truncation: TR_i = AP_i - 1
  #
  # Implements the negative Efron partial log-likelihood:
  #
  # sum_j sum_{r=0}^{O_j - 1}
  #   log( sum_{l in R_j} exp(phi_l)
  #        - r/O_j * sum_{s in O_j} exp(phi_s) )
  # - sum_j sum_{i in O_j} phi_i
  #
  # averaged by the number of events.

  log_h      <- log_h$view(c(-1))
  durations  <- durations$view(c(-1))
  events     <- events$view(c(-1))
  truncation <- truncation$view(c(-1))

  # These are labels only; they do not need gradients.
  durations_r <- as.numeric(durations)
  events_r    <- as.numeric(events)

  event_times <- sort(unique(durations_r[events_r == 1]))

  if (length(event_times) == 0L) {
    stop("No events found in cox_ph_loss_torch().")
  }

  exp_log_h <- torch::torch_exp(log_h)

  # Tensor-valued zero, preserving dtype/autograd graph.
  loss <- log_h$sum() * 0

  n_events <- sum(events_r == 1)

  for (t_j in event_times) {
    # Occurrence set O(t_j)
    event_ind <- ((durations == t_j) * (events == 1))$
      to(dtype = torch::torch_float32())

    # Risk/exposure set R(t_j), with left truncation.
    # Existing XGB code uses: starts_i < stop & stops_i >= stop.
    risk_ind <- ((durations >= t_j) * (truncation < t_j))$
      to(dtype = torch::torch_float32())

    d_j <- sum(durations_r == t_j & events_r == 1)

    if (d_j <= 0L) {
      next
    }

    risk_sum      <- torch::torch_sum(risk_ind * exp_log_h)
    event_sum_exp <- torch::torch_sum(event_ind * exp_log_h)
    event_sum_log <- torch::torch_sum(event_ind * log_h)

    # Efron correction:
    # denominator_r = risk_sum - (r / d_j) * event_sum_exp,
    # r = 0, ..., d_j - 1.
    efron_fraction <- torch::torch_tensor(
      seq(0, d_j - 1) / d_j,
      dtype = torch::torch_float32()
    )

    denominators <- risk_sum - efron_fraction * event_sum_exp

    # Do not silently repair invalid denominators. They indicate a coding or
    # data problem. This check is outside autograd, so it is only diagnostic.
    if (any(as.numeric(denominators) <= 0)) {
      stop(
        "Non-positive denominator in Efron partial likelihood. ",
        "Check risk sets, truncation times, and event times."
      )
    }

    loss <- loss +
      torch::torch_sum(torch::torch_log(denominators)) -
      event_sum_log
  }

  loss / n_events
}


pkg.env$nn_elastic_net_penalty_torch <- function(net,
                                                 xi,
                                                 eps) {
  # The paper writes:
  #   loss + rho * (epsilon * ||theta||_2^2
  #                 + (1 - epsilon) * ||theta||_1)
  #
  # Here I preserve the package's current convention:
  #   eps = penalty strength
  #   xi  = L2/L1 mixing parameter

  l2 <- torch::torch_tensor(0, dtype = torch::torch_float32())
  l1 <- torch::torch_tensor(0, dtype = torch::torch_float32())

  for (p in net$parameters) {
    # Penalize weight matrices, not biases.
    if (length(p$shape) >= 2L) {
      l2 <- l2 + torch::torch_sum(p^2)
      l1 <- l1 + torch::torch_sum(torch::torch_abs(p))
    }
  }

  eps * (xi * l2 + (1 - xi) * l1)
}

## Breslow baseline hazard estimator ----

pkg.env$compute_baseline_hazards_r <- function(net, input, df_target) {

  net$eval()

  x_tensor <- torch::torch_tensor(as.matrix(input), dtype = torch::torch_float32())

  predictions <- torch::with_no_grad({
    net(x_tensor)
  })

  pred <- exp(as.numeric(predictions))

  duration   <- df_target$duration
  event      <- df_target$event
  truncation <- df_target$truncation

  unique_times <- sort(unique(duration[event == 1]))

  bh <- numeric(length(unique_times))
  names(bh) <- as.character(unique_times)

  for (i in seq_along(unique_times)) {
    t_i <- unique_times[i]
    risk_mask <- (duration >= t_i) & (truncation < t_i)
    risk_sum  <- sum(pred[risk_mask])
    n_events  <- sum((duration == t_i) & (event == 1))
    bh[i] <- n_events / risk_sum
  }

  bh
}

## Training loop ----

pkg.env$train_deepsurv <- function(net,
                                   x_train,
                                   y_train,
                                   x_val,
                                   y_val,
                                   params,
                                   epochs,
                                   verbose,
                                   seed) {
  torch::torch_manual_seed(seed)

  optimizer <- switch(params$optim,
                      "Adam" = torch::optim_adam(net$parameters, lr = params$lr),
                      "SGD" = torch::optim_sgd(net$parameters, lr = params$lr),
                      "AdamW" = if (exists("optim_adamw", envir = asNamespace("torch"))) {
                        torch::optim_adamw(net$parameters, lr = params$lr)
                      } else {
                        torch::optim_adam(net$parameters, lr = params$lr)
                      },
                      stop("Unknown optimizer: ", params$optim)
  )

  x_train_t <- torch::torch_tensor(
    as.matrix(x_train),
    dtype = torch::torch_float32()
  )

  x_val_t <- torch::torch_tensor(
    as.matrix(x_val),
    dtype = torch::torch_float32()
  )

  y_train <- as.matrix(y_train)
  y_val   <- as.matrix(y_val)

  dur_train_t <- torch::torch_tensor(
    y_train[, 1],
    dtype = torch::torch_float32()
  )
  event_train_t <- torch::torch_tensor(
    y_train[, 2],
    dtype = torch::torch_float32()
  )
  trunc_train_t <- torch::torch_tensor(
    y_train[, 3],
    dtype = torch::torch_float32()
  )

  dur_val_t <- torch::torch_tensor(
    y_val[, 1],
    dtype = torch::torch_float32()
  )
  event_val_t <- torch::torch_tensor(
    y_val[, 2],
    dtype = torch::torch_float32()
  )
  trunc_val_t <- torch::torch_tensor(
    y_val[, 3],
    dtype = torch::torch_float32()
  )

  train_losses <- numeric(0)
  val_losses   <- numeric(0)

  best_val <- Inf
  wait     <- 0L

  for (epoch in seq_len(epochs)) {
    net$train()

    optimizer$zero_grad()

    # Full-batch NN output.
    # This is necessary because Cox risk sets are global.
    train_log_h <- net(x_train_t)$squeeze()

    cox_loss <- pkg.env$cox_ph_loss_torch(
      log_h      = train_log_h,
      durations  = dur_train_t,
      events     = event_train_t,
      truncation = trunc_train_t
    )

    reg <- pkg.env$nn_elastic_net_penalty_torch(
      net = net,
      xi  = params$xi,
      eps = params$eps
    )

    total_loss <- cox_loss + reg

    total_loss$backward()
    optimizer$step()

    # Epoch-level diagnostics.
    net$eval()

    torch::with_no_grad({
      train_log_h_eval <- net(x_train_t)$squeeze()
      t_loss <- pkg.env$cox_ph_loss_torch(
        log_h      = train_log_h_eval,
        durations  = dur_train_t,
        events     = event_train_t,
        truncation = trunc_train_t
      )

      val_log_h <- net(x_val_t)$squeeze()
      v_loss <- pkg.env$cox_ph_loss_torch(
        log_h      = val_log_h,
        durations  = dur_val_t,
        events     = event_val_t,
        truncation = trunc_val_t
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

        if (wait >= params$patience) {
          if (verbose) {
            message("Early stopping at epoch ", epoch)
          }
          break
        }
      }
    }

    if (verbose) {
      message(sprintf(
        "Epoch %d | train_loss: %.6f | val_loss: %.6f",
        epoch,
        t_loss_val,
        v_loss_val
      ))
    }
  }

  list(
    net = net,
    log = data.frame(
      train_loss = train_losses,
      val_loss   = val_losses
    )
  )
}
## Prediction ----

pkg.env$predict_deepsurv <- function(net, x) {

  net$eval()
  x_tensor <- torch::torch_tensor(as.matrix(x), dtype = torch::torch_float32())

  out <- torch::with_no_grad({
    net(x_tensor)
  })

  as.numeric(out)
}
