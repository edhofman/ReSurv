# R torch implementation of DeepSurv for ReSurv
#
# Replaces Python/reticulate NN backend with native R torch.
# Contains: network builder, Efron loss, training loop, prediction, baseline hazards.

## Network builder ----

pkg.env$build_deepsurv_net <- function(input_dim, params) {

  activation_map <- list(
    "ReLU"      = torch::nn_relu,
    "SELU"      = torch::nn_selu,
    "Tanh"      = torch::nn_tanh,
    "LeakyReLU" = torch::nn_leaky_relu
  )

  act_fn <- activation_map[[params$activation]]
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

pkg.env$cox_ph_loss_torch <- function(log_h, durations, events, truncation) {

  # 1. Sort by ascending duration
  sort_result <- torch::torch_sort(durations, descending = FALSE)
  idx <- sort_result[[2]]

  events     <- events[idx]
  durations  <- durations[idx]
  truncation <- truncation[idx]
  log_h      <- log_h[idx]$view(-1)$mul(events)  # zeros out censored

  # 2. Sorted labels as double
  labels <- torch::torch_sort(durations, descending = FALSE)[[1]]$to(torch::torch_double())

  unique_vals <- labels$unique()
  n_unique    <- unique_vals$shape[1]
  n           <- labels$shape[1]

  unique_list <- as.numeric(unique_vals)

  # 3. Build event_set, risk_set, efron_set matrices
  event_set  <- torch::torch_zeros(n_unique, n)
  risk_set   <- torch::torch_zeros(n_unique, n)
  efron_set  <- torch::torch_zeros(n_unique, n)

  for (i in seq_along(unique_list)) {
    u <- unique_list[i]
    event_set[i, ] <- (labels == u) * (events == 1)
    risk_set[i, ]  <- (durations >= u) * (truncation < u)

    d_i <- as.integer(torch::torch_sum(event_set[i, ])$item())
    if (d_i > 0) {
      tmp <- seq(0, d_i - 1)
      efron_set[i, tmp + 1L] <- torch::torch_tensor(as.numeric(tmp))
    }
  }

  # 4. Compute loss terms
  h_risk    <- risk_set$matmul(log_h$exp())
  h_events0 <- event_set$matmul(log_h$exp())

  efron_risk <- (1 / torch::torch_sum(event_set, 2L))$mul(h_events0)$reshape(c(-1, 1))$matmul(torch::torch_ones(1, n))
  efron_risk <- efron_risk$mul(efron_set)

  # Handle inf from division by zero (autograd-safe)
  h_risk <- torch::torch_where(h_risk$abs() == Inf, torch::torch_ones_like(h_risk), h_risk)

  events_sort <- torch::torch_sort(event_set, dim = 2, descending = TRUE)[[1]]

  efron_risk <- (h_risk$reshape(c(-1, 1))$sub(efron_risk))$log()$matmul(events_sort$t())$diagonal()
  efron_risk <- torch::torch_where(efron_risk$abs() == Inf, torch::torch_zeros_like(efron_risk), efron_risk)

  h_events <- event_set$matmul(log_h)

  log_like <- -h_events$sub(efron_risk)$sum()

  log_like$div(event_set$sum())
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

pkg.env$train_deepsurv <- function(net, x_train, y_train, x_val, y_val,
                                   params, epochs, verbose, seed) {

  torch::torch_manual_seed(seed)

  # Optimizer
  optimizer <- switch(params$optim,
    "Adam"  = torch::optim_adam(net$parameters, lr = params$lr),
    "SGD"   = torch::optim_sgd(net$parameters, lr = params$lr),
    "AdamW" = if (exists("optim_adamw", envir = asNamespace("torch"))) {
                torch::optim_adamw(net$parameters, lr = params$lr)
              } else {
                torch::optim_adam(net$parameters, lr = params$lr)
              },
    stop("Unknown optimizer: ", params$optim)
  )

  # Convert to tensors
  x_train_t <- torch::torch_tensor(as.matrix(x_train), dtype = torch::torch_float32())
  x_val_t   <- torch::torch_tensor(as.matrix(x_val),   dtype = torch::torch_float32())

  y_train <- as.matrix(y_train)
  y_val   <- as.matrix(y_val)

  dur_train_t   <- torch::torch_tensor(y_train[, 1], dtype = torch::torch_float32())
  event_train_t <- torch::torch_tensor(y_train[, 2], dtype = torch::torch_float32())
  trunc_train_t <- torch::torch_tensor(y_train[, 3], dtype = torch::torch_float32())

  dur_val_t   <- torch::torch_tensor(y_val[, 1], dtype = torch::torch_float32())
  event_val_t <- torch::torch_tensor(y_val[, 2], dtype = torch::torch_float32())
  trunc_val_t <- torch::torch_tensor(y_val[, 3], dtype = torch::torch_float32())

  n_train    <- nrow(y_train)
  batch_size <- as.integer(params$batch_size)

  train_losses <- numeric(0)
  val_losses   <- numeric(0)
  best_val     <- Inf
  wait         <- 0L

  for (epoch in seq_len(epochs)) {

    # --- Train ---
    net$train()
    perm <- sample.int(n_train)

    for (start in seq(1, n_train, by = batch_size)) {
      end <- min(start + batch_size - 1L, n_train)
      batch_idx <- perm[start:end]

      x_batch   <- x_train_t[batch_idx, , drop = FALSE]
      dur_b     <- dur_train_t[batch_idx]
      event_b   <- event_train_t[batch_idx]
      trunc_b   <- trunc_train_t[batch_idx]

      log_h <- net(x_batch)$squeeze(2)

      cox_loss <- pkg.env$cox_ph_loss_torch(log_h, dur_b, event_b, trunc_b)

      # Elastic net regularisation on weights (not biases)
      l2 <- torch::torch_tensor(0)
      l1 <- torch::torch_tensor(0)
      for (p in net$parameters) {
        if (length(p$shape) >= 2) {
          l2 <- l2 + torch::torch_sum(p^2)
          l1 <- l1 + torch::torch_sum(torch::torch_abs(p))
        }
      }
      reg <- params$eps * (params$xi * l2 + (1 - params$xi) * l1)
      total_loss <- cox_loss + reg

      optimizer$zero_grad()
      total_loss$backward()
      optimizer$step()
    }

    # --- Epoch-level losses (no grad, eval mode) ---
    net$eval()
    torch::with_no_grad({
      train_log_h <- net(x_train_t)$squeeze(2)
      t_loss <- pkg.env$cox_ph_loss_torch(train_log_h, dur_train_t, event_train_t, trunc_train_t)

      val_log_h <- net(x_val_t)$squeeze(2)
      v_loss <- pkg.env$cox_ph_loss_torch(val_log_h, dur_val_t, event_val_t, trunc_val_t)
    })

    t_loss_val <- t_loss$item()
    v_loss_val <- v_loss$item()

    train_losses <- c(train_losses, t_loss_val)
    val_losses   <- c(val_losses, v_loss_val)

    # --- Early stopping ---
    if (isTRUE(params$early_stopping)) {
      if (v_loss_val < best_val) {
        best_val <- v_loss_val
        wait <- 0L
      } else {
        wait <- wait + 1L
        if (wait >= params$patience) {
          if (verbose) message("Early stopping at epoch ", epoch)
          break
        }
      }
    }

    if (verbose) {
      message(sprintf("Epoch %d | train_loss: %.6f | val_loss: %.6f", epoch, t_loss_val, v_loss_val))
    }
  }

  list(
    net = net,
    log = data.frame(train_loss = train_losses, val_loss = val_losses)
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
