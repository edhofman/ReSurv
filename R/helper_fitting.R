# Model fitting helper functions
#
# Scalers, DeepSurv preprocessing, Cox and NN fitting routines.
#
## Deepsurv helpers ----

pkg.env$deep_surv_pp <- function(X,
                                 Y,
                                 training_test_split = 0.8,
                                 samples_TF = NULL) {

  tmp_order <- order(Y$DP_rev_i)

  X <- as.data.frame(X[tmp_order, , drop = FALSE])
  Y <- as.data.frame(Y[tmp_order, , drop = FALSE])

  tmp <- data.frame(id = seq_len(nrow(X)))

  if (is.null(samples_TF)) {
    if (!is.numeric(training_test_split) ||
        length(training_test_split) != 1L ||
        !is.finite(training_test_split) ||
        training_test_split <= 0 ||
        training_test_split > 1) {
      stop("`training_test_split` must be a number in (0, 1].",
           call. = FALSE)
    }

    if (training_test_split == 1) {
      sampled_id <- tmp$id
    } else {
      n_sample <- ceiling(nrow(tmp) * training_test_split)
      n_sample <- max(1L, min(nrow(tmp), n_sample))
      sampled_id <- sample(tmp$id, size = n_sample, replace = FALSE)
    }

    id_train <- tmp$id %in% sampled_id
  } else {
    if (length(samples_TF) != nrow(X)) {
      stop("`samples_TF` must have length equal to the number of rows in `X`.",
           call. = FALSE)
    }

    id_train <- as.logical(samples_TF)

    if (anyNA(id_train)) {
      stop("`samples_TF` must be coercible to TRUE/FALSE without NA values.",
           call. = FALSE)
    }
  }

  x_train <- as.matrix(X[id_train, , drop = FALSE])
  x_val   <- as.matrix(X[!id_train, , drop = FALSE])
  y_train <- as.matrix(Y[id_train, , drop = FALSE])
  y_val   <- as.matrix(Y[!id_train, , drop = FALSE])

  list(
    x_train = x_train,
    y_train = y_train,
    x_val   = x_val,
    y_val   = y_val,
    lkh_eval_data = list(
      data_train = X[id_train, , drop = FALSE],
      data_val   = X[!id_train, , drop = FALSE],
      y_train    = Y[id_train, , drop = FALSE],
      y_val      = Y[!id_train, , drop = FALSE]
    )
  )
}
## Fitting routines ----

pkg.env$fit_cox_model <- function(data,
                                  formula_ct,
                                  newdata){
  "This function is the fitting routine for the cox model."

  cox <- survival::coxph(formula_ct, data=data, ties="efron")
  cox_lp <- predict(cox,newdata=newdata,'lp',reference='zero')

  cox_training_lp <- predict(cox,newdata=data %>% dplyr::arrange(DP_rev_i) %>% as.data.frame(),'lp',reference='zero')

  out <- list(
    cox=cox,
    cox_lp=cox_lp,
    expg = exp(cox_lp),
    train_expg= cox_training_lp#exp(cox_training_lp)
  )

  return(out)
}


pkg.env$fit_deep_surv <- function(data,
                                  params,
                                  verbose,
                                  epochs,
                                  num_workers,
                                  seed = as.numeric(Sys.time()),
                                  network_structure=NULL,
                                  newdata){

  input_dim <- ncol(data$x_train)
  net <- pkg.env$build_deepsurv_net(input_dim, params)

  result <- pkg.env$train_deepsurv(
    net       = net,
    x_train   = data$x_train,
    y_train   = data$y_train,
    x_val     = data$x_val,
    y_val     = data$y_val,
    params    = params,
    epochs    = epochs,
    verbose   = verbose,
    seed      = seed
  )

  return(result)

}


