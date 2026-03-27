# Model fitting helper functions
#
# Scalers, DeepSurv preprocessing, Cox and NN fitting routines.
#
## Deepsurv helpers ----

pkg.env$deep_surv_pp <- function(X,
                                 Y,
                                 training_test_split,
                                 samples_TF=NULL){


  X <- cbind(X, DP_rev_i = Y$DP_rev_i) %>%
    arrange(DP_rev_i) %>%
    select(-DP_rev_i)

  Y <- Y %>%
    arrange(DP_rev_i) %>%
    as.data.frame()


  tmp <- as.data.frame(seq(1,dim(X)[1]))
  colnames(tmp) <- "id"

  if(is.null(samples_TF)){

    samples_cn <- tmp %>% sample_frac(size=training_test_split)
    id_train <- tmp$id %in% samples_cn$id

  }else{

    cond <- samples_TF
    samples_cn <- tmp %>% select(id) %>% filter(cond)
    id_train <- tmp$id %in% samples_cn$id
  }


  # Plain R matrices for native torch training
  x_train <- as.matrix(X[id_train, ])
  x_val   <- as.matrix(X[!id_train, ])
  y_train <- as.matrix(Y[id_train, ])   # columns: duration, event, truncation
  y_val   <- as.matrix(Y[!id_train, ])

  return(list(
    x_train = x_train,
    y_train = y_train,
    x_val   = x_val,
    y_val   = y_val,
    lkh_eval_data = list(data_train=X[id_train,],
                       data_val=X[!id_train,],
                       y_train=Y[id_train,],
                       y_val=Y[!id_train,])
  ))

}


## Fitting routines ----

pkg.env$fit_cox_model <- function(data,
                                  formula_ct,
                                  newdata){
  "This function is the fitting routine for the cox model."

  cox <- coxph(formula_ct, data=data, ties="efron")
  cox_lp <- predict(cox,newdata=newdata,'lp',reference='zero')

  cox_training_lp <- predict(cox,newdata=data %>% arrange(DP_rev_i) %>% as.data.frame(),'lp',reference='zero')

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


