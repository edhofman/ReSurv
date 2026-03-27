# Model fitting helper functions
#
# Scalers, DeepSurv preprocessing, Cox and NN fitting routines.
#
# @import reticulate
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


  #convert to array for later numpy transforamtion
  data_train <- as.array(as.matrix(X[id_train,]))
  data_val <- as.array(as.matrix(X[!id_train,]))
  y_train <- as.array(as.matrix(Y[id_train,]))
  y_val <- as.array(as.matrix(Y[!id_train,]))


  #create tuples holding target and validation values. Convert to same dtype to ensure safe pytorch handling.
  y_train <- reticulate::tuple(reticulate::np_array(y_train[,1], dtype = "float32"), #duration
                               reticulate::np_array(y_train[,2], dtype = "float32"), #event
                               reticulate::np_array(y_train[,3], dtype = "float32")) #truncation

  validation_data = reticulate::tuple(reticulate::np_array(data_val, dtype = "float32"),
                                      reticulate::tuple(reticulate::np_array(y_val[,1], dtype = "float32"), #duration
                                                        reticulate::np_array(y_val[,2], dtype = "float32"), #event
                                                        reticulate::np_array(y_val[,3], dtype = "float32"))) #truncation

  x_train = reticulate::np_array(data_train, dtype = "float32")


  return(list(
    x_train = x_train,
    y_train = y_train,
    validation_data = validation_data,
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


  # #Import python modules

  torchtuples <- reticulate::import("torchtuples")
  torch <- reticulate::import("torch")

  #Source python code for left truncated deepsurv
  reticulate::source_python(system.file("python", "coxnetwork_custom.py", package = "ReSurv"))

  torch$manual_seed(seed)

  net <- torch$nn$Sequential()
  input_shape =  data$x_train$shape[[1]]
  for( i in 1:(params$num_layers+1)){
    if( i > params$num_layers){
      net$add_module(paste0(i,"_l"),torch$nn$Linear(input_shape, as.integer(1), bias=FALSE))
    }
    else{
      net$add_module(paste0(i,"_l"),torch$nn$Linear(input_shape,as.integer(params[[paste0("node_",i)]] )))
      net$add_module(paste0(i,"_a"),torch$nn[[params$activation]]())
      input_shape = as.integer(params[[paste0("node_",i)]] )
    }
  }


  # Setup CoxPH model, as imported from python script. Seed is for weight initlization and comparability when doing cv.

  model <- CoxPH(
    net = net,
    optimizer = torchtuples$optim[[params$optim]](lr=params$lr),
    xi=params$xi,
    eps=params$eps,
    tie = params$tie
  )


  #If early stopping specified add to callbacks.
  if(params$early_stopping==TRUE){
    callbacks = list(torchtuples$callbacks$EarlyStopping(patience=as.integer(params$patience)))
  }else{
    callbacks = NULL
  }

  #fit model
  model$fit(
    input = data$x_train,
    target = data$y_train,
    batch_size = as.integer(params$batch_size),
    epochs = epochs,
    callbacks = r_to_py(callbacks),
    verbose = verbose,
    val_data=data$validation_data,
    val_batch_size=params$batch_size,
    num_workers=num_workers
  )

  return(model)

}


