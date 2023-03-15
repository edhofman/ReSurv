pkg.env <- new.env()

pkg.env$check.all.present <- function(x,check.on){

  "
  This function checks that you have all the periods in the data,
  from the minimum record to the maximum record.

  "

  tmp <- x

  v <- diff(as.integer(sort(unique(x))))

  if(sum(v>1)>0){

    warning(paste0("Some", check.on, "are missing in the data"))

  }

}

pkg.env$check.time.units <- function(input_time_unit,
                             output_time_unit){

  "
  This function checks that the input output time transformation is consistent.
  E.g. you cannot turn quarters into semesters, you can instead turn trimesters into semesters.

  "

  if((1/input_time_unit)%%(1/output_time_unit) != 0){

    stop('The provided time intervals are not subsettable.')

  }

}

pkg.env$check.traintestsplit <- function(x){

  "
  This function checks that the training test split is specified correctly

  "

  if(x>1 | x<0 ){
    warning(paste0("Traintestsplit has been put to ", x,". The value needs to be between 0 and 1, defaulting to 0.8."))
    return(.8)
    }else{return(x)}


}

pkg.env$encode.variables <- function(x){
  "
  This function encodes the periods.
  We impose that the indexization starts from 1.

  "
  seq <- min(x):max(x)

  dim = length(seq)

  v <- 1:length(seq)

  v[match(x,seq)]


}


pkg.env$formula.editor <- function(continuous_features,
                           categorical_features,
                           continuous_features_spline,
                           degrees_of_freedom,
                           input_output='i'){
  "
  This util edits creates the string that is used for model fitting in a compact way.
  @param continuous_features: string vector of continuous features to be included in the linear predictor.
  @param categorical_features: string vector of categorical features to be included in the linear predictor.
  @param continuous_features_spline: logical value, T if a spline is added to model the continuous features.
  @param degrees_of_freedom: degrees of freedom of the spline.
  @param input_output: set to input ('i') or output ('o') depending on the formula that we require.
  "


  tmp.cat <- switch(!is.null(categorical_features), paste(categorical_features, collapse='+'), NULL)
  tmp.cont <- switch(!is.null(continuous_features), paste(continuous_features, collapse='+'), NULL)
  tmp.splines <- switch((!is.null(continuous_features) & continuous_features_spline),paste0("pspline(",continuous_features, ",degree=3,df=",degrees_of_freedom,")"),NULL)

  string_formula<- paste(paste0("survival::Surv","(TR_",input_output,", DP_rev_",input_output,", I) ~ "),paste(c(tmp.cat,tmp.cont,tmp.splines), collapse='+'))
  string_formula


}


"This is a vectorized version of the grepl function.
See the grepl function documentation."
pkg.env$vgrepl <- Vectorize(grepl, vectorize.args = "pattern")


pkg.env$model.matrix.creator <- function(data,
                                 select_columns){
  "
  This function encodes the matrices that we need for model fitting.

  "

  #individual_data$training.data
  X <- data %>%
    fastDummies::dummy_cols(select_columns = select_columns, #individual_data$categorical_features
                            remove_selected_columns = T)

  tmp.cond=as.logical(apply(pkg.env$vgrepl(pattern=select_columns,
                                   x=colnames(X)), #individual_data$categorical_features
                            MARGIN=1,
                            sum))

  X <- X[,tmp.cond]

  return(X)

}


pkg.env$model.matrix.extract.hazard.names <- function(X,
                                                      string_formula,
                                                      data){

  formula_ct <- as.formula(string_formula)
  Y<-model.extract(model.frame(formula_ct, data=data),"response")

  enter <- Y[, 1]
  exit <- Y[, 2]
  event <- Y[, 3] != 0
  sco <- exp(rep(0, nrow(Y)))

  time <- sort(seq(1,max(exit[event]), by=1)) #might be times with no events

  X_unique <- unique(X)

  names_hazard <- (data.frame(X_unique) %>%
                     rowwise() %>%
                     mutate(name = paste0(names(.)[c_across() == 1], collapse = ',')))$name

  return(list(enter=enter,
              exit=exit,
              event=event,
              sco=sco,
              time=time,
              names_hazard=names_hazard))


}


pkg.env$MinMaxScaler <- function(x, na.rm = TRUE) {
  "MinMax Scaler"
  return((x- min(x)) /(max(x)-min(x)))
}

pkg.env$scaler <- function(continuous.features.scaling.method){
  ""
  if(continuous.features.scaling.method == "minmax" ){return(pkg.env$MinMaxScaler)}


}


pkg.env$deep_surv_pp <- function(X,
                    Y,
                    training_test_split){

  # data_transformed <- cbind(X, Y)

  id_train <- sample(c(TRUE,FALSE), nrow(X), replace=T, prob= c(training_test_split,1-training_test_split) )


  #convert to array for later numpy transforamtion
  data_train <- as.array(as.matrix(X[id_train,]))
  y_train <- as.array(as.matrix(Y[id_train,]))
  data_val <- as.array(as.matrix(X[!id_train,]))
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
    validation_data = validation_data
  ))

}


# Fitting routines

pkg.env$fit_cox_model <- function(data,
                          formula_ct,
                          X,
                          X_ams){
  "This function is the fitting routine for the cox model."

  cox <- coxph(formula_ct, data=data)
  beta_2 <- cox$coef
  beta<-c(0,beta_2)


  Xb <- as.matrix(X)%*%beta

  X_ams <- cbind(X_ams, Xb)

  beta_ams = unique(round(X_ams,10) )[,ncol(X_ams)] #if no round some systems has too high precision.


  list(
    cox=cox,
    beta=beta,
    beta_2=beta_2,
    Xb=Xb,
    X_ams=X_ams,
    beta_ams=beta_ams
  )
}



pkg.env$fit_deep_surv <- function(data,
                                  hparameters,
                                  network_structure=NULL){


  # #Import python modules
  torchtuples <- reticulate::import("torchtuples")
  torch <- reticulate::import("torch")

  #Source python code for left truncated deepsurv
  reticulate::source_python("./inst/python/coxnetwork_custom.py")

  # reticulate::py_run_file(system.file("python", "coxnetwork_custom.py", package = "ReSurv"))
  # print('test')
  #create network structure, if correct class specified in network_strucutre use this, otherwise use other variables.
  # if(!("torch.nn.modules.container.Sequential" %in% class(network_structure$net))){
  #   net <- torch$nn$Sequential()
  #   input_shape =  data$x_train$shape[[1]]
  #   for( i in 1:length(hparameters$num_nodes)+1){
  #     if( i > length(hparameters$num_nodes)){
  #       net$append(torch$nn$Linear(input_shape, as.integer(1)))
  #     }else{
  #       net$append(torch$nn$Linear(input_shape,as.integer(hparameters$num_nodes[i])))
  #       net$append(torch$nn[[activation]]() )
  #
  #       input_shape = as.integer(hparameters$num_nodes[i])
  #     }
  #   }
  # }

  if(!("torch.nn.modules.container.Sequential" %in% class(network_structure$net))){
    net <- torch$nn$Sequential()
    input_shape =  data$x_train$shape[[1]]
    for( i in 1:length(hparameters$num_nodes)+1){
      if( i > length(hparameters$num_nodes)){
        # net$append(torch$nn$Linear(input_shape, as.integer(1)))
        net$add_module('2',torch$nn$Linear(input_shape, as.integer(1)))
      }
      else{
        net$add_module('0',torch$nn$Linear(input_shape,as.integer(hparameters$num_nodes[i])))
        net$add_module('1',torch$nn[[hparameters$activation]]())
        input_shape = as.integer(hparameters$num_nodes[i])
      }
    }
  }


  #Setup batchsize, epochs and verbose settings
  batch_size = as.integer(hparameters$batch_size)

  epochs = as.integer(hparameters$epochs)

  verbose = T



  # Setup CoxPH model, as imported from python script.
  model <- CoxPH(
    net = net,
    optimizer = torchtuples$optim$Adam(lr = 0.005),
    xi=hparameters$xi,
    eps=hparameters$epsilon,
    tie = hparameters$tie
  )


  #If early stopping specified add to callbacks.
  if(hparameters$early_stopping==TRUE){
    callbacks = list(torchtuples$callbacks$EarlyStopping(patience=hparameters$patience))
  }else{
    callbacks = NULL
  }

  #fit model
  model$fit(
    input = data$x_train,
    target = data$y_train,
    batch_size = hparameters$batch_size,
    epochs = hparameters$epochs,
    callbacks = r_to_py(callbacks),
    verbose = verbose,
    val_data=data$validation_data,
    val_batch_size=hparameters$batch_size
  )

  return(model)

}

# Handling the baseline

pkg.env$hazard_baseline_model <- function(data,
                                  cox,
                                  hazard=NULL,
                                  baseline,
                                  conversion_factor,
                                  nk=50,
                                  nbin=48,
                                  phi=1){


  if(baseline == "breslow"){

    bs_hazard <- basehaz(cox, centered=FALSE) %>%
      mutate(hazard = hazard-lag(hazard,default=0))

    bs_hazard2 = tibble(DP_rev_i = bs_hazard$time,
                        hazard=bs_hazard$hazard) %>%
      mutate(hazard = hazard-lag(hazard, default=0))
  }

  if(baseline == "spline"){
    bs_hazard=bshazard::bshazard(pkg.env$formula.editor(continuous_features=NULL,
                                                categorical_features="1",
                                                continuous_features_spline=F),
                                 data=data[(data$AP_i-1)%%(conversion_factor^-1)==0 & data$claim_type==0,],
                                 nk=nk,
                                 nbin=nbin,
                                 phi=phi)
    bs_hazard2 <- tibble(DP_rev_i = bs_hazard$time,
                         hazard = bs_hazard$hazard)
  }

  return(list(bs_hazard=bs_hazard,
              bs_hazard2=bs_hazard2))

}

# Hazard computation

pkg.env$hazard_f<-function(i,
                   enter,
                   time,
                   exit,
                   event){
  "
  This function computes the hazard over a matrix.

  "
  rs <- (enter < time[i] & exit >= time[i])

  O<- sum(event[exit == time[i]])
  E <- sum(rs,-1/2*event[exit == time[i]] )

  c(O/E)}



pkg.env$m_to_q_hazard <- function(i,
                          hazard,
                          frame_tmp,
                          frame_tmp2){
  "
  Group monthly hazard to quarterly.

  "

  h_tmp <- data.frame(hazard[,(i*2+i-1):(i*2+i+1)]) %>%  mutate(DP_rev_i = row_number()) %>%
    filter(DP_rev_i<=max(frame_tmp$time_w)) %>%
    filter(DP_rev_i>=min(frame_tmp$time_w))


  h_tmp[4:5,1] <-0 ; h_tmp[c(1,5),2]<-0; h_tmp[c(1,2),3]<-0


  h_tmp[,1:3] <- cumsum(h_tmp[,1:3])

  h_tmp[4:5,1] <-0 ; h_tmp[c(1,5),2]<-0; h_tmp[c(1,2),3]<-0

  h_tmp2 <- h_tmp %>%  left_join(frame_tmp2, by=c("DP_rev_i"="time_w"))

  h_tmp2[,1:3] <- h_tmp2[,1:3] * h_tmp2[,5:7]

  hazard2<- sum(h_tmp2[,1:3], na.rm=T)/sum(h_tmp2[,5:7],na.rm=T)

  return(hazard2)
}










