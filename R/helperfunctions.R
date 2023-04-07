#' @importFrom fastDummies dummy_cols
#' @importFrom bshazard bshazard
#' @importFrom reshape2 melt
#' @import survival
#' @import forecast
#' @import reticulate
#' @import xgboost

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
                                         select_columns,
                                         remove_first_dummy = FALSE){
  "
  This function encodes the matrices that we need for model fitting.

  "

  #individual_data$training.data
  X <- data %>%
    dummy_cols(select_columns = select_columns, #individual_data$categorical_features
               remove_selected_columns = T,
               remove_first_dummy = remove_first_dummy)

  tmp.cond=as.logical(apply(pkg.env$vgrepl(pattern=select_columns,
                                           x=colnames(X)), #individual_data$categorical_features
                            MARGIN=1,
                            sum))

  X <- X %>%
    select(colnames(X)[tmp.cond] ) %>%
    as.data.frame()

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

pkg.env$scaler <- function(continuous_features_scaling_method){
  ""
  if(continuous_features_scaling_method == "minmax" ){return(pkg.env$MinMaxScaler)}


}


pkg.env$deep_surv_pp <- function(X,
                    Y,
                    training_test_split){

  # data_transformed <- cbind(X, Y)

  id_train <- sample(c(TRUE,FALSE), nrow(X), replace=T, prob= c(training_test_split,1-training_test_split) )


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
    validation_data = validation_data
  ))

}


# Fitting routines

pkg.env$fit_cox_model <- function(data,
                          formula_ct,
                          newdata){
  "This function is the fitting routine for the cox model."

  cox <- coxph(formula_ct, data=data)
  cox_lp <- predict(cox,newdata=newdata,'lp',reference='zero')

  # beta_2 <- cox$coef
  # beta<-c(0,beta_2)


  # Xb <- as.matrix(X)%*%beta

  # X_ams <- cbind(X_ams, Xb)

  # beta_ams = unique(round(X_ams,10) )[,ncol(X_ams)] #if no round some systems has too high precision.


  list(
    cox=cox,
    cox_lp=cox_lp,
    expg = exp(cox_lp)
    # beta=beta,
    # beta_2=beta_2,
    # Xb=Xb,
    # X_ams=X_ams,
    # beta_ams=beta_ams
  )
}



pkg.env$fit_deep_surv <- function(data,
                                  hparameters,
                                  network_structure=NULL,
                                  newdata){


  # #Import python modules

  torchtuples <- reticulate::import("torchtuples")
  torch <- reticulate::import("torch")

  #Source python code for left truncated deepsurv
  # reticulate::source_python(".\\inst\\python\\coxnetwork_custom.py")
  reticulate::source_python(system.file("python", "coxnetwork_custom.py", package = "ReSurv"))
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
  # batch_size = as.integer(hparameters$batch_size)
  #
  # epochs = as.integer(hparameters$epochs)


  # Setup CoxPH model, as imported from python script.
  model <- CoxPH(
    net = net,
    optimizer = torchtuples$optim$Adam(lr=hparameters$lr),
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
    verbose = hparameters$verbose,
    val_data=data$validation_data,
    val_batch_size=hparameters$batch_size,
    num_workers=hparameters$num_workers
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

    # bs_hazard2 = tibble(DP_rev_i = bs_hazard$time,
    #                     hazard=bs_hazard$hazard) %>%
    #   mutate(hazard = hazard-lag(hazard, default=0))
  }

  if(baseline == "spline"){
    bs_hazard=bshazard(pkg.env$formula.editor(continuous_features=NULL,
                                                categorical_features="1",
                                                continuous_features_spline=F),
                                 data=data[(data$AP_i-1)%%(conversion_factor^-1)==0 & data$claim_type==0,],
                                 nk=nk,
                                 nbin=nbin,
                                 phi=phi)
    bs_hazard <- tibble(time = bs_hazard$time,
                         hazard = bs_hazard$hazard)
  }

  return(list(bs_hazard=bs_hazard))

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

# pkg.env$dissect_hazard_name <- function(names_hazard, name = "AP"){
#
#   "
#   Get accident period and covariates to use for later grouping.
#
#   "
#
#   names_hazard <- as.character(names_hazard)
#   start_position <- gregexpr('_', names_hazard)[[1]][2] +1
#   end_position <- gregexpr(',', names_hazard)[[1]][1] -1
#
#   AP <- substr(names_hazard, start_position, end_position )
#
#   covariate <- substr(names_hazard, end_position+2, nchar(names_hazard) )
#   if(name == "AP"){
#     return(
#       as.numeric(AP)
#     )
#   }
#   if(name=="covariate"){
#     return(
#       covariate
#     )
#   }
# }

pkg.env$hazard_data_frame <- function(hazard,
                                      eta_old=1/2,
                                      categorical_features,
                                      continuous_features){

  "
  Convert hazard matrix to dataframe and add grouping variables.

  "

  #Need special handling if we ahve continuous variables
  if( (length(continuous_features)==1 & "AP_i" %in% continuous_features) | is.null(continuous_features)){
    #Make sure continuous is null, only applied when AP_i is only continuous feature
    continuous_features_group=NULL

    #Calculate input development factors and corresponding survival probabilities
    hazard_frame_tmp <- hazard %>%
      mutate(dev_f_i = (1+(1-eta_old)*hazard)/(1-eta_old*hazard) ) %>% #Follows from the assumption that claims are distributed evenly in the input period
      mutate(dev_f_i = ifelse(dev_f_i<0,1,dev_f_i)) %>%  #for initial development factor one can encounter negative values, we put to 0
      group_by(!!sym(categorical_features), AP_i) %>%
      arrange(DP_rev_i) %>%
      mutate(cum_dev_f_i = cumprod(dev_f_i)) %>%
      mutate(S_i = ifelse(cum_dev_f_i==0,0,1/cum_dev_f_i), # to handle the ifelse statement from above
             S_i_lead = lead(S_i, default = 0),
             S_i_lag = lag(S_i, default = 1)) %>%
      select(-c(expg, baseline, hazard))

     #we need the lead and lag values for later calcualtion of expected amounts, store in the dataset.
     hazard_frame <- hazard %>%
       left_join(hazard_frame_tmp, c(categorical_features,
                              "AP_i",
                              "DP_rev_i" )) %>%
       mutate(dev_f_i = coalesce(dev_f_i,1),
           S_i = coalesce(S_i,1),
           S_i_lead = coalesce(S_i_lead,1),
           S_i_lag = coalesce(S_i_lag, 1),
           cum_dev_f_i = coalesce(cum_dev_f_i,1))

    } else {
      #Equivalent to above expect we now also join by continuous features
      continuous_features_group=continuous_features[!("AP_i" %in% continuous_features)]
      hazard_frame_tmp <- hazard %>%
        mutate(dev_f_i = (1+(1-eta_old)*hazard)/(1-eta_old*hazard) ) %>%
        mutate(dev_f_i = ifelse(dev_f_i<0,0,dev_f_i)) %>%  #for initial development factor one can encounter negative values, we put to 0
        group_by(!!sym(categorical_features), AP_i) %>%
        arrange(DP_rev_i) %>%
        mutate(cum_dev_f_i = cumprod(dev_f_i)) %>%
        mutate(S_i = ifelse(cum_dev_f_i==0,0,1/cum_dev_f_i),
               S_i_lead = lead(S_i, default = 0),
               S_i_lag = lag(S_i, default = 1)) %>%
      select(-c(expg, baseline, hazard))

      hazard_frame <- hazard %>%
        left_join(hazard_frame_tmp, c(categorical_features,
                                      continuous_features_group,
                                      "AP_i",
                                      "DP_rev_i")) %>%
        mutate(dev_f_i = coalesce(dev_f_i,1),
               S_i = coalesce(S_i,1),
               S_i_lead = coalesce(S_i_lead,1),
               S_i_lag = coalesce(S_i_lag, 1),
               cum_dev_f_i = coalesce(cum_dev_f_i,1))
    }
  return(hazard_frame)
}

pkg.env$covariate_mapping <- function(hazard_frame,
                                      categorical_features,
                                      continuous_features,
                                      conversion_factor)
  {
  "
  Create a dimension table, that holds a link between inputted categorical features and the group, that is used for expected_values
  "

  #Need to handle Accident period effect seperatly
  if( (length(continuous_features)==1 & "AP_i" %in% continuous_features)){
    continuous_features_group = NULL
    }
  else{
    continuous_features_group=continuous_features[!("AP_i" %in% continuous_features)]
  }

  ## The next steps generate a grouping key, used for aggregating from input periods to output periods
  hazard_frame$covariate <- pkg.env$name_covariates(
    hazard_frame,
    categorical_features,
    continuous_features_group
  )

  #Group is pr. covariate, output accident period
  if("AP_i" %in% continuous_features){

    groups <- unique(data.frame(AP_i = hazard_frame$AP_i, covariate = hazard_frame$covariate)) %>%
      mutate(group_i = row_number())

    hazard_group <- hazard_frame %>%  left_join(groups, by=c("AP_i", "covariate"))

  }
  else{
    groups <- unique(data.frame(covariate = hazard_frame$covariate)) %>%
      mutate(group_i = row_number())

    expand.grid()

    hazard_group <- hazard_frame %>%  left_join(groups, by=c("covariate"))
  }

  #If we have to group for later output, add the relevant groups as well
  groups$group_o <- groups$group_i
  # The only time the groups will be different, is when we are including accident period as a covariate
  if(conversion_factor != 1 & "AP_i" %in% continuous_features){

      groups_o <- unique(data.frame(AP_o = ceiling(hazard_group$AP_i*conversion_factor), covariate = hazard_group$covariate)) %>%
        mutate(group_o = row_number())

      groups <- groups %>% select(-group_o) %>%
        left_join(groups_o, by=c("AP_o", "covariate"))

  }
  return(
    list(hazard_group=hazard_group, groups = groups)
  )



}


pkg.env$latest_observed_values_i <- function(data, groups,  categorical_features, continuous_features){
  "
  Retrieve total amount of observed claims

  "

  data_reserve <- data

  trunc = max(data_reserve$DP_i)

  max_observed_ap_dp <- data_reserve %>%
    group_by(AP_i) %>%
    summarize(max_DP_i = max(DP_i), .groups="drop")

  #create grid to hold observed values for all possible times (also where we have no observations)
  observed_grid <- expand.grid(AP_i = unique(data_reserve$AP_i),
                               DP_rev_i = unique(data_reserve$DP_rev_i),
                               group_i = groups$group_i ) %>%
    mutate(DP_i = trunc-DP_rev_i+1) %>%
    left_join(max_observed_ap_dp, by = "AP_i") %>%
    filter(DP_i <= max_DP_i) %>%
    select(-c(max_DP_i))

  #Max possible development time per accident period
  max_DP_i <- data_reserve %>% group_by(AP_i) %>%
    summarise(DP_max_rev =min(max(DP_rev_i)-DP_i)+1 ) %>%
    distinct()

  data_reserve2 <- data_reserve %>%
    select(AP_i, AP_o, DP_rev_i, DP_i, all_of(categorical_features), all_of(continuous_features), I) %>%
    mutate(AP_i = as.numeric(AP_i)) %>%
    left_join(max_DP_i, by="AP_i")

  #The reason for the if statement is due to the !!sym logic, because !!sym(NULL) is not valid
  if(is.null(continuous_features)){ #length(continuous_features) == 1 & "AP_i" %in% continuous_features
   observed_so_far <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                 DP_max_rev) %>%
    summarise(latest_I=sum(I), .groups = "drop")

  observed_dp_rev_i <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                 DP_rev_i, DP_i) %>%
    summarise(I=sum(I), .groups = "drop")

  #Combine covariate values into single variable
  observed_so_far$covariate <- pkg.env$name_covariates(
    observed_so_far,
    categorical_features,
    continuous_features
  )
  observed_dp_rev_i$covariate <- pkg.env$name_covariates(
    observed_dp_rev_i,
    categorical_features,
    continuous_features
  )

  # Latest cumulative
  observed_so_far_out <- observed_so_far %>%  left_join(groups, by=c("covariate")) %>%
    select(AP_i, group_i, DP_max_rev, latest_I)

  observed_dp_rev_i_tmp <- observed_dp_rev_i %>%  left_join(groups, by=c("covariate")) %>%
    select(AP_i, group_i, DP_rev_i, DP_i, I)

  #Observed pr. development period
  observed_dp_rev_i_out <- observed_grid %>%
    left_join(observed_dp_rev_i_tmp, by=c("AP_i", "DP_i", "DP_rev_i", "group_i"))



  } else{ #Very similiar to above, expect we now also take the continuous features into account

    if(length(continuous_features)==1 & "AP_i" %in% continuous_features){
      continuous_features<-NULL
    }
    else{
      continuous_features <- continuous_features[!("AP_i" %in% continuous_features)]
    }

    #If-statment due to grouping by continous vair
    if(is.null(continuous_features)){
    observed_so_far <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                   DP_max_rev) %>%
      summarise(latest_I=sum(I), .groups = "drop")
    observed_dp_rev_i <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                     DP_rev_i, DP_i) %>%
      summarise(I=sum(I), .groups = "drop")
    }
    else{
      observed_so_far <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                     !!sym(continuous_features),
                                                     DP_max_rev) %>%
        summarise(latest_I=sum(I), .groups = "drop")
      observed_dp_rev_i <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                       !!sym(continuous_features),
                                                       DP_rev_i, DP_i) %>%
        summarise(I=sum(I), .groups = "drop")
    }

    observed_so_far$covariate <- pkg.env$name_covariates(
      observed_so_far,
      categorical_features,
      continuous_features
    )

    observed_dp_rev_i$covariate <- pkg.env$name_covariates(
      observed_dp_rev_i,
      categorical_features,
      continuous_features
    )

    observed_so_far_out <- observed_so_far %>%  left_join(groups, by=c("AP_i", "covariate")) %>%
      select(AP_i, group_i, DP_max_rev,latest_I )

    observed_dp_rev_i_tmp <- observed_dp_rev_i %>%  left_join(groups, by=c("AP_i", "covariate")) %>%
      select(AP_i, group_i, DP_rev_i, DP_i, I)

    observed_dp_rev_i_out <- observed_grid %>%
      left_join(observed_dp_rev_i_tmp, by=c("AP_i", "DP_i", "DP_rev_i", "group_i")) %>%
      inner_join(groups[,c("AP_i", "group_i")], by =c("AP_i", "group_i")) #filter only relevant combinations

  }

  return(list(latest_cumulative = observed_so_far_out, observed_pr_dp = observed_dp_rev_i_out))

}

pkg.env$name_covariates <- function(data,
                                    categorical_features,
                                    continuous_features){
  "
  Create groups based upon combination of covariates.
  Here we craete the name of the group
  "
  if(is.null(continuous_features)){
    df <- data %>%  select(all_of(categorical_features))
    name_seperate <- suppressMessages(map2_dfc(colnames(df), df, paste, sep = '_'))
    name_combined <- apply( name_seperate , 1 , paste , collapse = ", " )
  return(name_combined)
  }
  else{
    df <- data %>%  select(all_of(categorical_features), all_of(continuous_features))
    name_seperate <- suppressMessages(map2_dfc(colnames(df), df, paste, sep = '_'))
    name_combined <- apply( name_seperate , 1 , paste , collapse = ", " )

  return(name_combined)
  }



}

pkg.env$predict_i <- function(hazard_data_frame,
                              latest_cumulative,
                              grouping_method
){
  "
  Calculate expected incremential claim number on input scale.
  Grouping is used when doing granularity increased development factors in i_to_o_development_factor
   "

  # #select relevant hazard values
  grouped_hazard_0 <- hazard_data_frame %>% #for the last development, if we included group '0', we would be extrapolating for half a parallelogram - doesn't make sense
    left_join(latest_cumulative, by=c("group_i", "AP_i"))

  # Predict expected numbers, this is also used grouping methodology
  expected <-  grouped_hazard_0 %>%
    select(DP_rev_i, AP_i, group_i, S_i, S_i_lag, DP_max_rev, latest_I ) %>%
    mutate(gm = grouping_method) %>%
    left_join(hazard_data_frame %>%
                mutate(DP_rev_i = DP_rev_i +1) %>%
                select(DP_rev_i, AP_i, group_i, S_i) %>%
                rename(S_ultimate_i = S_i), by=c("DP_max_rev"="DP_rev_i",
                                                 "AP_i" = "AP_i",
                                                 "group_i" = "group_i")) %>%
    mutate(U=case_when(
      grouping_method == "probability" ~ 1,
      S_ultimate_i ==0 ~ 0,
      AP_i != 1 ~ 1/S_ultimate_i * latest_I,
      TRUE ~ latest_I)) %>%
    mutate(I_expected = U*(S_i_lag-S_i)) %>%
    mutate(IBNR = ifelse(DP_rev_i < DP_max_rev, I_expected, NA)) %>%
    select(AP_i, group_i, DP_rev_i, I_expected, IBNR)

  return(expected)


}

pkg.env$retrieve_df_i <- function(hazard_data_frame,
                              groups
){
  "
  Return data frame only containing input development factors.
  "

  df_i <- hazard_data_frame %>%
    select(group_i, DP_rev_i, dev_f_i) %>%
    distinct() %>%
    reshape2::dcast(DP_rev_i ~group_i, value.var="dev_f_i") %>%
    select(-DP_rev_i)

  #We only have 5 columns in the case of AP being included as covariate
  if(ncol(groups) == 5){
    colnames(df_i) <- c(paste0("AP_i_",groups$AP_i,",", groups$covariate ))
  }
  else{
    colnames(df_i) <- c(groups$covariate )
  }
  #


  df_i <- as.data.frame(df_i[1:(nrow(df_i)-1),]) %>%
    map_df(rev) %>%
    mutate(DP_i=row_number())

  return(df_i)


}


pkg.env$input_hazard_frame <- function(
    hazard_frame,
    expected_i,
    categorical_features,
    continuous_features,
    df_i,
    groups)
{
  "
  Create a hazard frame with relevant input-period specific values for later output

  "
  if("AP_i" %in% continuous_features & length(continuous_features) == 1){
    continuous_features <- NULL
  }
  else{
    continuous_features <- continuous_features[!("AP_i" %in% continuous_features)]
  }
  hazard_frame_input_relevant <- hazard_frame %>%
    select(- c(cum_dev_f_i, S_i, S_i_lead, S_i_lag, covariate))

  #If AP is included as a grouping variable
  if(ncol(groups)==5){
    df_i_long <- df_i %>%
      reshape2::melt(id.vars="DP_i") %>%
      left_join(groups %>%
                  mutate(covariate = paste0("AP_i_", AP_i, ",", covariate) ) %>%
                           select(c(covariate, group_i)), by=c("variable" = "covariate"))

    colnames(df_i_long) <- c("DP_i", "covariate", "df_i", "group_i")

  }
  else{

    df_i_long <- df_i %>%
      reshape2::melt(id.vars="DP_i") %>%
      left_join(groups[,c("covariate", "group_i")], by=c("variable" = "covariate")) %>%
      mutate(DP_i = DP_i +1) #to get correct DP_i

    colnames(df_i_long) <- c("DP_i", "covariate", "df_i", "group_i")
  }


  max_DP_rev_i = max(expected_i$DP_rev_i)



  hazard_frame_input <- expected_i %>%
    mutate(DP_i = max_DP_rev_i-DP_rev_i +1) %>%
    left_join(hazard_frame_input_relevant, by =c("group_i", "AP_i", "DP_rev_i")) %>%
    left_join(df_i_long[, c("DP_i", "group_i", "df_i")], by = c("DP_i", "group_i")) %>%
    replace_na(list(df_i = 1))

  #Ordering
  hazard_frame_input <- hazard_frame_input[,c(categorical_features,
                                                continuous_features,
                                                "AP_i",
                                                "DP_rev_i",
                                                "expg",
                                                "baseline",
                                                "hazard",
                                                "df_i",
                                                "group_i",
                                                "I_expected",
                                                "IBNR")]
  return(hazard_frame_input)

}


pkg.env$predict_o <- function(
    expected_i,
    groups,
    conversion_factor
){
  "
  Calculate expected incremential claim number on output scale

   "

  # Predict expected numbers, this is also used grouping methodology
  expected <-  expected_i %>%
    left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    mutate(AP_o = ceiling(AP_i*conversion_factor),
           DP_rev_o = ceiling((DP_rev_i-(AP_i-1)%%(1/conversion_factor) ) *conversion_factor)) %>%
    filter(DP_rev_o >0) %>% #since for DP_rev_o = 0, we are working with half a parrallelogram in the end of the development time
    group_by(AP_o, DP_rev_o, group_o) %>%
    summarize(I_expected = sum(I_expected),
              IBNR = sum(IBNR), .groups="drop") %>%
    select(AP_o, group_o, DP_rev_o, I_expected, IBNR)

  return(expected)


}

pkg.env$i_to_o_development_factor <- function(i,
                                  hazard_data_frame,
                                  expected_i,
                                  dp_ranges,
                                  groups,
                                  observed_pr_dp,
                                  latest_cumulative,
                                  conversion_factor,
                                  grouping_method){
  "
  Group input development factor to output.

  "
  # Add output groupings to relevant frames
  hazard_data_frame <- hazard_data_frame %>%
    left_join(groups[,c("group_i", "group_o")], by =c("group_i"))

  observed_pr_dp_o  <- observed_pr_dp %>%
    left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    group_by(AP_i, group_o, DP_rev_i, DP_i) %>%
    summarize(I = sum(I, na.rm=T), .groups = "drop")

  latest_cumulative_o <- latest_cumulative %>%
    left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    group_by(AP_i, group_o, DP_max_rev) %>%
    summarize(latest_I = sum(latest_I, na.rm=T), .groups = "drop")

  #For probability approach to grouping method we assume equal exposure for each accident period
  if(grouping_method == "probability"){
  expected_i <-  pkg.env$predict_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    latest_cumulative = latest_observed$latest_cumulative,
    grouping_method = "probability"
  ) %>%
    left_join(groups[,c("group_i", "group_o")], by =c("group_i"))
  } else{
    expected_i <-  expected_i %>%
      left_join(groups[,c("group_i", "group_o")], by =c("group_i"))

  }


  # #select relevant hazard value group and add output variables, and other variables to help with grouping
   grouped_hazard_0 <- hazard_data_frame %>%
     filter(group_o==i) %>%
     mutate(DP_rev_o = ceiling((DP_rev_i-(AP_i-1)%%(1/conversion_factor) ) *conversion_factor)) %>%
     filter(DP_rev_o > 0) %>%  #for the last development, if we included group '0', we would be extrapolating for half a parallelogram - doesn't make sense
     left_join(dp_ranges, by=c("AP_i", "DP_rev_o")) %>%
     left_join(latest_cumulative_o, by=c("group_o", "AP_i")) %>%
     left_join(observed_pr_dp_o, by=c("group_o", "AP_i", "DP_rev_i"))

   # Create cumulative observed to find exposure for each period
   cumulative_observed <- observed_pr_dp_o %>%
     filter(group_o == i) %>%
     group_by(AP_i) %>%
     arrange(DP_i) %>%
     mutate(exposure = cumsum(ifelse(is.na(I),0,I) )) %>%
     mutate(DP_rev_i = DP_rev_i -1) %>%  #as we want this as exposure we join by the previous development period
     select(AP_i, group_o, DP_rev_i, exposure)

   exposures <- grouped_hazard_0 %>%
     group_by(AP_i, DP_rev_o, group_o) %>%
     filter(DP_rev_i == max(DP_rev_i)) %>%
     left_join(cumulative_observed, by=c("AP_i", "group_o",
                                        "max_dp"="DP_rev_i"))

   #Where we do not have any observed correct exposure we extrapolate based on fitted hazard
   no_exposure <-  exposures %>%
     select(DP_rev_i,  DP_rev_o, AP_i, group_o, S_i, DP_max_rev, latest_I ) %>%
     mutate(gm = grouping_method) %>%
     left_join(hazard_data_frame %>%
                 mutate(DP_rev_i = DP_rev_i +1) %>%
                 select(DP_rev_i, AP_i, group_o, S_i) %>%
                 rename(S_ultimate_i = S_i), by=c("DP_max_rev"="DP_rev_i",
                                                  "AP_i" = "AP_i",
                                                  "group_o" = "group_o")) %>%
     mutate(U=ifelse(
       S_ultimate_i ==0, 0,
       1/S_ultimate_i * latest_I) ) %>%
     mutate(U = ifelse(gm=="probability", 1 ,U)) %>%
     mutate(exposure_expected = U*(S_i)) %>%  #in theory one could say U*S_i- ifelse(DP_max_rev==DP_rev_i-1, latest_I, U*S_i_lead ), but this might lead to negative expected as we are not sure latest equal the same as distribution estimate
     select(AP_i, group_o, DP_rev_o, DP_rev_i, exposure_expected)

    #Take seen exposure if possible, otherwise extrapolated exposure
   exposures_combined <- exposures  %>%
     mutate(gm = grouping_method) %>%
     left_join(no_exposure, by  = c(   "AP_i",
                                    "DP_rev_o",
                                    "DP_rev_i",
                                    "group_o")) %>%
     mutate(exposure_combined = ifelse(gm == "probability",
                                       coalesce(exposure_expected,0),
                                       coalesce(exposure, exposure_expected))
            )

  #Take seen observed if possible otherwise extrapolated observed
   grouped_hazard_1 <- grouped_hazard_0 %>%
     mutate(gm = grouping_method) %>%
     left_join(expected_i, by  = c("AP_i",
                                        "group_o",
                                        "DP_rev_i")) %>%
     mutate(I_combined = ifelse(gm == "probability",
                                coalesce(I_expected,0),
                                coalesce(I, I_expected,0))
            )

   #group to output scale
   grouped_hazard_2 <- grouped_hazard_1 %>%
     group_by(AP_i, DP_rev_o, group_o) %>%
     summarize(observed = sum(I_combined), .groups="drop") %>%
     left_join(exposures_combined, by=c("AP_i", "group_o", "DP_rev_o"))

   output_dev_factor <- grouped_hazard_2 %>%
     group_by(DP_rev_o) %>%
     summarise(dev_f_o = (sum(observed)+  sum(exposure_combined))/sum(exposure_combined) )


   return(output_dev_factor$dev_f_o)

   }

pkg.env$output_hazard_frame <- function(
    hazard_frame_input,
    expected_o,
    categorical_features,
    continuous_features,
    df_o,
    groups)
  {
  "
  Create output hazard frame

  "
  if("AP_i" %in% continuous_features & length(continuous_features) == 1){
    continuous_features <- NULL
  }
  else{
    continuous_features <- continuous_features[!("AP_i" %in% continuous_features)]
  }

  #Relevant variables, we do not include hazard, baseline, expg as they currently only live on input-level
  hazard_frame_input_relevant <- hazard_frame_input %>%
    select(all_of(categorical_features), all_of(continuous_features), group_i) %>%
    left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    select(-c(group_i)) %>%
    distinct()

  #If AP is included as a grouping variable
  if(ncol(groups)==5){
    df_o_long <- df_o %>%
      reshape2::melt(id.vars="DP_o") %>%
      left_join(groups[,c("AP_o","covariate", "group_o")] %>%
                  mutate(covariate = paste0("AP_o_", AP_o, ", ", covariate)), by=c("variable" = "covariate"))

    colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")

  }
  else{

    df_o_long <- df_o %>%
      reshape2::melt(id.vars="DP_o") %>%
      left_join(groups[,c("covariate", "group_o")], by=c("variable" = "covariate")) %>%
      mutate(DP_o = DP_o +1) #to get correct

    colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")
  }


  max_DP_rev_o = max(expected_o$DP_rev_o)

  hazard_frame_output <- expected_o %>%
    mutate(DP_o = max_DP_rev_o-DP_rev_o +1) %>%
    left_join(hazard_frame_input_relevant, by =c("group_o")) %>%
    left_join(df_o_long[, c("DP_o", "group_o", "df_o")], by = c("DP_o", "group_o")) %>%
    replace_na(list(df_o = 1))

  hazard_frame_output <- hazard_frame_output[,c(categorical_features,
                         continuous_features,
                         "AP_o",
                         "DP_rev_o",
                         "df_o",
                         "group_o",
                         "I_expected",
                         "IBNR")]
  return(hazard_frame_output)

}

pkg.env$i_to_o_hazard<- function(i,
                                              hazard_data_frame,
                                              #frame_tmp,
                                              development_periods,
                                              observed_pr_dp,
                                              latest_cumulative,
                                              conversion_factor){
  "
  WIP
  Group input hazard to output.

  "

  " SHOULD be called from the following loop in the ReSurvIndividualData.R

    #group to quarters, this is relatively time consuming,
  #Note: the eta-approximation is not covariat dependent.
  for( i in 1:max_DP){ #Loop through each output period, to find weights
     # frame_tmp <- data.frame(IndividualData$training) %>% filter(TR_o<i) %>% #All claims that haven't been truncated at said reverse development
     #   filter(DP_rev_o>=i) %>% #Claims that still hasn't been reported
     #   mutate(time_w = round(ifelse(DP_rev_o==i, DP_rev_i, AP_i-1+1/(IndividualData$conversion_factor)*(i-AP_o+1) ),10) ) %>% #If a claim is reported in the corresponding development period save said reporting time, otherwise we need the corresponding limit for each acciedent period in the development period.
     #   #mutate(weight = ifelse(DP_rev_o==i, (DP_rev_i-1)%%(1/(IndividualData$conversion_factor))+1, 1/(IndividualData$conversion_factor) )) %>% #If reported in said period give weight corresponding to amount of input_time period spend in output time_period, otherwise give width length of output as weight.
     #   #mutate(weight_eta = ifelse(DP_rev_o==i, (DP_rev_i-1)%%(1/(IndividualData$conversion_factor))+1, 0 )) %>% #If reported in said period give weight corresponding to amount of input_time period spend in output time_period, otherwise give width length of output as weight.
     #   mutate(observed = ifelse(DP_rev_o==i, 1, 0 )) %>%
     #   mutate(exposure = ifelse(DP_rev_o==i, 0, 1 )) %>%
     #   mutate(p_month = (AP_i-1)%%(1/(IndividualData$conversion_factor))+1) %>% #Entering month in development period
     #   group_by(p_month, AP_i, time_w) %>%
     #   dplyr::summarise(observed=sum(observed),
     #                    exposure = sum(exposure), .groups='drop )


  Grouoped hazard called from something of the like:
      development_factor_o[i,] <- mapply(pkg.env$i_to_o_development_factor,
                           1:max(hazard_data_frame$groups$group),
                           MoreArgs=list(hazard_data_frame=hazard_data_frame$hazard_group,
                                         development_periods = development_periods,
                                         observed_pr_dp = latest_observed$observed_pr_dp,
                                         latest_cumulative = latest_observed$latest_cumulative,
                                         conversion_factor = IndividualData$conversion_factor))


    #hazard_o[i,] <- unlist(grouped_hazard[1,])
    #eta_o[i] <- unlist(grouped_hazard[2,1]) #since not accident-period dependent, eta is the same for every output period

  }




  "


 #Get correct values and convert to ultimates before calculating weights.
  input_period_weights <- grouped_hazard_0 %>%
    group_by(AP_i, DP_max_rev, group ) %>%
    summarize(I=mean(I), .groups = "drop") %>% #we have duplicates pr. AP_i but they will always be same pr. AP_i
    left_join(hazard_data_frame[,c("AP_i", "group", "DP_rev_i", "S_i")], by=c("AP_i" = "AP_i",
                                                                              "group" = "group",
                                                                     "DP_max_rev" = "DP_rev_i")) %>%
    mutate(U=I*1/S_i) %>%
    mutate(i_period_weight = U/sum(U))

  #Join the weights, which also give an indication of which development periods, pr. accident period is relevant in the parallelogram
   grouped_hazard_1 <- grouped_hazard_0 %>%   inner_join(frame_tmp, by = c("p_month" = "p_month", "DP_rev_i" = "time_w", "AP_i" = "AP_i")) %>%
    mutate(spend = (DP_rev_i-AP_i)%%(1/conversion_factor)+1) %>% #assuming average reporting times of 0.5 pr. development period in
    #filter(!is.na(weight)) %>%
    group_by(AP_i) %>%
    mutate(#value_cum = cumsum(value),
           #weight_new = ifelse(spend==1/conversion_factor,1/2*(S_i - S_i_lead)+(S_i_lead),1/2*(S_i - S_i_lead)),
           weight_new_2 = 1/2*(S_i + S_i_lead)) %>% #probability weighted approach, for each parellelogram we assume reporting half way through
    #arrange(AP_i, desc(DP_rev_i)) %>%  #here one weights by cumulative remaining hazard in parallelogram
    #mutate(value_cum_rev = cumsum(value)) %>%
    ungroup()

 #If we have no rows in grouped_hazard_1 it is because we are in the lower triangle, here we allow for taking the average of the observed during other periods, or the probability approach
  if(nrow(grouped_hazard_1 == 0)){
    frame_tmp_2 <- frame_tmp %>%  group_by(p_month, time_w) %>%
      dplyr::summarise(weight_eta=sum(weight_eta), .groups="drop" )

    grouped_hazard_1 <- grouped_hazard_0 %>%   inner_join(frame_tmp_2, by = c("p_month" = "p_month", "DP_rev_i" = "time_w")) %>%
      mutate(spend = (DP_rev_i-AP_i)%%(1/conversion_factor)+1) %>% #assuming average reporting times of 0.5 pr. development period in
      group_by(AP_i) %>%
      mutate(#value_cum = cumsum(value),
        #weight_new = ifelse(spend==1/conversion_factor,1/2*(S_i - S_i_lead)+(S_i_lead),1/2*(S_i - S_i_lead)),
        weight_new_2 = 1/2*(S_i + S_i_lead)) %>% #probability weighted approach, for each parellelogram we assume reporting half way through
      #arrange(AP_i, desc(DP_rev_i)) %>%  #here one weights by cumulative remaining hazard in parallelogram
      #mutate(value_cum_rev = cumsum(value)) %>%
      ungroup()

  }




    #hazard2 <- sum(h_tmp_3$value_cum*h_tmp_3$weight) / sum(h_tmp_3$weight) #weighing by observed claims in each area in the parallelogram
    #hazard2_2 <- sum(h_tmp_3$value_cum*(1-h_tmp_3$value_cum_s_i)) / sum(1-h_tmp_3$value_cum_s_i)
    hazard_o <- sum(sum(grouped_hazard_1$value*grouped_hazard_1$weight_new_2) / sum(grouped_hazard_1$weight_new_2)) #weighting by probability
    eta_o <- 1/ hazard_o - sum(grouped_hazard_1$exposure)/sum(grouped_hazard_1$observed)

    #(1+(1-eta_o)*hazard_o)/(1-eta_o*hazard_o)

    return(list(hazard_o, eta_o) )


}


pkg.env$spline_hp <- function(hparameters,IndividualData){
  "
  Returns spline hyperparameters in case they are not provided from the user.

  "

  tmp <- list()

  tmp$nk <- ifelse(is.null(hparameters$nk),nrow(IndividualData$training.data)/4,hparameters$nk)

  tmp$nbin <- ifelse(is.null(hparameters$nbin),NULL,hparameters$nbin)

  tmp$phi <- ifelse(is.null(hparameters$phi),NULL,hparameters$phi)

  return(tmp)

}


pkg.env$create.df.2.fcst <- function(IndividualData){

  l1 <- lapply(IndividualData$training.data %>% select(IndividualData$categorical_features), levels)
  l2 <- lapply(IndividualData$training.data %>% select(IndividualData$continuous_features), unique)
  l3 <- list()
  l4 <- list()

  if(!("AP_i"%in%c(IndividualData$categorical_features,IndividualData$continuous_features))){
    l3$AP_i <- unique(IndividualData$full.data[,'AP_i'])
  }else{
    l3 <- NULL
  }

  l4$DP_rev_i <- min(IndividualData$training.data[,'DP_i']):max(IndividualData$training.data[,'DP_i'])

  cross_df(c(l1,l2,l3,l4)) %>%
    as.data.frame()

}



pkg.env$df.2.fcst.nn.pp <- function(data,
                                    newdata,
                                    continuous_features,
                                    categorical_features){

  tmp <- newdata[continuous_features]

  for(cft in continuous_features){

    mnv <- min(data[cft])
    mxv <- max(data[cft])

    tmp[,cft] <-(tmp[,cft]-mnv)/(mxv-mnv)

  }

  Xc=as.matrix.data.frame(tmp)

  X=pkg.env$model.matrix.creator(data= newdata,
                                 select_columns = categorical_features)

  return(cbind(X,Xc))

}


pkg.env$df.2.fcst.xgboost.pp <- function(data,
                                    newdata,
                                    continuous_features,
                                    categorical_features){
  tmp <- Xc <- NULL

  if(!is.null(continuous_features)){
    tmp <- newdata[continuous_features]

    for(cft in continuous_features){

      mnv <- min(data[cft])
      mxv <- max(data[cft])

      tmp[,cft] <-(tmp[,cft]-mnv)/(mxv-mnv)

    }
    Xc=as.matrix.data.frame(tmp)

    }



  X=pkg.env$model.matrix.creator(data= newdata,
                                 select_columns = categorical_features,
                                 remove_first_dummy = T)

  if(!is.null(Xc)){X <- cbind(X,Xc)}
  # browser()
  ds_train_fcst <- xgboost::xgb.DMatrix(as.matrix.data.frame(X), label=rep(1, dim(X)[1]))

  return(ds_train_fcst)

}

## xgboost ----

pkg.env$xgboost_pp <-function(X,Y, training_test_split){

  tmp=cbind(X,Y) %>%
    arrange(DP_rev_i) %>%
    as.data.frame()

  tmp[,'id'] = seq(1,dim(tmp)[1])

  samples_cn <- tmp %>% select(id) %>% sample_frac(size=training_test_split)

  tmp_train <- tmp %>%
    semi_join(samples_cn)%>%
    arrange(DP_rev_i) %>%
    group_by(DP_rev_i) %>%
    mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame()

  tmp_test <- tmp %>%
    anti_join(samples_cn)%>%
    arrange(DP_rev_i) %>%
    group_by(DP_rev_i) %>%
    mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame()


  # ds_all_m <- xgboost::xgb.DMatrix( as.matrix(tmp,ncol=1),
  #                          label=tmp$I)
  ds_train_m <- xgboost::xgb.DMatrix( as.matrix.data.frame(tmp_train %>% select(colnames(X))), label=tmp_train$I)
  ds_test_m <- xgboost::xgb.DMatrix( as.matrix.data.frame(tmp_test %>% select(colnames(X)), label=tmp_test$I))

  attr(ds_train_m, 'truncation') <- tmp_train$TR_i
  attr(ds_test_m, 'truncation') <- tmp_test$TR_i
  attr(ds_train_m, 'claim_arrival') <- tmp_train$DP_rev_i
  attr(ds_test_m, 'claim_arrival') <- tmp_test$DP_rev_i


  attr(ds_train_m, 'risk_sets') <- risks_in_the_tie(starts_i=tmp_train$TR_i,
                                                    stops_i=tmp_train$DP_rev_i,
                                                    stops = unique(tmp_train$DP_rev_i))

  attr(ds_test_m, 'risk_sets') <- risks_in_the_tie(starts_i=tmp_test$TR_i,
                                                   stops_i=tmp_test$DP_rev_i,
                                                   stops = unique(tmp_test$DP_rev_i))

  #

  attr(ds_train_m, 'event_sets') <- events_in_the_tie(starts_i=tmp_train$TR_i,
                                                      stops_i=tmp_train$DP_rev_i,
                                                      stops = unique(tmp_train$DP_rev_i))
  attr(ds_test_m, 'event_sets') <- events_in_the_tie(starts_i=tmp_test$TR_i,
                                                     stops_i=tmp_test$DP_rev_i,
                                                     stops = unique(tmp_test$DP_rev_i))

  #
  attr(ds_train_m, 'efron_c') <- tmp_train$efron_c
  attr(ds_test_m, 'efron_c') <- tmp_test$efron_c

  attr(ds_train_m, 'tieid') <- unname(table(tmp_train$DP_rev_i))
  attr(ds_test_m, 'tieid') <- unname(table(tmp_test$DP_rev_i))


  attr(ds_train_m, 'groups') <- rep( as.integer(names(table(tmp_train$end_time))),
                                     attr(ds_train_m, 'tieid'))
  attr(ds_test_m, 'groups') <- rep( as.integer(names(table(tmp_test$end_time))),
                                    attr(ds_test_m, 'tieid'))

  return(list(ds_train_m=ds_train_m,
              ds_test_m=ds_test_m))

}


pkg.env$fit_xgboost <- function(datads_pp,
                                hparameters){


  out <- xgboost::xgb.train(params = hparameters$params,
                            data =datads_pp$ds_train_m,
                        obj=cox_loss_objective2,
                        nrounds = hparameters$nrounds,
                        feval= cox_evaluation_metrix,
                        watchlist = list(train=datads_pp$ds_train_m,
                                         eval=datads_pp$ds_test_m),
                        print_every_n = hparameters$print_every_n,
                        early_stopping_rounds = hparameters$early_stopping_rounds,
                        maximize = F)

  return(out)



}

pkg.env$baseline.xgboost <- function(model.out, datads_pp){

  risk_sets <- attr(datads_pp$ds_train_m, 'risk_sets')
  preds <- predict(model.out, datads_pp$ds_train_m)
  tieid<- attr(datads_pp$ds_train_m, 'tieid')
  exp_p_sum <- sapply(risk_sets,FUN=exp_sum_computer, ypred=preds)

  tieid/(exp_p_sum)

}






