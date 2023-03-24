#' @importFrom fastDummies dummy_cols
#' @importFrom bshazard bshazard
#' @importFrom reshape2 melt
#' @import survival
#' @import forecast

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
    dummy_cols(select_columns = select_columns, #individual_data$categorical_features
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
    bs_hazard <- tibble(DP_rev_i = bs_hazard$time,
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

pkg.env$dissect_hazard_name <- function(names_hazard, name = "AP"){

  "
  Get accident period and covariates to use for later grouping.

  "

  names_hazard <- as.character(names_hazard)
  start_position <- gregexpr('_', names_hazard)[[1]][2] +1
  end_position <- gregexpr(',', names_hazard)[[1]][1] -1

  AP <- substr(names_hazard, start_position, end_position )

  covariate <- substr(names_hazard, end_position+2, nchar(names_hazard) )
  if(name == "AP"){
    return(
      as.numeric(AP)
    )
  }
  if(name=="covariate"){
    return(
      covariate
    )
  }
}

pkg.env$hazard_data_frame <- function(hazard, conversion_factor,
                                      eta_old=1/2,
                                      continuous_features){

  "
  Convert hazard matrix to dataframe and add grouping variables.

  "
  dev_f <- (1+(1-eta_old)*hazard)/(1-eta_old*hazard)

  dev_f <- dev_f[1:(nrow(dev_f)-1),]

  dev_f <- rbind(dev_f, c(rep(1, ncol(dev_f))))

  dev_f <- dev_f[c(nrow(dev_f),1:(nrow(dev_f)-1)),]

  cum_dev_factor <- cumprod(dev_f)

  s_curve <- 1/cum_dev_factor

  hazard <- as.data.frame(hazard)

  hazard$DP_rev_i <- 1:nrow(hazard)
  s_curve$DP_rev_i <- 1:nrow(hazard)

  hazard_t <- reshape2::melt(hazard, id.vars = "DP_rev_i")
  cum_dev_factor_t <- reshape2::melt(s_curve, id.vars = "DP_rev_i") %>%
    group_by(variable) %>%
    mutate(value = ifelse(value <0,0,value)) %>%
    mutate(S_i_lead = lead(value, default = 0)) #for the initial development the factor doesn't make sense, and could become negatie, we set to 1

  hazard_t$AP_i <- sapply(hazard_t$variable, pkg.env$dissect_hazard_name)
  hazard_t$covariate <- mapply(pkg.env$dissect_hazard_name, hazard_t$variable,  MoreArgs = list(name= "covariate"))
  hazard_t$AP_o  <- ceiling(hazard_t$AP_i*conversion_factor)

  hazard_t$S_i = cum_dev_factor_t$value
  hazard_t$S_i_lead = cum_dev_factor_t$S_i_lead

  #Group is pr. covariate, output accident period
  if("AP_i" %in% continuous_features){
  groups <- unique(data.frame(AP_o = hazard_t$AP_o, covariate = hazard_t$covariate)) %>%
    mutate(group = row_number())
  hazard_group <- hazard_t %>%  left_join(groups, by=c("AP_o", "covariate")) %>%
    select(DP_rev_i, AP_i, value, group, S_i, S_i_lead)

  }
  else{
    groups <- unique(data.frame(covariate = hazard_t$covariate)) %>%
      mutate(group = row_number())
    hazard_group <- hazard_t %>%  left_join(groups, by=c("covariate")) %>%
      select(DP_rev_i, AP_i, value, group, S_i, S_i_lead)
  }


  return(
    list(hazard_group=hazard_group, groups = groups)
  )
}

pkg.env$latest_observed_values <- function(data, groups, time_scale = "i", categorical_features, continuous_features){
  "
  Retrieve total amount of observed claims

  "

  data_reserve <- data %>%
    filter(DP_rev_i>TR_i)

  if(time_scale == "i"){

  max_DP_i <- data_reserve %>% group_by(AP_i) %>%
    summarise(DP_max_rev =min(max(DP_rev_i)-DP_i)+1 ) %>%
    distinct()

  data_reserve2 <- data_reserve %>%
    select(AP_i, AP_o, DP_rev_i, DP_i, categorical_features, continuous_features, I) %>%
    mutate(AP_i = as.numeric(AP_i)) %>%
    left_join(max_DP_i, by="AP_i")

  if(is.null(continuous_features) | length(continuous_features) == 1 & "AP_i" %in% continuous_features){
  observed_so_far <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                 DP_max_rev) %>%
    summarise(latest_I=sum(I), .groups = "drop")

  observed_dp_rev_i <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                 DP_rev_i, DP_i) %>%
    summarise(I=sum(I), .groups = "drop")


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

  observed_so_far_out <- observed_so_far %>%  left_join(groups, by=c("covariate")) %>%
    select(AP_i, group, DP_max_rev, latest_I)

  observed_dp_rev_i_out <- observed_dp_rev_i %>%  left_join(groups, by=c("covariate")) %>%
    select(AP_i, group, DP_rev_i, DP_i, I)



  } else{

    continuous_features <- continuous_features[!("AP_i" %in% continuous_features)]

    observed_so_far <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                   !!sym(categorical_features),
                                                   DP_max_rev) %>%
      summarise(latest_I=sum(I), .groups = "drop")
    observed_dp_rev_i <- data_reserve2 %>%  group_by(AP_i, AP_o, !!sym(categorical_features),
                                                     !!sym(categorical_features),
                                                     DP_rev_i, DP_i) %>%
      summarise(I=sum(I), .groups = "drop")

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

    observed_so_far_out <- observed_so_far %>%  left_join(groups, by=c("AP_o", "covariate")) %>%
      select(AP_i, group, DP_max_rev,latest_I )

    observed_dp_rev_i_out <- observed_dp_rev_i %>%  left_join(groups, by=c("AP_o", "covariate")) %>%
      select(AP_i, group, DP_rev_i, DP_i, I)


  }

  }
  else{

  ' implement for continuous feateus'
  }

  return(list(latest_cumulative = observed_so_far_out, observed_pr_dp = observed_dp_rev_i_out))

}

pkg.env$name_covariates <- function(data,
                                    categorical_features,
                                    continuous_features){
  if(is.null(continuous_features)){
    df <- data[,categorical_features]
    name_seperate <- suppressMessages(map2_dfc(colnames(df), df, paste, sep = '_'))
    name_combined <- apply( name_seperate , 1 , paste , collapse = ", " )
  return(name_combined)
  }
  else{
    df <- observed_so_far[,c(IndividualDatacategorical_features,continuous_features)]
    name_seperate <- suppressMessages(map2_dfc(colnames(df), df, paste, sep = '_'))
    name_combined <- apply( name_seperate , 1 , paste , collapse = ", " )

  return(name_combined)
  }


}

pkg.env$i_to_o_development_factor <- function(i,
                                  hazard_data_frame,
                                  #frame_tmp,
                                  development_periods,
                                  observed_pr_dp,
                                  latest_cumulative,
                                  conversion_factor){
  "
  now: Group input development factor to output.
  OLD: Group input hazard to output.

  "


  # #select relevant hazard values
   grouped_hazard_0 <- hazard_data_frame %>%  filter(group==i) %>%
      left_join(development_periods, by=c("AP_i")) %>%
     filter(DP_rev_i >= min_dp,
            DP_rev_i<=max_dp) %>%
     left_join(latest_cumulative, by=c("group", "AP_i")) %>%
     left_join(observed_pr_dp, by=c("group", "AP_i", "DP_rev_i"))

   #Where we do not have any observed we extrapolate based on fitted hazard
   no_observations <-  grouped_hazard_0 %>%
     filter(is.na(I)) %>%
     select(DP_rev_i, AP_i, group, S_i, S_i_lead, DP_max_rev, latest_I ) %>%
     left_join(hazard_data_frame %>%
                 select(DP_rev_i, AP_i, group, S_i) %>%
                 rename(S_ultimate_i = S_i), by=c("DP_max_rev"="DP_rev_i",
                                       "AP_i" = "AP_i",
                                       "group" = "group")) %>%
     mutate(U=1/S_ultimate_i * latest_I) %>%
     mutate(I_expected = U*(S_i-S_i_lead)) %>%  #in theory one could say U*S_i- ifelse(DP_max_rev==DP_rev_i-1, latest_I, U*S_i_lead ), but this might lead to negative expected as we are not sure latest equal the same as distribution estimate
     select(AP_i, group, DP_rev_i, I_expected)

   cumulative_observed <- observed_pr_dp %>%
     group_by(AP_i, group) %>%
     arrange(DP_i) %>%
     mutate(exposure = cumsum(I)) %>%
     mutate(DP_rev_i = DP_rev_i -1) %>%  #as we want this as exposure we join by the previous development period
     select(AP_i, group, DP_rev_i, exposure)

   exposures <- grouped_hazard_0 %>%
     group_by(AP_i, group) %>%
     filter(DP_rev_i == max(DP_rev_i)) %>%
     left_join(cumulative_observed, by=c("AP_i", "group", "DP_rev_i"))

   #Where we do not have any observed correct exposure we extrapolate based on fitted hazard
   no_exposure <-  exposures %>%
     filter(is.na(exposure)) %>%
     select(DP_rev_i, AP_i, group, S_i, S_i_lead, DP_max_rev, latest_I ) %>%
     left_join(hazard_data_frame %>%
                 select(DP_rev_i, AP_i, group, S_i) %>%
                 rename(S_ultimate_i = S_i), by=c("DP_max_rev"="DP_rev_i",
                                                  "AP_i" = "AP_i",
                                                  "group" = "group")) %>%
     mutate(U=1/S_ultimate_i * latest_I) %>%
     mutate(exposure_expected = U*(S_i_lead)) %>%  #in theory one could say U*S_i- ifelse(DP_max_rev==DP_rev_i-1, latest_I, U*S_i_lead ), but this might lead to negative expected as we are not sure latest equal the same as distribution estimate
     select(AP_i, group, exposure_expected)


   exposures_combined <- exposures  %>%
     left_join(no_exposure, by  = c("AP_i",
                                        "group")) %>%
     mutate(exposure_combined = coalesce(exposure, exposure_expected))


   grouped_hazard_1 <- grouped_hazard_0 %>%
     left_join(no_observations, by  = c("AP_i",
                                        "group",
                                        "DP_rev_i")) %>%
     mutate(I_combined = coalesce(I, I_expected))

   grouped_hazard_2 <- grouped_hazard_1 %>%
     group_by(AP_i, group) %>%
     summarize(observed = sum(I_combined), .groups="drop") %>%
     left_join(exposures_combined, by=c("AP_i", "group"))


   # OLD #Get correct values and convert to ultimates before calculating weights.
   # input_period_weights <- grouped_hazard_0 %>%
   #   group_by(AP_i, DP_max_rev, group ) %>%
   #   summarize(I=mean(I), .groups = "drop") %>% #we have duplicates pr. AP_i but they will always be same pr. AP_i
   #   left_join(hazard_data_frame[,c("AP_i", "group", "DP_rev_i", "S_i")], by=c("AP_i" = "AP_i",
   #                                                                             "group" = "group",
   #                                                                    "DP_max_rev" = "DP_rev_i")) %>%
   #   mutate(U=I*1/S_i) %>%
   #   mutate(i_period_weight = U/sum(U))

  # OLD: #Join the weights, which also give an indication of which development periods, pr. accident period is relevant in the parallelogram
  #  grouped_hazard_1 <- grouped_hazard_0 %>%   inner_join(frame_tmp, by = c("p_month" = "p_month", "DP_rev_i" = "time_w", "AP_i" = "AP_i")) %>%
  #   mutate(spend = (DP_rev_i-AP_i)%%(1/conversion_factor)+1) %>% #assuming average reporting times of 0.5 pr. development period in
  #   #filter(!is.na(weight)) %>%
  #   group_by(AP_i) %>%
  #   mutate(#value_cum = cumsum(value),
  #          #weight_new = ifelse(spend==1/conversion_factor,1/2*(S_i - S_i_lead)+(S_i_lead),1/2*(S_i - S_i_lead)),
  #          weight_new_2 = 1/2*(S_i + S_i_lead)) %>% #probability weighted approach, for each parellelogram we assume reporting half way through
  #   #arrange(AP_i, desc(DP_rev_i)) %>%  #here one weights by cumulative remaining hazard in parallelogram
  #   #mutate(value_cum_rev = cumsum(value)) %>%
  #   ungroup()

  #  old #If we have no rows in grouped_hazard_1 it is because we are in the lower triangle, here we allow for taking the average of the observed during other periods, or the probability approach
  # if(nrow(grouped_hazard_1 == 0)){
  #   frame_tmp_2 <- frame_tmp %>%  group_by(p_month, time_w) %>%
  #     dplyr::summarise(weight_eta=sum(weight_eta), .groups="drop" )
  #
  #   grouped_hazard_1 <- grouped_hazard_0 %>%   inner_join(frame_tmp_2, by = c("p_month" = "p_month", "DP_rev_i" = "time_w")) %>%
  #     mutate(spend = (DP_rev_i-AP_i)%%(1/conversion_factor)+1) %>% #assuming average reporting times of 0.5 pr. development period in
  #     group_by(AP_i) %>%
  #     mutate(#value_cum = cumsum(value),
  #       #weight_new = ifelse(spend==1/conversion_factor,1/2*(S_i - S_i_lead)+(S_i_lead),1/2*(S_i - S_i_lead)),
  #       weight_new_2 = 1/2*(S_i + S_i_lead)) %>% #probability weighted approach, for each parellelogram we assume reporting half way through
  #     #arrange(AP_i, desc(DP_rev_i)) %>%  #here one weights by cumulative remaining hazard in parallelogram
  #     #mutate(value_cum_rev = cumsum(value)) %>%
  #     ungroup()
  #
  # }



#   Old
#   #hazard2 <- sum(h_tmp_3$value_cum*h_tmp_3$weight) / sum(h_tmp_3$weight) #weighing by observed claims in each area in the parallelogram
#   #hazard2_2 <- sum(h_tmp_3$value_cum*(1-h_tmp_3$value_cum_s_i)) / sum(1-h_tmp_3$value_cum_s_i)
#   hazard_o <- sum(sum(grouped_hazard_1$value*grouped_hazard_1$weight_new_2) / sum(grouped_hazard_1$weight_new_2)) #weighting by probability
#   eta_o <- 1/ hazard_o - sum(grouped_hazard_1$exposure)/sum(grouped_hazard_1$observed)
#
#   #(1+(1-eta_o)*hazard_o)/(1-eta_o*hazard_o)

#  return(list(hazard_o, eta_o) )

   development_factor = sum(grouped_hazard_2$observed, grouped_hazard_2$exposure_combined)/sum(grouped_hazard_2$exposure_combined)
   return(development_factor)

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






