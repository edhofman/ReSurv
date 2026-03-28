# Baseline and data handling helper functions
#
# @importFrom bshazard bshazard
## Baseline calculation ----

pkg.env$benchmark_id <- function(X,
                                 Y,
                                 newdata.mx,
                                 remove_first_dummy=FALSE
){
  "
  Find benchmark value used in baseline calculation.
  "

  # benchmark <- cbind(X,DP_rev_i = Y$DP_rev_i) %>%
  #   arrange(DP_rev_i) %>%
  #   first() %>%
  #   select(-DP_rev_i) %>%
  #   as.vector() %>%
  #   unlist() %>%
  #   unname()
  # 

  #new fast code
  DT <- cbind(X, DP_rev_i = Y$DP_rev_i)
  benchmark <- DT[order(DP_rev_i)][1, .SD, .SDcols = !'DP_rev_i']

  ## probably not useful rows below
  # benchmark <- unname(unlist(res, use.names = FALSE))
  # benchmark <- as.list(benchmark)
  # benchmark <- as.data.table(benchmark)
  # setnames(benchmark, names(newdata.mx))

  if(remove_first_dummy==TRUE){
    #old code
    # newdata.mx <- data.frame(newdata.mx[,colnames(newdata.mx) %in% names(X)])
    newdata.mx<- newdata.mx[,.SD,.SDcols = colnames(newdata.mx) %in% names(X)]
  }

  #old code
  # benchmark_id <- which(apply(newdata.mx, 1, function(x) sum(benchmark == x) == length(benchmark) ))[1]

  benchmark_id <- newdata.mx[benchmark, on = names(newdata.mx), which = TRUE][1]



  return(benchmark_id)

}

#Note that we for all methods apply xgboost naming convention

pkg.env$baseline.efron <- function(preds, dtrain){

  risk_sets <- attr(dtrain, 'risk_sets')
  event_sets <- attr(dtrain, 'event_sets')
  # efron_c<-attr(dtrain, 'efron_c')
  tieid<- attr(dtrain, 'tieid')

  exp_p_sum <- sapply(risk_sets,FUN=exp_sum_computer, ypred=preds)
  exp_p_tie <- sapply(event_sets,FUN=exp_sum_computer, ypred=preds)

  exp_p_sum <- rep(sapply(risk_sets,FUN=exp_sum_computer, ypred=preds), tieid)
  exp_p_tie <-  rep(sapply(event_sets,FUN=exp_sum_computer, ypred=preds), tieid)

  # alpha_i <- 1/(exp_p_sum-efron_c*exp_p_tie)

  alpha_i <- 1/(exp_p_sum-.5*exp_p_tie)

  baseline <- sapply(event_sets, FUN = function(x,values){sum(values[x]) }, values=alpha_i)

  baseline

}

pkg.env$baseline.calc <- function(hazard_model,
                                  model.out,
                                  X,
                                  Y,
                                  training_df = NULL){


  #for baseline need full training data
  datads_pp <- pkg.env$xgboost_pp(X,Y, training_test_split = 1)

  if(hazard_model=="COX"){

    predict_bsln <- model.out$train_expg

  }

  if(hazard_model=="NN"){
    datads_pp_nn = pkg.env$deep_surv_pp(X=X,
                                        Y=Y,
                                        training_test_split = 1)

    predict_bsln <- pkg.env$predict_deepsurv(model.out$net, datads_pp_nn$x_train)

  }

  if(hazard_model == "XGB"){
    predict_bsln <- predict(model.out,datads_pp$ds_train_m)
  }


  predict_bsln <- predict_bsln - predict_bsln[1] #make relative to initial value, same approach as cox
  bsln <- pkg.env$baseline.efron(predict_bsln,
                                 datads_pp$ds_train_m)

  bsln

}



## Data handling ----

pkg.env$fix.double.ap<-function(features,accident_period){
  if(is.null(features)){
    return(NULL)
  }
  features[features==accident_period] <- "AP_i"

  return(features)

}

pkg.env$create.om.df<-function(training.data,
                               input_time_granularity,
                               years){

  tmp <- training.data %>%
    group_by(DP_rev_i) %>%
    summarise(Om= sum(I))

  tmp.v <- tmp$DP_rev_i
  sequ.v <- seq(1,pkg.env$maximum.time(years,input_time_granularity))

  cond <- !sequ.v %in% tmp.v

  if(sum(cond) > 0){

    tmp2 <- data.frame(DP_rev_i=sequ.v[cond],
                       Om=0)

    tmp <- bind_rows(tmp,tmp2)

  }

  tmp <- tmp %>% as.data.frame()


  return(tmp)

}

pkg.env$simplified_fill_data_frame<-function(data,
                                  continuous_features,
                                  categorical_features,
                                  years,
                                  input_time_granularity,
                                  conversion_factor){


  # browser()
  #Take the features unique values
  tmp.ls <- data %>%
    filter((pkg.env$maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1))

  setDT(tmp.ls)

  cols <- c(categorical_features,
            continuous_features)


  tmp.ls <- tmp.ls[,.(.N),by=cols][,.(DP_i=1:pkg.env$maximum.time(years,input_time_granularity)),by=cols] #


  #Take only the training data
  tmp.existing <- data %>%
    filter((pkg.env$maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1)) %>%
    select(all_of(continuous_features),
           all_of(categorical_features),
           AP_i,
           DP_i) %>%
    unique() %>%
    as.data.frame()


  tmp.missing <- dplyr::setdiff(x=tmp.ls,y=tmp.existing)

  if(dim(tmp.missing)[1]==0){
    return(NULL)
  }else{

    max_dp_i = pkg.env$maximum.time(years,input_time_granularity)
    tmp.missing<- tmp.missing %>%
      mutate(DP_rev_i = pkg.env$maximum.time(years,input_time_granularity) - DP_i+1,
             TR_i = AP_i-1, #just setting truncation to max year simulated. and accounting for
             I=0)%>%
      filter(DP_rev_i > TR_i) %>%
      mutate(
        DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
        AP_o = ceiling(AP_i*conversion_factor)
      ) %>%
      mutate(TR_o= AP_o-1) %>%
      mutate(across(all_of(categorical_features),
                    as.factor)) %>%
      select(all_of(categorical_features),
             all_of(continuous_features),
             AP_i,
             AP_o,
             DP_i,
             DP_rev_i,
             DP_rev_o,
             TR_i,
             TR_o,
             I) %>%
      as.data.frame()

    return(tmp.missing)

  }}

pkg.env$fill_data_frame<-function(data,
                                  continuous_features,
                                  categorical_features,
                                  years,
                                  input_time_granularity,
                                  conversion_factor){


  #Take the features unique values
  tmp.ls <- data %>%
    select(all_of(continuous_features),
           all_of(categorical_features)) %>%
    as.data.frame() %>%
    lapply(FUN=unique)


  #Take only the training data
  tmp.existing <- data %>%
    filter((pkg.env$maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1)) %>%
    select(all_of(continuous_features),
           all_of(categorical_features),
           AP_i,
           DP_i) %>%
    unique() %>%
    as.data.frame()

  # accidents <- sort(unique(data$AP_i))
  # developments <- sort(unique(data$DP_i))

  # v1 <- diff(as.integer(accidents))
  # v2 <- diff(as.integer(sort(unique(developments))))

  #Take the complete sequence
  tmp1 <- min(data$AP_i):max(data$AP_i)
  tmp2 <- 1:max(data$DP_i)

  tmp.ls$AP_i <- tmp1
  tmp.ls$DP_i <- tmp2

  tmp.full <- expand.grid(tmp.ls) %>%
    as.data.frame() %>%
    filter((pkg.env$maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1))

  tmp.missing <- dplyr::setdiff(x=tmp.full,y=tmp.existing)

  if(dim(tmp.missing)[1]==0){
    return(NULL)
  }else{

    max_dp_i = pkg.env$maximum.time(years,input_time_granularity)
    tmp.missing<- tmp.missing %>%
      mutate(DP_rev_i = pkg.env$maximum.time(years,input_time_granularity) - DP_i+1,
             TR_i = AP_i-1, #just setting truncation to max year simulated. and accounting for
             I=0)%>%
      filter(DP_rev_i > TR_i) %>%
      mutate(
        DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
        AP_o = ceiling(AP_i*conversion_factor)
      ) %>%
      mutate(TR_o= AP_o-1) %>%
      mutate(across(all_of(categorical_features),
                    as.factor)) %>%
      select(all_of(categorical_features),
             all_of(continuous_features),
             AP_i,
             AP_o,
             DP_i,
             DP_rev_i,
             DP_rev_o,
             TR_i,
             TR_o,
             I) %>%
      as.data.frame()

    return(tmp.missing)

  }}




