# Baseline and data handling helper functions
#
# @importFrom bshazard bshazard
## Baseline calculation ----

benchmark_id <- function(X,
                                 Y,
                                 newdata.mx,
                                 remove_first_dummy=FALSE
){
  "
  Find benchmark value used in baseline calculation.
  "

  # benchmark <- cbind(X,DP_rev_i = Y$DP_rev_i) %>%
  #   dplyr::arrange(DP_rev_i) %>%
  #   first() %>%
  #   dplyr::select(-DP_rev_i) %>%
  #   as.vector() %>%
  #   unlist() %>%
  #   unname()
  #

  #new fast code
  DT <- data.table::as.data.table(cbind(X, DP_rev_i = Y$DP_rev_i))
  benchmark <- DT[order(DP_rev_i)][1, .SD, .SDcols = !'DP_rev_i']
  newdata.mx <- data.table::as.data.table(newdata.mx)

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

baseline.efron <- function(preds,
                                   dtrain,
                                   eta = 0.5) {
  eta <- validate_eta(eta)

  risk_sets  <- attr(dtrain, "risk_sets")
  event_sets <- attr(dtrain, "event_sets")

  risk_sum <- vapply(
    risk_sets,
    FUN = exp_sum_computer,
    ypred = preds,
    FUN.VALUE = numeric(1)
  )

  event_sum <- vapply(
    event_sets,
    FUN = exp_sum_computer,
    ypred = preds,
    FUN.VALUE = numeric(1)
  )

  n_events <- lengths(event_sets)

  denom <- risk_sum - eta * event_sum

  if (any(!is.finite(denom)) || any(denom <= 0)) {
    stop(
      "Non-positive denominator in baseline hazard calculation. ",
      "Check `eta`, fitted risk scores, and event/risk sets.",
      call. = FALSE
    )
  }

  baseline <- n_events / denom

  baseline
}

baseline.calc <- function(hazard_model,
                                  model.out,
                                  X,
                                  Y,
                                  training_df = NULL,
                                  eta = 0.5) {
  eta <- validate_eta(eta)


  #for baseline need full training data
  datads_pp <- xgboost_pp(X,Y, training_test_split = 1)

  if(hazard_model=="COX"){

    predict_bsln <- model.out$train_expg

  }

  if(hazard_model=="NN"){
    datads_pp_nn = deep_surv_pp(X=X,
                                        Y=Y,
                                        training_test_split = 1)

    predict_bsln <- predict_deepsurv(model.out$net, datads_pp_nn$x_train)

  }

  if(hazard_model == "XGB"){
    predict_bsln <- predict(model.out,datads_pp$ds_train_m)
  }


  predict_bsln <- predict_bsln - predict_bsln[1] #make relative to initial value, same approach as cox
  bsln <- baseline.efron(
    preds  = predict_bsln,
    dtrain = datads_pp$ds_train_m,
    eta    = eta
  )

  bsln

}



## Data handling ----

fix.double.ap<-function(features,accident_period){
  if(is.null(features)){
    return(NULL)
  }
  features[features==accident_period] <- "AP_i"

  return(features)

}

create.om.df<-function(training.data,
                               input_time_granularity,
                               years){

  tmp <- training.data %>%
    dplyr::group_by(DP_rev_i) %>%
    dplyr::summarise(Om= sum(I))

  tmp.v <- tmp$DP_rev_i
  sequ.v <- seq(1,maximum.time(years,input_time_granularity))

  cond <- !sequ.v %in% tmp.v

  if(sum(cond) > 0){

    tmp2 <- data.frame(DP_rev_i=sequ.v[cond],
                       Om=0)

    tmp <- bind_rows(tmp,tmp2)

  }

  tmp <- tmp %>% as.data.frame()


  return(tmp)

}

simplified_fill_data_frame<-function(data,
                                  continuous_features,
                                  categorical_features,
                                  years,
                                  input_time_granularity,
                                  conversion_factor){


  # browser()
  #Take the features unique values
  tmp.ls <- data %>%
    dplyr::filter((maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1))

  setDT(tmp.ls)

  cols <- c(categorical_features,
            continuous_features)


  tmp.ls <- tmp.ls[,.(.N),by=cols][,.(DP_i=1:maximum.time(years,input_time_granularity)),by=cols] #


  #Take only the training data
  tmp.existing <- data %>%
    dplyr::filter((maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1)) %>%
    dplyr::select(dplyr::all_of(continuous_features),
           dplyr::all_of(categorical_features),
           AP_i,
           DP_i) %>%
    unique() %>%
    as.data.frame()


  tmp.missing <- dplyr::setdiff(x=tmp.ls,y=tmp.existing)

  if(dim(tmp.missing)[1]==0){
    return(NULL)
  }else{

    max_dp_i = maximum.time(years,input_time_granularity)
    tmp.missing<- tmp.missing %>%
      dplyr::mutate(DP_rev_i = maximum.time(years,input_time_granularity) - DP_i+1,
             TR_i = AP_i-1, #just setting truncation to max year simulated. and accounting for
             I=0)%>%
      dplyr::filter(DP_rev_i > TR_i) %>%
      dplyr::mutate(
        DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
        AP_o = ceiling(AP_i*conversion_factor)
      ) %>%
      dplyr::mutate(TR_o= AP_o-1) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(categorical_features),
                    as.factor)) %>%
      dplyr::select(dplyr::all_of(categorical_features),
             dplyr::all_of(continuous_features),
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

fill_data_frame<-function(data,
                                  continuous_features,
                                  categorical_features,
                                  years,
                                  input_time_granularity,
                                  conversion_factor){


  #Take the features unique values
  tmp.ls <- data %>%
    dplyr::select(dplyr::all_of(continuous_features),
           dplyr::all_of(categorical_features)) %>%
    as.data.frame() %>%
    lapply(FUN=unique)


  #Take only the training data
  tmp.existing <- data %>%
    dplyr::filter((maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1)) %>%
    dplyr::select(dplyr::all_of(continuous_features),
           dplyr::all_of(categorical_features),
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
    dplyr::filter((maximum.time(years,input_time_granularity) - DP_i+1) > (AP_i-1))

  tmp.missing <- dplyr::setdiff(x=tmp.full,y=tmp.existing)

  if(dim(tmp.missing)[1]==0){
    return(NULL)
  }else{

    max_dp_i = maximum.time(years,input_time_granularity)
    tmp.missing<- tmp.missing %>%
      dplyr::mutate(DP_rev_i = maximum.time(years,input_time_granularity) - DP_i+1,
             TR_i = AP_i-1, #just setting truncation to max year simulated. and accounting for
             I=0)%>%
      dplyr::filter(DP_rev_i > TR_i) %>%
      dplyr::mutate(
        DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
        AP_o = ceiling(AP_i*conversion_factor)
      ) %>%
      dplyr::mutate(TR_o= AP_o-1) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(categorical_features),
                    as.factor)) %>%
      dplyr::select(dplyr::all_of(categorical_features),
             dplyr::all_of(continuous_features),
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




