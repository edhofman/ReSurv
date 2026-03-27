# Evaluation and prediction helper functions
#
# Likelihood evaluation, triangle completion, prediction adjustment, and survival CRPS helpers.
#
# @importFrom dplyr reframe lag full_join rename
# @importFrom tidyr replace_na
## Evaluation metrics ----

pkg.env$evaluate_lkh_nn <-function(X_train,
                                   Y_train,
                                   model){


  # data_transformed <- cbind(X, Y)
  data_train <- cbind(X_train, DP_rev_i = Y_train$DP_rev_i) %>%
    arrange(DP_rev_i) %>%
    select(-DP_rev_i) %>%
    as.matrix()

  preds <- pkg.env$predict_deepsurv(model$net, data_train)
  preds <-preds-preds[1]


  xy_tr=cbind(X_train,Y_train)

  tmp_tr=xy_tr %>%
    arrange(DP_rev_i) %>%
    as.data.frame()


  tmp_train <- tmp_tr %>%
    arrange(DP_rev_i) %>%
    group_by(DP_rev_i) %>%
    mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame()


  ds_train_m <- tmp_train %>%
    arrange(DP_rev_i) %>%
    group_by(DP_rev_i) %>%
    mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame()


  attr(ds_train_m, 'truncation') <- tmp_train$TR_i
  attr(ds_train_m, 'claim_arrival') <- tmp_train$DP_rev_i


  attr(ds_train_m, 'risk_sets') <- risks_in_the_tie(starts_i=tmp_train$TR_i,
                                                    stops_i=tmp_train$DP_rev_i,
                                                    stops = unique(tmp_train$DP_rev_i))

  attr(ds_train_m, 'event_sets') <- events_in_the_tie(starts_i=tmp_train$TR_i,
                                                      stops_i=tmp_train$DP_rev_i,
                                                      stops = unique(tmp_train$DP_rev_i))

  attr(ds_train_m, 'efron_c') <- tmp_train$efron_c

  attr(ds_train_m, 'tieid') <- unname(table(tmp_train$DP_rev_i))

  attr(ds_train_m, 'groups') <- rep( as.integer(names(table(tmp_train$end_time))),
                                     attr(ds_train_m, 'tieid'))







  train_lkh=cox_evaluation_metrics(dtrain=ds_train_m,
                                  preds=as.vector(preds))


  return(train_lkh)


}

pkg.env$evaluate_lkh_xgb <-function(X_train,
                                    Y_train,
                                    dset,
                                    samples_cn,
                                    model){

  xy_tr=cbind(X_train,Y_train) %>%
    arrange(DP_rev_i) %>%
    as.data.frame()

  id <- seq(1, dim(X_train)[1])
  cond <- id %in% samples_cn$id

  if(dset=='os'){cond <- !cond}

  tmp_tr=xy_tr[cond,] %>%
    arrange(DP_rev_i) %>%
    as.data.frame()

  tmp_train <- tmp_tr %>%
    arrange(DP_rev_i) %>%
    group_by(DP_rev_i) %>%
    mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame()



  ds_train_m <- xgboost::xgb.DMatrix( as.matrix.data.frame(tmp_train %>% select(colnames(X_train))),
                                      label=tmp_train$I)

  attr(ds_train_m, 'truncation') <- tmp_train$TR_i
  attr(ds_train_m, 'claim_arrival') <- tmp_train$DP_rev_i


  attr(ds_train_m, 'risk_sets') <- risks_in_the_tie(starts_i=tmp_train$TR_i,
                                                    stops_i=tmp_train$DP_rev_i,
                                                    stops = unique(tmp_train$DP_rev_i))

  attr(ds_train_m, 'event_sets') <- events_in_the_tie(starts_i=tmp_train$TR_i,
                                                      stops_i=tmp_train$DP_rev_i,
                                                      stops = unique(tmp_train$DP_rev_i))

  attr(ds_train_m, 'efron_c') <- tmp_train$efron_c

  attr(ds_train_m, 'tieid') <- unname(table(tmp_train$DP_rev_i))

  attr(ds_train_m, 'groups') <- rep( as.integer(names(table(tmp_train$end_time))),
                                     attr(ds_train_m, 'tieid'))




  # if(hazard_model == "XGB"){
  preds_tr <- predict(model,ds_train_m)
  preds_tr <- preds_tr - preds_tr[1]



  train_lkh=cox_evaluation_metrics(dtrain=ds_train_m,
                                  preds=preds_tr)


  return(train_lkh)


}


pkg.env$evaluate_lkh_cox <-function(X_train,
                                    Y_train,
                                    model){

  xy_tr=cbind(X_train,Y_train)


  tmp_tr=xy_tr %>%
    arrange(DP_rev_i) %>%
    as.data.frame()


  tmp_train <- tmp_tr %>%
    arrange(DP_rev_i) %>%
    group_by(DP_rev_i) %>%
    mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame()


  ds_train_m <- X_train
  # if(hazard_model == "XGB"){
  #   ds_train_m <- xgboost::xgb.DMatrix( as.matrix.data.frame(tmp_train %>% select(colnames(X_train))),
  #                                       label=tmp_train$I)}

  attr(ds_train_m, 'truncation') <- tmp_train$TR_i
  attr(ds_train_m, 'claim_arrival') <- tmp_train$DP_rev_i


  attr(ds_train_m, 'risk_sets') <- risks_in_the_tie(starts_i=tmp_train$TR_i,
                                                    stops_i=tmp_train$DP_rev_i,
                                                    stops = unique(tmp_train$DP_rev_i))

  attr(ds_train_m, 'event_sets') <- events_in_the_tie(starts_i=tmp_train$TR_i,
                                                      stops_i=tmp_train$DP_rev_i,
                                                      stops = unique(tmp_train$DP_rev_i))

  attr(ds_train_m, 'efron_c') <- tmp_train$efron_c

  attr(ds_train_m, 'tieid') <- unname(table(tmp_train$DP_rev_i))

  attr(ds_train_m, 'groups') <- rep( as.integer(names(table(tmp_train$end_time))),
                                     attr(ds_train_m, 'tieid'))


  preds_tr <- predict(model$cox,ds_train_m)

  train_lkh=cox_evaluation_metrics(dtrain=ds_train_m,
                                  preds=preds_tr)


  return(train_lkh)


}


adjust.predictions <- function(ResurvFit,
                               hazard_model,
                               idata){

  formula_ct <- idata$string_formula_i

  newdata <- create.df.2.fcst(IndividualDataPP=idata,
                              hazard_model=hazard_model)

  # create data frame of occurrencies to weight development factors
  Om.df <-   ResurvFit$Om.df


  if(hazard_model=="COX"){

    data=idata$training.data
    X=data %>%
      select(c(idata$continuous_features,idata$categorical_features))

    Y=data[,c("DP_rev_i", "I", "TR_i")]

    model.out <- ResurvFit$model.out$model.out
    coxlp <-  predict(model.out$cox,
              newdata=newdata,
              'lp',
              reference='zero')

    expg <- exp(coxlp)

    bs_hazard <- basehaz(model.out$cox,
                         newdata=newdata, # here the baseline is refitted
                         centered=FALSE) %>%
      mutate(hazard = hazard-lag(hazard,default=0))

    bsln <- data.frame(baseline=bs_hazard$hazard,
                       DP_rev_i=ceiling(bs_hazard$time))  #$hazard

    hazard_frame <- cbind(newdata, expg)
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"

  }

  if(hazard_model=="NN"){

    X <- pkg.env$model.matrix.creator(data= idata$training.data,
                                      select_columns = idata$categorical_features)

    scaler <- pkg.env$scaler(continuous_features_scaling_method='minmax')

    Xc <- idata$training.data %>%
      reframe(across(all_of(idata$continuous_features),
                     scaler))


    X = cbind(X,Xc)

    Y=idata$training.data[,c("DP_rev_i", "I", "TR_i")]

    datads_pp = pkg.env$deep_surv_pp(X=X,
                                     Y=Y,
                                     training_test_split = 1)

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = ResurvFit$model.out$model.out,
                                  X=X,
                                  Y=Y)

    newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=idata$training.data,
                                          newdata=newdata,
                                          continuous_features=idata$continuous_features,
                                          categorical_features=idata$categorical_features)



    x_fc = as.matrix(newdata.mx)



    beta_ams <- pkg.env$predict_deepsurv(ResurvFit$model.out$model.out$net, x_fc)

    #make to hazard relative to initial model, to have similiar interpretation as standard cox

    benchmark_id <- pkg.env$benchmark_id(X = X,
                                         Y =Y ,
                                         newdata.mx = newdata.mx
    )

    pred_relative <- beta_ams - beta_ams[benchmark_id]

    expg <- exp(pred_relative)
    hazard_frame <- cbind(newdata,expg)
    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(idata$training.data$DP_rev_i))))


  }


  hazard_frame <- hazard_frame %>%
    full_join(bsln,
              by="DP_rev_i") %>%
    as.data.frame() %>%
    replace_na(list(baseline=0))

  hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']

  #Add development and relevant survival values to the hazard_frame
  hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                    Om.df=Om.df,
                                                    categorical_features = idata$categorical_features,
                                                    continuous_features = idata$continuous_features,
                                                    calendar_period_extrapolation = idata$calendar_period_extrapolation)


  return(hazard_frame_updated)

}

# survival crps ----



survival_information<-function(x,
                               group,
                               hazard_list){

  tmp <-hazard_list[[group]]

  x.vals = tmp$x.vals
  cdf2_i = tmp$cdf2_i
  DP_rev_i = tmp$DP_rev_i
  S2_i = tmp$S2_i

  crps=sum((x.vals*cdf2_i)[DP_rev_i<=x])+sum((x.vals*S2_i)[DP_rev_i>x])

  return(crps)

}


pkg.env$complete_lt_predictions_i <- function(dt,max_dp){

  "
  Add the missing combinations of AP_i and DP_i to the long format output long_tr_input to create a triangle data.frame.
  "

  seq1_main <- unique(dt$AP_i)
  seq2_main <- unique(dt$DP_i)

  complete_seq <- 1:max_dp

  diff1<-setdiff(complete_seq, seq1_main)
  diff2<-setdiff(complete_seq, seq2_main)


  if(length(diff1)==0){

    if(length(diff2)==0){

      return(NULL)

    }else{

      return(CJ(seq1_main,diff2))

    }

  }else{

    if(length(diff2)==0){

      return(CJ(diff1,seq2_main))

    }else{

      return(CJ(diff1,diff2))

    }


  }


}


pkg.env$complete_lt_predictions_o <- function(dt,max_dp){

  "
  Add the missing combinations of AP_o and DP_o to the long format output long_tr_output to create a triangle data.frame.
  "

  seq1_main <- unique(dt$AP_o)
  seq2_main <- unique(dt$DP_o)

  complete_seq <- 1:max_dp

  diff1<-setdiff(complete_seq, seq1_main)
  diff2<-setdiff(complete_seq, seq2_main)


  if(length(diff1)==0){

    if(length(diff2)==0){

      return(NULL)

    }else{

      return(CJ(seq1_main,diff2))

    }

  }else{

    if(length(diff2)==0){

      return(CJ(diff1,seq2_main))

    }else{

      return(CJ(diff1,diff2))

    }


  }


}

pkg.env$find_lt_input <- function(dt,max_dp){

  "
  Return the lower triangular output in a data.frame format (input granularity).
  "


  dt <- as.data.table(dt)

  dt<-dt[,.(value=sum(IBNR,na.rm=TRUE)),by=.(AP_i,DP_i)]

  add_up <- pkg.env$complete_lt_predictions_i(dt,max_dp)

  if(!is.null(add_up)){

    colnames(add_up) <- c("AP_i","DP_i")

    add_up[["value"]] <- 0

    dt <- rbind(dt,add_up)

  }

  dt.w<-dcast(dt, AP_i ~ DP_i , value.var = "value") %>%
    as.data.frame()

  rownames(dt.w) <- dt.w$AP_i
  dt.w <- dt.w[,-1]

  for(i in 1:max_dp){

    for(j in 1:max_dp){

      if((i+j-1) <= (max_dp)){

        dt.w[i,j]<-NA

      }

    }

  }

  return(dt.w)

}


pkg.env$find_lt_output <- function(dt,
                                   max_dp,
                                   cut_point){


  "
  Return the lower triangular output in a data.frame format (output granularity).
  "

  dt <- as.data.table(dt)

  dt<-dt[,.(value=sum(IBNR,na.rm=TRUE)),by=.(AP_o,DP_o)]

  add_up <- pkg.env$complete_lt_predictions_o(dt,max_dp)

  if(!is.null(add_up)){

    colnames(add_up) <- c("AP_o","DP_o")

    add_up[["value"]] <- 0

    dt <- rbind(dt,add_up)

  }

  dt.w<-dcast(dt, AP_o ~ DP_o , value.var = "value") %>%
    as.data.frame()

  rownames(dt.w) <- dt.w$AP_o
  dt.w <- dt.w[,-1]



  for(i in 1:max_dp){

    for(j in 1:max_dp){

      if((i+j-1) <= (max_dp)){

        if(is.na(dt.w[i,j]) | dt.w[i,j]==0){dt.w[i,j]<-NA}

      }

    }

  }

  return(dt.w[1:cut_point,])

}



manually_extract_info_for_scoring_cont <- function(ReSurvFit,
                                                   hazard_model,
                                                   IndividualDataPP,
                                                   tie = "efron",
                                                   baseline = "spline",
                                                   continuous_features_scaling_method = "minmax",
                                                   random_seed = 1,
                                                   hparameters = list(),
                                                   percentage_data_training = .8,
                                                   grouping_method = "exposure",
                                                   check_value = 1.85,
                                                   eta=0.5,
                                                   simplifier=TRUE){


  set.seed(random_seed)

  formula_ct <- as.formula(IndividualDataPP$string_formula_i)


    cont_f <- IndividualDataPP$continuous_features
    cat_f <- IndividualDataPP$categorical_features
    columns_for_grouping <- unique(c(cont_f,cat_f,"AP_i"))

    tmp <- as.data.table(IndividualDataPP$full.data)

    out <- tmp[,.(.N),by=columns_for_grouping][,..columns_for_grouping]

    l4 <- list()

    l4$DP_rev_i <- min(IndividualDataPP$full.data[,'DP_rev_i']):max(IndividualDataPP$full.data[,'DP_rev_i'])

    l4<-do.call(CJ, c(l4, sorted = FALSE))

    newdata<-as.data.frame(setkey(out[,c(k=1,.SD)],k)[l4[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL])



  # logical: check if we work with a baseline model
  is_baseline_model = is.null(c(IndividualDataPP$categorical_features,
                                IndividualDataPP$continuous_features))


  if(hazard_model=="COX"){
    #
    data=IndividualDataPP$training.data

    X=data %>%
      select(c(IndividualDataPP$continuous_features,IndividualDataPP$categorical_features))

    Y=IndividualDataPP$full.data[,c("DP_rev_i", "I", "TR_i")]

    cox <- coxph(formula_ct, data=data, ties="efron")
    cox_lp <- predict(cox,newdata=newdata,'lp',reference='zero')

    cox_training_lp <- predict(cox,newdata=newdata %>% arrange(DP_rev_i) %>% as.data.frame(),'lp',reference='zero')


    model.out <- list(cox=cox,
                      cox_lp=cox_lp,
                      expg = exp(cox_lp))

    ## NEW BASELINE COMPUTATION (RESURV)

      scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

      Xc_tmp_bsln <- IndividualDataPP$full.data %>%
        reframe(across(all_of(IndividualDataPP$continuous_features),
                       scaler))


      if(!is.null(IndividualDataPP$categorical_features)){


        X_tmp_bsln <- pkg.env$model.matrix.creator(data= IndividualDataPP$full.data,
                                                   select_columns = IndividualDataPP$categorical_features,
                                                   remove_first_dummy=T)


        X_tmp_bsln=cbind(X_tmp_bsln,Xc_tmp_bsln)

      }else{

        X_tmp_bsln= Xc_tmp_bsln

      }



    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X_tmp_bsln,
                                  Y=Y)

    #

    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualDataPP$training.data$DP_rev_i))))

    ### make it relative


      newdata.bs <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$full.data,
                                            newdata=newdata,
                                            continuous_features=IndividualDataPP$continuous_features,
                                            categorical_features=IndividualDataPP$categorical_features)

      benchmark_id <- pkg.env$benchmark_id(X = X_tmp_bsln,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=T)








    pred_relative <- cox_lp-cox_lp[benchmark_id]

    ###

    hazard_frame <- cbind(newdata, exp(pred_relative))
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"



    is_lkh <- pkg.env$evaluate_lkh_cox(X_train=X,
                                       Y_train=Y,
                                       model=model.out)


    os_lkh <- NULL


  }

  if(hazard_model=="NN"){


    Y=IndividualDataPP$training.data[,c("DP_rev_i", "I", "TR_i")]

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    if(is_baseline_model){

      X <- data.frame(intercept_1 = rep(1,dim(Y)[1]))

    }else{

      scaler <- pkg.env$scaler(continuous_features_scaling_method=continuous_features_scaling_method)

      Xc <- IndividualDataPP$training.data %>%
        reframe(across(all_of(IndividualDataPP$continuous_features),
                       scaler))

      if(!is.null(IndividualDataPP$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                          select_columns = IndividualDataPP$categorical_features)

        X = cbind(X,Xc)


      }else{


        X <- Xc

      }

    }

    datads_pp = pkg.env$deep_surv_pp(X=X,
                                     Y=Y,
                                     training_test_split = training_test_split)

    hparameters <- pkg.env$nn_hparameter_nodes_grid(hparameters)

    hparameters <- list(params=as.list.data.frame(hparameters),
                        verbose=hparameters$verbose,
                        epochs = hparameters$epochs,
                        num_workers = hparameters$num_workers)


    model.out <- pkg.env$fit_deep_surv(datads_pp,
                                       params=hparameters$params,
                                       verbose = hparameters$verbose,
                                       epochs = hparameters$epochs,
                                       num_workers = hparameters$num_workers,
                                       seed = random_seed)


    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X,
                                  Y=Y)

    if(is_baseline_model){

      newdata.mx <- data.frame(intercept_1= rep(1,dim(newdata)[1]))

    }else{

      newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$training.data,
                                            newdata=newdata,
                                            continuous_features=IndividualDataPP$continuous_features,
                                            categorical_features=IndividualDataPP$categorical_features)}



    x_fc = as.matrix(newdata.mx)


    beta_ams <- pkg.env$predict_deepsurv(model.out$net, x_fc)

    #make to hazard relative to initial model, to have similiar interpretation as standard cox

    benchmark_id <- pkg.env$benchmark_id(X = X,
                                         Y =Y ,
                                         newdata.mx = newdata.mx
    )

    pred_relative <- beta_ams - beta_ams[benchmark_id]

    expg <- exp(pred_relative)
    hazard_frame <- cbind(newdata,expg)
    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualDataPP$training.data$DP_rev_i))))



    if(!inherits(datads_pp$lkh_eval_data$data_train,"data.frame")){


      is_lkh <- pkg.env$evaluate_lkh_nn(X_train=as.data.frame(datads_pp$lkh_eval_data$data_train),
                                        Y_train=datads_pp$lkh_eval_data$y_train,
                                        model=model.out)

      os_lkh <- pkg.env$evaluate_lkh_nn(X_train=as.data.frame(datads_pp$lkh_eval_data$data_val),
                                        Y_train=datads_pp$lkh_eval_data$y_val,
                                        model=model.out)


    }else{

      is_lkh <- pkg.env$evaluate_lkh_nn(X_train=datads_pp$lkh_eval_data$data_train,
                                        Y_train=datads_pp$lkh_eval_data$y_train,
                                        model=model.out)

      os_lkh <- pkg.env$evaluate_lkh_nn(X_train=datads_pp$lkh_eval_data$data_val,
                                        Y_train=datads_pp$lkh_eval_data$y_val,
                                        model=model.out)

    }




  }

  if(hazard_model == "XGB"){

    Y=IndividualDataPP$training.data[,c("DP_rev_i", "I", "TR_i")]

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    if(is_baseline_model){

      X= data.frame(intercept_1 = rep(1,dim(Y)[1]))

    }else{

      scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

      Xc <- IndividualDataPP$training.data %>%
        reframe(across(all_of(IndividualDataPP$continuous_features),
                       scaler))


      if(!is.null(IndividualDataPP$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                          select_columns = IndividualDataPP$categorical_features,
                                          remove_first_dummy=T)

        X=cbind(X,Xc)
      }else{


        X <- Xc


      }

    }



    datads_pp <- pkg.env$xgboost_pp(X=X,
                                    Y=Y,
                                    training_test_split=training_test_split)

    model.out <- pkg.env$fit_xgboost(datads_pp,
                                     hparameters=hparameters)

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X,
                                  Y=Y)


    if(is_baseline_model){

      newdata.mx <- xgboost::xgb.DMatrix(as.matrix(rep(1, dim(newdata)[1])))

    }else{

      newdata.mx <- pkg.env$df.2.fcst.xgboost.pp(data=IndividualDataPP$training.data,
                                                 newdata=newdata,
                                                 continuous_features=IndividualDataPP$continuous_features,
                                                 categorical_features=IndividualDataPP$categorical_features)
    }

    pred <- predict(model.out,newdata.mx)



    if(is_baseline_model){

      newdata.bs <- data.frame(intercept_1 = rep(1, dim(newdata)[1]))

      benchmark_id <- pkg.env$benchmark_id(X = X,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=F)

    }else{
      #make to hazard relative to initial model, to have similiar interpretation as standard cox
      newdata.bs <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$training.data,
                                            newdata=newdata,
                                            continuous_features=IndividualDataPP$continuous_features,
                                            categorical_features=IndividualDataPP$categorical_features)

      benchmark_id <- pkg.env$benchmark_id(X = X,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=T)}


    pred_relative <- pred - pred[benchmark_id]

    expg <- exp(pred_relative)

    hazard_frame <- cbind(newdata,expg)

    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualDataPP$training.data$DP_rev_i))))

    # compute the likelihood of the fitted model (upper triangle)
    is_lkh <- pkg.env$evaluate_lkh_xgb(X_train=X,
                                       Y_train=Y,
                                       dset='is',
                                       samples_cn=datads_pp$samples_cn,
                                       model=model.out)

    os_lkh <- pkg.env$evaluate_lkh_xgb(X_train=X,
                                       Y_train=Y,
                                       dset='os',
                                       samples_cn=datads_pp$samples_cn,
                                       model=model.out)

  }

  ##################################################################################


  hazard_frame <- hazard_frame %>%
    full_join(bsln,
              by="DP_rev_i") %>%
    as.data.frame() %>%
    replace_na(list(baseline=0))

  hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']




  #Add development and relevant survival values to the hazard_frame
  hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                    # Om.df=Om.df,
                                                    eta_old=eta,
                                                    categorical_features = IndividualDataPP$categorical_features,
                                                    continuous_features = IndividualDataPP$continuous_features,
                                                    calendar_period_extrapolation = IndividualDataPP$calendar_period_extrapolation)


  out_hz_frame <-  hazard_frame_updated %>%
    mutate(DP_i=pkg.env$maximum.time(IndividualDataPP$years, IndividualDataPP$input_time_granularity)-DP_rev_i+1) %>%
    relocate(DP_i, .after =  AP_i) %>%
    rename(f_i=dev_f_i,
           cum_f_i=cum_dev_f_i)

  out=list(model.out=list(data=X,
                          model.out=model.out),
           simplifier=simplifier,
           is_lkh=is_lkh,
           os_lkh=os_lkh,
           hazard_frame = out_hz_frame,
           hazard_model = hazard_model,
           IndividualDataPP = IndividualDataPP)

  class(out) <- c('ReSurvFit')

  return(out)












}











