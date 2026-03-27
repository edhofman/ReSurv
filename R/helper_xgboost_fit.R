# XGBoost helper functions
#
# XGBoost data preprocessing and model fitting.
#
# @import xgboost
## xgboost ----

pkg.env$xgboost_pp <-function(X,
                              Y,
                              samples_TF=NULL,
                              training_test_split=.1){

  if(!is.null(samples_TF)){
    xy=cbind(X,Y,samples_TF)
  }else{
    xy= cbind(X,Y)
  }

  tmp=xy %>%
    arrange(DP_rev_i) %>%
    as.data.frame()

  tmp[,'id'] = seq(1,dim(tmp)[1])

  if(is.null(samples_TF)){

    samples_cn <- tmp %>% select(id) %>% sample_frac(size=training_test_split)

  }else{

    cond <- tmp$samples_TF
    samples_cn <- tmp %>% select(id) %>% filter(cond)
    tmp <- tmp %>% select(-samples_TF)
  }

  suppressMessages(
    tmp_train <- tmp %>%
      semi_join(samples_cn)%>%
      arrange(DP_rev_i) %>%
      group_by(DP_rev_i) %>%
      mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame())

  ds_train_m <- xgboost::xgb.DMatrix( as.matrix.data.frame(tmp_train %>% select(colnames(X))), label=tmp_train$I)
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

  if(training_test_split<1){

    suppressMessages(
      tmp_test <- tmp %>%
        anti_join(samples_cn)%>%
        arrange(DP_rev_i) %>%
        group_by(DP_rev_i) %>%
        mutate(efron_c=(1:length(DP_rev_i)-1)/length(DP_rev_i))%>% as.data.frame())


    # ds_all_m <- xgboost::xgb.DMatrix( as.matrix(tmp,ncol=1),
    #                          label=tmp$I)
    ds_test_m <- xgboost::xgb.DMatrix( as.matrix.data.frame(tmp_test %>% select(colnames(X)), label=tmp_test$I))


    attr(ds_test_m, 'truncation') <- tmp_test$TR_i
    attr(ds_test_m, 'claim_arrival') <- tmp_test$DP_rev_i

    attr(ds_test_m, 'risk_sets') <- risks_in_the_tie(starts_i=tmp_test$TR_i,
                                                     stops_i=tmp_test$DP_rev_i,
                                                     stops = unique(tmp_test$DP_rev_i))

    #
    attr(ds_test_m, 'event_sets') <- events_in_the_tie(starts_i=tmp_test$TR_i,
                                                       stops_i=tmp_test$DP_rev_i,
                                                       stops = unique(tmp_test$DP_rev_i))
    #
    attr(ds_test_m, 'efron_c') <- tmp_test$efron_c

    attr(ds_test_m, 'tieid') <- unname(table(tmp_test$DP_rev_i))

    attr(ds_test_m, 'groups') <- rep( as.integer(names(table(tmp_test$end_time))),
                                      attr(ds_test_m, 'tieid'))

    return(list(ds_train_m=ds_train_m,
                ds_test_m=ds_test_m,
                samples_cn=samples_cn))
  }
  else{
    return(list(ds_train_m=ds_train_m,
                ds_test_m=NULL,
                samples_cn=samples_cn))
  }
}


pkg.env$fit_xgboost <- function(datads_pp,
                                hparameters=list(params=list(booster="gbtree",
                                                             eta=.01,
                                                             subsample=.5,
                                                             alpha=1,
                                                             lambda=1,
                                                             min_child_weight=.2),
                                                 print_every_n = NULL,
                                                 nrounds=10,
                                                 verbose=FALSE,
                                                 early_stopping_rounds = 500)){


  out <- xgboost::xgb.train(params = hparameters$params,
                            data =datads_pp$ds_train_m,
                            obj=cox_loss_objective,
                            nrounds = hparameters$nrounds,
                            feval= cox_evaluation_metrics,
                            watchlist = list(train=datads_pp$ds_train_m,
                                             eval=datads_pp$ds_test_m),
                            verbose= hparameters$verbose,
                            print_every_n = hparameters$print_every_n,
                            early_stopping_rounds = hparameters$early_stopping_rounds,
                            maximize = FALSE)

  return(out)



}



# Cross-validation

