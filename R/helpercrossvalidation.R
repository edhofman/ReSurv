cv_deep_surv <- function(hp,
                         IndividualDataPP,
                         continuous_features_scaling_method,
                         folds,
                         kfolds,
                         random_seed,
                         verbose,
                         epochs,
                         num_workers,
                         hparameters.f){

  start <- Sys.time()
  hparameters <- list(params=as.list.data.frame(hparameters.f[hp,]),
                      verbose=verbose,
                      epochs = epochs,
                      num_workers = num_workers)

  tmp.train.lkh <- vector("numeric",
                          length=folds)
  tmp.test.lkh <- vector("numeric",
                         length=folds)

  X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                    select_columns = IndividualDataPP$categorical_features)

  scaler <- pkg.env$scaler(continuous_features_scaling_method=continuous_features_scaling_method)

  Xc <- IndividualDataPP$training.data %>%
    summarize(across(all_of(IndividualDataPP$continuous_features),
                     scaler))
  X=cbind(X,Xc)

  Y=IndividualDataPP$training.data[,c("DP_rev_i", "I", "TR_i")]
  for(i in c(1:folds)){

    datads_pp = pkg.env$deep_surv_pp(X=X,
                                     Y=Y,
                                     samples_TF= c(kfolds!=i))



    model.out.k <- do.call(pkg.env$fit_deep_surv, list(data=datads_pp,
                                                       params=hparameters$params,
                                                       verbose = hparameters$verbose,
                                                       epochs = hparameters$epochs,
                                                       num_workers = hparameters$num_workers,
                                                       seed= random_seed))


    best.it <- model.out.k$log$to_pandas()[,1] == min(model.out.k$log$to_pandas()[,1])

    if(any(is.na(best.it))){
      best.it<-1
    }
    tmp.train.lkh[i] <- unname(unlist(model.out.k$log$to_pandas()[best.it,]['train_loss']))
    tmp.test.lkh[i] <- unname(unlist(model.out.k$log$to_pandas()[best.it,]['val_loss']))

  }

  time <- as.numeric(difftime(Sys.time(), start, units='mins'))

  c(mean(tmp.train.lkh),mean(tmp.test.lkh),time)

}

cv_xgboost <- function( hp,
                        IndividualDataPP,
                        folds,
                        kfolds,
                        print_every_n ,
                        nrounds,
                        verbose,
                        early_stopping_rounds ,
                        hparameters.f){
    start <- Sys.time()
    hparameters <- list(params=as.list.data.frame(hparameters.f[hp,]),
                        print_every_n=print_every_n,
                        nrounds=nrounds,
                        verbose=verbose,
                        early_stopping_rounds=early_stopping_rounds)

    tmp.train.lkh <- vector("numeric",
                            length=folds)
    tmp.test.lkh <- vector("numeric",
                           length=folds)

    for(i in c(1:folds)){

      X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                        select_columns = IndividualDataPP$categorical_features,
                                        remove_first_dummy=T)

      scaler <- pkg.env$scaler(continuous_features_scaling_method = "minmax")

      Xc <- IndividualDataPP$training.data %>%
        summarize(across(all_of(IndividualDataPP$continuous_features),
                         scaler))

      X=cbind(X,Xc)

      Y=IndividualDataPP$training.data[,c("DP_rev_i", "I", "TR_i")]

      datads_pp =  pkg.env$xgboost_pp(X,
                                      Y,
                                      samples_TF= c(kfolds!=i))



      model.out.k <- do.call(pkg.env$fit_xgboost, list(datads_pp=datads_pp,
                                                       hparameters=hparameters))


      best.it <- model.out.k$best_iteration
      tmp.train.lkh[i] <- model.out.k$evaluation_log$`train_log_partial likelihood`[best.it]
      tmp.test.lkh[i] <- model.out.k$evaluation_log$`eval_log_partial likelihood`[best.it]

    }
    time <- as.numeric(difftime(Sys.time(), start, units='mins'))

    c(mean(tmp.train.lkh),mean(tmp.test.lkh),time)
}
