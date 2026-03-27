pkg.env$xgboost_cv <- function(IndividualDataPP,
                               folds,
                               kfolds,
                               print_every_n = 1L,
                               nrounds= NULL,
                               verbose=1,
                               early_stopping_rounds = NULL,
                               hparameters.f,
                               out,
                               verbose.cv=FALSE,
                               parallel=FALSE,
                               ncores=1,
                               random_seed){

  "Function to perform K-fold cross-validation with xgboost"
  if(parallel == TRUE){
    # handle UNIx-operated systems seperatly?.Platform$OS.type
    require(parallel)

    cl <- makeCluster(ncores)

    objects_export <- list(
      "random_seed"
    )
    clusterExport(cl, objects_export, envir = environment())

    clusterEvalQ(cl, {library("ReSurv")
      set.seed(random_seed)} )

    out[,c("train.lkh","test.lkh", "time")] <- t(parSapply(cl, 1:dim(hparameters.f)[1],  FUN =cv_xgboost,
                                                           IndividualDataPP=IndividualDataPP,
                                                           folds=folds,
                                                           kfolds=kfolds,
                                                           print_every_n=print_every_n,
                                                           nrounds=nrounds,
                                                           verbose=FALSE,
                                                           early_stopping_rounds=early_stopping_rounds,
                                                           hparameters.f=hparameters.f))
    stopCluster(cl)
  }
  else{
    for(hp in 1:dim(hparameters.f)[1]){

      if(verbose.cv){cat(as.character(Sys.time()),
                         "Testing hyperparameters combination",
                         hp,
                         "out of",
                         dim(hparameters.f)[1], "\n")}


      out[hp,c("train.lkh","test.lkh", "time")] <- cv_xgboost(hp,
                                                              IndividualDataPP,
                                                              folds,
                                                              kfolds,
                                                              print_every_n,
                                                              nrounds,
                                                              verbose,
                                                              early_stopping_rounds,
                                                              hparameters.f)
    }
  }
  return(out)

}


# nn cv -----
pkg.env$nn_hparameter_nodes_grid <- function(hparameters, cv = FALSE){
  "
  Expand hyperparameter grid for network structure
  "
  if("num_layers" %in% names(hparameters)){
    if(cv == TRUE){
      # for ( i in 1:max(hparameters$num_layers)){
      #   hparameters[[paste0("node_",i)]] <- hparameters$num_nodes
      #  }
      names <- sapply(1:max(hparameters$num_layers), function(x){paste0("node_",x)})

      suppressWarnings (
        hparameters <-hparameters %>% rowwise() %>%
          mutate(new = paste(paste(rep(num_nodes, num_layers),collapse=","),
                             paste(rep("NA", max(hparameters$num_layers) - num_layers), collapse = ","),
                             sep =",")) %>%
          separate(new, into = names, sep=",") %>%  ungroup() %>%
          mutate(across(starts_with("node_"), as.integer))
      )

    }
    else{
      if (hparameters$num_layers == length(hparameters$num_nodes)){
        for ( i in 1:hparameters$num_layers){
          hparameters[[paste0("node_",i)]] <- hparameters$num_nodes[i]
        }
      } else if(length(hparameters$num_nodes) == 1) {
        for ( i in 1:hparameters$num_layers){
          hparameters[[paste0("node_",i)]] <- hparameters$num_nodes
        }
      } else{
        warning(paste0("Num_nodes hyperparameter not inputted correctly.
                      Please either input one number, which will be used for all layer, or the same amount of nodes as layers.
                       Defaulting to first element in Num_nodes list for all layers."))
        for ( i in 1:hparameters$num_layers){
          hparameters[[paste0("node_",i)]] <- hparameters$num_nodes[1]
        }
      }

    }
    hparameters[["num_nodes"]] <- NULL }
  return(hparameters)
}

pkg.env$deep_surv_cv <- function(IndividualDataPP,
                                 continuous_features_scaling_method,
                                 folds,
                                 kfolds,
                                 random_seed,
                                 verbose=0,
                                 epochs,
                                 num_workers,
                                 hparameters.f,
                                 out,
                                 parallel,
                                 ncores,
                                 verbose.cv=FALSE){

  "Function to perform K-fold cross-validation with xgboost"

  if(parallel == TRUE){
    # handle UNIx-operated systems seperatly?.Platform$OS.type
    require(parallel)

    cl <- makeCluster(ncores)

    objects_export <- list(
      "random_seed"
    )
    clusterExport(cl, objects_export, envir = environment())

    clusterEvalQ(cl, {library("ReSurv")
      library("fastDummies")
      library("reticulate")
      set.seed(random_seed)} )

    out[,c("train.lkh","test.lkh", "time")] <- t(parSapply(cl, 1:dim(hparameters.f)[1],  FUN =cv_deep_surv,
                                                           IndividualDataPP = IndividualDataPP,
                                                           continuous_features_scaling_method = continuous_features_scaling_method,
                                                           folds= folds,
                                                           kfolds =kfolds,
                                                           random_seed=random_seed,
                                                           verbose=verbose,
                                                           epochs=epochs,
                                                           num_workers=num_workers,
                                                           hparameters.f=hparameters.f))
    stopCluster(cl)
  }
  else{
    for(hp in 1:dim(hparameters.f)[1]){
      if(verbose.cv){cat(as.character(Sys.time()),
                         "Testing hyperparameters combination",
                         hp,
                         "out of",
                         dim(hparameters.f)[1], "\n")}


      out[hp,c("train.lkh","test.lkh", "time")] = cv_deep_surv(hp,
                                                               IndividualDataPP,
                                                               continuous_features_scaling_method,
                                                               folds,
                                                               kfolds,
                                                               random_seed,
                                                               verbose=verbose,
                                                               epochs,
                                                               num_workers,
                                                               hparameters.f)

    }
  }
  return(out)

}


