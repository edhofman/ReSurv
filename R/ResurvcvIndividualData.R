#' K fold cross-validation of ReSurv model.
#'
#' This function computes a K fold cross-validation of a pre-specified ReSurv model for a given grid of parameters.
#'
#' @param IndividualData IndividualData object to use for the ReSurv fit cross-validation.
#' @param hparameters_grid grid of the hyperparameters to cross-validate.
#' @param folds number of folds
#'
#' @return Best ReSurv model fit.
#'
#' @import reticulate
#' @import tidyverse
#' @import xgboost
#'
#' @export
ReSurvCV <- function(IndividualData,
                     model,
                     hparameters_grid,
                     folds,
                     random_seed,
                     print_every_n = NULL,
                     nrounds= NULL,
                     early_stopping_rounds = NULL){

  UseMethod("ReSurvCV")

}

#' K fold cross-validation of ReSurv model.
#'
#' This function computes a K fold cross-validation of a pre-specified ReSurv model for a given grid of parameters.
#'
#' @param IndividualData IndividualData object to use for the ReSurv fit cross-validation.
#' @param hparameters_grid grid of the hyperparameters to cross-validate.
#' @param folds number of folds
#'
#' @return Best ReSurv model fit.
#'
#' @export
ReSurvCV.default <- function(IndividualData,
                             model,
                             hparameters_grid,
                             folds,
                             random_seed,
                             print_every_n = NULL,
                             nrounds= NULL,
                             early_stopping_rounds = NULL){

  message('The object provided must be of class IndividualData')

}

#' K fold cross-validation of ReSurv model.
#'
#' This function computes a K fold cross-validation of a pre-specified ReSurv model for a given grid of parameters.
#'
#' @param IndividualData IndividualData object to use for the ReSurv fit cross-validation.
#' @param hparameters_grid grid of the hyperparameters to cross-validate.
#' @param folds number of folds.
#' @param random_seed random seed to make the results replicable.
#'
#' @return Best ReSurv model fit.
#'
#' @export
ReSurvCV.IndividualData <- function(IndividualData,
                                  model,
                                  hparameters_grid,
                                  folds,
                                  random_seed,
                                  print_every_n = NULL,
                                  nrounds= NULL,
                                  early_stopping_rounds = NULL){


  set.seed(random_seed)

  kfolds <- sample(1:folds,size=nrow(IndividualData$training.data),
                   replace=TRUE,
                   prob=rep(1/folds,folds))


  hparameters.f <- expand.grid(hparameters_grid,
                               KEEP.OUT.ATTRS = FALSE)

  train.lkh <- vector("numeric",
                      length=dim(hparameters.f)[1])

  test.lkh <- vector("numeric",
                     length=dim(hparameters.f)[1])

  out <- cbind(hparameters.f,
               train.lkh,
               test.lkh)


  if(model == "xgboost"){

    out.cv <- pkg.env$xgboost_cv(IndividualData,
                              folds,
                              kfolds,
                              print_every_n = print_every_n,
                              nrounds= nrounds,
                              early_stopping_rounds = early_stopping_rounds,
                              hparameters.f,
                              out)



  }

  # for(hp in 1:dim(hparameters.f)[1]){
  #
  #   hparameters <- list(params=as.list.data.frame(hparameters.f[hp,]),
  #                       print_every_n=print_every_n,
  #                       nrounds=nrounds,
  #                       early_stopping_rounds=early_stopping_rounds)
  #
  #   tmp.train.lkh <- vector("numeric",length=folds)
  #   tmp.test.lkh <- vector("numeric",length=folds)
  #
  #   for(i in c(1:folds)){
  #
  #     X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
  #                                       select_columns = IndividualData$categorical_features,
  #                                       remove_first_dummy=T)
  #
  #     scaler <- pkg.env$scaler(continuous_features_scaling_method = "minmax")
  #
  #     Xc <- IndividualData$training.data %>%
  #       summarize(across(all_of(IndividualData$continuous_features),
  #                        scaler))
  #
  #     X=cbind(X,Xc)
  #
  #     Y=individual_data$training.data[,c("DP_rev_i", "I", "TR_i")]
  #
  #     datads_pp =  pkg.env$xgboost_pp(X,
  #                                     Y,
  #                                     samples_TF= c(kfolds!=i))
  #
  #
  #
  #     model.out.k <- do.call(pkg.env$fit_xgboost, list(datads_pp=datads_pp,
  #                                                      hparameters=hparameters))
  #
  #
  #     best.it <- model.out.k$best_iteration
  #     tmp.train.lkh[i] <- model.out.k$evaluation_log$`train_log_partial likelihood`[best.it]
  #     tmp.test.lkh[i] <- model.out.k$evaluation_log$`eval_log_partial likelihood`[best.it]
  #
  #
  #
  #   }
  #
  #
  #   out[hp,c("train.lkh","test.lkh")] = c(mean(tmp.train.lkh),mean(tmp.test.lkh))
  #
  #   }
  #

  # Take the best result oos
  out.best.oos <- out.cv %>%
    filter(test.lkh==min(test.lkh)) %>%
    as.data.frame()

  # List the output of the cv and the best result OOS
  out <- list(
    out.cv = out.cv,
    out.cv.best.oos = out.best.oos
  )

  class(out) <- c('ReSurvCV')

  return(out)
}




