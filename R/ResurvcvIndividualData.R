#' K fold cross-validation of a \code{ReSurv} model.
#'
#' This function computes a K fold cross-validation of a pre-specified machine learning model supported from the \code{ReSurv} package for a given grid of hyperparameters.
#' The hyperparameters to be tested are provided in a list, namely \code{hparameters_grid}.
#' Conversely, the parameters for the models run are provided separately as arguments and they are specific for each machine learning model support from.
#'
#'
#' @param IndividualData \code{IndividualData} object to use for the \code{ReSurv} fit cross-validation.
#' @param hparameters_grid \code{list}, grid of the hyperparameters to cross-validate.
#' @param folds \code{integer}, number of folds (i.e. K).
#' @param random_seed \code{integer}, random seed for making the code reproducible.
#' @param continuous_features_scaling_method \code{character}, method for scaling continuous features.
#' @param print_every_n \code{integer}, specific to the \code{xgboost} approach, see \code{xgboost::xgb.train} documentation.
#' @param early_stopping_rounds \code{integer}, specific to the \code{xgboost} approach, see \code{xgboost::xgb.train} documentation.
#' @param epochs \code{integer}, specific to the \code{deepsurv} approach, epochs to be checked.
#' @param parallel \code{logical}, specific to the \code{deepsurv} approach, whether to use parallel computing.
#' @param num_workers \code{numeric}, number of workers for the \code{deepsurv} approach, multi-process data loading with the specified number of loader worker processes.
#' @param verbose \code{logical}, whether messages from the machine learning models must be printed.
#' @param verbose.cv \code{logical}, whether messages from cross-validation must be printed.
#' @return Best \code{ReSurv} model fit. The output is different depending on the machine learning approach that is required for cross-validation.
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
                     continuous_features_scaling_method="minmax",
                     print_every_n = 1L,
                     nrounds= NULL,
                     early_stopping_rounds = NULL,
                     epochs=1,
                     parallel=F,
                     ncores = 1,
                     num_workers  =0,
                     verbose.cv=F,
                     verbose=F){

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
                             continuous_features_scaling_method="minmax",
                             print_every_n = 1L,
                             nrounds= NULL,
                             early_stopping_rounds = NULL,
                             epochs=1,
                             parallel=F,
                             ncores = 1,
                             num_workers  =0,
                             verbose = F,
                             verbose.cv){

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
                                  continuous_features_scaling_method="minmax",
                                  print_every_n = 1L,
                                  nrounds= NULL,
                                  epochs=NULL,
                                  num_workers = 0,
                                  parallel=F,
                                  ncores = 1,
                                  verbose.cv=F,
                                  verbose = F,
                                  early_stopping_rounds = NULL){


  set.seed(random_seed)

  kfolds <- sample(1:folds,size=nrow(IndividualData$training.data),
                   replace=TRUE,
                   prob=rep(1/folds,folds))

  #Allow for different number og nodes pr. layer - in blind gridsearch this is too computational expensive.
  #hparameters_grid <- pkg.env$nn_hparameter_nodes_grid(hparameters_grid, cv=T)

  hparameters.f <- expand.grid(hparameters_grid,
                               KEEP.OUT.ATTRS = FALSE,
                               stringsAsFactors = FALSE)

  hparameters.f <-  pkg.env$nn_hparameter_nodes_grid(hparameters.f, cv=T)

  if(model == "LTRCtrees"){

    formula_ct <- as.formula(IndividualData$string_formula_i)

    out.cv <- pkg.env$ltrcart_cv(IndividualData=IndividualData,
                                 folds=folds,
                                 formula_ct=formula_ct,
                                 hparameters.f=hparameters.f,
                                 verbose.cv=verbose.cv)

    # Take the best result oos
    out.best.oos <- out.cv %>%
      filter(is_lkh==min(is_lkh)) %>%
      as.data.frame()

    # List the output of the cv and the best result OOS
    out <- list(
      out.cv = out.cv,
      out.cv.best.oos = out.best.oos
    )

    class(out) <- c('ReSurvCV')

    return(out)

  }

  train.lkh <- vector("numeric",
                      length=dim(hparameters.f)[1])

  test.lkh <- vector("numeric",
                     length=dim(hparameters.f)[1])
  time <- vector("numeric",
                 length=dim(hparameters.f)[1])

  out <- cbind(hparameters.f,
               train.lkh,
               test.lkh,
               time)


  if(model == "xgboost"){

    out.cv <- pkg.env$xgboost_cv(IndividualData,
                              folds,
                              kfolds,
                              random_seed = random_seed,
                              print_every_n = print_every_n,
                              nrounds= nrounds,
                              early_stopping_rounds = early_stopping_rounds,
                              hparameters.f,
                              out,
                              parallel=parallel,
                              ncores=ncores,
                              verbose=verbose,
                              verbose.cv=verbose.cv)



  }
  if(model == "deepsurv"){
    out.cv <- pkg.env$deep_surv_cv(IndividualData,
                                 continuous_features_scaling_method = continuous_features_scaling_method,
                                 folds,
                                 kfolds,
                                 random_seed = random_seed,
                                 hparameters.f,
                                 epochs = epochs,
                                 num_workers = num_workers,
                                 out,
                                 parallel=parallel,
                                 ncores=ncores,
                                 verbose=verbose,
                                 verbose.cv=verbose.cv)

  }

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




