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
ReSurv.cv <- function(IndividualData,
                      hparameters_grid,
                      folds){

  UseMethod("ReSurv")

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
ReSurv.default <- function(IndividualData,
                           hparameters_grid,
                           folds){

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
ReSurv.IndividualData <- function(IndividualData,
                                  model,
                                  hparameters_grid,
                                  folds,
                                  random_seed){


  set.seed(random_seed)

  class(out) <- c('ReSurvCV')

  return(out)
}




