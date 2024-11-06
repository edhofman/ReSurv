#' K fold cross-validation of a \code{ReSurv} model.
#'
#' This function computes a K fold cross-validation of a pre-specified machine learning model supported from the \code{ReSurv} package for a given grid of hyperparameters.
#' The hyperparameters to be tested are provided in a list, namely \code{hparameters_grid}.
#' Conversely, the parameters for the models run are provided separately as arguments and they are specific for each machine learning model support from.
#'
#'
#' @param IndividualDataPP \code{IndividualDataPP} object to use for the \code{ReSurv} fit cross-validation.
#' @param model \code{character}, machine learning for cross validation.
#' @param hparameters_grid \code{list}, grid of the hyperparameters to cross-validate.
#' @param folds \code{integer}, number of folds (i.e. K).
#' @param random_seed \code{integer}, random seed for making the code reproducible.
#' @param continuous_features_scaling_method \code{character}, method for scaling continuous features.
#' @param print_every_n \code{integer}, specific to the \code{XGB} approach, see \code{xgboost::xgb.train} documentation.
#' @param early_stopping_rounds \code{integer}, specific to the \code{XGB} approach, see \code{xgboost::xgb.train} documentation.
#' @param epochs \code{integer}, specific to the \code{NN} approach, epochs to be checked.
#' @param parallel \code{logical}, specific to the \code{NN} approach, whether to use parallel computing.
#' @param num_workers \code{numeric}, number of workers for the \code{NN} approach, multi-process data loading with the specified number of loader worker processes.
#' @param verbose \code{logical}, whether messages from the machine learning models must be printed.
#' @param verbose.cv \code{logical}, whether messages from cross-validation must be printed.
#' @param nrounds \code{integer}, specific to \code{XGB}, max number of boosting iterations.
#' @param ncores \code{integer}, specific to \code{NN}, max number of cores used.
#'
#' @return Best \code{ReSurv} model fit. The output is different depending on the machine learning approach that is required for cross-validation. A list containing:
#'  \itemize{
#' \item{\code{out.cv}: \code{data.frame}, total output of the cross-validation (all the input parameters combinations). }
#' \item{\code{out.cv.best.oos}:  \code{data.frame}, combination with the best out of sample likelihood. }
#' }
#'
#' For XGB the columns in \code{out.cv} and \code{out.cv.best.oos} are the hyperparameters \code{booster}, \code{eta}, \code{max_depth}, \code{subsample}, \code{alpha}, \code{lambda}, \code{min_child_weight}. They also contain the metrics \code{train.lkh}, \code{test.lkh}, and the computational time \code{time}. For NN the columns in \code{out.cv} and \code{out.cv.best.oos} are the hyperparameters \code{num_layers}, \code{optim}, \code{activation}, \code{lr}, \code{xi}, \code{eps}, \code{tie}, \code{batch_size}, \code{early_stopping}, \code{patience}, \code{node} train.lkh test.lkh. They also contain the metrics \code{train.lkh}, \code{test.lkh}, and the computational time \code{time}.
#'
#' @import reticulate
#' @import tidyverse
#' @import xgboost
#'
#' @examples
#' ## Not run
#' input_data <- data_generator(random_seed = 1964)
#'
#' individual_data <- IndividualDataPP(input_data,
#'                                   id="claim_number",
#'                                   continuous_features=NULL,
#'                                   categorical_features="claim_type",
#'                                   accident_period="AP",
#'                                   calendar_period="RP",
#'                                   input_time_granularity = "months",
#'                                   output_time_granularity = "quarters",
#'                                   years=4,
#'                                   continuous_features_spline=NULL,
#'                                   calendar_period_extrapolation=FALSE)
#'
#' resurv.cv.xgboost <- ReSurvCV(IndividualDataPP=individual_data,
#'                               model="XGB",
#'                               hparameters_grid=list(booster="gbtree",
#'                               eta=c(.001,.01,.2,.3),
#'                               max_depth=c(3,6,8),
#'                               subsample=c(1),
#'                               alpha=c(0,.2,1),
#'                               lambda=c(0,.2,1),
#'                                min_child_weight=c(.5,1)),
#'                                print_every_n = 1L,
#'                                nrounds=1L, ##set to one to run quickly
#'                                verbose=FALSE,
#'                                verbose.cv=TRUE,
#'                                early_stopping_rounds = 100L,
#'                                folds=5L,
#'                                parallel=TRUE,
#'                                ncores=2L,
#'                                random_seed=1L)
#'
#'
#' @references
#' Munir, H., Emil, H., & Gabriele, P. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' @export
ReSurvCV <- function(IndividualDataPP,
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
                     verbose=F,
                     verbose.cv=F){

  UseMethod("ReSurvCV")

}

#' K fold cross-validation of ReSurv model.
#'
#' This function computes a K fold cross-validation of a pre-specified ReSurv model for a given grid of parameters.
#'
#' @param IndividualDataPP \code{IndividualDataPP} object to use for the \code{ReSurv} fit cross-validation.
#' @param model \code{character}, machine learning for cross validation.
#' @param hparameters_grid \code{list}, grid of the hyperparameters to cross-validate.
#' @param folds \code{integer}, number of folds (i.e. K).
#' @param random_seed \code{integer}, random seed for making the code reproducible.
#' @param continuous_features_scaling_method \code{character}, method for scaling continuous features.
#' @param print_every_n \code{integer}, specific to the \code{XGB} approach, see \code{xgboost::xgb.train} documentation.
#' @param early_stopping_rounds \code{integer}, specific to the \code{XGB} approach, see \code{xgboost::xgb.train} documentation.
#' @param epochs \code{integer}, specific to the \code{NN} approach, epochs to be checked.
#' @param parallel \code{logical}, specific to the \code{NN} approach, whether to use parallel computing.
#' @param num_workers \code{numeric}, number of workers for the \code{NN} approach, multi-process data loading with the specified number of loader worker processes.
#' @param verbose \code{logical}, whether messages from the machine learning models must be printed.
#' @param verbose.cv \code{logical}, whether messages from cross-validation must be printed.
#' @param nrounds \code{integer}, specific to \code{XGB}, max number of boosting iterations.
#' @param ncores \code{integer}, specific to \code{NN}, max number of cores used.
#'
#' @return Best \code{ReSurv} model fit. The output is different depending on the machine learning approach that is required for cross-validation. A list containing:
#'  \itemize{
#' \item{\code{out.cv}: \code{data.frame}, total output of the cross-validation (all the input parameters combinations). }
#' \item{\code{out.cv.best.oos}:  \code{data.frame}, combination with the best out of sample likelihood. }
#' }
#'
#' For XGB the columns in \code{out.cv} and \code{out.cv.best.oos} are the hyperparameters \code{booster}, \code{eta}, \code{max_depth}, \code{subsample}, \code{alpha}, \code{lambda}, \code{min_child_weight}. They also contain the metrics \code{train.lkh}, \code{test.lkh}, and the computational time \code{time}. For NN the columns in \code{out.cv} and \code{out.cv.best.oos} are the hyperparameters \code{num_layers}, \code{optim}, \code{activation}, \code{lr}, \code{xi}, \code{eps}, \code{tie}, \code{batch_size}, \code{early_stopping}, \code{patience}, \code{node} train.lkh test.lkh. They also contain the metrics \code{train.lkh}, \code{test.lkh}, and the computational time \code{time}.
#'
#' @references
#' Munir, H., Emil, H., & Gabriele, P. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' @export
ReSurvCV.default <- function(IndividualDataPP,
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
                             verbose.cv=F){

  message('The object provided must be of class IndividualDataPP')

}

#' K fold cross-validation of ReSurv model.
#'
#' This function computes a K fold cross-validation of a pre-specified ReSurv model for a given grid of parameters.
#'
#' @param IndividualDataPP \code{IndividualDataPP} object to use for the \code{ReSurv} fit cross-validation.
#' @param model \code{character}, machine learning for cross validation.
#' @param hparameters_grid \code{list}, grid of the hyperparameters to cross-validate.
#' @param folds \code{integer}, number of folds (i.e. K).
#' @param random_seed \code{integer}, random seed for making the code reproducible.
#' @param continuous_features_scaling_method \code{character}, method for scaling continuous features.
#' @param print_every_n \code{integer}, specific to the \code{XGB} approach, see \code{xgboost::xgb.train} documentation.
#' @param early_stopping_rounds \code{integer}, specific to the \code{XGB} approach, see \code{xgboost::xgb.train} documentation.
#' @param epochs \code{integer}, specific to the \code{NN} approach, epochs to be checked.
#' @param parallel \code{logical}, specific to the \code{NN} approach, whether to use parallel computing.
#' @param num_workers \code{numeric}, number of workers for the \code{NN} approach, multi-process data loading with the specified number of loader worker processes.
#' @param verbose \code{logical}, whether messages from the machine learning models must be printed.
#' @param verbose.cv \code{logical}, whether messages from cross-validation must be printed.
#' @param nrounds \code{integer}, specific to \code{XGB}, max number of boosting iterations.
#' @param ncores \code{integer}, specific to \code{NN}, max number of cores used.
#'
#' @return Best \code{ReSurv} model fit. The output is different depending on the machine learning approach that is required for cross-validation. A list containing:
#'  \itemize{
#' \item{\code{out.cv}: \code{data.frame}, total output of the cross-validation (all the input parameters combinations). }
#' \item{\code{out.cv.best.oos}:  \code{data.frame}, combination with the best out of sample likelihood. }
#' }
#'
#' For XGB the columns in \code{out.cv} and \code{out.cv.best.oos} are the hyperparameters \code{booster}, \code{eta}, \code{max_depth}, \code{subsample}, \code{alpha}, \code{lambda}, \code{min_child_weight}. They also contain the metrics \code{train.lkh}, \code{test.lkh}, and the computational time \code{time}. For NN the columns in \code{out.cv} and \code{out.cv.best.oos} are the hyperparameters \code{num_layers}, \code{optim}, \code{activation}, \code{lr}, \code{xi}, \code{eps}, \code{tie}, \code{batch_size}, \code{early_stopping}, \code{patience}, \code{node} train.lkh test.lkh. They also contain the metrics \code{train.lkh}, \code{test.lkh}, and the computational time \code{time}.
#'
#'
#' @export
ReSurvCV.IndividualDataPP <- function(IndividualDataPP,
                                  model,
                                  hparameters_grid,
                                  folds,
                                  random_seed,
                                  continuous_features_scaling_method="minmax",
                                  print_every_n = 1L,
                                  nrounds= NULL,
                                  early_stopping_rounds = NULL,
                                  epochs=NULL,
                                  parallel=F,
                                  ncores = 1,
                                  num_workers = 0,
                                  verbose = F,
                                  verbose.cv=F){


  set.seed(random_seed)

  kfolds <- sample(1:folds,size=nrow(IndividualDataPP$training.data),
                   replace=TRUE,
                   prob=rep(1/folds,folds))

  #Allow for different number og nodes pr. layer - in blind gridsearch this is too computational expensive.
  #hparameters_grid <- pkg.env$nn_hparameter_nodes_grid(hparameters_grid, cv=T)

  hparameters.f <- expand.grid(hparameters_grid,
                               KEEP.OUT.ATTRS = FALSE,
                               stringsAsFactors = FALSE)

  hparameters.f <-  pkg.env$nn_hparameter_nodes_grid(hparameters.f, cv=T)

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


  if(model == "XGB"){

    out.cv <- pkg.env$xgboost_cv(IndividualDataPP,
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
  if(model == "NN"){
    out.cv <- pkg.env$deep_surv_cv(IndividualDataPP,
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




