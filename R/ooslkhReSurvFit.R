#' Compute the out-of-the sample likelihood
#'
#' When the lower triangle data are available, this method computes the likelihood on the lower triangle.
#'
#' @param object \code{ReSurvFit} object.
#' @param ... Other arguments to pass to ooslkh.
#'
#' @return \code{numeric}, out of the sample likelihood.
#'
#' @export
ooslkh <- function(object,
                    ...){

  UseMethod("ooslkh")

}

#' Compute the out-of-the sample likelihood
#'
#' When the lower triangle data are available, this method computes the likelihood on the lower triangle.
#'
#' @param object \code{ReSurvFit} object.
#' @param ... Other arguments to pass to ooslkh.
#'
#' @return \code{numeric}, out of the sample likelihood.
#'
#' @export
ooslkh.default <- function(object,
                            ...){

  message('The object provided must be of class ReSurvFit')

}

#' Compute the out-of-the sample likelihood
#'
#' When the lower triangle data are available, this method computes the likelihood on the lower triangle.
#'
#' @param object \code{ReSurvFit} object.
#' @param ... Other arguments to pass to ooslkh.
#'
#' @return \code{numeric}, out of the sample likelihood.
#'
#' @export
ooslkh.ReSurvFit <- function(object,
                              ...){

  # Extract quantities that you need
  starting.data <- object$IndividualDataPP$full.data
  fitted.model <- object$model.out
  hazard_model<- object$hazard_model
  categorical_features <- object$IndividualDataPP$categorical_features
  continuous_features <- object$IndividualDataPP$continuous_features


  # Perform the computations
  test.data <- starting.data %>%
      filter(DP_rev_i <= TR_i) %>%
    mutate(across(all_of(categorical_features),
                  as.factor)) %>%
    mutate(TR_i=0)


  if(hazard_model=="COX"){

    X=test.data %>%
      select(c(continuous_features,categorical_features))

    Y=test.data[,c("DP_rev_i", "I", "TR_i")]

    lkh <- pkg.env$evaluate_lkh_cox(X_train=X,
                                    Y_train=Y,
                                    model=fitted.model$model.out)


  }

  if(hazard_model=="NN"){

    X <- pkg.env$model.matrix.creator(data= test.data,
                                      select_columns = categorical_features)

    scaler <- pkg.env$scaler(continuous_features_scaling_method='minmax')

    Xc <- test.data %>%
      reframe(across(all_of(continuous_features),
                       scaler))

    X = cbind(X,Xc)

    Y=test.data[,c("DP_rev_i", "I", "TR_i")]

    lkh <- pkg.env$evaluate_lkh_nn(X_train=X,
                                   Y_train=Y,
                                   model=fitted.model$model.out)


  }

  if(hazard_model=="XGB"){

    X <- pkg.env$model.matrix.creator(data= test.data,
                                      select_columns = categorical_features,
                                      remove_first_dummy=T)

    scaler <- pkg.env$scaler(continuous_features_scaling_method='minmax')

    Xc <- test.data %>%
      reframe(across(all_of(continuous_features),
                     scaler))


    X = cbind(X,Xc)

    Y=test.data[,c("DP_rev_i", "I", "TR_i")]

    lkh <- pkg.env$evaluate_lkh_xgb(X_train=X,
                                    Y_train=Y,
                                    dset='is',
                                    samples_cn=data.frame(id=seq(1,dim(X)[1])),
                                    model=fitted.model$model.out)


  }

  if(hazard_model=="LTRCtrees"){

    Y=test.data[,c("DP_rev_i", "I", "TR_i")]
    lkh <- pkg.env$evaluate_lkh_LTRCtrees(X_train=test.data %>% select(c(categorical_features,
                                                                         continuous_features)),
                                             Y_train=Y,
                                             model=fitted.model$model.out)

  }

  return(lkh)

}









