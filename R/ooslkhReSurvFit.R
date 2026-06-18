#' Compute the out-of-sample likelihood
#'
#' When the lower triangle data are available, this method computes the likelihood on the lower triangle.
#'
#' @param object \code{ReSurvFit} object.
#' @param ... Other arguments to pass to ooslkh.
#'
#' @return \code{numeric}, out-of-sample likelihood.
#'
#' @importFrom dplyr reframe
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
#' @return \code{numeric}, out-of-sample likelihood.
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
#' @return \code{numeric}, out-of-sample likelihood.
#'
#' @export
ooslkh.ReSurvFit <- function(object,
                              ...){

  # Extract quantities that you need
  if (is.null(object$data_information$data_for_reserving)) {
    stop(
      "Out-of-sample likelihood requires `object$data_information$data_for_reserving`, which is not available in this ReSurvFit object.",
      call. = FALSE
    )
  }

  data_for_reserving <- object$data_information$data_for_reserving
  fitted.model <- object$model.out
  hazard_model <- object$hazard_model
  if (is.null(hazard_model)) {
    hazard_model <- object$fit_information$hazard_model
  }
  if (is.null(hazard_model)) {
    stop("Out-of-sample likelihood requires `object$fit_information$hazard_model`.", call. = FALSE)
  }
  categorical_features <- object$data_information$categorical_features
  continuous_features <- object$data_information$continuous_features


  # Perform the computations
  test.data <- data_for_reserving %>%
      dplyr::filter(DP_rev_i <= TR_i) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(categorical_features),
                  as.factor)) %>%
    dplyr::mutate(TR_i=0)

  if (nrow(test.data) == 0L) {
    stop(
      "Out-of-sample likelihood requires lower-triangle observations in `object$data_information$data_for_reserving`.",
      call. = FALSE
    )
  }


  if(hazard_model=="COX"){

    X=test.data %>%
      dplyr::select(c(continuous_features,categorical_features))

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
      dplyr::reframe(dplyr::across(dplyr::all_of(continuous_features),
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
      dplyr::reframe(dplyr::across(dplyr::all_of(continuous_features),
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

    stop("Out-of-sample likelihood for LTRCtrees is not yet implemented.")

  }

  if(!hazard_model %in% c("COX", "NN", "XGB", "LTRCtrees")){
    stop("Unsupported hazard model for out-of-sample likelihood: ", hazard_model, call. = FALSE)
  }

  return(lkh)

}









