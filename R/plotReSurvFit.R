#' Plot for machine learning models
#'
#' This function plots the mean absolute SHAP values for the ReSurv fits of machine learning models.
#'
#' @param x \code{ReSurvFit} x.
#' @param nsamples \code{integer}, number of observations to sample for neural networks features importance plot.
#' @param ... Other arguments to be passed to plot.
#'
#' @return \code{ggplot2} of the SHAP values for an \code{"XGB"} model or a \code{"NN"} model.
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#'
#' @export
#' @method plot ReSurvFit
plot.ReSurvFit <- function(x,
                           nsamples=NULL,
                           ...
                           ){

  hazard_model <- x$hazard_model
  if (is.null(hazard_model)) {
    hazard_model <- x$fit_information$hazard_model
  }
  if (is.null(hazard_model)) {
    stop("Cannot plot this ReSurvFit object: missing `fit_information$hazard_model`.", call. = FALSE)
  }

  output.fit <- x$model.out

  if(hazard_model=="XGB"){

    if(!requireNamespace("SHAPforxgboost", quietly = TRUE))
      stop("Package 'SHAPforxgboost' is required for XGB feature importance plots. Install it with install.packages('SHAPforxgboost').")

    #we need the following
    shap_values <- SHAPforxgboost::shap.values(xgb_model = output.fit$model.out,
                               X_train = as.matrix(output.fit$data))

    df.2.plot <- apply(abs(shap_values$shap_score),2,mean)
    plot.color <- "royalblue"

  }

  if(hazard_model!="XGB"){
    stop(
      "Feature importance plots are currently supported only for ReSurvFit objects fitted with hazard_model = 'XGB'.",
      call. = FALSE
    )
  }


  data.frame(value = df.2.plot) %>%
    tibble::rownames_to_column(var = "feature") %>%
  ggplot2::ggplot(ggplot2::aes(x=feature, y=value)) +
    ggplot2::geom_bar(stat = "identity", fill=plot.color) +
    ggplot2::coord_flip() +
    ggplot2::labs(title=" ",
         x="",
         y="mean(|SHAP|)") +
    ggplot2::theme_bw()



}
