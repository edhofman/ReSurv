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

  if(hazard_model=="NN"){
    stop("Feature importance plots for NN models require the 'shap' Python package via reticulate, which is no longer the default backend. Use ReSurv with hazard_model='XGB' for SHAP-based feature importance.")
  }


  data.frame(value = df.2.plot) %>%
    rownames_to_column(var = "feature") %>%
  ggplot(aes(x=feature, y=value)) +
    geom_bar(stat = "identity", fill=plot.color) +
    coord_flip() +
    labs(title=" ",
         x="",
         y="mean(|SHAP|)") +
    theme_bw()



}

#
# shap <- reticulate::import("shap")
#
#

# Kernel


# compute SHAP values
# explainer = shap$DeepExplainer(output.fit$model.out$predict,
#                                x_fc2)
# shap_values = explainer.shap_values(x_fc)

#
# x_fc2 <- shap$sample(x_fc,as.integer(5))
#
# x_fc2 <- reticulate::np_array(as.matrix(x_fc2), dtype = "float32")
#
#
# shap_values = explainer$shap_values(x_fc2)
#
# shap$summary_plot(shap_values[[1]],x_fc2)
#
# resurv.fit.deepsurv$model.out$model.out



