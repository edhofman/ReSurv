#' Plot for machine learning models
#'
#' This function plots the mean absolute SHAP values for the ReSurv fits of machine learning models.
#'
#' @param object \code{ReSurvFit} object.
#' @param nsamples \code{integer}, number of observations to sample for neural networks features importance plot.
#'
#' @return predictions
#'
#' @import SHAPforxgboost
#'
#' @export
#' @method plot ReSurvFit
plot.ReSurvFit <- function(object,
                           nsamples=NULL){

  hazard_model <- object$hazard_model
  output.fit <- object$model.out

  if(hazard_model=="xgboost"){

    #we need the following
    shap_values <- shap.values(xgb_model = output.fit$model.out,
                               X_train = as.matrix(output.fit$data))

    df.2.plot <- apply(abs(shap_values$shap_score),2,mean)


  }

  if(hazard_model=="deep_surv"){

    shap <- reticulate::import("shap")

    x_fc= reticulate::np_array(as.matrix( output.fit$data), dtype = "float32")

    explainer = shap$KernelExplainer(output.fit$model.out$predict,
                                     x_fc)

    if(!is.null(nsamples)){
      x_fc <- shap$sample(x_fc,as.integer(nsamples))
      x_fc <- reticulate::np_array(as.matrix(x_fc), dtype = "float32")

      }

    shap_values = explainer$shap_values(x_fc)[[1]]
    colnames(shap_values) <- colnames(output.fit$data)

    df.2.plot <- apply(abs(shap_values),MARGIN = 2,mean)

  }


  df.2.plot %>%
    reshape2::melt(df.2.plot, na.rm = FALSE, value.name = "value", id = NULL) %>%
    rownames_to_column(var = "feature") %>%
  ggplot(aes(x=feature, y=value)) +
    geom_bar(stat = "identity", fill="navy") +
    coord_flip() +
    labs(title=" ",
         x="",
         y="mean(|SHAP|)") +
    theme_bw()



}

#
# library(SHAPforxgboost)
# dataX=resurv.fit.xgb$model.out$data
# #
#
# library(xgboost)
# featImp_RBNS <- xgb.importance(model=resurv.fit.xgb$model.out$model.out)
# xgb.plot.importance(featImp_RBNS, main="Feature Importance - RBNS")
#
# shap_values <- shap.values(xgb_model = resurv.fit.xgb$model.out$model.out, X_train = dataX)
# shap_values <- shap.values(xgb_model = resurv.fit.xgb$model.out$model.out, X_train = as.matrix(dataX))
# #
# shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)
# shap_long <- shap.prep(xgb_model = resurv.fit.xgb$model.out$model.out, X_train = as.matrix(resurv.fit.xgb$model.out$data))
#
# #
# # Return the SHAP values and ranked features by mean|SHAP|
# shap_values <- shap.values(xgb_model = xgb_RBNS_Fit, X_train = as.matrix(df.RBNS_train))
#
# # Prepare the long-format data:
# shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train =  as.matrix(df.RBNS_train))
#
# # **SHAP summary plot**
# shap.plot.summary(shap_long)

#
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




