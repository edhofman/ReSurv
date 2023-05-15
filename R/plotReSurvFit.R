#' Plot for machine learning models
#'
#' This function plots the results from the ReSurv fits.
#'
#' @param object "ResurvFit" object specifying start time, end time and status.
#'
#' @return predictions
#'
#' @import SHAPforxgboost
#'
#' @export
#' @method plot ReSurvFit
plot.ReSurvFit <- function(object){

  hazard_model <- object$hazard_model
  output.fit <- object$model.out

  if(hazard_model=="xgboost"){
    shap_long <- shap.prep(xgb_model = output.fit$model.out, X_train = as.matrix(output.fit$data))
    shap.plot.summary(shap_long)
  }

  # if(hazard_model=="deep_surv"){
  #   shap <- reticulate::import("shap")

    # shap_long <- shap.prep(xgb_model = output.fit$model.out, X_train = as.matrix(output.fit$data))
    # shap.plot.summary(shap_long)
  # }

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
# x_fc= reticulate::np_array(as.matrix( resurv.fit.deepsurv$model.out$data), dtype = "float32")
#
# explainer = shap$KernelExplainer(resurv.fit.deepsurv$model.out$model.out$predict,
#                                  x_fc)
#
#
# x_fc2 <- shap$sample(x_fc,as.integer(100))
#
# x_fc2 <- reticulate::np_array(as.matrix(x_fc2), dtype = "float32")
#
#
# shap_values = explainer$shap_values(x_fc2)
#
# shap$summary_plot(shap_values[[1]],x_fc2)
#
# resurv.fit.deepsurv$model.out$model.out




