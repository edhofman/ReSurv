#' Plot for machine learning models
#'
#' This function plots the results from the ReSurv fits.
#'
#' @param object "ResurvFit" object specifying start time, end time and status.
#'
#' @return predictions
#'
#' @export
#' @method plot ReSurvFit
plot.ReSurvFit <- function(object){





}

#
# library(SHAPforxgboost)
# dataX=resurv.fit.xgb$data
# #
#
# featImp_RBNS <- xgb.importance(model=resurv.fit.xgb$model.out)
# xgb.plot.importance(featImp_RBNS, main="Feature Importance - RBNS")
#
# shap_values <- shap.values(xgb_model = resurv.fit.xgb$model.out, X_train = dataX)
# shap_values <- shap.values(xgb_model = resurv.fit.xgb$model.out, X_train = as.matrix(dataX))
# #
# shap_long <- shap.prep(xgb_model = resurv.fit.xgb$model.out, X_train = as.matrix.data.frame(dataX))
# shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)
# #
# # Return the SHAP values and ranked features by mean|SHAP|
# shap_values <- shap.values(xgb_model = xgb_RBNS_Fit, X_train = as.matrix(df.RBNS_train))
#
# # Prepare the long-format data:
# shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train =  as.matrix(df.RBNS_train))
#
# # **SHAP summary plot**
# shap.plot.summary(shap_long, dilute = nrow(df.RBNS_train)/10000)




