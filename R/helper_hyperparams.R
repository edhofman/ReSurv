# Internal hyperparameter and forecast-frame helpers
#
# These functions normalize fitting hyperparameters and prepare forecast design
# frames for the Cox, XGBoost, and NN model paths.

## hyperparameters and prepare data for fitting ----

spline_hp <- function(hparameters,IndividualDataPP){
  "
  Returns spline hyperparameters in case they are not provided from the user.

  "
  if(length(hparameters)>0){
    tmp <- list()

    tmp$nk <- ifelse(is.null(hparameters$nk),nrow(IndividualDataPP$training.data)/4,hparameters$nk)

    tmp$nbin <- ifelse(is.null(hparameters$nbin),NULL,hparameters$nbin)

    tmp$phi <- ifelse(is.null(hparameters$phi),NULL,hparameters$phi)

    return(tmp)
  }
}


simplified_df_2_fcst<- function(IndividualDataPP,
                                hazard_model){

  if (!is.null(IndividualDataPP$data_information)) {
    cont_f <- IndividualDataPP$data_information$continuous_features
    cat_f <- IndividualDataPP$data_information$categorical_features
    years <- IndividualDataPP$data_information$years
    input_time_granularity <- IndividualDataPP$data_information$input_time_granularity
    calendar_period_extrapolation <- IndividualDataPP$data_information$calendar_period_extrapolation
  } else {
    cont_f <- IndividualDataPP$continuous_features
    cat_f <- IndividualDataPP$categorical_features
    years <- IndividualDataPP$years
    input_time_granularity <- IndividualDataPP$input_time_granularity
    calendar_period_extrapolation <- IndividualDataPP$calendar_period_extrapolation
  }

  time_features <- c("DP_i", "DP_rev_i", "RP_i")
  columns_for_grouping <- unique(c(
    setdiff(cont_f, time_features),
    setdiff(cat_f, time_features),
    "AP_i"
  ))

  tmp <- as.data.table(IndividualDataPP$training.data)

  out <- tmp[,.(.N),by=columns_for_grouping][,..columns_for_grouping]

  l4 <- list()

  l4$DP_rev_i <- min(IndividualDataPP$training.data[['DP_rev_i']]):max(IndividualDataPP$training.data[['DP_rev_i']])

  l4<-do.call(CJ, c(l4, sorted = FALSE))

  out <- setkey(out[, c(k = 1, .SD)], k)[
    l4[, c(k = 1, .SD)],
    allow.cartesian = TRUE
  ][, k := NULL]

  max_dp_i <- maximum.time(years, input_time_granularity)
  out[, DP_i := max_dp_i - DP_rev_i + 1L]

  if (isTRUE(calendar_period_extrapolation) || "RP_i" %in% c(cont_f, cat_f)) {
    out[, RP_i := AP_i + DP_i - 1L]
  }

  if (!is.null(cat_f)) {
    time_cat_f <- intersect(cat_f, time_features)
    for (cft in time_cat_f) {
      out[[cft]] <- factor(out[[cft]], levels = levels(IndividualDataPP$training.data[[cft]]))
    }
  }

  return(as.data.frame(out))


}


create.df.2.fcst <- function(IndividualDataPP,
                             hazard_model){

  simplified_df_2_fcst(
    IndividualDataPP = IndividualDataPP,
    hazard_model = hazard_model
  )

}



df.2.fcst.nn.pp <- function(data,
                                    newdata,
                                    continuous_features,
                                    categorical_features){

  tmp <- newdata[continuous_features]

  setDT(tmp)
  setDT(data)

  for(cft in continuous_features){

    mnv <- min(data[[cft]])
    mxv <- max(data[[cft]])

    tmp[[cft]] <-2*(tmp[[cft]]-mnv)/(mxv-mnv)-1

  }

  Xc=as.matrix(tmp)

  if(!is.null(categorical_features)){

    X=model.matrix.creator(data= newdata,
                                   select_columns = categorical_features)

    out <- cbind(X,Xc)

    }else{

    out <- Xc

    }

  return(out)

}


df.2.fcst.xgboost.pp <- function(data,
                                         newdata,
                                         continuous_features,
                                         categorical_features){
  tmp <- Xc <- NULL

  if(!is.null(continuous_features)){
    tmp <- newdata[continuous_features]

    for(cft in continuous_features){

      mnv <- min(data[[cft]])
      mxv <- max(data[[cft]])

      tmp[,cft] <-2*(tmp[,cft]-mnv)/(mxv-mnv)-1

    }
    Xc=as.matrix(tmp)

  }


  if(!is.null(categorical_features)){

    X=model.matrix.creator(data= newdata,
                                   select_columns = categorical_features,
                                   remove_first_dummy = TRUE)
  }


  if(!is.null(Xc)){

    if(!is.null(categorical_features)){

    X <- cbind(X,Xc)}else{

      X <- Xc

    }}

  ds_train_fcst <- xgboost::xgb.DMatrix(as.matrix(X), label=rep(1, dim(X)[1]))

  return(ds_train_fcst)

}
