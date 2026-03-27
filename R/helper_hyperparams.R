## hyperparameters and prepare data for fitting ----

pkg.env$spline_hp <- function(hparameters,IndividualDataPP){
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

  cont_f <- IndividualDataPP$data_information$continuous_features
  cat_f <- IndividualDataPP$data_information$categorical_features
  columns_for_grouping <- unique(c(cont_f,cat_f,"AP_i"))

  tmp <- as.data.table(IndividualDataPP$training.data)

  out <- tmp[,.(.N),by=columns_for_grouping][,..columns_for_grouping]

  l4 <- list()

  l4$DP_rev_i <- min(IndividualDataPP$training.data[,'DP_rev_i']):max(IndividualDataPP$training.data[,'DP_rev_i'])

  l4<-do.call(CJ, c(l4, sorted = FALSE))

  out<-as.data.frame(setkey(out[,c(k=1,.SD)],k)[l4[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL])



  return(out)


}


create.df.2.fcst <- function(IndividualDataPP,
                             hazard_model){

  l1 <- lapply(IndividualDataPP$training.data %>% select(IndividualDataPP$data_information$categorical_features), levels)
  l2 <- lapply(IndividualDataPP$training.data %>% select(IndividualDataPP$data_information$continuous_features), unique)
  l3 <- list()
  l4 <- list()
  l5 <- list()

  if(!('AP_i'%in%c(IndividualDataPP$data_information$categorical_features,IndividualDataPP$data_information$continuous_features))){
    l3$AP_i <- unique(IndividualDataPP$training.data[,'AP_i'])
  }else{
    l3 <- NULL
  }

  l4$DP_rev_i <- min(IndividualDataPP$training.data[,'DP_rev_i']):max(IndividualDataPP$training.data[,'DP_rev_i'])

  # OLD
  # l1 <- as.data.table(cross_df(l1))
  # l2 <- as.data.table(cross_df(l2))
  # data.table alternative
  l1<-do.call(CJ, c(l1, sorted = FALSE))
  l2<-do.call(CJ, c(l2, sorted = FALSE))

  tmp<-setkey(l1[,c(k=1,.SD)],k)[l2[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL]

  if(!is.null(l3)){
    # OLD
    # l3 <- as.data.table(cross_df(l3))
    l3<-do.call(CJ, c(l3, sorted = FALSE))
    tmp<-setkey(tmp[,c(k=1,.SD)],k)[l3[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL]}
  # OLD
  # l4 <- as.data.table(cross_df(l4))
  l4<-do.call(CJ, c(l4, sorted = FALSE))
  tmp<-setkey(tmp[,c(k=1,.SD)],k)[l4[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL]
  tmp <- tmp %>%
    as.data.frame()
  # e2 <- Sys.time()
  # Time difference of 0.6656282 secs


  if(IndividualDataPP$data_information$calendar_period_extrapolation & (hazard_model=='COX')){
    tmp$RP_i <- tmp$AP_i+tmp$DP_rev_i-1
  }else{
    if(IndividualDataPP$data_information$calendar_period_extrapolation){
      warning("The calendar year component extrapolation is disregarded.
             The current implementation supports this feature only for the Cox model")}

  }

  return(tmp)

}



pkg.env$df.2.fcst.nn.pp <- function(data,
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

  Xc=as.matrix.data.frame(tmp)

  if(!is.null(categorical_features)){

    X=pkg.env$model.matrix.creator(data= newdata,
                                   select_columns = categorical_features)

    out <- cbind(X,Xc)

    }else{

    out <- Xc

    }

  return(out)

}


pkg.env$df.2.fcst.xgboost.pp <- function(data,
                                         newdata,
                                         continuous_features,
                                         categorical_features){
  tmp <- Xc <- NULL

  if(!is.null(continuous_features)){
    tmp <- newdata[continuous_features]

    for(cft in continuous_features){

      mnv <- min(data[cft])
      mxv <- max(data[cft])

      tmp[,cft] <-2*(tmp[,cft]-mnv)/(mxv-mnv)-1

    }
    Xc=as.matrix.data.frame(tmp)

  }


  if(!is.null(categorical_features)){

    X=pkg.env$model.matrix.creator(data= newdata,
                                   select_columns = categorical_features,
                                   remove_first_dummy = TRUE)
  }


  if(!is.null(Xc)){

    if(!is.null(categorical_features)){

    X <- cbind(X,Xc)}else{

      X <- Xc

    }}

  ds_train_fcst <- xgboost::xgb.DMatrix(as.matrix.data.frame(X), label=rep(1, dim(X)[1]))

  return(ds_train_fcst)

}
