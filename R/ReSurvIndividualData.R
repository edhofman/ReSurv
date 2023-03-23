#' Fit individual chain ladder plus models.
#'
#' This function fits and computes the reserves for the ReSurv models.
#'
#' @param IndividualData IndividualData object to use for the ReSurv fit.
#' @param hazard_model hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{'cox': Standard Cox model for the hazard.}
#' \item{'deep_surv': Deep Survival Neural Network.}
#' \item{'xgboost': Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed random seed set for reproducibility
#' @param hparameters list of hyperparameters for 'deep_surv' and 'xgboost'. It will be disregarded for 'cox'.
#' @param percentage_data_training percentage of data used for training on the upper triangle.
#'
#'
#' @return ReSurv fit.
#'
#' @import reticulate
#' @import tidyverse
#'
#' @export
ReSurv <- function(IndividualData,
                   hazard_model="cox",
                   tie='efron',
                   baseline="spline",
                   continuous_features_scaling_method="minmax",
                   random_seed=1,
                   hparameters=list(),
                   percentage_data_training=.8
){

  UseMethod("ReSurv")

}

#' Fit individual chain ladder plus models.
#'
#' This function fits and computes the reserves for the ReSurv models.
#'
#' @param IndividualData IndividualData object to use for the ReSurv fit.
#' @param hazard_model hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{'cox': Standard Cox model for the hazard.}
#' \item{'deep_surv': Deep Survival Neural Network.}
#' \item{'xgboost': Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed random seed set for reproducibility
#' @param hparameters list of hyperparameters for 'deep_surv' and 'xgboost'. It will be disregarded for 'cox'.
#' @param percentage_data_training percentage of data used for training on the upper triangle.
#'
#' @return ReSurv fit.
#' @export
ReSurv.default <- function(IndividualData,
                           hazard_model="cox",
                           tie='efron',
                           baseline="spline",
                           continuous_features_scaling_method="minmax",
                           random_seed=1,
                           hparameters=list(),
                           percentage_data_training=.8){

  message('The object provided must be of class IndividualData')

}



#' Fit chain-ladder+ to reverse time triangles.
#'
#' This function fits and computes the reserves for the ReSurv models.
#'
#' @param IndividualData IndividualData object to use for the ReSurv fit.
#' @param hazard_model hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{'cox': Standard Cox model for the hazard.}
#' \item{'deep_surv': Deep Survival Neural Network.}
#' \item{'xgboost': Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed random seed set for reproducibility
#' @param hparameters list of hyperparameters for 'deep_surv' and 'xgboost'. It will be disregarded for 'cox'.
#' @param percentage_data_training percentage of data used for training on the upper triangle.
#'
#' @return ReSurv fit.
#' @export
ReSurv.IndividualData <- function(IndividualData,
                               hazard_model="cox",
                               tie='efron',
                               baseline="spline",
                               continuous_features_scaling_method="minmax",
                               random_seed=1,
                               hparameters=list(),
                               percentage_data_training=.8
                               ){


  set.seed(random_seed)

  formula_ct <- as.formula(IndividualData$string_formula_i)

  # X_i <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
  #                                       select_columns = c('AP_i',IndividualData$categorical_features))

  # hz_names_i = pkg.env$model.matrix.extract.hazard.names(X=X_i,
  #                                                        string_formula=IndividualData$string_formula_i,
  #                                                        data=IndividualData$training.data)

  newdata <- pkg.env$create.df.2.fcst(IndividualData)

  if(hazard_model=="cox"){

    model.out <- pkg.env$fit_cox_model(data=IndividualData$training.data,
                                       formula_ct=formula_ct,
                                       newdata=newdata)

    tmp <- pkg.env$spline_hp(hparameters,IndividualData)

    baseline_out <- pkg.env$hazard_baseline_model(data=IndividualData$training.data,
                                                  cox=model.out$cox,
                                                  hazard=NULL,
                                                  baseline=baseline,
                                                  conversion_factor=IndividualData$conversion_factor,
                                                  nk=tmp$nk,
                                                  nbin=tmp$nbin,
                                                  phi=tmp$phi)

    bsln <- data.frame(baseline=baseline_out$bs_hazard$hazard,
                       DP_rev_i=ceiling(baseline_out$bs_hazard$time))  #$hazard
    # colnames(bsln)=c('baseline','DP_rev_i')

    # expg <- exp(model.out$beta_ams)

    hazard_frame <- cbind(newdata, model.out$expg)
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"

    # hazard_frame <- hazard_frame %>%
    #   full_join(bsln,
    #             by="DP_rev_i") %>%
    #   as.data.frame()
    #
    # hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']
    #
    # return(hazard_frame)

    }

  if(hazard_model=="deepsurv"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features)

    scaler <- pkg.env$scaler(continuous_features_scaling_method=continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      summarize(across(all_of(IndividualData$continuous_features),
                       scaler))

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    datads_pp = pkg.env$deep_surv_pp(X=cbind(X,Xc),
                           Y=IndividualData$training.data[,c("DP_rev_i", "I", "TR_i")],
                           training_test_split = training_test_split)


    model.out <- pkg.env$fit_deep_surv(datads_pp,
                                       hparameters=hparameters)



    bsln <- model.out$compute_baseline_hazards(
      input = datads_pp$x_train,
      target = datads_pp$y_train,
      batch_size = hparameters$batch_size)


    newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=IndividualData$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualData$continuous_features,
                                          categorical_features=IndividualData$categorical_features)

    x_fc= reticulate::np_array(as.matrix(newdata.mx), dtype = "float32")

    beta_ams <- model.out$predict(input=x_fc,
                                  batch_size=hparameters$batch_size,
                                  num_workers=hparameters$num_workers)

    # X_ams <- cbind(X_i, Xb)
    #
    # beta_ams = unique(round(X_ams,10) )[,ncol(X_ams)] #if no round some systems has too high precision.
    #

    expg <- exp(beta_ams)

    hazard_frame <- cbind(newdata,expg)
    bsln <- data.frame(baseline=bsln, DP_rev_i=as.integer(names(bsln)))

  }


  ##################################################################################


  hazard_frame <- hazard_frame %>%
    full_join(bsln,
              by="DP_rev_i") %>%
    as.data.frame()

  hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']

  return(hazard_frame)

  #need placeholder for latest i mirror cl behaviour
  tmp= sapply(1:l,function(x) bsln*expg[x])

  hazard_cl <- (sapply(seq_along(hz_names_i$time),
                       pkg.env$hazard_f,
                       enter= hz_names_i$enter,
                       time=hz_names_i$time,
                       exit=hz_names_i$exit,
                       event=hz_names_i$event))

  hazard = as.matrix(cbind(hazard_cl,
                           tmp ))

  colnames(hazard) <- c("CL",hz_names_i$names_hazard )

  max_DP <- max(IndividualData$training$DP_rev_o)

  ############################################################
  #check
  hazard_q <- matrix(nrow=max_DP, ncol=(ncol(hazard)-1)*IndividualData$conversion_factor)
  ############################################################


  #group to quarters, this is relatively time consuming
  for( i in 1:max_DP){ #Loop through each output period, to find weights
    frame_tmp <- data.table(IndividualData$training) %>% filter(TR_o<i) %>% #All claims that haven't been truncated at said reverse development
      filter(DP_rev_o>=i) %>% #Claims that still hasn't been reported
      mutate(time_w = round(ifelse(DP_rev_o==i, DP_rev_i, AP_i-1+1/(IndividualData$conversion_factor)*(i-AP_o+1) ),10) ) %>% #If a claim is reported in the corresponding development period save said reporting time, otherwise we need the corresponding limit for each acciedent period in the development period.
      mutate(weight = ifelse(DP_rev_o==i, (DP_rev_i-1)%%(1/(IndividualData$conversion_factor))+1, 1/(IndividualData$conversion_factor) )) %>% #If reported in said period give weight corresponding to amount of input_time period spend in output time_period, otherwise give width length of output as weight.
      mutate(p_month = (AP_i-1)%%(1/(IndividualData$conversion_factor))+1) %>% #Entering month in development period
      group_by(p_month, time_w) %>%
      summarize(weight=sum(weight*I) )

    #frame_tmp2 <- frame_tmp %>%  reshape2::dcast(time_w ~p_month, value.var="weight") #create parrelogram structure for output development period, that hold weights for each input development period

    hazard_data_frame <- pkg.env$hazard_data_frame(hazard=hazard,
                                                   conversion_factor = IndividualData$conversion_factor)

    # hazard_q[i,] <- mapply(pkg.env$m_to_q_hazard,
    #                        1:ncol(hazard_q),
    #                        MoreArgs=list(hazard=hazard,
    #                        frame_tmp=frame_tmp,
    #                        frame_tmp2=frame_tmp2,
    #                        conversion_factor = IndividualData$conversion_factor))

    hazard_q[i,] <- mapply(pkg.env$m_to_q_hazard_2,
                           2:max(hazard_data_frame$group), #start from 2 since group 1 is chain ladder
                           MoreArgs=list(hazard_data_frame=hazard_data_frame,
                                         frame_tmp=frame_tmp,
                                         conversion_factor = IndividualData$conversion_factor))
    }

  #create monthly development factors, though not outputted later on.


  df_monthly <- (2+hazard)/(2-hazard)

  #Remove last row,since doesn't make sense
  df_monthly <- as.data.frame(df_monthly[1:(nrow(df_monthly)-1),]) %>%
    map_df(rev) %>%
    mutate(DM=row_number())

  X_o <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                        select_columns = c('AP_o',IndividualData$categorical_features))

  hz_names_o = pkg.env$model.matrix.extract.hazard.names(X=X_o,
                                                         string_formula=IndividualData$string_formula_o,
                                                         data=IndividualData$training.data)

hazard_cl <- (sapply(seq_along(hz_names_o$time),
                     pkg.env$hazard_f,
                     enter= hz_names_o$enter,
                     time= hz_names_o$time,
                     exit= hz_names_o$exit,
                     event= hz_names_o$event))

q_hazard = as.matrix(cbind(hazard_cl, hazard_q )   )

colnames(q_hazard) <- c("CL",hz_names_o$names_hazard )


df_quarterly <- (2+q_hazard)/(2-q_hazard)


df_quarterly <- as.data.frame(df_quarterly[1:(nrow(df_quarterly)-1),]) %>%
  map_df(rev) %>%
  mutate(DP_i=row_number())

out=list(df = df_quarterly,
         hazard= q_hazard,
         IndividualData=IndividualData)

class(out) <- c('ReSurvFit')

out

}




