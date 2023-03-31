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
#' @param grouping_method Use probability or exposure approach to group from input to output development factors.
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
                   percentage_data_training=.8,
                   grouping_method = "exposure"
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
#' @param grouping_method Use probability or exposure approach to group from input to output development factors.
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
                           percentage_data_training=.8,
                           grouping_method = "exposure"){

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
#' @param grouping_method Use probability or exposure approach to group from input to output development factors.
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
                               percentage_data_training=.8,
                               grouping_method = "exposure"
                               ){


  set.seed(random_seed)

  formula_ct <- as.formula(IndividualData$string_formula_i)

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

    hazard_frame <- cbind(newdata, model.out$expg)
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"



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

  #return(hazard_frame)

  #need placeholder for latest i mirror cl behaviour

  # hazard_cl <- (sapply(seq_along(hz_names_i$time),
  #                      pkg.env$hazard_f,
  #                      enter= hz_names_i$enter,
  #                      time=hz_names_i$time,
  #                      exit=hz_names_i$exit,
  #                      event=hz_names_i$event))



  max_DP <- max(IndividualData$training$DP_rev_o)

  ############################################################
  #check
  #hazard_q <- matrix(nrow=max_DP, ncol=(ncol(hazard)-1)*IndividualData$conversion_factor)
  #eta_o <- c()


  ############################################################

  #Add development and relevant survival values to the hazard_frame
  hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                 categorical_features = IndividualData$categorical_features,
                                                 continuous_features = IndividualData$continuous_features)

  hazard_frame_grouped <- pkg.env$covariate_mapping(
    hazard_frame = hazard_frame_updated,
    categorical_features = IndividualData$categorical_features,
    continuous_features = IndividualData$continuous_features,
    conversion_factor = IndividualData$conversion_factor
    )


  latest_observed <- pkg.env$latest_observed_values_i(
    data=IndividualData$training.data,
    groups = hazard_frame_grouped$groups,
    categorical_features = IndividualData$categorical_features,
    continuous_features = IndividualData$continuous_features
  )

  expected_i <- pkg.env$predict_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    latest_cumulative = latest_observed$latest_cumulative,
    grouping_method = "exposure"
  )

  df_i <- pkg.env$retrieve_df_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    groups = hazard_frame_grouped$groups
  )

  hazard_frame_input <- pkg.env$input_hazard_frame(
    hazard_frame = hazard_frame_grouped$hazard_group,
    expected_i = expected_i ,
    categorical_features = IndividualData$categorical_features,
    continuous_features = IndividualData$continuous_features,
    df_i = df_i,
    groups = hazard_frame_grouped$groups)

  ############################################################
  # If input dimension not equal to output dimension, we perform a grouping #

  if(IndividualData$conversion_factor != 1){

  development_factor_o <- matrix(nrow=max_DP, ncol=(nrow(hazard_frame_grouped$groups)) )

  development_periods <- distinct(select(data.frame(IndividualData$training), AP_i, AP_o))

  # Calculate the minimum and maximum development periods for each row in development_periods for each DP_rev_o
  dp_ranges <- t(lapply(1:max_DP, function(DP_rev_o) {
    cbind(DP_rev_o, development_periods,
          min_dp = with(development_periods, AP_i+1/(IndividualData$conversion_factor)*(DP_rev_o-AP_o)),
          max_dp = with(development_periods, AP_i-1+1/(IndividualData$conversion_factor)*(DP_rev_o-AP_o+1)))
  }
  ))

  dp_ranges <- do.call(rbind, dp_ranges)

  expected_o <-pkg.env$predict_o(expected_i = expected_i,
                                 groups = hazard_frame_grouped$groups,
                                 conversion_factor = IndividualData$conversion_factor)

  development_factor_o <- mapply(pkg.env$i_to_o_development_factor,
                           1:max(hazard_frame_grouped$groups$group_o),
                           MoreArgs=list(hazard_data_frame=hazard_frame_grouped$hazard_group,
                                         expected_i = expected_i,
                                         dp_ranges = dp_ranges,
                                         groups = hazard_frame_grouped$groups,
                                         observed_pr_dp = latest_observed$observed_pr_dp,
                                         latest_cumulative = latest_observed$latest_cumulative,
                                         conversion_factor = IndividualData$conversion_factor,
                                         grouping_method = grouping_method))



#We only have 3 groups if AP is included as a covariate
  if(ncol(hazard_frame_grouped$groups) == 5){
    colnames(development_factor_o) <- unique(c(paste0("AP_o_",hazard_frame_grouped$groups$AP_o,",", hazard_frame_grouped$groups$covariate )))
  }
  else{
    colnames(development_factor_o) <- c(hazard_frame_grouped$groups$covariate )
  }

  df_o <- as.data.frame(development_factor_o[1:(nrow(development_factor_o)-1),]) %>%
    map_df(rev) %>%
    mutate(DP_o=row_number())


  hazard_frame_output <- pkg.env$output_hazard_frame(
    hazard_frame_input=hazard_frame_input,
    expected_o=expected_o,
    categorical_features=IndividualData$categorical_features,
    continuous_features=IndividualData$continuous_features,
    df_o=df_o,
    groups = hazard_frame_grouped$groups
  )

  out=list(df_output = df_o,
           df_input = df_i,
           hazard_frame_input = hazard_frame_input,
           hazard_frame_output = hazard_frame_output,
           IndividualData=IndividualData)

  class(out) <- c('ReSurvFit')

  return(out)
  }

  out=list(df_input = df_i,
           hazard_frame_input = hazard_frame_input,
           IndividualData=IndividualData)

  class(out) <- c('ReSurvFit')

  return(out)
}




