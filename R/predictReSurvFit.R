#' Predict IBNR frequency
#'
#' This function predicts the results from the ReSurv fits.
#'
#' @param object \code{ResurvFit} object specifying start time, end time and status.
#' @param newdata \code{IndividualData} object that contains new data to predict.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#'
#' @return predictions
#'
#' @export
#' @method predict ReSurvFit
predict.ReSurvFit <- function(object,
                              newdata=NULL,
                              grouping_method = "probability",
                              check_value = 1.85){


  if(!is.null(newdata)){
    pkg.env$check.newdata(newdata=newdata,
                          pastdata=object$IndividualData)

    idata <- newdata
    # hazard_frame <- adjust.predictions(ResurvFit=object,
    #                                    hazard_model=object$hazard_model,
    #                                    idata=idata)


  }else{

    idata <- object$IndividualData


    }

  hazard_frame <-object$hazard_frame
  # browser()
  hazard_frame_grouped <- pkg.env$covariate_mapping(
    hazard_frame = hazard_frame,
    categorical_features = idata$categorical_features,
    continuous_features = idata$continuous_features,
    conversion_factor = idata$conversion_factor,
    calendar_period_extrapolation = idata$calendar_period_extrapolation
  )

  # browser()
  missing.obsevations <- pkg.env$fill_data_frame(data=idata$full.data,
                                                 continuous_features=idata$continuous_features,
                                                 categorical_features=idata$categorical_features,
                                                 years=idata$years,
                                                 input_time_granularity=idata$input_time_granularity,
                                                 conversion_factor=idata$conversion_factor)


  latest_observed <- pkg.env$latest_observed_values_i(
    data_reserve= bind_rows(idata$training.data, missing.obsevations),
    groups = hazard_frame_grouped$groups,
    categorical_features = idata$categorical_features,
    continuous_features = idata$continuous_features,
    calendar_period_extrapolation = idata$calendar_period_extrapolation
  )

  max_DP <- max(bind_rows(idata$training.data, missing.obsevations)$DP_rev_o)


  expected_i <- pkg.env$predict_i(
    hazard_data_frame = lazy_dt(hazard_frame_grouped$hazard_group),
    latest_cumulative = latest_observed$latest_cumulative,
    grouping_method = "exposure",
    min_DP_rev_i = min(hazard_frame_grouped$hazard_group$DP_rev_i)
  )

  df_i <- pkg.env$retrieve_df_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    groups = hazard_frame_grouped$groups
  )

  hazard_frame_input <- pkg.env$input_hazard_frame(
    hazard_frame = hazard_frame_grouped$hazard_group,
    expected_i = expected_i ,
    categorical_features = idata$categorical_features,
    continuous_features = idata$continuous_features,
    df_i = df_i,
    groups = hazard_frame_grouped$groups)

  # browser()
  if(idata$conversion_factor != 1){

    development_periods <- distinct(select(data.frame(idata$training), AP_i, AP_o))

    # Calculate the minimum and maximum development periods for each row in development_periods for each DP_rev_o
    dp_ranges <- t(lapply(1:max_DP, function(DP_rev_o) {
      cbind(DP_rev_o, development_periods,
            min_dp = with(development_periods, AP_i+1/(idata$conversion_factor)*(DP_rev_o-AP_o)),
            max_dp = with(development_periods, AP_i-1+1/(idata$conversion_factor)*(DP_rev_o-AP_o+1)))
    }
    ))

    dp_ranges <- do.call(rbind, dp_ranges)


    #check_input_hazard <- pkg.env$check_input_hazard(hazard_frame_input,
    #                           check_value=check_value)

    #If we exeed the check value, we calculate on ouput granularity, predict on output granularity, and distribute evenly in the relevant input-periods.
    #From here we do simple chain-ladder to calculate new development factor.
    if(#check_input_hazard
      FALSE # I guess we  can remove this part if not used?
    ){
      development_factor_o <- pkg.env$i_to_o_development_factor(
        hazard_data_frame=hazard_frame_grouped$hazard_group,
        expected_i = expected_i,
        dp_ranges = dp_ranges,
        groups = hazard_frame_grouped$groups,
        observed_pr_dp = latest_observed$observed_pr_dp,
        latest_cumulative = latest_observed$latest_cumulative,
        conversion_factor = idata$conversion_factor,
        grouping_method = "probability",
        min_DP_rev_i = min(hazard_frame_grouped$hazard_group$DP_rev_i)
      )

      if(ncol(hazard_frame_grouped$groups) == 5){
        colnames(development_factor_o) <- unique(c(paste0("AP_o_",hazard_frame_grouped$groups$AP_o,",", hazard_frame_grouped$groups$covariate )))
      }
      else{
        colnames(development_factor_o) <- c(hazard_frame_grouped$groups$covariate )
      }

      df_o <- as.data.frame(development_factor_o[1:(nrow(development_factor_o)-1),]) %>%
        map_df(rev) %>%
        mutate(DP_o=row_number())

      #We only update for relevant periods, hence for example for accident periods, where we have already seen the development, we just put to 1.
      hazard_frame_grouped$hazard_group <- pkg.env$update_hazard_frame(
        hazard_frame_input=hazard_frame_input,
        hazard_frame_grouped=hazard_frame_grouped$hazard_group,
        df_o=df_o,
        latest_observed_i = latest_observed$observed_pr_dp,
        groups = hazard_frame_grouped$groups,
        categorical_features=idata$categorical_features,
        continuous_features = idata$continuous_features,
        conversion_factor = idata$conversion_factor,
        check_value = check_value
      )

      expected_i <- pkg.env$predict_i(
        hazard_data_frame = hazard_frame_grouped$hazard_group,
        latest_cumulative = latest_observed$latest_cumulative,
        grouping_method = "exposure"
      )

      df_i <- pkg.env$retrieve_df_i(
        hazard_data_frame = hazard_frame_grouped$hazard_group,
        groups = hazard_frame_grouped$groups,
        adjusted=T
      )

      hazard_frame_input <- pkg.env$input_hazard_frame(
        hazard_frame = hazard_frame_grouped$hazard_group,
        expected_i = expected_i ,
        categorical_features = idata$categorical_features,
        continuous_features = idata$continuous_features,
        df_i = df_i,
        groups = hazard_frame_grouped$groups,
        adjusted=T)



    }

    expected_o <-pkg.env$predict_o(expected_i = expected_i,
                                   groups = hazard_frame_grouped$groups,
                                   conversion_factor = idata$conversion_factor,
                                   years = object$IndividualData$years,
                                   input_time_granularity = object$IndividualData$input_time_granularity)



    development_factor_o <- pkg.env$i_to_o_development_factor(
      hazard_data_frame=hazard_frame_grouped$hazard_group,
      expected_i = expected_i,
      dp_ranges = dp_ranges,
      groups = hazard_frame_grouped$groups,
      observed_pr_dp = latest_observed$observed_pr_dp,
      latest_cumulative = latest_observed$latest_cumulative,
      conversion_factor = idata$conversion_factor,
      grouping_method = grouping_method,
      min_DP_rev_i = min(hazard_frame_grouped$hazard_group$DP_rev_i),
      years = object$IndividualData$years,
      input_time_granularity = object$IndividualData$input_time_granularity
    )



    #We only have 5 groups if AP is included as a covariate
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
      categorical_features=idata$categorical_features,
      continuous_features=idata$continuous_features,
      df_o=df_o,
      groups = hazard_frame_grouped$groups
    )

    out=list(ReSurvFit = object,
             df_output = df_o,
             df_input = df_i,
             hazard_frame_input = hazard_frame_input,
             hazard_frame_output = hazard_frame_output,
             grouping_method = grouping_method)

    class(out) <- c('ReSurvPredict')

    return(out)
  }

  out=list(ReSurvFit = object,
           df_input = df_i,
           hazard_frame_input = hazard_frame_input,
           grouping_method = grouping_method)

  class(out) <- c('ReSurvPredict')

  return(out)


}




