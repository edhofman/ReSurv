#' Individual data set
#'
#' This function pre-processes the data for the application of a \code{ReSurv} model.
#'
#' @param data \code{data.frame}, for the individual reserving. The number of development periods can be larger than the number of accident periods.
#' @param id \code{character}, \code{data} column that contains the policy identifier.
#' @param continuous_features \code{character}, continuous features columns to be scaled.
#' @param categorical_features \code{character}, categorical features columns to be one-hot encoded.
#' @param accident_period \code{character}, it contains the name of the column in data corresponding to the accident period.
#' @param calendar_period \code{character}, it contains the name of the column in data corresponding to the calendar period.
#' @param calendar_period_extrapolation \code{character}, whether a spline for calendar extrapolation should be considered in the cox model fit.
#'                                       Default is `FALSE`.
#' @param input_time_granularity \code{character}, time unit of the input data. Granularity supported:
#' \itemize{
#' \item{\code{"months"}: the input data are monthly.}
#' \item{\code{"quarters"}: the input data are quarterly}
#' \item{\code{"years"}: the input data are yearly.}
#' }
#' Default to \code{months}.
#'
#' @param output_time_granularity \code{character}, time unit of the output data. The granularity supported is the same as for the input data:
#'  \itemize{
#' \item{\code{"months"}: the input data are monthly.}
#' \item{\code{"quarters"}: the input data are quarterly}
#' \item{\code{"years"}: the input data are yearly.}
#' }
#' The output granularity must be bigger than the input granularity.
#' Also, the output granularity must be consistent with the input granularity, meaning that the time conversion must be possible.
#' E.g., it is possible to group quarters to years. It is not possible to group quarters to semesters.
#' Default to \code{quarters}.
#'
#' @param years \code{numeric}, number of development years in the study.
#' @param continuous_features_spline \code{logical}, weather a spline for smoothing continuous features should be added.
#' @param degrees_cf \code{numeric}, degrees of the spline for smoothing continuous features.
#' @param degrees_of_freedom_cf \code{numeric}, degrees of freedom of the splines for smoothing continuous features.
#' @param degrees_cp \code{numeric}, degrees of the spline for smoothing the calendar period effect.
#' @param degrees_of_freedom_cp \code{numeric}, degrees of freedom of the splines for smoothing the calendar period effect.
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr row_number
#' @importFrom dplyr summarize
#' @importFrom purrr map_df
#' @importFrom stats as.formula
#' @importFrom data.table data.table
#'
#'
#' @return \code{IndividualData} object. A list containing:
#'  \itemize{
#' \item{\code{full.data}: the input data after pre-processing.}
#' \item{\code{starting.data}: the input data as they were provided from the user.}
#' \item{\code{training.data}: the input data pre-processed for training.}
#' \item{\code{conversion_factor}: the conversion factor for going from input granularity to output granularity. E.g, the conversion factor for going from months to quarters is 1/3.}
#' \item{\code{string_formula_i}: string of the \code{survival} formula to model the data in input granularity.}
#' \item{\code{string_formula_o}: string of the \code{survival} formula to model the in data output granularity.}
#' \item{\code{continuous_features}: the continuous features names as provided from the user.}
#' \item{\code{categorical_features}: the categorical features names as provided from the user.}
#' \item{\code{calendar_period_extrapolation}: the \code{logical} value specifying wheter a calendar period component is extrapolated.}
#'
#' }
#'
#' After pre-processing, we provide a standard encoding for the time components. This regards the output in \code{training.data} and \code{full.data}.
#' In the \code{ReSurv} notation:
#'\itemize{
#'\item{\code{AP_i}: input accident period.}
#'\item{\code{AP_o}: output accident period.}
#'\item{\code{DP_i}: input development period in forward time.}
#'\item{\code{DP_rev_i}: input development period in reverse time.}
#'\item{\code{DP_rev_o}: output development period in reverse time.}
#'\item{\code{TR_i}: input truncation time.}
#'\item{\code{TR_o}: output truncation time.}
#'\item{\code{I}: event indicator, under this framework is equal to one for each entry. }
#'}
#'
#'
#'
#' @references
#' Munir, H., Emil, H., & Gabriele, P. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' @export
IndividualData <- function(data,
                           id=NULL,
                           continuous_features,
                           categorical_features,
                           accident_period,
                           calendar_period,
                           input_time_granularity="months",
                           output_time_granularity="quarters",
                           years=4,
                           calendar_period_extrapolation=FALSE,
                           continuous_features_spline=NULL,
                           degrees_cf=3,
                           degrees_of_freedom_cf=4,
                           degrees_cp=3,
                           degrees_of_freedom_cp=4){



  # Work on a copy of the input data
  tmp <- as.data.frame(data)
  # browser()

  # Accident periods encoding
  x.ap <- pkg.env$check.dates.consistency(tmp[,accident_period],
                                          input_time_granularity=input_time_granularity,
                                          ap1=min(tmp[,accident_period]))
  tmp.ap <- pkg.env$encode.variables(x.ap)

  # Calendar periods encoding
  x.cp <- pkg.env$check.dates.consistency(tmp[,calendar_period],
                                          input_time_granularity=input_time_granularity,
                                          ap1=min(tmp[,accident_period]))
  tmp.cp <- pkg.env$encode.variables.cp(x.cp,
                                        ap1=min(x.ap))
  # Development periods encoding
  # browser()
  tmp.dp <- tmp.cp-tmp.ap+1

  # Check the ap among features
  continuous_features<-pkg.env$fix.double.ap(features=continuous_features,accident_period=accident_period)
  categorical_features<-pkg.env$fix.double.ap(features=categorical_features,accident_period=accident_period)

  # The following checks warn you if there is a missing accident period or reporting period
  # in the data. They do not interrupt the code.
  pkg.env$check.all.present(tmp.ap, check.on='accident periods')
  pkg.env$check.all.present(tmp.cp, check.on='calendar periods')

  dim1=max(tmp.ap)
  dim2=max(tmp.dp)

  # You need at least the number of development periods on the rows.
  # if(dim1>dim2){
  #
  #   stop("The number of accident periods is bigger than the number of development periods")
  #
  # }


  # We need a conversion factor from input_time_granularity to output_time_granularity

  # conversion_factor <- input_time_granularity*(1/output_time_granularity)

  conversion_factor <- pkg.env$conversion.factor.of.time.units(input_time_granularity,
                                                               output_time_granularity)



  # Build the variables you need
  tmp = tmp %>%
    mutate(AP_i=tmp.ap,
           DP_i=tmp.dp,
           RP_i=tmp.cp,
           DP_rev_i = pkg.env$maximum.time(years,input_time_granularity) - DP_i+1,
           TR_i = AP_i-1, #just setting truncation to max year simulated. and accounting for
           I=1) %>%
    as.data.frame()

  max_dp_i =  pkg.env$maximum.time(years,input_time_granularity)

  # Take the training data (upper triangle) and convert it from input_time_granularitys to output_time_granularitys
  train= tmp %>%
    filter(DP_rev_i > TR_i) %>%
    mutate(
      DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
      AP_o = ceiling(AP_i*conversion_factor)
    ) %>%
    mutate(TR_o= AP_o-1) %>%
    mutate(across(all_of(categorical_features),
                  as.factor)) %>%
    select(id,
           all_of(categorical_features),
           all_of(continuous_features),
           all_of(switch(calendar_period_extrapolation, 'RP_i', NULL)),
           AP_i,
           AP_o,
           DP_i,
           DP_rev_i,
           DP_rev_o,
           TR_i,
           TR_o,
           I) %>%
    as.data.frame()



  if(calendar_period_extrapolation){
    train= train %>%
      mutate(
        RP_o=ceiling(RP_i*conversion_factor)) %>%
      as.data.frame()

  }

  string_formula_i <- pkg.env$formula.editor(continuous_features=continuous_features,
                                             categorical_features=categorical_features,
                                             continuous_features_spline=continuous_features_spline,
                                             degree_cf=degrees_cf,
                                             degrees_of_freedom_cf=degrees_of_freedom_cf,
                                             calendar_period="RP_i",
                                             calendar_period_extrapolation=calendar_period_extrapolation,
                                             degree_cp=degree_cp,
                                             degrees_of_freedom_cp=degrees_of_freedom_cp,
                                             input_output='i')

  string_formula_o <- pkg.env$formula.editor(continuous_features=continuous_features,
                                             categorical_features=categorical_features,
                                             continuous_features_spline=continuous_features_spline,
                                             degree_cf=degrees_cf,
                                             degrees_of_freedom_cf=degrees_of_freedom_cf,
                                             calendar_period="RP_o",
                                             calendar_period_extrapolation=calendar_period_extrapolation,
                                             degree_cp=degree_cp,
                                             degrees_of_freedom_cp=degrees_of_freedom_cp,
                                             input_output='o')


  # Create and organize the output
  out <- list(full.data = tmp,
              starting.data = data,
              training.data = train,
              conversion_factor=conversion_factor,
              string_formula_i=string_formula_i,
              string_formula_o=string_formula_o,
              continuous_features=continuous_features,
              categorical_features=categorical_features,
              calendar_period_extrapolation=calendar_period_extrapolation,
              years=years,
              accident_period=accident_period,
              calendar_period=calendar_period,
              input_time_granularity=input_time_granularity,
              output_time_granularity=output_time_granularity)

  # Return the correct output
  class(out) <- "IndividualData"

  out

  }








