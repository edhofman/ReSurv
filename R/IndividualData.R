#' Individual data set
#'
#' This function inputs and pre-processes the data for the application.
#'
#' @param data data.frame for the individual reserving. The number of development periods can be bigger than the number of accident periods
#' @param id column that contains the policy id.
#' @param continuous_features continuous features columns to be scaled.
#' @param categorical_features categorical features columns to be one-hot encoded.
#' @param accident_period string that contains the name of the column in data corresponding to accident_period
#' @param calendar_period string that contains the name of the column in data corresponding to the calendar_period.
#' @param calendar_period_extrapolation logical, whether a spline for calendar extrapolation should be considered in the cox model fit. Default is FALSE.
#' @param input_time_granularity time unit of the input data with reference to a year
#' @param output_time_granularity time unit of the output data with reference to a year
#' @param continuous_features_spline weather a spline for smoothing continuous features should be added.
#' @param degree_cf degrees of the spline for smoothing continuous features.
#' @param degrees_of_freedom_cf degrees of freedom of the splines for smoothing continuous features.
#' @param degree_cp degrees of the spline for smoothing the calendar period effect.
#' @param degrees_of_freedom_cp degrees of freedom of the splines for smoothing the calendar period effect.
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
#' @return Pre-processed data ready for individual reserving.
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
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


  # Check the conversion is possible
  # pkg.env$check.time.units(input_time_granularity,
                           # output_time_granularity)

  # Work on a copy of the input data
  tmp <- as.data.frame(data)

  # Accident periods encoding
  tmp.ap <- pkg.env$encode.variables(tmp[,accident_period])
  # Calendar periods encoding
  tmp.cp <- pkg.env$encode.variables(tmp[,calendar_period])
  # Development periods encoding
  tmp.dp <- tmp.cp-tmp.ap+1

  # The following checks warn you if there is a missing accident period or reporting period
  # in the data. They do not interrupt the code.
  pkg.env$check.all.present(tmp.ap, check.on='accident periods')
  pkg.env$check.all.present(tmp.cp, check.on='calendar periods')

  dim1=max(tmp.ap)
  dim2=max(tmp.dp)

  # You need at least the number of development periods on the rows.
  if(dim1>dim2){

    stop("The number of accident periods is bigger than the number of development periods")

  }


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



  # Take the training data (upper triangle) and convert it from input_time_granularitys to output_time_granularitys
  train= tmp %>%
    filter(DP_rev_i > TR_i) %>%
    mutate(
      DP_rev_o = floor(max(DP_i)*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
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
              calendar_period_extrapolation=calendar_period_extrapolation)

  # Return the correct output
  class(out) <- "IndividualData"

  out

  }








