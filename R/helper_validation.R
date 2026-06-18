## Checks ----
pkg.env$validate_eta <- function(eta) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta)) {
    stop("`eta` must be a single finite numeric value.", call. = FALSE)
  }

  if (eta < 0 || eta > 1) {
    stop("`eta` must lie in [0, 1].", call. = FALSE)
  }

  eta
}

pkg.env$check.all.present <- function(x,check.on){

  "
  This function checks that you have all the periods in the data,
  from the minimum record to the maximum record.

  x: integer or numeric, input data to check.
  check.on: character, it specifies the variable to check. I.e., accident period and calendar period.

  "

  tmp <- x

  v <- diff(as.integer(sort(unique(x))))

  if(sum(v>1)>0){

    warning(paste("Some", check.on, "are missing in the data"))

  }

}

pkg.env$check.time.units <- function(input_time_unit,
                                     output_time_unit){

  "
  This function checks that the input output time transformation is consistent.
  E.g. you cannot turn quarters into semesters, you can instead turn trimesters into semesters.

  input_time_unit: numeric, input time unit with respect to one year. E.g., 1/12 for months.
  output_time_unit: numeric, output time unit with respect to one year. E.g., 1/4 for quarters.
  "

  if((1/input_time_unit)%%(1/output_time_unit) != 0){

    stop('The provided time granularities are not subsettable.')

  }

}


pkg.env$maximum.time <- function(years,
                                 input_time_granularity){

  "
  This function returns the triangle width.

  years: numeric, number of years in the triangle.
  input_time_granularity: numeric, input data granularity with respect to the one year reference. E.g., 1/12 for months.

  "

  time_unit_string <- c('days','months','quarters', 'semesters', 'years')
  time_unit_numeric <- c(1/360, 1/12, 1/4, 1/2, 1)

  input.pos <- which(time_unit_string%in%intersect(input_time_granularity,time_unit_string))

  years/time_unit_numeric[input.pos]

}

pkg.env$conversion.factor.of.time.units <- function(input_time_unit,
                                                    output_time_unit){

  "
  This function computes the conversion factor of the time units.
  Given an input granularity and an output granularity, it returns the numeric conversion factor.
  E.g., the conversion factor is 1/3 to go from months to quarters.

  input_time_unit: character, input time granularity.
  output_time_unit: character, output time granularity.

  returns: numeric, conversion factor.

  "
  time_unit_string <- c('days', 'months', 'quarters', 'semesters', 'years')
  time_unit_numeric <- c(1/360, 1/12, 1/4, 1/2, 1)

  input.pos <- which(time_unit_string%in%intersect(input_time_unit,time_unit_string))
  output.pos <- which(time_unit_string%in%intersect(output_time_unit,time_unit_string))

  input_numeric <- time_unit_numeric[input.pos]
  output_numeric <- time_unit_numeric[output.pos]

  pkg.env$check.time.units(input_numeric,
                           output_numeric)



  conversion_factor <- input_numeric*(1/output_numeric)
  conversion_factor

}


pkg.env$total.years.in.the.data <- function(input_time_unit,
                                            development_period){

  "
  This function computes the total number of years in the data, if not provided by the user.

  input_time_unit: character, input time granularity.
  development_period: numeric, vector of dp_i.

  returns: numeric, number of years in the data.

  "
  time_unit_string <- c('days', 'months', 'quarters', 'semesters', 'years')
  time_unit_numeric <- c(1/360, 1/12, 1/4, 1/2, 1)

  input.pos <- which(time_unit_string%in%intersect(input_time_unit,time_unit_string))

  input_numeric <- time_unit_numeric[input.pos]
  output_numeric <- 1

  conversion_factor <- input_numeric*(1/output_numeric)
  return(ceiling(max(development_period)*conversion_factor))

}

pkg.env$check.traintestsplit <- function(x){

  "
  This function checks that the training test split is specified correctly.

  x: numeric, training test split.

  returns: numeric, the default split of eighty percent when the specified training test split is not between zero and one.

  "

  if(x>1 | x<0 ){
    warning(paste0("Traintestsplit has been put to ", x,". The value needs to be between 0 and 1, defaulting to 0.8."))
    return(.8)
  }else{return(x)}


}


pkg.env$check_input_hazard <- function(hazard_frame_input, check_value=1.9){
  check <- hazard_frame_input %>%  dplyr::filter(hazard > check_value & DP_rev_i < max(DP_rev_i))

  if(nrow(check)>0){
    warning(paste0("Hazard value on input granularity exceeds ", check_value,
                   " for reverse development periods ", unique(check$DP_rev_i),". This is most likely due to low exposure, we calculate 'probability'-grouped output development factors, and from here adjust input development factor. ")
    )
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}


pkg.env$check.newdata <- function(newdata,
                                  pastdata){


  # cf <- pkg.env$conversion.factor.of.time.units(pastdata$input_time_granularity,
  #                                               newdata$output_time_granularity)

  if(!identical(pastdata$input_time_unit,newdata$data_information$input_time_unit)){

    stop('newdata must have the same input granularity as pastdata.')

  }

  # old:class(newdata) != "IndividualDataPP"
  if(!inherits(newdata, "IndividualDataPP")){

    stop('newdata must be an IndividualDataPP object.')

  }

  newfeatures <- c(newdata$data_information$categorical_features, newdata$data_information$continuous_features)
  pastfeatures <- c(pastdata$categorical_features, pastdata$continuous_features)

  if(!identical(newfeatures,pastfeatures)){

    stop('newdata must have the same features as pastdata.')

  }



}

