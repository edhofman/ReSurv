# Encoding, formula, model matrix, and scaler helper functions
#
# @importFrom fastDummies dummy_cols
## Encoding and formula ----

pkg.env$check.dates.consistency <- function(x,
                                            input_time_granularity,
                                            ap1){
"
This function checks weather the accident date and the reporting date are of 'Date' class.
In case they are, it transforms them into numeric.
"

  if(inherits(x, "Date")){

    if(input_time_granularity %in% c('quarters','semesters')){
      time_unit_string <- c('quarters', 'semesters', 'years')
      # BE CAREFUL: different from other codes, here we will bring everything to months and divide by six or four. Simpler.
      time_unit_numeric <- c(1/4, 1/6)
      input.pos <- which(time_unit_string%in%intersect(input_time_granularity,time_unit_string))
      divide.by <- time_unit_numeric[input.pos]
      diff.operator <- 'months'

    }else{

      divide.by <- 1
      diff.operator <- input_time_granularity

    }

    out <- floor(time_length(x-ap1,diff.operator)*divide.by)

    return(out)

  }else{

    return(x)

  }


}


pkg.env$encode.variables <- function(x,ap1){
  "
  This function encodes the periods.
  We impose that the indexization starts from 1.

  "

  seq <- min(c(x,ap1)):max(x)

  dim = length(seq)

  v <- 1:length(seq)

  v[match(x,seq)]


}

pkg.env$encode.variables.cp <- function(x,ap1){
  "
  This function encodes the periods.
  We impose that the indexization starts from 1.

  "
  seq <- ap1:max(x)

  v <- 1:length(seq)

  v[match(x,seq)]


}


pkg.env$formula.editor <- function(continuous_features,
                                   categorical_features,
                                   continuous_features_spline,
                                   degree_cf,
                                   degrees_of_freedom_cf,
                                   calendar_period,
                                   calendar_period_extrapolation,
                                   degree_cp,
                                   degrees_of_freedom_cp,
                                   input_output='i'){
  "
  This util edits creates the string that is used for model fitting in a compact way.
  continuous_features: character, vector of continuous features to be included in the linear predictor.
  categorical_features: character, vector of categorical features to be included in the linear predictor.
  continuous_features_spline: logical, T if a spline is added to model the continuous features.
  degree_cf: numeric, degrees of the spline for continuous features.
  degrees_of_freedom_cf: numeric, degrees of freedom of the spline for continuous features.
  degree_cp: numeric, degrees of the spline for calendar period.
  degrees_of_freedom_cp: numeric, degrees of freedom of the spline for calendar period features.
  input_output: character, set to input ('i') or output ('o') depending on the formula that we require.


  returns: the character that can be converted to a survival package formula object for the fit.

  "

  tmp.cat <- switch(!is.null(categorical_features), paste(categorical_features, collapse='+'), NULL)
  tmp.spline.pos <- which(continuous_features%in%intersect(continuous_features,continuous_features_spline))
  tmp.cont.pos <- which(!(continuous_features%in%intersect(continuous_features,continuous_features_spline)))
  tmp.cont <- switch(!is.null(continuous_features[tmp.cont.pos]) & length(continuous_features[tmp.cont.pos])>0, paste(continuous_features[tmp.cont.pos], collapse='+'), NULL)
  tmp.splines <- switch((!is.null(continuous_features[tmp.spline.pos]) & !is.null(continuous_features_spline)),paste0("pspline(",continuous_features[tmp.spline.pos], ",degree=",degree_cf,",df=",degrees_of_freedom_cf,")"),NULL)
  tmp.calendar <- switch(calendar_period_extrapolation,paste0("pspline(",calendar_period, ",degree=",degree_cf,",df=",degrees_of_freedom_cp,")"),NULL)

  tmp.all <- c(tmp.cat,tmp.cont,tmp.splines,tmp.calendar)

  if(is.null(tmp.all)){

    string_formula<- paste(paste0("survival::Surv","(TR_",input_output,", DP_rev_",input_output,", I) ~ "), "1")

  }else{

    string_formula<- paste(paste0("survival::Surv","(TR_",input_output,", DP_rev_",input_output,", I) ~ "),paste(tmp.all, collapse='+'))

  }


  string_formula


}




"This is a vectorized version of the grepl function.
See the grepl function documentation."
pkg.env$vgrepl <- Vectorize(grepl, vectorize.args = "pattern")

## Model Matrix helpers ----

pkg.env$model.matrix.creator <- function(data,
                                         select_columns,
                                         remove_first_dummy = FALSE){
  "
  This function encodes the matrices that we need for model fitting.

  "
  #
  #individual_data$training.data
  # X <- data %>%
  #   dummy_cols(select_columns = select_columns, #individual_data$categorical_features
  #              remove_selected_columns = TRUE,
  #              remove_first_dummy = remove_first_dummy)

  X <- dummy_cols(data,
                  select_columns = select_columns, #individual_data$categorical_features
                  remove_selected_columns = TRUE,
                  remove_first_dummy = remove_first_dummy)

  tmp.cond=as.logical(apply(pkg.env$vgrepl(pattern=select_columns,
                                           x=colnames(X)), #individual_data$categorical_features
                            MARGIN=1,
                            sum))

  setDT(X)

  X <- X[,.SD,.SDcols = colnames(X)[tmp.cond]]

  # X <- X %>%
  #   select(colnames(X)[tmp.cond] ) %>%
  #   as.data.frame()

  return(X)

}


pkg.env$model.matrix.extract.hazard.names <- function(X,
                                                      string_formula,
                                                      data){

  formula_ct <- as.formula(string_formula)
  Y<-model.extract(model.frame(formula_ct, data=data),"response")

  enter <- Y[, 1]
  exit <- Y[, 2]
  event <- Y[, 3] != 0
  sco <- exp(rep(0, nrow(Y)))

  time <- sort(seq(1,max(exit[event]), by=1)) #might be times with no events

  X_unique <- unique(X)

  names_hazard <- (data.frame(X_unique) %>%
                     rowwise() %>%
                     mutate(name = paste0(names(.)[c_across() == 1], collapse = ',')))$name

  return(list(enter=enter,
              exit=exit,
              event=event,
              sco=sco,
              time=time,
              names_hazard=names_hazard))


}

## Scalers ----

pkg.env$MinMaxScaler <- function(x, na.rm = TRUE) {
  "MinMax Scaler"
  return(2*(x- min(x)) /(max(x)-min(x))-1)
}
pkg.env$StandardScaler <- function(x, na.rm = TRUE) {
  "Standard Scaler"
  return( (x-mean(x))/sd(x) )
}

pkg.env$scaler <- function(continuous_features_scaling_method){
  "Apply the scaling method"
  if(continuous_features_scaling_method == "minmax" ){return(pkg.env$MinMaxScaler)}
  if(continuous_features_scaling_method == "standard" ){return(pkg.env$StandardScaler)}


}
