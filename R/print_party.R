#' @export
print.IndividualDataPP <- function(x, ...) {

  data_information <- x$data_information

  n_training <- tryCatch(
    nrow(x$training.data),
    error = function(e) NA_integer_
  )

  categorical_features <- data_information$categorical_features
  continuous_features <- data_information$continuous_features

  if (is.null(categorical_features) || length(categorical_features) == 0L) {
    categorical_features <- "none"
  } else {
    categorical_features <- paste(as.character(categorical_features), collapse = ", ")
  }

  if (is.null(continuous_features) || length(continuous_features) == 0L) {
    continuous_features <- "none"
  } else {
    continuous_features <- paste(as.character(continuous_features), collapse = ", ")
  }

  cat("\nIndividualDataPP object\n")
  cat("-----------------------\n")

  cat("Rows in training data: ",
      ifelse(is.na(n_training), "not available", n_training),
      "\n",
      sep = "")

  cat("Categorical features: ", categorical_features, "\n", sep = "")
  cat("Continuous features: ", continuous_features, "\n", sep = "")

  cat("Input time granularity: ",
      ifelse(is.null(data_information$input_time_granularity),
             "not available",
             data_information$input_time_granularity),
      "\n",
      sep = "")

  cat("Output time granularity: ",
      ifelse(is.null(data_information$output_time_granularity),
             "not available",
             data_information$output_time_granularity),
      "\n",
      sep = "")

  cat("Years: ",
      ifelse(is.null(data_information$years),
             "not available",
             data_information$years),
      "\n",
      sep = "")

  cat("Conversion factor: ",
      ifelse(is.null(data_information$conversion_factor),
             "not available",
             data_information$conversion_factor),
      "\n",
      sep = "")

  invisible(x)
}
## ------------------------------------------------------------------
## print method for ReSurvPredict objects
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## print method for ReSurvFit objects
## ------------------------------------------------------------------

#' @export
print.ReSurvFit <- function(x, ...) {

  data_information <- x$data_information
  fit_information <- x$fit_information

  hazard_model <- NULL

  if (!is.null(x$hazard_model)) {
    hazard_model <- x$hazard_model
  }

  if (is.null(hazard_model) && !is.null(fit_information$hazard_model)) {
    hazard_model <- fit_information$hazard_model
  }

  if (is.null(hazard_model) && !is.null(x$model.out$hazard_model)) {
    hazard_model <- x$model.out$hazard_model
  }

  if (is.null(hazard_model) && !is.null(x$model.out$model.out$hazard_model)) {
    hazard_model <- x$model.out$model.out$hazard_model
  }

  if (is.null(hazard_model) || length(hazard_model) == 0L) {
    hazard_model <- "not available"
  }

  categorical_features <- data_information$categorical_features
  continuous_features <- data_information$continuous_features

  if (is.null(categorical_features) || length(categorical_features) == 0L) {
    categorical_features <- "none"
  } else {
    categorical_features <- paste(as.character(categorical_features), collapse = ", ")
  }

  if (is.null(continuous_features) || length(continuous_features) == 0L) {
    continuous_features <- "none"
  } else {
    continuous_features <- paste(as.character(continuous_features), collapse = ", ")
  }

  n_hazard <- tryCatch(
    nrow(x$hazard_frame),
    error = function(e) NA_integer_
  )

  cat("\nReSurvFit object\n")
  cat("----------------\n")

  cat("Hazard model: ", hazard_model, "\n", sep = "")
  cat("Categorical features: ", categorical_features, "\n", sep = "")
  cat("Continuous features: ", continuous_features, "\n", sep = "")

  cat("Rows in hazard frame: ",
      ifelse(is.na(n_hazard), "not available", n_hazard),
      "\n",
      sep = "")

  cat("Input time granularity: ",
      ifelse(is.null(data_information$input_time_granularity),
             "not available",
             data_information$input_time_granularity),
      "\n",
      sep = "")

  cat("Output time granularity: ",
      ifelse(is.null(data_information$output_time_granularity),
             "not available",
             data_information$output_time_granularity),
      "\n",
      sep = "")

  cat("Years: ",
      ifelse(is.null(data_information$years),
             "not available",
             data_information$years),
      "\n",
      sep = "")

  if (!is.null(fit_information$LKH)) {
    cat("Likelihood: ", fit_information$LKH, "\n", sep = "")
  }

  invisible(x)
}
#' @export
print.ReSurvPredict <- function(x, ...) {

  fit <- x$ReSurvFit
  data_information <- fit$data_information
  fit_information <- fit$fit_information

  hazard_model <- NULL

  if (!is.null(fit$hazard_model)) {
    hazard_model <- fit$hazard_model
  }

  if (is.null(hazard_model) && !is.null(fit_information$hazard_model)) {
    hazard_model <- fit_information$hazard_model
  }

  if (is.null(hazard_model) && !is.null(fit$model.out$hazard_model)) {
    hazard_model <- fit$model.out$hazard_model
  }

  if (is.null(hazard_model) && !is.null(fit$model.out$model.out$hazard_model)) {
    hazard_model <- fit$model.out$model.out$hazard_model
  }

  if (is.null(hazard_model) || length(hazard_model) == 0L) {
    hazard_model <- "not available"
  }

  categorical_features <- data_information$categorical_features
  continuous_features <- data_information$continuous_features

  if (is.null(categorical_features) || length(categorical_features) == 0L) {
    categorical_features <- "none"
  } else {
    categorical_features <- paste(as.character(categorical_features), collapse = ", ")
  }

  if (is.null(continuous_features) || length(continuous_features) == 0L) {
    continuous_features <- "none"
  } else {
    continuous_features <- paste(as.character(continuous_features), collapse = ", ")
  }

  input_rows <- tryCatch(
    nrow(x$long_triangle_format_out$input_granularity),
    error = function(e) NA_integer_
  )

  output_rows <- tryCatch(
    nrow(x$long_triangle_format_out$output_granularity),
    error = function(e) NA_integer_
  )

  has_output <- !is.null(x$long_triangle_format_out$output_granularity)
  has_lower_triangle <- !is.null(x$lower_triangle)

  predicted_counts <- x$predicted_counts

  if (is.null(predicted_counts) || length(predicted_counts) == 0L) {
    predicted_counts <- "not available"
  }

  cat("\nReSurvPredict object\n")
  cat("--------------------\n")

  cat("Hazard model: ", hazard_model, "\n", sep = "")
  cat("Categorical features: ", categorical_features, "\n", sep = "")
  cat("Continuous features: ", continuous_features, "\n", sep = "")
  cat("Predicted IBNR count: ", predicted_counts, "\n", sep = "")

  cat("Input-granularity rows: ",
      ifelse(is.na(input_rows), "not available", input_rows),
      "\n",
      sep = "")

  cat("Output granularity available: ", has_output, "\n", sep = "")

  if (has_output) {
    cat("Output-granularity rows: ",
        ifelse(is.na(output_rows), "not available", output_rows),
        "\n",
        sep = "")
  }

  cat("Lower triangle available: ", has_lower_triangle, "\n", sep = "")

  invisible(x)
}
