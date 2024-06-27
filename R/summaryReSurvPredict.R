#' Summary of IBNR predictions
#'
#' Gives overview of IBNR predictions
#'
#' @param object "ReSurvPredict" object specifying hazard and development factors.
#' @param granularity \code{character}, specify if which granularity the summary should be on.
#' \itemize{
#' \item{\code{"input"}}
#' \item{\code{"output"}}
#' }
#' Default is \code{"input"}.
#' @return summary of predictions
#'
#' @export
#' @method summary ReSurvPredict

summary.ReSurvPredict <- function(object, granularity = "input")
{
  handle <- match(granularity, c("input","output"))

  IBNR_AP <- switch(handle,
                    data.table(object$hazard_frame_input)[, .(IBNR=sum(IBNR, na.rm=T)), by = AP_i],
                    data.table(object$hazard_frame_output)[, .(IBNR=sum(IBNR, na.rm=T)), by = AP_o]
  )

  development_factor = switch(handle,
                              object$df_input,
                              object$df_output)

  keep = c(
    "model.out",
    "IndividualDataPP",
    "hazard_model"
  )

  summary <- list(
    IBNR_AP = IBNR_AP,
    total_IBNR = sum(IBNR_AP$IBNR),
    development_factor=development_factor,
    grouping_method = object$grouping_method,
    granularity = granularity,
    ReSurvFit = object$ReSurvFit[keep]
  )

  class(summary) <- "summaryReSurvPredict"
  return(summary)
}

#' Print summary of IBNR predictions
#'
#' Gives overview of IBNr predictions
#'
#' @param x "ReSurvPredict" object specifying hazard and development factors.
#' @param digits \code{numeric}, number of digits to print for IBNR and Likelihood.
#' @return print of summary of predictions
#'
#' @export
#' @method print summaryReSurvPredict

print.summaryReSurvPredict <-
  function (x, digits = max(3L, getOption("digits") - 3L))
  {
    cat("\n Hazard model:\n",
        paste(deparse(x$ReSurvFit$hazard_model), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    #cat("Likelihood: \n")
    #  xx <- x$ReSurvFit$model.out$likelihood
    #print.default(xx, digits = digits, na.print = "", print.gap = 2L)

    if(is.null(x$ReSurvFit$IndividualDataPP$categorical_features) & is.null(x$ReSurvFit$IndividualDataPP$continuous_features)) {
      cat("\nNo Features \n")
    } else {
      categorical_features<-NULL
      continuous_features <- NULL
      if(!is.null(x$ReSurvFit$IndividualDataPP$categorical_features)){
        categorical_features <- sprintf("\nCategorical Features:\n%s",
                                      paste(x$ReSurvFit$IndividualDataPP$categorical_features ,  collapse="\n"))
      }
      if(!is.null(x$ReSurvFit$IndividualDataPP$continuous_features)){
        continuous_features <- sprintf("\nContinuous Features:\n%s",
                                        paste(x$ReSurvFit$IndividualDataPP$continuous_features ,  collapse="\n"))
      }
      cat(categorical_features, continuous_features)

    }
    ##
    cat("\nTotal IBNR level: \n")
    print.default(x$total_IBNR, digits = digits, na.print = "", print.gap = 2L)

    # handle <- match(x$granularity, c("input","output"))
    #
    # df_string <- paste0("\n Development factors for ", x$granularity, " granularity.",
    #                     switch(
    #                       handle,
    #                       "",
    #                       paste0(" Estimated by ", x$ReSurvPredict$grouping_method, " priciple.")
    #                     ))
    #
    # cat(df_string
    # )
    # print(x$development_factor)


    invisible(x)
  }
