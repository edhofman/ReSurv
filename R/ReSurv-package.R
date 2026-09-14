#' ReSurv: Individual claim-count reserving
#'
#' Fit reverse-time hazard models with \code{\link{ReSurv}}, following
#' preprocessing with \code{\link{IndividualDataPP}}. Use
#' \code{\link{ReSurvCV}} for hyperparameter selection,
#' \code{\link{predictReserve}} for count forecasts, and
#' \code{\link{Score_Reserving}} to evaluate realized claims.
#'
#' @importFrom stats ave model.extract model.frame pnorm predict rexp rlnorm runif sd
#' @importFrom lubridate time_length
#' @importFrom actuar rtrgamma
#' @importFrom dplyr pick rows_update rowwise ungroup
#' @importFrom ggplot2 ggplot aes rel geom_line ylim labs scale_x_continuous theme_bw theme element_text
#' @keywords internal
"_PACKAGE"
