#' Individual data generator
#'
#' This function generates the monthly individual claims data in the accompanying methodological paper using the \code{SynthETIC} package.
#' This simple function allows to simulate from a sand-box to test out the \code{ReSurv} approach.
#' Some parameters of the simulation can be changed.
#'
#' @param ref_claim \code{integer}, reference claim size.
#' @param time_unit \code{numeric}, output time unit.
#' @param years \code{integer}, number of years to be simulated.
#' @param random_seed \code{integer}, random seed for replicable code.
#' @param period_exposure \code{integer}, volume (number of policies) underwritten each period.
#' @param period_frequency \code{numeric}, expected frequency in each period.
#' @param scenario \code{character} or \code{numeric}, one of the scenarios described in the accompanying manuscript. Possible choices are
#'                  'alpha' (0), 'beta' (1), 'gamma'(2), 'delta'(3),'epsilon'(4). Our simulated data are constituted of a mix of short tail claims (\code{claim_type 0}) and claims with longer resolution (\code{claim_type 1}).
#'                  We chose the parameter of the simulator to resemble a mix of property damage (\code{claim_type 0}) and bodily injuries (\code{claim_type 1}). each scenario has distinctive characteristics.
#'                  Scenario Alpha is a mix of \code{claim_type 0} and  \code{claim_type 1} with same number of claims volume at each accident period.
#'                  Differently from scenario Alpha, in scenario Beta the volumes of \code{claim_type 1} are decreasing in the most recent accident periods.
#'                  In scenario Gamma we add an interaction between \code{claim_type 1} and accident period: in a real world setting this can be motivated by a change in consumer behavior or company policies resulted in different reporting patterns over time.
#'                  In scenario Delta, we introduce a seasonality effect dependent on the accident period for \code{claim_type 0} and \code{claim_type 1}.
#'                  In the real word, scenario Delta resembles seasonal changes in the workforce composition. Scenario Epsilon does not satisfy the proportionality assumption.
#'
#' @import SynthETIC
#'
#' @examples
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 2,
#' period_exposure = 100)
#'
#'
#'
#'
#' @return Individual claims data. It contains the following columns:
#' \itemize{
#' \item{\code{claim_number}: Policy ID.}
#' \item{\code{claim_type}: Type of claim. It can be either 0 or 1.}
#' \item{\code{AP}: Accident period}
#' \item{\code{RP}: Reporting period.}
#' }
#'
#' @references
#' Avanzi, B., Taylor, G., Wang, M., & Wong, B. (2021). SynthETIC: an individual insurance claim simulator with feature control. Insurance: Mathematics and Economics, 100, 296-308.
#'
#' Hiabu, M., Hofman, E., & Pittarello, G. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' @export
data_generator <- function(ref_claim = 200000,
                           time_unit = 1/360,
                           years = 4,
                           random_seed = 1964,
                           period_exposure = 200,
                           period_frequency = 0.2,
                           scenario = 1){

  set.seed(random_seed)

  scenario <- pkg.env$check_scenario(scenario)

  parameters <- list(ref_claim=ref_claim,
                     time_unit=time_unit,
                     years=years,
                     yearly_exposure=period_exposure,
                     yearly_frequency=period_frequency)



  if(scenario == 0){ return(do.call(pkg.env$scenario0_simulator,parameters)
                            )}

  if(scenario == 1){ return(do.call(pkg.env$scenario1_simulator,parameters)
                            )}

  if(scenario == 2){ return(do.call(pkg.env$scenario2_simulator,parameters)
                            )}

  if(scenario == 3){ return(do.call(pkg.env$scenario3_simulator,parameters)
                            )}

  if(scenario == 4){ return(do.call(pkg.env$scenario4_simulator,parameters)
  )}

  if(scenario == 5){ return(do.call(pkg.env$scenario5_simulator,parameters)
  )}

  if(scenario == 6){ return(do.call(pkg.env$scenario6_simulator,parameters)
  )}

}









