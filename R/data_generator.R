#' Individual data generator
#'
#' This function generates the monthly individual claims data in the accompanying methodological paper using the \code{SynthETIC} package.
#' This simple function allows to simulate from a sand-box to test out the \code{ReSurv} approach.
#' Some parameters of the simulation can be changed.
#'
#' @param ref_claim \code{integer}, reference claim size.
#' @param time_unit \code{numeric}, output time unit.
#' @param years \code{integer}, number of years to be simulated.
#' @param yearly_exposure \code{integer}, volume underwritten each year.
#' @param yearly_frequency \code{numeric}, yearly frequency.
#' @param scenario \code{numeric}, one of the scenarios shown in the accompanying paper.
#'
#' @import SynthETIC
#'
#' @examples
#' ## Not run
#' input_data <- data_generator(random_seed = 1964)
#'
#'
#'
#' @return Individual claims data. It contains the following columns:
#' \itemize{
#' \item{\code{AT}: Accident month in continuous time.}
#' \item{\code{AM}: Accident month.}
#' \item{\code{RT}: Reporting month in continuous time.}
#' \item{\code{RM}: Reporting month.}
#' \item{\code{DT}: Development month in continuous time.}
#' \item{\code{DM}: Development month.}
#' \item{\code{DT_rev}:Development month in continuous reverse time.}
#' \item{\code{DM_rev}: Development month in reverse time. }
#' \item{\code{TR}: Truncation time. }
#' \item{\code{I}: Event indicator.  }
#' }
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
#'
#' @export
data_generator <- function(ref_claim = 200000,
                           time_unit = 1/12,
                           years = 4,
                           random_seed=1964,
                           yearly_exposure= 6000,
                           yearly_frequency=0.2,
                           scenario=1){

  set.seed(random_seed)

  parameters <- list(ref_claim=ref_claim,
                     time_unit=time_unit,
                     years=years,
                     yearly_exposure=yearly_exposure,
                     yearly_frequency=yearly_frequency)


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
}









