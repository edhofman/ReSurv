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
                           time_unit = 1 / 360,
                           years = 4,
                           random_seed = 1964,
                           period_exposure = 200,
                           period_frequency = 0.2,
                           scenario = 1) {

  set.seed(random_seed)

  ## ------------------------------------------------------------------
  ## Inline check_scenario()
  ## ------------------------------------------------------------------

  available_scenarios <- c(0, 1, 2, 3, 4, 5, 6)
  available_scenario_char <- c(
    "alpha",
    "beta",
    "gamma",
    "delta",
    "epsilon",
    "zeta",
    "eta"
  )

  if (is.numeric(scenario)) {

    tmp <- scenario %in% available_scenarios

    if (!tmp) {
      stop(
        "Scenario must be one of 'alpha','beta','gamma','delta','epsilon', 'zeta','eta'."
      )
    }
  }

  if (is.character(scenario)) {

    tmp <- scenario %in% available_scenario_char

    if (!tmp) {
      stop(
        "Scenario must be one of 'alpha','beta','gamma','delta','epsilon', 'zeta','eta'."
      )
    }

    input.pos <- which(scenario == available_scenario_char)

    scenario <- available_scenarios[input.pos]
  }

  yearly_exposure <- period_exposure
  yearly_frequency <- period_frequency

  ## ------------------------------------------------------------------
  ## Local utility functions
  ## ------------------------------------------------------------------

  S_df_local <- function(s) {

    if (s < 30) {
      return(0)
    } else {
      p_trun <- stats::pnorm(s^0.2, 9.5, 3) -
        stats::pnorm(30^0.2, 9.5, 3)

      p_rescaled <- p_trun / (1 - stats::pnorm(30^0.2, 9.5, 3))

      return(p_rescaled)
    }
  }

  RTFWD_inverse_local <- function(n, alpha, beta, lambda, k, b) {
    U <- stats::runif(n)

    (
      -log(U) /
        (beta^alpha * lambda^(alpha * k)) +
        b^(-alpha * k)
    )^(1 / (-alpha * k))
  }

  period_function_local <- function(x) {

    tmp <- floor((x - 1) / 30)

    if ((tmp %% 12) %in% c(2, 3, 4)) {
      return(-0.3)
    }

    if ((tmp %% 12) %in% c(5, 6, 7)) {
      return(0.4)
    }

    if ((tmp %% 12) %in% c(8, 9, 10)) {
      return(-0.7)
    }

    if ((tmp %% 12) %in% c(11, 0, 1)) {
      return(0.1)
    }
  }

  notidel_param_0_local <- function(claim_size,
                                    occurrence_period,
                                    scenario,
                                    years,
                                    time_unit) {

    if (scenario %in% c(0, 1, 2)) {
      return(c(
        alpha = 0.5,
        beta = 2 * 30,
        lambda = 0.1 * exp(1.15129)^(1 / 0.5),
        k = 1,
        b = years / time_unit
      ))
    }

    if (scenario == 3) {
      return(c(
        alpha = 0.5,
        beta = 2 * 30,
        lambda = 0.1 *
          exp(1.15129 + period_function_local(ceiling(occurrence_period)))^(1 / 0.5),
        k = 1,
        b = years / time_unit
      ))
    }

    if (scenario == 4) {
      return(c(
        alpha = 0.5,
        beta = (2 + 0.5 * 1.15129) * 30,
        lambda = 0.1 * exp(1.15129)^(1 / 0.5) + 0.5 * 1.15129,
        k = 1,
        b = years / time_unit
      ))
    }
  }

  notidel_param_1_local <- function(claim_size,
                                    occurrence_period,
                                    scenario,
                                    years,
                                    time_unit) {

    if (scenario %in% c(0, 1)) {
      return(c(
        alpha = 0.5,
        beta = 2 * 30,
        lambda = 0.1 * exp(1.95601)^(1 / 0.5),
        k = 1,
        b = years / time_unit
      ))
    }

    if (scenario == 2) {
      return(c(
        alpha = 0.5,
        beta = 2 * 30,
        lambda = 0.1 *
          exp(1.95601 - 0.02120623 * sqrt(ceiling(occurrence_period)))^(1 / 0.5),
        k = 1,
        b = years / time_unit
      ))
    }

    if (scenario == 3) {
      return(c(
        alpha = 0.5,
        beta = 2 * 30,
        lambda = 0.1 *
          exp(1.95601 + period_function_local(ceiling(occurrence_period)))^(1 / 0.5),
        k = 1,
        b = years / time_unit
      ))
    }

    if (scenario == 4) {
      return(c(
        alpha = 0.5,
        beta = (2 + 0.5 * 1.95601) * 30,
        lambda = 0.1 * exp(1.95601)^(1 / 0.5) + 0.5 * 1.95601,
        k = 1,
        b = years / time_unit
      ))
    }
  }

  notification_delay_scenario5_local <- function(x) {

    pv <- as.numeric(x[["property_value"]]) / 10
    bu <- x[["business_use"]]

    a <- 1 + 1 / pv
    d <- 1 - (1 + (bu == "Y")) / 10

    target_mean <- 850

    num <- gamma(a + 1 / d)
    denom <- target_mean * gamma(a)

    rate <- (num / denom)^d

    actuar::rtrgamma(
      1,
      shape1 = a,
      shape2 = d,
      rate = rate
    )
  }

  notification_delay_scenario6_local <- function(x) {

    ap <- as.numeric(x[["AP"]])
    bu <- x[["business_use"]]

    actuar::rtrgamma(
      1,
      shape1 = 1 + 1 / ap + sin(ap) / (ap^2) + cos(ap) / (ap^3),
      shape2 = 1 - (1 + (bu == "Y")) / 10 + (1 - (bu == "N")) / 100,
      rate = 0.2
    )
  }

  finalize_two_claim_types <- function(occurrence_times_0,
                                       occurrence_times_1,
                                       notidel_claim_type_0,
                                       notidel_claim_type_1,
                                       years,
                                       time_unit) {

    ct0 <- data.table::data.table(
      AT = unlist(occurrence_times_0),
      RT = unlist(occurrence_times_0) + unlist(notidel_claim_type_0),
      claim_type = 0
    )

    ct1 <- data.table::data.table(
      AT = unlist(occurrence_times_1),
      RT = unlist(occurrence_times_1) + unlist(notidel_claim_type_1),
      claim_type = 1
    )

    simulated_dataframe_RM_CT <- data.table::rbindlist(list(ct0, ct1))

    simulated_dataframe_RM_CT[
      ,
      `:=`(
        claim_number = seq_len(.N),
        AP = ceiling(AT),
        RP = ceiling(RT),
        DT = RT - AT,
        DP = ceiling(RT) - ceiling(AT) + 1,
        DP_rev = years / time_unit -
          (ceiling(RT) - ceiling(AT) + 1) + 1,
        DT_rev = years / time_unit - (RT - AT),
        TR = ceiling(AT) - 1,
        I = 1
      )
    ]

    simulated_dataframe_RM_CT[
      ,
      .(claim_number, claim_type, AP, RP)
    ]
  }

  simulate_two_claim_type_scenario <- function(scenario_id) {

    I <- years / time_unit
    E <- c(rep(yearly_exposure, I))
    lambda <- c(rep(yearly_frequency, I))

    ## Important: the original scenario 0, 2, 3, and 4 simulators
    ## perform these two random operations even though the objects are
    ## not subsequently used. They must remain here to preserve the
    ## exact RNG stream and therefore the exact simulated output.

    if (scenario_id %in% c(0, 2, 3, 4)) {
      n_vector <- claim_frequency(I = I, E = E, freq = lambda)
      occurrence_times <- claim_occurrence(frequency_vector = n_vector)
    }

    if (scenario_id == 1) {

      E_1 <- c(rep(yearly_exposure, I)) +
        round(seq(from = 0, by = -0.1, length = I))

      n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
      n_vector_1 <- claim_frequency(I = I, E = E_1, freq = lambda)

    } else {

      n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
      n_vector_1 <- claim_frequency(I = I, E = E, freq = lambda)
    }

    occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
    occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

    claim_sizes_0 <- claim_size(
      frequency_vector = n_vector_0,
      simfun = S_df_local,
      type = "p",
      range = c(0, 1e24)
    )

    claim_sizes_1 <- claim_size(
      frequency_vector = n_vector_1,
      simfun = S_df_local,
      type = "p",
      range = c(0, 1e24)
    )

    notidel_claim_type_0 <- claim_notification(
      n_vector_0,
      claim_sizes_0,
      rfun = RTFWD_inverse_local,
      paramfun = notidel_param_0_local,
      scenario = scenario_id,
      years = years,
      time_unit = time_unit
    )

    notidel_claim_type_1 <- claim_notification(
      n_vector_1,
      claim_sizes_1,
      rfun = RTFWD_inverse_local,
      paramfun = notidel_param_1_local,
      scenario = scenario_id,
      years = years,
      time_unit = time_unit
    )

    finalize_two_claim_types(
      occurrence_times_0 = occurrence_times_0,
      occurrence_times_1 = occurrence_times_1,
      notidel_claim_type_0 = notidel_claim_type_0,
      notidel_claim_type_1 = notidel_claim_type_1,
      years = years,
      time_unit = time_unit
    )
  }

  simulate_gamma_scenario <- function() {

    max_day <- as.integer(years / time_unit)
    break_day <- floor(max_day / 2)

    beta0 <- 1.15129
    beta1 <- 1.95601
    beta_break <- -1.2

    age_mean <- 45
    age_sd <- 4
    age_min <- 32
    age_max <- 58
    age_cut1 <- 42
    age_cut2 <- 49
    age_low <- -0.10
    age_mid <- 0.00
    age_high <- 0.25

    lambda_n <- as.numeric(period_exposure) * as.numeric(period_frequency)

    grid <- data.table::data.table(
      AP = rep(seq_len(max_day), times = 2L),
      claim_type = rep(c(0L, 1L), each = max_day)
    )

    grid[, n_claims := stats::rpois(.N, lambda = lambda_n)]

    claims <- grid[n_claims > 0L, .(
      AP = rep(AP, n_claims),
      claim_type = rep(claim_type, n_claims)
    )]

    if (nrow(claims) == 0L) {
      stop("No claims simulated. Increase period_exposure or period_frequency.")
    }

    claims[, AT := AP - 1 + stats::runif(.N)]

    claims[
      ,
      `:=`(
        age_cont = stats::rnorm(.N, mean = age_mean, sd = age_sd)
      )
    ]
    claims[, age := round(pmin(pmax(age_cont, age_min), age_max))]

    claims[
      ,
      age_effect := data.table::fcase(
        age < age_cut1, age_low,
        age < age_cut2, age_mid,
        default = age_high
      )
    ]

    claims[
      ,
      phi := beta0 * as.numeric(claim_type == 0L) +
        beta1 * as.numeric(claim_type == 1L) +
        beta_break * as.numeric(claim_type == 1L & AP > break_day) +
        age_effect
    ]

    alpha <- 0.5
    beta <- 2 * 30
    k <- 1
    b <- max_day

    lambda_delay <- 0.1 * exp(claims$phi)^(1 / alpha)

    claims[
      ,
      delay := RTFWD_inverse_local(
        n = .N,
        alpha = alpha,
        beta = beta,
        lambda = lambda_delay,
        k = k,
        b = b
      )
    ]

    claims[
      ,
      `:=`(
        RT = AT + delay,
        claim_number = seq_len(.N)
      )
    ]
    claims[, RP := as.integer(ceiling(RT))]

    claims[
      ,
      .(
        claim_number = claim_number,
        claim_type = claim_type,
        age = age,
        AP = as.integer(AP),
        RP = RP
      )
    ]
  }

  ## ------------------------------------------------------------------
  ## Scenarios 0--4
  ## ------------------------------------------------------------------

  if (scenario == 2) {
    return(simulate_gamma_scenario())
  }

  if (scenario %in% c(0, 1, 3, 4)) {
    return(simulate_two_claim_type_scenario(scenario))
  }

  ## ------------------------------------------------------------------
  ## Scenario 5
  ## ------------------------------------------------------------------

  if (scenario == 5) {

    I <- years / time_unit
    E <- c(rep(floor(yearly_exposure), I))
    lambda <- c(rep(yearly_frequency, I))

    E_1 <- c(rep(floor(yearly_exposure), I)) +
      round(seq(from = 0, by = -0.1, length = I))

    n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
    n_vector_1 <- claim_frequency(I = I, E = E_1, freq = lambda)

    occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
    occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

    claim_sizes <- claim_size(frequency_vector = c(n_vector_0, n_vector_1))
    n_of_claims <- length(unlist(claim_sizes))

    bu_covariates_dataset <- data.table::data.table(
      claim_number = seq_len(n_of_claims),
      business_use = c(
        rep("Y", sum(n_vector_0)),
        rep("N", sum(n_vector_1))
      )
    )

    age_range <- 50:55
    probabilties_age <- rep(0.01, length(age_range))
    probabilties_age <- probabilties_age / sum(probabilties_age)

    covariates_dataset <- data.table::data.table(
      claim_number = seq_len(n_of_claims),
      age = sample(
        age_range,
        n_of_claims,
        replace = TRUE,
        prob = probabilties_age
      ),
      property_value = stats::rlnorm(
        n_of_claims,
        meanlog = 3.034513,
        sdlog = 0.4087569
      )
    )

    covariates_dataset <- covariates_dataset[
      bu_covariates_dataset,
      on = "claim_number"
    ]

    ## Keep apply() here to preserve the exact old row-wise RNG behaviour.
    ## This is not a dplyr dependency.

    rdelay <- apply(
      X = covariates_dataset,
      MARGIN = 1,
      FUN = notification_delay_scenario5_local
    )

    rdelay <- pmin(rdelay, years / time_unit)

    dt_dates <- data.table::data.table(
      claim_number = seq_len(n_of_claims),
      AP = ceiling(c(
        unlist(occurrence_times_0),
        unlist(occurrence_times_1)
      )),
      RP = ceiling(c(
        unlist(occurrence_times_0),
        unlist(occurrence_times_1)
      ) + rdelay)
    )

    dt <- dt_dates[
      covariates_dataset,
      on = "claim_number"
    ]

    return(dt)
  }

  ## ------------------------------------------------------------------
  ## Scenario 6
  ## ------------------------------------------------------------------

  if (scenario == 6) {

    I <- years / time_unit
    E <- c(rep(yearly_exposure, I))
    lambda <- c(rep(yearly_frequency, I))

    n_vector <- claim_frequency(I = I, E = E, freq = lambda)
    occurrence_times <- claim_occurrence(frequency_vector = n_vector)
    claim_sizes <- claim_size(frequency_vector = n_vector)

    n_of_claims <- length(unlist(claim_sizes))

    covariates_dataset <- data.table::data.table(
      claim_number = seq_len(n_of_claims),
      AP = ceiling(unlist(occurrence_times)),
      business_use = sample(c("Y", "N"), n_of_claims, replace = TRUE)
    )

    ## Keep apply() here to preserve the exact old row-wise RNG behaviour.
    ## This is not a dplyr dependency.

    rdelay <- apply(
      X = covariates_dataset,
      MARGIN = 1,
      FUN = notification_delay_scenario6_local
    )

    rdelay <- pmin(rdelay, years / time_unit)

    dt_dates <- data.table::data.table(
      claim_number = seq_len(n_of_claims),
      RP = ceiling(unlist(occurrence_times) + rdelay)
    )

    dt <- dt_dates[
      covariates_dataset,
      on = "claim_number"
    ]

    return(dt)
  }
}





