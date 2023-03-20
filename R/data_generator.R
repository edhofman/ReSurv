#' Individual data generator
#'
#' This function generates individual claim data.
#'
#' @param ref_claim reference claim size
#' @param time_unit output time unit
#' @param years simulated number of years
#' @param yearly_exposure volume underwritten
#' @param yearly_frequency yearly frequency
#'
#' @import SynthETIC
#'
#' @return individual claim data
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
                           yearly_frequency=0.2){

  set.seed(random_seed)
  # Letting Claimtype 1 have decreasing frequency ----------------------------------------------
  # set_parameters(ref_claim = 200000, time_unit = 1/4)
  # ref_claim = 200000
  # time_unit <- 1/12
  # years <- 4
  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))

  #Decreasing the exposure, and hence lowering the claims occurred
  E_1 <- c(rep(yearly_exposure, I)) + seq(from = 0, by = -100, length = I)
  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E_1*2, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

  #Claim sizes
  #Not needed for reporting analysis
  S_df <- function(s) {
    # truncate and rescale
    if (s < 30) {
      return(0)
    } else {
      p_trun <- pnorm(s^0.2, 9.5, 3) - pnorm(30^0.2, 9.5, 3)
      p_rescaled <- p_trun/(1 - pnorm(30^0.2, 9.5, 3))
      return(p_rescaled)
    }
  }
  claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
                              simfun = S_df, type = "p", range = c(0, 1e24))
  claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
                              simfun = S_df, type = "p", range = c(0, 1e24))

  RTFWD_inverse <- function(n, alpha, beta, lambda, k,b){
    U<-runif(n)
    (-log(U)/(beta^alpha*lambda^(alpha*k))+b^(-alpha*k))^(1/(-alpha*k))
  }

  notidel_param_0 <- function(claim_size, occurrence_period) {


    #c(scale = 5/(exp(-25.7845)^(1/12.15343)),
    #  shape = 12.15343) #cv 0.1
    c(alpha=0.5,
      beta=2,
      lambda=0.1*exp(1.15129)^(1/0.5),
      k=1,
      b=years / time_unit)

  }

  notidel_param_1 <- function(claim_size, occurrence_period) {

    c(alpha=0.5,
      beta=2,
      lambda=0.1*exp(1.95601)^(1/0.5),
      k=1,
      b=years / time_unit)

  }


  ## output
  # simulate notification delays from the transformed gamma
  notidel_claim_type_0 <- claim_notification(n_vector_0, claim_sizes_0,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_0)

  notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_1)


  ct0 <- tibble(AT = unlist(occurrence_times_0),
                RT = unlist(occurrence_times_0) + unlist(notidel_claim_type_0),
                claim_type = 0)

  ct1 <- tibble(AT = unlist(occurrence_times_1),
                RT = unlist(occurrence_times_1) + unlist(notidel_claim_type_1),
                claim_type = 1)

  simulated_dataframe_RM_CT <- ct0 %>%  bind_rows(ct1) %>%
    mutate(
      claim_number = row_number(),
    )  %>%  mutate(
      AM = ceiling(AT),
      RM = ceiling(RT),
      DT = RT-AT,
      DM = RM-AM+1,
      DM_rev = years/time_unit - DM+1,
      DT_rev = years/time_unit - DT,
      TR = AM-1, #just setting truncation to max year simulated. and accounting for
      I=1
    ) %>%
    select(claim_number, AT, RT, claim_type, AM, RM, DT, DM, DM_rev, DT_rev, TR, I)

  simulated_dataframe_RM_CT

}
