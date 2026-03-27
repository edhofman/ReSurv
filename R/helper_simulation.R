# Simulation helper functions
#
# Utility functions and scenario simulators for data generation.
#
# @importFrom actuar rztpois rtrgamma
# @importFrom stats runif pnorm
# Utils individual claims generators

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

RTFWD_inverse <- function(n, alpha, beta, lambda, k,b){
  U<-runif(n)
  (-log(U)/(beta^alpha*lambda^(alpha*k))+b^(-alpha*k))^(1/(-alpha*k))
}

period_function <-function(x){

  "
  Add monthly seasonal effect starting from daily input.

  "

  tmp <- floor((x-1)/30)

  if((tmp%%12) %in% (c(2,3,4))){
    return(-0.3)
  }
  if((tmp%%12) %in% (c(5,6,7))){
    return(0.4)
  }
  if((tmp%%12) %in% (c(8,9,10))){
    return(-0.7)
  }
  if((tmp%%12) %in% (c(11,0,1))){ #0 instead of 12
    return(0.1)
  }
}

notidel_param_0 <- function(claim_size,
                            occurrence_period,
                            scenario,
                            years,
                            time_unit) {

  if(scenario %in% c(0,1,2)){

    return(c(alpha=0.5,
             beta=2*30,
             lambda=0.1*exp(1.15129)^(1/0.5),
             k=1,
             b=years / time_unit))

  }

  if(scenario==3){

    return(c(alpha=0.5,
             beta=2*30,
             lambda=0.1*exp(1.15129+period_function(ceiling(occurrence_period)))^(1/0.5),
             k=1,
             b=years / time_unit))

  }

  if(scenario==4){

    return(c(alpha=0.5,
             beta=(2+0.5*1.15129)*30,
             lambda=0.1*exp(1.15129)^(1/0.5)+0.5*1.15129,
             k=1,
             b=years / time_unit))

  }



}



notidel_param_1 <- function(claim_size,
                            occurrence_period,
                            scenario,
                            years,
                            time_unit) {

  if(scenario%in%c(0,1)){
    return(c(alpha=0.5,
             beta=2*30,
             lambda=0.1*exp(1.95601)^(1/0.5),
             k=1,
             b=years / time_unit))}

  if(scenario==2){
    return(c(alpha=0.5,
             beta=2*30,
             lambda=0.1*exp(1.95601-0.02120623*sqrt(ceiling(occurrence_period)))^(1/0.5),
             k=1,
             b=years / time_unit))}

  if(scenario==3){
    return(c(alpha=0.5,
             beta=2*30,
             lambda=0.1*exp(1.95601+period_function(ceiling(occurrence_period)) )^(1/0.5),
             k=1,
             b=years / time_unit))}

  if(scenario==4){
    return(c(alpha=0.5,
             beta=(2+0.5*1.95601)*30,
             lambda=0.1*exp(1.95601)^(1/0.5)+0.5*1.95601,
             k=1,
             b=years / time_unit))}



}


## Internal functions for W18 comparison

notification_delay_scenario5 <- function(x) {
  pv <- as.numeric(x[['property_value']]) / 10
  bu <- x[['business_use']] # "Y" or "N"

  a <- 1 + 1 / pv
  d <- 1 - (1 + (bu == "Y")) / 10

  # Target mean
  target_mean <- 850

  # Compute required rate
  num <- gamma(a + 1 / d)
  denom <- target_mean * gamma(a)
  rate <- (num / denom)^d

  out <- rtrgamma(1, shape1 = a, shape2 = d, rate = rate)
  return(out)
}


notification_delay_scenario6 <- function(x) {
  ap <- as.numeric(x[['AP']])
  bu <- x[['business_use']]



  # out <- rweibull(1,shape = 1 + 1 / pv, scale=1 - (1 + (bu == "Y")) / 10)
  out <- rtrgamma(
    1,
    shape1 = 1 + 1 / ap + sin(ap)/(ap^2)+cos(ap)/(ap^3),
    shape2 = 1 - (1 + (bu == "Y")) / 10 + (1-(bu == "N"))/100,
    rate = .2
  )

  return(out)
}


## Data generator ----

pkg.env$check_scenario <- function(scenario){

  available_scenarios <- c(0,1,2,3,4,5,6)
  available_scenario_char <- c('alpha','beta','gamma','delta','epsilon','zeta','eta')

  if(is.numeric(scenario)){

    tmp <- scenario %in% available_scenarios

    if(!tmp){stop("Scenario must be one of 'alpha','beta','gamma','delta','epsilon', 'zeta','eta'.")}
  }


  if(is.character(scenario)){

    tmp <- scenario %in% available_scenario_char

    if(!tmp){stop("Scenario must be one of 'alpha','beta','gamma','delta','epsilon', 'zeta','eta'.")}

    input.pos <- which(scenario==available_scenario_char)

    scenario <- available_scenarios[input.pos]

  }

  return(scenario)

}


# Superimposed inflation:
# 1) With respect to occurrence "time" (continuous scale)
SI_occurrence <- function(occurrence_time, claim_size) {
  if (occurrence_time <= 20 / 4 / time_unit) {1}
  else {1 - 0.4*max(0, 1 - claim_size/(0.25 * ref_claim))}
}
# 2) With respect to payment "time" (continuous scale)
# -> compounding by user-defined time unit
SI_payment <- function(payment_time, claim_size) {
  period_rate <- (1 + 0.30)^(time_unit) - 1
  beta <- period_rate * max(0, 1 - claim_size/ref_claim)
  (1 + beta)^payment_time
}


pkg.env$scenario0_simulator <- function(ref_claim,
                                        time_unit,
                                        years,
                                        random_seed,
                                        yearly_exposure,
                                        yearly_frequency){

  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=0


  #Frequency simulation
  n_vector <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times <- claim_occurrence(frequency_vector = n_vector)


  #idea for simulating dependency with claim types
  #claim_types <- c(round(runif(sum(n_vector),0,1) ) )


  # No difference, same simulations with different reporting mean and number numbers ------------------------


  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

  #Doesn't matter since only looking at observed numbers, but needed for the other functions
  claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
                              simfun = S_df, type = "p", range = c(0, 1e24))
  claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
                              simfun = S_df, type = "p", range = c(0, 1e24))


  ## output
  # simulate notification delays from the transformed gamma
  notidel_claim_type_0 <- claim_notification(n_vector_0,
                                             claim_sizes_0,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_0,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

  notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_1,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

  # graphically compare the result with the default Weibull distribution

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
      AP = ceiling(AT),
      RP = ceiling(RT),
      DT = RT-AT,
      DP = RP-AP+1,
      DP_rev = years/time_unit - DP+1,
      DT_rev = years/time_unit - DT,
      TR = AP-1, #just setting truncation to max year simulated. and accounting for
      I=1
    ) %>%
    select(claim_number,
           # AT,
           # RT,
           claim_type,
           AP,
           RP#,
           # DT,
           # DP,
           # DP_rev,
           # DT_rev,
           # TR,
           # I
           ) %>%
    as.data.frame()

  return(simulated_dataframe_RM_CT)

}

pkg.env$scenario1_simulator <- function(ref_claim,
                                        time_unit,
                                        years,
                                        random_seed,
                                        yearly_exposure,
                                        yearly_frequency){


  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=1

  #Decreasing the exposure, and hence lowering the claims occurred
  E_1 <- c(rep(yearly_exposure, I)) + round(seq(from = 0, by = -.1, length = I))# now adjusted for days: monthly code was seq(from = 0, by = -100, length = I)
  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E_1, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)




  claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
                              simfun = S_df, type = "p", range = c(0, 1e24))
  claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
                              simfun = S_df, type = "p", range = c(0, 1e24))



  ## output
  # simulate notification delays from the transformed gamma
  notidel_claim_type_0 <- claim_notification(n_vector_0, claim_sizes_0,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_0,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

  notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_1,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)


  ct0 <- tibble(AT = unlist(occurrence_times_0),
                RT = unlist(occurrence_times_0) + unlist(notidel_claim_type_0),
                claim_type = 0)

  ct1 <- tibble(AT = unlist(occurrence_times_1),
                RT = unlist(occurrence_times_1) + unlist(notidel_claim_type_1),
                claim_type = 1)

  simulated_dataframe_RM_CT <- ct0 %>%
    bind_rows(ct1) %>%
    mutate(
      claim_number = row_number(),
    )  %>%  mutate(
      AP = ceiling(AT),
      RP = ceiling(RT),
      DT = RT-AT,
      DP = RP-AP+1,
      DP_rev = years/time_unit - DP+1,
      DT_rev = years/time_unit - DT,
      TR = AP-1, #just setting truncation to max year simulated. and accounting for
      I=1
    ) %>%
    select(claim_number,
           # AT,
           # RT,
           claim_type,
           AP,
           RP#,
           # DT,
           # RP,
           # DP_rev,
           # DT_rev,
           # TR,
           # I
           ) %>% as.data.frame()

  #simulated_dataframe_RM_CT

  # setDT(simulated_dataframe_RM_CT)

  return(simulated_dataframe_RM_CT)

}


pkg.env$scenario2_simulator <- function(ref_claim,
                                        time_unit,
                                        years,
                                        random_seed,
                                        yearly_exposure,
                                        yearly_frequency){


  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=2

  #Frequency simulation
  n_vector <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times <- claim_occurrence(frequency_vector = n_vector)



  #idea for simulating dependency with claim types
  #claim_types <- c(round(runif(sum(n_vector),0,1) ) )


  # No difference, same simulations with different reporting mean and number numbers ------------------------

  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

  #Doesn't matter since only looking at observed numbers, but needed for the other functions
  claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
                              simfun = S_df, type = "p", range = c(0, 1e24))
  claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
                              simfun = S_df, type = "p", range = c(0, 1e24))





  ## output
  # simulate notification delays from the transformed gamma
  notidel_claim_type_0 <- claim_notification(n_vector_0, claim_sizes_0,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_0,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

  notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_1,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)


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
      AP = ceiling(AT),
      RP = ceiling(RT),
      DT = RT-AT,
      DP = RP-AP+1,
      DP_rev = years/time_unit - DP+1,
      DT_rev = years/time_unit - DT,
      TR = AP-1, #just setting truncation to max year simulated. and accounting for
      I=1
    ) %>%
    select(claim_number,
           # AT,
           # RT,
           claim_type,
           AP,
           RP)#,
           # DT,
           # DP,
           # DP_rev, DT_rev, TR, I)

  # simulated_dataframe_RM_CT
  # setDT(simulated_dataframe_RM_CT)

  return(simulated_dataframe_RM_CT)

  }

pkg.env$scenario3_simulator <- function(ref_claim,
                                        time_unit,
                                        years,
                                        random_seed,
                                        yearly_exposure,
                                        yearly_frequency){


  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=3

  #Frequency simulation
  n_vector <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times <- claim_occurrence(frequency_vector = n_vector)


  #idea for simulating dependency with claim types
  #claim_types <- c(round(runif(sum(n_vector),0,1) ) )


  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

  #Doesn't matter since only looking at observed numbers, but needed for the other functions
  claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
                              simfun = S_df, type = "p", range = c(0, 1e24))
  claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
                              simfun = S_df, type = "p", range = c(0, 1e24))


  # simulate notification delays from the transformed gamma
  notidel_claim_type_0 <- claim_notification(n_vector_0, claim_sizes_0,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_0,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

  notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_1,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)


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
      AP = ceiling(AT),
      RP = ceiling(RT),
      DT = RT-AT,
      DP = RP-AP+1,
      DP_rev = years/time_unit - DP+1,
      DT_rev = years/time_unit - DT,
      TR = AP-1, #just setting truncation to max year simulated. and accounting for
      I=1
    ) %>%
    select(claim_number,
           #AT,
           #RT,
           claim_type,
           AP,
           RP) %>% as.data.frame()#,
           #DT, DP, DP_rev, DT_rev, TR, I)

  # simulated_dataframe_RM_CT

  # setDT(simulated_dataframe_RM_CT)

  return(simulated_dataframe_RM_CT)



}

pkg.env$scenario4_simulator <- function(ref_claim,
                                        time_unit,
                                        years,
                                        random_seed,
                                        yearly_exposure,
                                        yearly_frequency){



  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=4

  #Frequency simulation
  n_vector <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times <- claim_occurrence(frequency_vector = n_vector)

  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)

  #Doesn't matter since only looking at observed numbers, but needed for the other functions
  claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
                              simfun = S_df, type = "p", range = c(0, 1e24))
  claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
                              simfun = S_df, type = "p", range = c(0, 1e24))


  #sum(unlist(notidel_claim_type_0)>48)/sum(unlist(notidel_claim_type_0)>0)
  ## output
  # simulate notification delays from the transformed gamma
  notidel_claim_type_0 <- claim_notification(n_vector_0, claim_sizes_0,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_0,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

  notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
                                             rfun = RTFWD_inverse,
                                             paramfun = notidel_param_1,
                                             scenario=scenario,
                                             years=years,
                                             time_unit=time_unit)

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
      AP = ceiling(AT),
      RP = ceiling(RT),
      DT = RT-AT,
      DP = RP-AP+1,
      DP_rev = years/time_unit - DP+1,
      DT_rev = years/time_unit - DT,
      TR = AP-1, #just setting truncation to max year simulated. and accounting for
      I=1
    ) %>%
    select(claim_number, #AT, RT,
           claim_type,
           AP,
           RP#,
           #DT, DP, DP_rev, DT_rev, TR, I


           ) %>% as.data.frame()

  # simulated_dataframe_RM_CT

  # setDT(simulated_dataframe_RM_CT)

  return(simulated_dataframe_RM_CT)




}

generate_proportions <- function(I) {
       raw <- rexp(I, rate = runif(I, min = 0.1, max = 2))  # exponential with varying rates ÔåÆ irregular
       proportions <- raw / sum(raw)
       return(proportions)
   }

## Wuethrich 18 comparisons

pkg.env$scenario5_simulator <- function(ref_claim=200000,
                                        time_unit=1/4,
                                        years=10,
                                        random_seed,
                                        yearly_exposure=120000,
                                        yearly_frequency=0.08){


  I <- years / time_unit
  E <- c(rep(floor(yearly_exposure), I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=5

  #Frequency simulation -- business use 0
  # n_vector <- claim_frequency(I = I, E = E, freq = lambda)
  # occurrence_times <- claim_occurrence(frequency_vector = n_vector)


  #Decreasing the exposure, and hence lowering the claims occurred -- business use 1
  E_1 <- c(rep(floor(yearly_exposure), I)) + round(seq(from = 0, by = -.1, length = I))# now adjusted for days: monthly code was seq(from = 0, by = -100, length = I)
  #Frequency simulation
  n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
  n_vector_1 <- claim_frequency(I = I, E = E_1, freq = lambda)
  occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
  occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)


  claim_sizes <- claim_size(frequency_vector = c(n_vector_0,n_vector_1))
  n_of_claims <- length(unlist(claim_sizes))

  bu_covariates_dataset <- data.frame(
    "claim_number"=1:n_of_claims,
    "business_use" = c(rep("Y",sum(n_vector_0)),
                       rep("N",sum(n_vector_1)))

    )


  age_range <- 50:55
  probabilties_age <- rep(.01,length(age_range))

  probabilties_age <- probabilties_age/sum(probabilties_age)

  covariates_dataset <- data.frame(
    "claim_number"=1:n_of_claims,
    "age" = sample(age_range,n_of_claims,replace=TRUE,prob=probabilties_age),
    "property_value"= rlnorm(n_of_claims, meanlog = 3.034513, sdlog = 0.4087569)#,
    # "business_use" = sample(c("Y","N"),n_of_claims,replace = TRUE)
  )

  covariates_dataset <-merge(covariates_dataset,
                             bu_covariates_dataset,
                             by="claim_number")

  rdelay = apply(FUN = notification_delay_scenario5 ,
                 covariates_dataset,
                 MARGIN = 1)



  rdelay = pmin(rdelay, years / time_unit)


  dt_dates <- data.frame(
    claim_number=1:n_of_claims,
    AP=ceiling(c(unlist(occurrence_times_0),
                 unlist(occurrence_times_1))),
    RP=ceiling(c(unlist(occurrence_times_0),
                 unlist(occurrence_times_1))+rdelay))


  dt <- merge(dt_dates,covariates_dataset,by.x="claim_number",by.y="claim_number",all=TRUE)



  return(dt)




}

pkg.env$scenario6_simulator <- function(ref_claim=200000,
                                        time_unit=1/4,
                                        years=10,
                                        random_seed,
                                        yearly_exposure=120000,
                                        yearly_frequency=0.08){



  I <- years / time_unit
  E <- c(rep(yearly_exposure, I))
  lambda <- c(rep(yearly_frequency, I))
  scenario=5

  #Frequency simulation
  n_vector <- claim_frequency(I = I, E = E, freq = lambda)
  occurrence_times <- claim_occurrence(frequency_vector = n_vector)
  claim_sizes <- claim_size(frequency_vector = n_vector)

  n_of_claims <- length(unlist(claim_sizes))

  age_range <- 50:55
  probabilties_age <- rep(.01,length(age_range))
  # probabilties_age[age_range >= 40 & age_range <= 45] <- .2
  # probabilties_age[age_range >= 20 & age_range <= 30] <- .1
  # probabilties_age[age_range >= 30 & age_range <= 39] <- .15
  # probabilties_age[age_range > 45 & age_range <= 55] <- .3

  probabilties_age <- probabilties_age/sum(probabilties_age)

  covariates_dataset <- data.frame(
    "claim_number"=1:n_of_claims,
    "AP" = ceiling(unlist(occurrence_times)),
    "business_use" = sample(c("Y","N"),n_of_claims,replace = TRUE)
  )

  rdelay = apply(FUN = notification_delay_scenario6 ,
                 covariates_dataset,
                 MARGIN = 1)

  rdelay = pmin(rdelay, years / time_unit)


  dt_dates <- data.frame(
    claim_number=1:n_of_claims,
    # AP=ceiling(unlist(occurrence_times)), (we have it already in the covariates.)
    RP=ceiling(unlist(occurrence_times)+rdelay))


  dt <- merge(dt_dates,covariates_dataset,by.x="claim_number",by.y="claim_number",all=TRUE)



  return(dt)




}






