library(ReSurv)

input_data <- data_generator(random_seed = 1,
                             scenario = "zeta",
                             time_unit = 1 / 360,
                             years = 4,
                             period_exposure = 200)

individual_data <- IndividualDataPP(data = input_data,
                                    categorical_features = "business_use",
                                    continuous_features = c("age","property_value","AP"),
                                    accident_period = "AP",
                                    calendar_period = "RP",
                                    input_time_granularity = "days",
                                    output_time_granularity = "quarters",
                                    years = 4)


#when simplifier is set to T, I ease the computations both for fitting and predicting.
#  aka, there is an if cycle in ReSurv and another one in ReSurv.predict.
resurv_fit_cox <- ReSurv(individual_data,
                         hazard_model = "COX",
                         grouping_method="probability",
                         simplifier = TRUE)

resurv_fit_cox2 <- ReSurv(individual_data,
                         hazard_model = "COX",
                         grouping_method="probability",
                         simplifier = FALSE)


pr1 <- predict(resurv_fit_cox)
pr2 <- predict(resurv_fit_cox2)


## we get the same ibnr level

print(summary(pr1))
print(summary(pr2))


## below you compute the same aretot for both but different are cal.

##simplifier T


conversion_factor <- individual_data$conversion_factor

max_dp_i <- 1440

true_output <- resurv_fit_cox$IndividualDataPP$full.data %>%
  mutate(
    DP_rev_o = floor(max_dp_i * conversion_factor) -
      ceiling(
        DP_i * conversion_factor +
          ((AP_i - 1) %% (
            1 / conversion_factor
          )) * conversion_factor) + 1,
    AP_o = ceiling(AP_i * conversion_factor),
    TR_o = AP_o - 1
  ) %>%
  filter(DP_rev_o <= TR_o) %>%
  group_by(claim_type, AP_o, DP_rev_o) %>%
  mutate(claim_type = as.character(claim_type)) %>%
  summarize(I = sum(I), .groups = "drop") %>%
  filter(DP_rev_o > 0)

out_list <- pr1$long_triangle_format_out
out <- out_list$output_granularity

#Total output
score_total <- out[, c("claim_type", "AP_o", "DP_o", "expected_counts")] %>%
  mutate(DP_rev_o = 16 - DP_o + 1) %>%
  inner_join(true_output, by = c("claim_type", "AP_o", "DP_rev_o")) %>%
  mutate(ave = I - expected_counts, abs_ave = abs(ave)) %>%
  # from here it is reformulated for the are tot
  ungroup() %>%
  group_by(AP_o, DP_rev_o) %>%
  reframe(abs_ave = abs(sum(ave)), I = sum(I))

are_tot <- sum(score_total$abs_ave) / sum(score_total$I)


dfs_output <- out %>%
  mutate(DP_rev_o = 16 - DP_o + 1) %>%
  select(AP_o, claim_type, DP_rev_o, f_o) %>%
  mutate(DP_rev_o = DP_rev_o) %>%
  distinct()

#Cashflow on output scale.Etc quarterly cashflow development
score_diagonal <- resurv_fit_cox$IndividualDataPP$full.data  %>%
  mutate(
    DP_rev_o = floor(max_dp_i * conversion_factor) -
      ceiling(
        DP_i * conversion_factor +
          ((AP_i - 1) %% (
            1 / conversion_factor
          )) * conversion_factor) + 1,
    AP_o = ceiling(AP_i * conversion_factor)
  ) %>%
  group_by(claim_type, AP_o, DP_rev_o) %>%
  mutate(claim_type = as.character(claim_type)) %>%
  summarize(I = sum(I), .groups = "drop") %>%
  group_by(claim_type, AP_o) %>%
  arrange(desc(DP_rev_o)) %>%
  mutate(I_cum = cumsum(I)) %>%
  mutate(I_cum_lag = lag(I_cum, default = 0)) %>%
  left_join(dfs_output, by = c("AP_o", "claim_type", "DP_rev_o")) %>%
  mutate(I_cum_hat =  I_cum_lag * f_o,
         RP_o = max(DP_rev_o) - DP_rev_o + AP_o) %>%
  inner_join(true_output[, c("AP_o", "DP_rev_o")] %>%  distinct()
             , by = c("AP_o", "DP_rev_o")) %>%
  group_by(AP_o, DP_rev_o) %>%
  reframe(abs_ave2_diag = abs(sum(I_cum_hat) - sum(I_cum)), I = sum(I))

are_cal_q <- sum(score_diagonal$abs_ave2_diag) / sum(score_diagonal$I)


##simplifier F


conversion_factor <- individual_data$conversion_factor

max_dp_i <- 1440

true_output <- resurv_fit_cox2$IndividualDataPP$full.data %>%
  mutate(
    DP_rev_o = floor(max_dp_i * conversion_factor) -
      ceiling(
        DP_i * conversion_factor +
          ((AP_i - 1) %% (
            1 / conversion_factor
          )) * conversion_factor) + 1,
    AP_o = ceiling(AP_i * conversion_factor),
    TR_o = AP_o - 1
  ) %>%
  filter(DP_rev_o <= TR_o) %>%
  group_by(claim_type, AP_o, DP_rev_o) %>%
  mutate(claim_type = as.character(claim_type)) %>%
  summarize(I = sum(I), .groups = "drop") %>%
  filter(DP_rev_o > 0)

out_list <- pr2$long_triangle_format_out
out <- out_list$output_granularity

#Total output
score_total <- out[, c("claim_type", "AP_o", "DP_o", "expected_counts")] %>%
  mutate(DP_rev_o = 16 - DP_o + 1) %>%
  inner_join(true_output, by = c("claim_type", "AP_o", "DP_rev_o")) %>%
  mutate(ave = I - expected_counts, abs_ave = abs(ave)) %>%
  # from here it is reformulated for the are tot
  ungroup() %>%
  group_by(AP_o, DP_rev_o) %>%
  reframe(abs_ave = abs(sum(ave)), I = sum(I))

are_tot <- sum(score_total$abs_ave) / sum(score_total$I)


dfs_output <- out %>%
  mutate(DP_rev_o = 16 - DP_o + 1) %>%
  select(AP_o, claim_type, DP_rev_o, f_o) %>%
  mutate(DP_rev_o = DP_rev_o) %>%
  distinct()

#Cashflow on output scale.Etc quarterly cashflow development
score_diagonal <- resurv_fit_cox2$IndividualDataPP$full.data  %>%
  mutate(
    DP_rev_o = floor(max_dp_i * conversion_factor) -
      ceiling(
        DP_i * conversion_factor +
          ((AP_i - 1) %% (
            1 / conversion_factor
          )) * conversion_factor) + 1,
    AP_o = ceiling(AP_i * conversion_factor)
  ) %>%
  group_by(claim_type, AP_o, DP_rev_o) %>%
  mutate(claim_type = as.character(claim_type)) %>%
  summarize(I = sum(I), .groups = "drop") %>%
  group_by(claim_type, AP_o) %>%
  arrange(desc(DP_rev_o)) %>%
  mutate(I_cum = cumsum(I)) %>%
  mutate(I_cum_lag = lag(I_cum, default = 0)) %>%
  left_join(dfs_output, by = c("AP_o", "claim_type", "DP_rev_o")) %>%
  mutate(I_cum_hat =  I_cum_lag * f_o,
         RP_o = max(DP_rev_o) - DP_rev_o + AP_o) %>%
  inner_join(true_output[, c("AP_o", "DP_rev_o")] %>%  distinct()
             , by = c("AP_o", "DP_rev_o")) %>%
  group_by(AP_o, DP_rev_o) %>%
  reframe(abs_ave2_diag = abs(sum(I_cum_hat) - sum(I_cum)), I = sum(I))

are_cal_q2 <- sum(score_diagonal$abs_ave2_diag) / sum(score_diagonal$I)


