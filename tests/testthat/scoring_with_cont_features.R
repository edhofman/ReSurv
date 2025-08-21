library(ReSurv)
library(data.table)
library(dplyr)
require(parallel)


seed=1
scenario="zeta"

input_data <- data_generator(random_seed = seed,
                             scenario=scenario,
                             time_unit =  1 / 360,
                             years = 4,
                             period_exposure  = 200)

categorical_features = c("business_use")
continuous_features = c("AP","age","property_value")

individual_data <- IndividualDataPP(
  input_data,
  categorical_features = categorical_features,
  continuous_features = continuous_features,
  accident_period = "AP",
  calendar_period = "RP",
  input_time_granularity = "days",
  output_time_granularity = "quarters",
  years = 4
)

individual_data_y <- IndividualDataPP(
  input_data,
  categorical_features = categorical_features,
  continuous_features = continuous_features,
  accident_period = "AP",
  calendar_period = "RP",
  input_time_granularity = "days",
  output_time_granularity = "years",
  years = 4
)

resurv.fit <- ReSurv(individual_data,
                     hazard_model = "COX",
                     grouping_method="probability",
                     simplifier=TRUE)

resurv.fit.predict <- predict(resurv.fit)



## Error one
### here it happens sometimes (?)

sum(resurv.fit.predict$long_triangle_format_out$output_granularity$IBNR)
sum(resurv.fit.predict$long_triangle_format_out$input_granularity$IBNR,na.rm=T)


### Below I wrote the code to fix the are tot / are cal.. We have the problem that
### we do not save the continuous features in resurv.fit.predict$long_triangle_format_out$xxxxxxxxx_granularity.
### I then added an encoder output:

# cfs <- c("age","property_value","group_o")
# encoder <- resurv.fit.predict$groups_encoding[,..cfs]

### One should think a bit how to score.. It can be tricky.

conversion_factor <- individual_data$conversion_factor

max_dp_i <- 1440

true_output <- resurv.fit$IndividualDataPP$full.data %>%
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
  filter(DP_rev_o <= (TR_o)) %>%
  group_by(business_use, age, property_value, AP_o, DP_rev_o) %>%
  mutate(business_use = as.character(business_use)) %>%
  summarize(I = sum(I), .groups = "drop") %>%
  filter(DP_rev_o > 0)

out_list <- resurv.fit.predict$long_triangle_format_out
out <- out_list$output_granularity


cfs <- c("age","property_value","group_o")
encoder <- resurv.fit.predict$groups_encoding[,..cfs]

setDT(out)

out <- merge(out,encoder,
             by=c("group_o"))

cs <- c(categorical_features, setdiff(continuous_features,"AP"), "AP_o", "DP_o", "expected_counts")

out <-  out[, ..cs]
out[,business_use:=as.character(business_use)]

out[,c("TR_o",
       "DP_rev_o"):=list(AP_o - 1,
                         16 - DP_o + 1)]

setDT(true_output)


# why no matching? ----
tmp <- merge(true_output,
             out,
             by = c(categorical_features, setdiff(continuous_features,"AP"), "AP_o", "DP_rev_o"))


#Total output
score_total <- out%>%
  as.data.frame() %>%
  inner_join(true_output, by = c(categorical_features, setdiff(continuous_features,"AP"), "AP_o", "DP_rev_o")) %>%
  mutate(ave = I - expected_counts, abs_ave = abs(ave)) %>%
  # from here it is reformulated for the are tot
  ungroup() %>%
  group_by(AP_o, DP_rev_o) %>%
  reframe(abs_ave = abs(sum(ave)), I = sum(I))

are_tot <- sum(score_total$abs_ave) / sum(score_total$I)
ei_r <- sum(resurv.fit.predict$predicted_counts) / sum(score_total$I)


dfs_output <- out %>%
  mutate(DP_rev_o = 16 - DP_o + 1) %>%
  select(AP_o, claim_type, DP_rev_o, f_o) %>%
  mutate(DP_rev_o = DP_rev_o) %>%
  distinct()

#Cashflow on output scale.Etc quarterly cashflow development
# are cal quarterly ----
score_diagonal <- resurv.fit$IndividualDataPP$full.data  %>%
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



# are cal yearly ----

conversion_factor <- individual_data_y$conversion_factor

max_dp_i <- 1440

true_output <- individual_data_y$full.data %>%
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


resurv.fit.predict.y <- predict(resurv.fit,
                                newdata = individual_data_y,
                                grouping_method = "probability")

out_list <- resurv.fit.predict.y$long_triangle_format_out
out <- out_list$output_granularity


score_diagonal <- resurv.fit$IndividualDataPP$full.data  %>%
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

are_cal_y <- sum(score_diagonal$abs_ave2_diag) / sum(score_diagonal$I)







