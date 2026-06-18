test_that("simplified forecast frame derives input time features", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = NULL,
    continuous_features = "DP_i",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )

  newdata <- simplified_df_2_fcst(idata, hazard_model = "XGB")
  max_dp_i <- pkg.env$maximum.time(
    idata$data_information$years,
    idata$data_information$input_time_granularity
  )

  expect_true(all(c("AP_i", "DP_rev_i", "DP_i") %in% names(newdata)))
  expect_equal(newdata$DP_i, max_dp_i - newdata$DP_rev_i + 1L)
  observed_dp_rev_i <- seq(
    min(idata$training.data$DP_rev_i),
    max(idata$training.data$DP_rev_i)
  )
  expect_equal(nrow(newdata), length(unique(idata$training.data$AP_i)) * length(observed_dp_rev_i))
})
