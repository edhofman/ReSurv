test_that("summary, plot, and ooslkh use current ReSurvFit structure", {
  fit <- list(
    model.out = list(data = data.frame(x = c(0, 1)), model.out = NULL),
    hazard_frame = data.frame(
      AP_i = c(1, 1, 1),
      DP_i = c(1, 2, 3),
      DP_rev_i = c(3, 2, 1),
      expg = 1,
      baseline = 0.01,
      hazard = 0.01,
      f_i = 1.01,
      cum_f_i = c(1.01, 1.01^2, 1.01^3),
      S_i = 1 / c(1.01, 1.01^2, 1.01^3),
      S_i_lag = c(1, 1 / 1.01, 1 / 1.01^2),
      S_i_lead = c(1 / 1.01^2, 1 / 1.01^3, 0)
    ),
    data_information = list(
      continuous_features = "x",
      categorical_features = NULL,
      data_for_reserving = data.frame(
        AP_i = c(1, 2),
        DP_rev_i = c(1, 1),
        TR_i = c(1, 1),
        I = c(0, 0),
        x = c(0, 1)
      )
    ),
    fit_information = list(hazard_model = "COX")
  )
  class(fit) <- "ReSurvFit"

  pred <- list(
    long_triangle_format_out = list(
      input_granularity = data.frame(
        AP_i = c(1, 1, 1),
        DP_i = c(1, 2, 3),
        group_i = 1,
        f_i = c(1, 1.05, 1.1),
        IBNR = c(NA, 1, 2)
      ),
      output_granularity = data.frame(
        AP_o = c(1, 1, 1),
        DP_o = c(1, 2, 3),
        group_o = 1,
        f_o = c(1, 1.05, 1.1),
        IBNR = c(NA, 1, 2)
      )
    ),
    grouping_method = "exposure",
    ReSurvFit = fit
  )
  class(pred) <- "ReSurvPredict"

  smry <- summary(pred)
  expect_s3_class(smry, "summaryReSurvPredict")
  expect_output(print(smry), "Continuous Features")

  expect_s3_class(plot(pred), "ggplot")
  expect_error(
    plot(fit),
    "Feature importance plots are currently supported only"
  )

  fit_missing_data <- fit
  fit_missing_data$data_information$data_for_reserving <- NULL
  expect_error(
    ooslkh(fit_missing_data),
    "data_information\\$data_for_reserving"
  )
})
