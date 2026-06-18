test_that("hazard-to-development-factor formula keeps legacy negative-factor repair", {
  hazard <- 0
  eta <- 0.5
  expect_equal((1 + (1 - eta) * hazard) / (1 - eta * hazard), 1)

  hazard <- 0.1
  eta <- 0.5
  expect_equal((1 + (1 - eta) * hazard) / (1 - eta * hazard), 1.05 / 0.95)

  hazard <- 0.1
  eta <- 0.25
  expect_equal((1 + (1 - eta) * hazard) / (1 - eta * hazard), 1.075 / 0.975)

  hazard <- (1 / eta) - 1e-8
  expect_true(is.finite((1 + (1 - eta) * hazard) / (1 - eta * hazard)))
  expect_gt((1 + (1 - eta) * hazard) / (1 - eta * hazard), 1e6)

  hazard_frame <- data.frame(
    AP_i = 1,
    DP_rev_i = 1,
    expg = 1,
    baseline = 1,
    hazard = (1 / 0.5) + 0.1
  )

  expect_equal(
    pkg.env$hazard_data_frame(
      hazard = hazard_frame,
      eta = 0.5,
      categorical_features = NULL,
      continuous_features = NULL,
      calendar_period_extrapolation = FALSE
    )$dev_f_i,
    1
  )
})
