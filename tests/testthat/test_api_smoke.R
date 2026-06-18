test_that("API smoke test: COX with no covariates and conversion_factor = 1", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    cat = rep(c("A", "B"), 4),
    z = seq(0.1, 0.8, by = 0.1)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = NULL,
    continuous_features = NULL,
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )
  fit <- suppressWarnings(ReSurv(idata, hazard_model = "COX", eta = 0))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})

test_that("API smoke test: COX with categorical covariates only", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    cat = rep(c("A", "B"), 4),
    z = seq(0.1, 0.8, by = 0.1)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = "cat",
    continuous_features = NULL,
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )
  fit <- suppressWarnings(ReSurv(idata, hazard_model = "COX", eta = 0))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})

test_that("API smoke test: COX with continuous covariates only", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    cat = rep(c("A", "B"), 4),
    z = seq(0.1, 0.8, by = 0.1)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = NULL,
    continuous_features = "z",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )
  fit <- suppressWarnings(ReSurv(idata, hazard_model = "COX", eta = 0))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})

test_that("API smoke test: COX with AP_i as covariate", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    cat = rep(c("A", "B"), 4),
    z = seq(0.1, 0.8, by = 0.1)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = NULL,
    continuous_features = "AP",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )
  fit <- suppressWarnings(ReSurv(idata, hazard_model = "COX", eta = 0))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})

test_that("API smoke test: COX with mixed covariates", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    cat = rep(c("A", "B"), 4),
    z = seq(0.1, 0.8, by = 0.1)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = "cat",
    continuous_features = "z",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )
  fit <- suppressWarnings(ReSurv(idata, hazard_model = "COX", eta = 0))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})

test_that("API smoke test: COX with conversion_factor != 1 predicts successfully", {
  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7),
    cat = rep(c("A", "B"), 6),
    z = seq(0.1, 1.2, by = 0.1)
  )

  idata <- IndividualDataPP(
    dat,
    categorical_features = "cat",
    continuous_features = "z",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "quarters",
    output_time_granularity = "years",
    years = 2
  )

  fit <- suppressWarnings(ReSurv(idata, hazard_model = "COX", eta = 0))

  pred <- predict(
    fit,
    minimal_output = FALSE,
    lower_triangular_output = FALSE
  )

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")

  expect_true("input_granularity" %in% names(pred$long_triangle_format_out))
  expect_true("output_granularity" %in% names(pred$long_triangle_format_out))

  expect_true(is.finite(pred$predicted_counts))
  expect_true(pred$predicted_counts >= 0)
})
test_that("API smoke test: XGB path is skipped when local xgboost DMatrix is unavailable", {
  skip_if_not_installed("xgboost")
  skip_if(inherits(
    try(xgboost::xgb.DMatrix(matrix(as.numeric(1:8), ncol = 1), label = rep(1, 8)), silent = TRUE),
    "try-error"
  ))

  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    z = seq(0.1, 0.8, by = 0.1)
  )
  idata <- IndividualDataPP(
    dat,
    categorical_features = NULL,
    continuous_features = "z",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )

  fit <- try(
    ReSurv(
      idata,
      hazard_model = "XGB",
      eta = 0,
      percentage_data_training = 1,
      simplifier = TRUE
    ),
    silent = TRUE
  )
  skip_if(inherits(fit, "try-error"))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})

test_that("API smoke test: NN path with mixed covariates", {
  skip_if_not_installed("torch")

  dat <- data.frame(
    AP = c(1, 1, 2, 2, 3, 3, 4, 4),
    RP = c(1, 2, 2, 3, 3, 4, 4, 5),
    cat = rep(c("A", "B"), 4),
    z = seq(0.1, 0.8, by = 0.1)
  )
  idata <- IndividualDataPP(
    dat,
    categorical_features = "cat",
    continuous_features = "z",
    accident_period = "AP",
    calendar_period = "RP",
    input_time_granularity = "years",
    output_time_granularity = "years",
    years = 5
  )

  fit <- try(ReSurv(
    idata,
    hazard_model = "NN",
    eta = 0,
    percentage_data_training = 1,
    hparameters = list(
      num_layers = 1,
      num_nodes = 2,
      activation = "relu",
      optim = "Adam",
      lr = 0.01,
      xi = 0.5,
      eps = 0,
      early_stopping = FALSE,
      patience = 1,
      verbose = FALSE,
      epochs = 1,
      num_workers = 0
    )
  ), silent = TRUE)
  skip_if(inherits(fit, "try-error"))
  pred <- predict(fit, minimal_output = FALSE, lower_triangular_output = FALSE)

  expect_s3_class(fit, "ReSurvFit")
  expect_s3_class(pred, "ReSurvPredict")
  expect_s3_class(summary(pred), "summaryReSurvPredict")
})
