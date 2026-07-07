make_clmplus_score_fixture <- function() {
  data_env <- new.env(parent = emptyenv())
  utils::data("sifa.mtpl", package = "clmplus", envir = data_env)
  triangle <- as.matrix(data_env$sifa.mtpl)
  n <- nrow(triangle)

  incremental <- triangle
  incremental[, -1L] <-
    triangle[, -1L, drop = FALSE] - triangle[, -ncol(triangle), drop = FALSE]
  observed <- which(!is.na(incremental), arr.ind = TRUE)

  fit <- structure(
    list(data_information = list(
      conversion_factor = 1,
      years = n,
      input_time_granularity = "years",
      output_time_granularity = "years",
      data_for_reserving = data.frame(
        AP_i = observed[, 1L],
        DP_i = observed[, 2L],
        I = incremental[observed]
      )
    )),
    class = c("clmplus_mock_resurv", "ReSurvFit")
  )

  factors <- vapply(
    seq_len(n - 1L),
    function(jj) {
      rows <- seq_len(n - jj)
      sum(triangle[rows, jj + 1L]) / sum(triangle[rows, jj])
    },
    numeric(1L)
  )
  projected <- triangle
  for (aa in 2:n) {
    latest_dp <- n - aa + 1L
    for (dd in (latest_dp + 1L):n) {
      projected[aa, dd] <- projected[aa, dd - 1L] * factors[dd - 1L]
    }
  }
  projected_incremental <- projected
  projected_incremental[, -1L] <-
    projected[, -1L, drop = FALSE] - projected[, -ncol(projected), drop = FALSE]
  lower <- which(row(projected) + col(projected) > n + 1L, arr.ind = TRUE)

  actual <- data.frame(
    AP = lower[, 1L],
    DP = lower[, 2L],
    CP = lower[, 1L] + lower[, 2L] - 1L,
    actual = projected_incremental[lower]
  )

  upper_incremental <- data.table::data.table(
    AP = observed[, 1L],
    DP = observed[, 2L],
    CP = observed[, 1L] + observed[, 2L] - 1L,
    I = incremental[observed]
  )
  data.table::setorder(upper_incremental, AP, DP)
  upper_incremental[, C := cumsum(I), by = AP]

  list(
    fit = fit,
    actual = actual,
    upper_incremental = upper_incremental,
    max_dp = n
  )
}

predictReserve.clmplus_mock_resurv <- function(object, ...) {
  data.table::data.table(
    AP = integer(),
    DP = integer(),
    CP = integer(),
    IBNR = numeric()
  )
}

register_clmplus_mock_method <- function() {
  registerS3method(
    "predictReserve",
    "clmplus_mock_resurv",
    predictReserve.clmplus_mock_resurv,
    envir = asNamespace("ReSurv")
  )
}

test_that("clmplus a benchmark replicates the chain-ladder benchmark", {
  skip_if_not_installed("clmplus")
  register_clmplus_mock_method()
  fixture <- make_clmplus_score_fixture()

  score <- Score_Reserving(
    fixture$fit,
    fixture$actual,
    scoring_metrics = c("EI", "R-cell-wise")
  )
  clmplus_a <- suppressWarnings(ReSurv:::.clmplus_benchmark_predictions(
    fixture$upper_incremental,
    fixture$max_dp,
    "a"
  ))[[1L]]

  ei <- score[["EI"]]
  cell_score <- score[["R-cell-wise"]]
  expect_equal(ei[model == "CL", score], 1, tolerance = 1e-10)
  expect_equal(cell_score[model == "CL", score], 0, tolerance = 1e-10)
  comparison <- merge(
    clmplus_a,
    fixture$actual,
    by = c("AP", "DP", "CP")
  )
  expect_equal(comparison$IBNR, comparison$actual, tolerance = 1e-10)
})

test_that("clmplus ac and apc benchmarks are added when requested", {
  skip_if_not_installed("clmplus")
  register_clmplus_mock_method()
  fixture <- make_clmplus_score_fixture()

  score <- suppressWarnings(Score_Reserving(
    fixture$fit,
    fixture$actual,
    scoring_metrics = "EI",
    chain_ladder = FALSE,
    clmplus_benchmark = c("ac", "apc")
  ))

  ei <- score[["EI"]]
  expect_setequal(ei$model, c("CLMplus-ac", "CLMplus-apc", "Model"))
  expect_true(all(is.finite(ei[model != "Model", score])))
  expect_false("CL" %in% ei$model)
})

test_that("clmplus benchmark input is validated", {
  expect_error(
    Score_Reserving(list(), data.frame(), clmplus_benchmark = "unknown"),
    "Unsupported `clmplus_benchmark`"
  )
})
