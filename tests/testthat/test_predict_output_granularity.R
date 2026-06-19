test_that("output development factors are coherent with input development factors", {
  years <- 1
  input_time_granularity <- "months"
  output_time_granularity <- "quarters"
  conversion_factor <- 1 / 3

  q <- 1.05

  max_dp_i <- maximum.time(
    years,
    input_time_granularity
  )

  max_dp_o <- maximum.time(
    years,
    output_time_granularity
  )

  hazard_frame <- data.table::CJ(
    AP_i = seq_len(max_dp_i),
    DP_i = seq_len(max_dp_i),
    z = 1
  )

  hazard_frame[
    ,
    DP_rev_i := max_dp_i - DP_i + 1
  ]

  hazard_frame <- hazard_frame[
    DP_rev_i > AP_i - 1
  ]

  hazard_frame[
    ,
    `:=`(
      f_i       = q,
      cum_f_i   = q^DP_i,
      S_i       = 1 / q^DP_i,
      S_i_lag   = 1 / q^pmax(DP_i - 1, 0),
      expg      = 1,
      baseline  = 0.01,
      hazard    = 0.01
    )
  ]

  data_reserve <- data.table::copy(hazard_frame)[
    ,
    .(AP_i, DP_i, DP_rev_i, z)
  ]

  data_reserve[
    ,
    `:=`(
      AP_o = ceiling(AP_i * conversion_factor),
      I    = as.numeric(DP_i == 1)
    )
  ]

  development_periods <- unique(data_reserve[, .(AP_i, AP_o)])

  dp_ranges <- development_periods[
    ,
    .(DP_rev_o = seq_len(max_dp_o)),
    by = .(AP_i, AP_o)
  ][
    ,
    `:=`(
      min_dp = AP_i + (DP_rev_o - AP_o) / conversion_factor,
      max_dp = AP_i - 1 + (DP_rev_o - AP_o + 1) / conversion_factor
    )
  ]

  object <- list(
    hazard_frame = hazard_frame,
    data_information = list(
      conversion_factor = conversion_factor,
      string_formula_i = NULL,
      string_formula_o = NULL,
      continuous_features = c("AP_i", "z"),
      categorical_features = NULL,
      calendar_period_extrapolation = FALSE,
      years = years,
      accident_period = NULL,
      calendar_period = NULL,
      input_time_granularity = input_time_granularity,
      output_time_granularity = output_time_granularity,
      data_for_reserving = data_reserve,
      dp_ranges = dp_ranges
    )
  )

  class(object) <- "ReSurvFit"

  pred <- predict(
    object,
    minimal_output = FALSE,
    lower_triangular_output = FALSE
  )

  out_i <- data.table::as.data.table(
    pred$long_triangle_format_out$input_granularity
  )

  out_o <- data.table::as.data.table(
    pred$long_triangle_format_out$output_granularity
  )

  testthat::expect_true("f_i" %in% names(out_i))
  testthat::expect_true("f_o" %in% names(out_o))

  testthat::expect_false(anyNA(out_o$f_o))
  testthat::expect_true(all(is.finite(out_o$f_o)))
  testthat::expect_true(all(out_o$f_o > 0))

  non_trivial_factors <- unique(out_o[
    f_o != 1,
    .(group_o, DP_o, f_o)
  ])

  testthat::expect_gt(nrow(non_trivial_factors), 0)

  testthat::expect_gt(nrow(non_trivial_factors), 0)

  testthat::expect_true(
    all(non_trivial_factors$f_o > 1)
  )

  testthat::expect_true(
    all(non_trivial_factors$f_o <= q^3 + 1e-8)
  )

  testthat::expect_equal(
    sum(out_o$expected_counts, na.rm = TRUE),
    sum(out_i$expected_counts, na.rm = TRUE),
    tolerance = 1e-8
  )

  testthat::expect_equal(
    sum(out_o$IBNR, na.rm = TRUE),
    sum(out_i$IBNR, na.rm = TRUE),
    tolerance = 1e-8
  )
})


test_that("output development factors increase when input development factors increase", {
  make_prediction <- function(q) {
    years <- 1
    input_time_granularity <- "months"
    output_time_granularity <- "quarters"
    conversion_factor <- 1 / 3

    max_dp_i <- maximum.time(
      years,
      input_time_granularity
    )

    max_dp_o <- maximum.time(
      years,
      output_time_granularity
    )

    hazard_frame <- data.table::CJ(
      AP_i = seq_len(max_dp_i),
      DP_i = seq_len(max_dp_i),
      z = 1
    )

    hazard_frame[
      ,
      DP_rev_i := max_dp_i - DP_i + 1
    ]

    hazard_frame <- hazard_frame[
      DP_rev_i > AP_i - 1
    ]

    hazard_frame[
      ,
      `:=`(
        f_i       = q,
        cum_f_i   = q^DP_i,
        S_i       = 1 / q^DP_i,
        S_i_lag   = 1 / q^pmax(DP_i - 1, 0),
        expg      = 1,
        baseline  = 0.01,
        hazard    = 0.01
      )
    ]

    data_reserve <- data.table::copy(hazard_frame)[
      ,
      .(AP_i, DP_i, DP_rev_i, z)
    ]

    data_reserve[
      ,
      `:=`(
        AP_o = ceiling(AP_i * conversion_factor),
        I    = as.numeric(DP_i == 1)
      )
    ]

    development_periods <- unique(data_reserve[, .(AP_i, AP_o)])

    dp_ranges <- development_periods[
      ,
      .(DP_rev_o = seq_len(max_dp_o)),
      by = .(AP_i, AP_o)
    ][
      ,
      `:=`(
        min_dp = AP_i + (DP_rev_o - AP_o) / conversion_factor,
        max_dp = AP_i - 1 + (DP_rev_o - AP_o + 1) / conversion_factor
      )
    ]

    object <- list(
      hazard_frame = hazard_frame,
      data_information = list(
        conversion_factor = conversion_factor,
        string_formula_i = NULL,
        string_formula_o = NULL,
        continuous_features = c("AP_i", "z"),
        categorical_features = NULL,
        calendar_period_extrapolation = FALSE,
        years = years,
        accident_period = NULL,
        calendar_period = NULL,
        input_time_granularity = input_time_granularity,
        output_time_granularity = output_time_granularity,
        data_for_reserving = data_reserve,
        dp_ranges = dp_ranges
      )
    )

    class(object) <- "ReSurvFit"

    pred <- predict(
      object,
      minimal_output = FALSE,
      lower_triangular_output = FALSE
    )

    data.table::as.data.table(
      pred$long_triangle_format_out$output_granularity
    )[
      f_o != 1,
      .(AP_o, group_o, DP_o, f_o)
    ]
  }

  low_q <- make_prediction(1.02)
  high_q <- make_prediction(1.08)

  cmp <- merge(
    low_q,
    high_q,
    by = c("AP_o", "group_o", "DP_o"),
    suffixes = c("_low", "_high")
  )

  testthat::expect_gt(nrow(cmp), 0)

  testthat::expect_true(
    all(cmp$f_o_high >= cmp$f_o_low - 1e-8)
  )

  testthat::expect_true(
    any(cmp$f_o_high > cmp$f_o_low + 1e-8)
  )
})
