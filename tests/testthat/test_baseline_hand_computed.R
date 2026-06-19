test_that("hand-computed Efron baseline matches package baseline", {
  skip_if_not_installed("xgboost")

  toy <- data.frame(
    DP_rev_i = c(2, 2, 3, 4, 4),
    I        = c(1, 1, 1, 1, 1),
    TR_i     = c(0, 0, 1, 2, 0),
    group    = c("A", "B", "A", "B", "A"),
    phi      = c(0.10, -0.20, 0.30, 0.05, -0.10)
  )

  event_times <- sort(unique(toy$DP_rev_i[toy$I == 1]))
  dtrain <- xgboost::xgb.DMatrix(
    as.matrix(toy["phi"]),
    label = toy$I
  )
  attr(dtrain, "risk_sets") <- risks_in_the_tie(
    starts_i = toy$TR_i,
    stops_i = toy$DP_rev_i,
    stops = event_times
  )
  attr(dtrain, "event_sets") <- events_in_the_tie(
    starts_i = toy$TR_i,
    stops_i = toy$DP_rev_i,
    stops = event_times
  )

  risk_sets_manual <- vector("list", length(event_times))
  event_sets_manual <- vector("list", length(event_times))
  manual_baseline_eta_05 <- numeric(length(event_times))
  manual_baseline_eta_025 <- numeric(length(event_times))

  for (j in seq_along(event_times)) {
    t_j <- event_times[j]
    event_set <- which(toy$DP_rev_i == t_j & toy$I == 1)
    risk_set  <- which(toy$DP_rev_i >= t_j & toy$TR_i < t_j)
    d_j       <- length(event_set)
    risk_sum  <- sum(exp(toy$phi[risk_set]))
    event_sum <- sum(exp(toy$phi[event_set]))

    risk_sets_manual[[j]] <- risk_set
    event_sets_manual[[j]] <- event_set
    manual_baseline_eta_05[j] <- d_j / (risk_sum - 0.5 * event_sum)
    manual_baseline_eta_025[j] <- d_j / (risk_sum - 0.25 * event_sum)
  }

  expect_equal(attr(dtrain, "risk_sets"), risk_sets_manual)
  expect_equal(attr(dtrain, "event_sets"), event_sets_manual)

  package_baseline_eta_05 <- baseline.efron(
    preds = toy$phi,
    dtrain = dtrain,
    eta = 0.5
  )
  package_baseline_eta_025 <- baseline.efron(
    preds = toy$phi,
    dtrain = dtrain,
    eta = 0.25
  )

  expect_equal(package_baseline_eta_05, manual_baseline_eta_05, tolerance = 1e-10)
  expect_equal(package_baseline_eta_025, manual_baseline_eta_025, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(package_baseline_eta_05, package_baseline_eta_025)))
})
