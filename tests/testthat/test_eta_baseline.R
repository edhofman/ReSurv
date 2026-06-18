test_that("baseline.efron uses eta consistently", {
  skip_if_not_installed("xgboost")

  X <- data.frame(x = c(0, 1, 0, 1))
  Y <- data.frame(
    DP_rev_i = c(2, 2, 3, 4),
    I        = c(1, 1, 1, 1),
    TR_i     = c(0, 0, 1, 2)
  )

  dtrain <- pkg.env$xgboost_pp(
    X = X,
    Y = Y,
    training_test_split = 1
  )$ds_train_m

  preds <- c(0.1, -0.2, 0.3, 0.0)

  b_eta_05 <- pkg.env$baseline.efron(
    preds  = preds,
    dtrain = dtrain,
    eta    = 0.5
  )

  b_eta_00 <- pkg.env$baseline.efron(
    preds  = preds,
    dtrain = dtrain,
    eta    = 0
  )

  risk_sets  <- attr(dtrain, "risk_sets")
  event_sets <- attr(dtrain, "event_sets")

  manual <- function(eta) {
    vapply(seq_along(event_sets), function(j) {
      risk_sum  <- sum(exp(preds[risk_sets[[j]]]))
      event_sum <- sum(exp(preds[event_sets[[j]]]))
      d_j       <- length(event_sets[[j]])

      d_j / (risk_sum - eta * event_sum)
    }, numeric(1))
  }

  expect_equal(b_eta_05, manual(0.5), tolerance = 1e-10)
  expect_equal(b_eta_00, manual(0), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(b_eta_05, b_eta_00)))
})
