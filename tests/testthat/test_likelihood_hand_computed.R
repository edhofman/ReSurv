test_that("hand-computed Efron likelihood matches package likelihood", {
  toy <- data.frame(
    DP_rev_i = c(2, 2, 3, 4, 4),
    I        = c(1, 1, 1, 1, 1),
    TR_i     = c(0, 0, 1, 2, 0),
    group    = c("A", "B", "A", "B", "A"),
    phi      = c(0.10, -0.20, 0.30, 0.05, -0.10)
  )

  event_times <- sort(unique(toy$DP_rev_i[toy$I == 1]))

  manual_loss <- 0

  for (t_j in event_times) {
    event_set <- which(toy$DP_rev_i == t_j & toy$I == 1)
    risk_set  <- which(toy$DP_rev_i >= t_j & toy$TR_i < t_j)
    d_j       <- length(event_set)

    manual_loss <- manual_loss +
      sum(log(
        sum(exp(toy$phi[risk_set])) -
          (seq(0, d_j - 1) / d_j) * sum(exp(toy$phi[event_set]))
      )) -
      sum(toy$phi[event_set])
  }

  manual_avg_loss <- manual_loss / sum(toy$I == 1)

  skip_if_not_installed("xgboost")
  dtrain <- xgboost::xgb.DMatrix(
    as.matrix(toy["phi"]),
    label = toy$I
  )
  attr(dtrain, "truncation") <- toy$TR_i
  attr(dtrain, "claim_arrival") <- toy$DP_rev_i
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
  attr(dtrain, "tieid") <- as.integer(table(toy$DP_rev_i))
  attr(dtrain, "efron_c") <- unlist(lapply(
    as.integer(table(toy$DP_rev_i)),
    function(d_j) seq(0, d_j - 1) / d_j
  ))

  expect_equal(
    cox_evaluation_metrics(preds = toy$phi, dtrain = dtrain)$value,
    manual_avg_loss,
    tolerance = 1e-10
  )

  cox <- survival::coxph(
    survival::Surv(TR_i, DP_rev_i, I) ~ offset(phi),
    data = toy,
    ties = "efron"
  )
  cox_eval <- pkg.env$evaluate_lkh_cox(
    X_train = toy["phi"],
    Y_train = toy[c("DP_rev_i", "I", "TR_i")],
    model = list(cox = cox)
  )
  expect_equal(cox_eval$value, manual_avg_loss, tolerance = 1e-10)

  skip_if_not_installed("torch")
  torch_loss <- pkg.env$cox_ph_loss_torch(
    log_h = torch::torch_tensor(toy$phi, dtype = torch::torch_float32()),
    durations = torch::torch_tensor(toy$DP_rev_i, dtype = torch::torch_float32()),
    events = torch::torch_tensor(toy$I, dtype = torch::torch_float32()),
    truncation = torch::torch_tensor(toy$TR_i, dtype = torch::torch_float32())
  )
  expect_equal(torch_loss$item(), manual_avg_loss, tolerance = 1e-6)
})
