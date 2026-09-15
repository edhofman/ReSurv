# CI must fail on backend errors instead of silently skipping model tests.
library(ReSurv)
models <- c("COX", "XGB")
if (!"--skip-torch" %in% commandArgs(trailingOnly = TRUE)) {
  stopifnot(torch::torch_is_installed())
  torch::torch_set_num_threads(2)
  models <- c(models, "NN")
}
claims <- data_generator(
  random_seed = 1964, scenario = "alpha", time_unit = 1,
  years = 4, period_exposure = 100
)
individual <- IndividualDataPP(
  claims, categorical_features = "claim_type",
  accident_period = "AP", calendar_period = "RP",
  input_time_granularity = "years", output_time_granularity = "years",
  years = 4
)
for (model in models) {
  hp <- switch(model,
    COX = list(),
    XGB = list(params = list(max_depth = 1, eta = 0.1, nthread = 2),
               nrounds = 2, verbose = 0),
    NN = list(num_layers = 1, num_nodes = 2, activation = "relu",
              optim = "Adam", lr = 0.01, xi = 0.5, eps = 0,
              early_stopping = FALSE, patience = 1, epochs = 2,
              verbose = FALSE, num_workers = 0)
  )
  fit <- ReSurv(individual, hazard_model = model, eta = 0, hparameters = hp)
  reserve <- predictReserve(fit)
  stopifnot(inherits(fit, "ReSurvFit"), nrow(reserve) > 0,
            all(is.finite(reserve$IBNR)), all(reserve$IBNR >= 0))
  message(model, " backend passed")
}
