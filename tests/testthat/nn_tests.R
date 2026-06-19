test_that("torch NN Cox-Efron loss matches R reference", {
  skip_if_not_installed("torch")

  log_h <- c(0.10, -0.20, 0.30, 0.05, -0.10)
  durations <- c(2, 2, 3, 4, 4)
  events <- c(1, 1, 1, 1, 1)
  truncation <- c(0, 0, 1, 2, 0)

  ref <- cox_ph_loss_reference_r(
    log_h      = log_h,
    durations  = durations,
    events     = events,
    truncation = truncation
  )

  torch_loss <- cox_ph_loss_torch(
    log_h = torch::torch_tensor(log_h, dtype = torch::torch_float32()),
    durations = torch::torch_tensor(durations, dtype = torch::torch_float32()),
    events = torch::torch_tensor(events, dtype = torch::torch_float32()),
    truncation = torch::torch_tensor(truncation, dtype = torch::torch_float32())
  )$item()

  expect_equal(torch_loss, ref, tolerance = 1e-5)
})


test_that("torch NN Cox-Efron loss matches R reference", {
  skip_if_not_installed("torch")

  log_h <- c(0.10, -0.20, 0.30, 0.05, -0.10)
  durations <- c(2, 2, 3, 4, 4)
  events <- c(1, 1, 1, 1, 1)
  truncation <- c(0, 0, 1, 2, 0)

  ref <- cox_ph_loss_reference_r(
    log_h      = log_h,
    durations  = durations,
    events     = events,
    truncation = truncation
  )

  torch_loss <- cox_ph_loss_torch(
    log_h = torch::torch_tensor(log_h, dtype = torch::torch_float32()),
    durations = torch::torch_tensor(durations, dtype = torch::torch_float32()),
    events = torch::torch_tensor(events, dtype = torch::torch_float32()),
    truncation = torch::torch_tensor(truncation, dtype = torch::torch_float32())
  )$item()

  expect_equal(torch_loss, ref, tolerance = 1e-5)
})
