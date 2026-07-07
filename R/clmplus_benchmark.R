.clmplus_benchmark_predictions <- function(upper_incremental,
                                           max_dp,
                                           benchmark_models) {
  cumulative_triangle <- matrix(
    NA_real_,
    nrow = max_dp,
    ncol = max_dp
  )
  cumulative_triangle[
    cbind(upper_incremental$AP, upper_incremental$DP)
  ] <- upper_incremental$C

  aggregate_data <- clmplus::AggregateDataPP(
    cumulative.payments.triangle = cumulative_triangle,
    eta = 1 / 2
  )

  predictions <- vector("list", length(benchmark_models))
  names(predictions) <- paste0("CLMplus-", benchmark_models)

  for (ii in seq_along(benchmark_models)) {
    model_name <- benchmark_models[[ii]]

    result <- tryCatch(
      {
        invisible(utils::capture.output(
          fit <- clmplus::clmplus(
            aggregate_data,
            hazard.model = model_name
          )
        ))
        invisible(utils::capture.output(
          prediction <- stats::predict(fit)
        ))
        prediction
      },
      error = function(e) {
        stop(
          "The clmplus '", model_name, "' benchmark failed: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    full_triangle <- as.matrix(result$full_triangle)

    if (!identical(dim(full_triangle), c(max_dp, max_dp)) ||
        any(!is.finite(full_triangle[!is.na(full_triangle)]))) {
      stop(
        "The clmplus '", model_name,
        "' benchmark returned an invalid full triangle.",
        call. = FALSE
      )
    }

    incremental_triangle <- full_triangle
    if (max_dp > 1L) {
      incremental_triangle[, 2:max_dp] <-
        full_triangle[, 2:max_dp, drop = FALSE] -
        full_triangle[, 1:(max_dp - 1L), drop = FALSE]
    }

    lower_cells <- which(
      row(incremental_triangle) + col(incremental_triangle) > max_dp + 1L,
      arr.ind = TRUE
    )

    ap <- as.integer(lower_cells[, 1L])
    dp <- as.integer(lower_cells[, 2L])
    pred <- data.table::data.table(
      AP = ap,
      DP = dp,
      CP = ap + dp - 1L,
      IBNR = as.numeric(incremental_triangle[lower_cells])
    )

    if (any(!is.finite(pred$IBNR))) {
      stop(
        "The clmplus '", model_name,
        "' benchmark returned non-finite reserve predictions.",
        call. = FALSE
      )
    }

    predictions[[ii]] <- pred
  }

  predictions
}
