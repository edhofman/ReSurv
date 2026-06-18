#' Predict deterministic reserve table
#'
#' @export
predictReserve <- function(object, ...) {
  UseMethod("predictReserve")
}


#' Predict deterministic reserve table from a ReSurvFit object
#'
#' @param object A ReSurvFit object.
#' @param granularity Character. Either "output" or "input".
#' @param ... Additional arguments passed to predict.ReSurvFit().
#'
#' @return A data.table with columns AP, DP, CP, IBNR.
#'
#' @export
predictReserve.ReSurvFit <- function(object,
                                     granularity = c("output", "input"),
                                     ...) {

  granularity <- match.arg(granularity)

  dot_args <- list(...)

  dot_args$object <- object
  dot_args$minimal_output <- TRUE
  dot_args$lower_triangular_output <- FALSE

  pred <- do.call(predict, dot_args)

  if (is.null(pred$long_triangle_format_out$input_granularity)) {
    stop(
      "`predict(object, minimal_output = TRUE)` did not return input-granularity predictions.",
      call. = FALSE
    )
  }

  dt <- data.table::copy(
    data.table::as.data.table(
      pred$long_triangle_format_out$input_granularity
    )
  )

  if (!all(c("AP_i", "DP_i", "IBNR") %in% names(dt))) {
    stop(
      "Input-granularity predictions must contain `AP_i`, `DP_i`, and `IBNR`.",
      call. = FALSE
    )
  }

  dt <- dt[
    !is.na(IBNR)
  ]

  if (granularity == "input") {

    out <- dt[
      ,
      .(IBNR = sum(IBNR, na.rm = TRUE)),
      by = .(
        AP = as.integer(AP_i),
        DP = as.integer(DP_i)
      )
    ]

    out[
      ,
      CP := AP + DP - 1L
    ]

  } else {

    data_information <- object$data_information

    conversion_factor <- data_information$conversion_factor
    years <- data_information$years
    input_time_granularity <- data_information$input_time_granularity
    output_time_granularity <- data_information$output_time_granularity

    time_unit_string <- c("days", "months", "quarters", "semesters", "years")
    time_unit_numeric <- c(1 / 360, 1 / 12, 1 / 4, 1 / 2, 1)

    input_pos <- match(input_time_granularity, time_unit_string)
    output_pos <- match(output_time_granularity, time_unit_string)

    if (is.na(input_pos)) {
      stop(
        "`input_time_granularity` must be one of: ",
        paste(time_unit_string, collapse = ", "),
        call. = FALSE
      )
    }

    if (is.na(output_pos)) {
      stop(
        "`output_time_granularity` must be one of: ",
        paste(time_unit_string, collapse = ", "),
        call. = FALSE
      )
    }

    if (is.null(conversion_factor) || length(conversion_factor) != 1L) {
      conversion_factor <- time_unit_numeric[input_pos] /
        time_unit_numeric[output_pos]
    }

    dt[
      ,
      `:=`(
        AP = as.integer(ceiling(AP_i * conversion_factor)),
        CP = as.integer(ceiling((AP_i + DP_i - 1L) * conversion_factor))
      )
    ]

    dt[
      ,
      DP := CP - AP + 1L
    ]

    out <- dt[
      ,
      .(IBNR = sum(IBNR, na.rm = TRUE)),
      by = .(AP, DP, CP)
    ]
  }

  data.table::setorder(out, AP, DP, CP)

  data.table::setcolorder(
    out,
    c("AP", "DP", "CP", "IBNR")
  )

  out[]
}
