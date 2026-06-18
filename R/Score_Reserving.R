#' Score reserving predictions
#'
#' @export
Score_Reserving <- function(models,
                            newdata,
                            scoring_metrics = c(
                              "EI",
                              "R-tot",
                              "R-cell-wise",
                              "R-cal-wise"
                            ),
                            granularity = c("output", "input"),
                            chain_ladder = TRUE,
                            ...) {

  granularity <- match.arg(granularity)

  allowed_metrics <- c(
    "EI",
    "R-tot",
    "R-cell-wise",
    "R-cal-wise",
    "CRPS"
  )

  if (any(!scoring_metrics %in% allowed_metrics)) {
    stop(
      "`scoring_metrics` must be a subset of: ",
      paste(allowed_metrics, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.list(models) || inherits(models, "ReSurvFit")) {
    models <- list(Model = models)
  }

  if (is.null(names(models)) || any(names(models) == "")) {
    names(models) <- paste0("Model", seq_along(models))
  }

  if ("CRPS" %in% scoring_metrics) {

    is_resurv_for_crps <- vapply(
      models,
      function(z) inherits(z, "ReSurvFit"),
      logical(1L)
    )

    if (!any(is_resurv_for_crps)) {
      stop(
        "Score_Reserving computes CRPS only for objects of class 'ReSurvFit'.",
        call. = FALSE
      )
    }

    if (any(!is_resurv_for_crps)) {
      warning(
        "Score_Reserving computes CRPS only for ReSurvFit models. ",
        "Non-ReSurvFit models are omitted from the CRPS table.",
        call. = FALSE
      )
    }
  }

  first_resurv <- NULL

  for (ii in seq_along(models)) {
    if (inherits(models[[ii]], "ReSurvFit")) {
      first_resurv <- models[[ii]]
      break
    }
  }

  if (is.null(first_resurv)) {
    stop(
      "At least one object in `models` must inherit from class 'ReSurvFit'. ",
      "This is needed to recover the reserving time scale and, if requested, the chain-ladder benchmark.",
      call. = FALSE
    )
  }

  data_information <- first_resurv$data_information

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

  max_dp_i <- as.integer(years / time_unit_numeric[input_pos])
  max_dp_o <- as.integer(years / time_unit_numeric[output_pos])

  max_dp <- if (granularity == "output") {
    max_dp_o
  } else {
    max_dp_i
  }

  ## ------------------------------------------------------------------
  ## Realized lower triangle: standardize newdata to AP, DP, CP, actual
  ## ------------------------------------------------------------------

  actual_dt <- data.table::copy(data.table::as.data.table(newdata))

  if (all(c("AP", "DP", "CP") %in% names(actual_dt))) {

    actual_value_col <- NULL

    if ("actual" %in% names(actual_dt)) {
      actual_value_col <- "actual"
    }

    if (is.null(actual_value_col) && "I" %in% names(actual_dt)) {
      actual_value_col <- "I"
    }

    if (is.null(actual_value_col) && "IBNR" %in% names(actual_dt)) {
      actual_value_col <- "IBNR"
    }

    if (is.null(actual_value_col)) {
      actual_dt[, actual := 1]
    } else {
      actual_dt[, actual := get(actual_value_col)]
    }

    actual_dt <- actual_dt[
      ,
      .(actual = sum(actual, na.rm = TRUE)),
      by = .(
        AP = as.integer(AP),
        DP = as.integer(DP),
        CP = as.integer(CP)
      )
    ]

  } else {

    actual_value_col <- NULL

    if ("I" %in% names(actual_dt)) {
      actual_value_col <- "I"
    }

    if (is.null(actual_value_col) && "actual" %in% names(actual_dt)) {
      actual_value_col <- "actual"
    }

    if (is.null(actual_value_col) && "IBNR" %in% names(actual_dt)) {
      actual_value_col <- "IBNR"
    }

    if (is.null(actual_value_col)) {
      actual_dt[, actual := 1]
    } else {
      actual_dt[, actual := get(actual_value_col)]
    }

    ap_col <- NULL

    if ("AP_i" %in% names(actual_dt)) {
      ap_col <- "AP_i"
    }

    if (is.null(ap_col) && "AP" %in% names(actual_dt)) {
      ap_col <- "AP"
    }

    if (is.null(ap_col)) {
      stop(
        "`newdata` must contain either `AP`, `AP_i`, or standardized columns `AP`, `DP`, `CP`.",
        call. = FALSE
      )
    }

    rp_col <- NULL

    if ("RP_i" %in% names(actual_dt)) {
      rp_col <- "RP_i"
    }

    if (is.null(rp_col) && "RP" %in% names(actual_dt)) {
      rp_col <- "RP"
    }

    dp_col <- NULL

    if ("DP_i" %in% names(actual_dt)) {
      dp_col <- "DP_i"
    }

    if (is.null(dp_col) && "DP" %in% names(actual_dt)) {
      dp_col <- "DP"
    }

    if (is.null(rp_col) && is.null(dp_col)) {
      stop(
        "`newdata` must contain either `RP`/`RP_i` or `DP`/`DP_i`.",
        call. = FALSE
      )
    }

    actual_dt[, AP_i_tmp := as.integer(get(ap_col))]

    if (!is.null(rp_col)) {
      actual_dt[, RP_i_tmp := as.integer(get(rp_col))]
      actual_dt[, DP_i_tmp := RP_i_tmp - AP_i_tmp + 1L]
    } else {
      actual_dt[, DP_i_tmp := as.integer(get(dp_col))]
      actual_dt[, RP_i_tmp := AP_i_tmp + DP_i_tmp - 1L]
    }

    if (granularity == "input") {

      actual_dt[
        ,
        `:=`(
          AP = AP_i_tmp,
          DP = DP_i_tmp,
          CP = RP_i_tmp
        )
      ]

    } else {

      actual_dt[
        ,
        `:=`(
          AP = as.integer(ceiling(AP_i_tmp * conversion_factor)),
          CP = as.integer(ceiling(RP_i_tmp * conversion_factor))
        )
      ]

      actual_dt[
        ,
        DP := CP - AP + 1L
      ]
    }

    actual_dt <- actual_dt[
      ,
      .(actual = sum(actual, na.rm = TRUE)),
      by = .(AP, DP, CP)
    ]
  }

  ## Keep the same reserving horizon as the old scoring code.
  ## The old syntax used DP_rev > 0, which is equivalent to DP <= max_dp.
  ## It also scores the lower triangle, CP > max_dp.

  actual_dt <- actual_dt[
    !is.na(AP) &
      !is.na(DP) &
      !is.na(CP) &
      AP >= 1L &
      AP <= max_dp &
      DP >= 1L &
      DP <= max_dp &
      CP > max_dp
  ]

  ## ------------------------------------------------------------------
  ## Model predictions through predictReserve()
  ## ------------------------------------------------------------------

  prediction_list <- list()

  for (mm in names(models)) {

    object <- models[[mm]]

    if (inherits(object, "data.table") || inherits(object, "data.frame")) {

      pred_dt <- data.table::copy(data.table::as.data.table(object))

      if (!all(c("AP", "DP", "CP", "IBNR") %in% names(pred_dt))) {
        stop(
          "Prediction tables supplied in `models` must contain columns `AP`, `DP`, `CP`, and `IBNR`.",
          call. = FALSE
        )
      }

      pred_dt <- pred_dt[
        ,
        .(IBNR = sum(IBNR, na.rm = TRUE)),
        by = .(
          AP = as.integer(AP),
          DP = as.integer(DP),
          CP = as.integer(CP)
        )
      ]

    } else {

      pred_dt <- predictReserve(
        object,
        granularity = granularity,
        ...
      )

      pred_dt <- data.table::as.data.table(pred_dt)

      if (!all(c("AP", "DP", "CP", "IBNR") %in% names(pred_dt))) {
        stop(
          "`predictReserve()` must return columns `AP`, `DP`, `CP`, and `IBNR`.",
          call. = FALSE
        )
      }

      pred_dt <- pred_dt[
        ,
        .(IBNR = sum(IBNR, na.rm = TRUE)),
        by = .(
          AP = as.integer(AP),
          DP = as.integer(DP),
          CP = as.integer(CP)
        )
      ]
    }

    prediction_list[[mm]] <- pred_dt
  }

  ## ------------------------------------------------------------------
  ## Chain-ladder benchmark
  ## ------------------------------------------------------------------

  ## ------------------------------------------------------------------
  ## Chain-ladder benchmark
  ## ------------------------------------------------------------------

  if (isTRUE(chain_ladder)) {

    cl_raw <- data.table::copy(
      data.table::as.data.table(
        first_resurv$data_information$data_for_reserving
      )
    )

    if (!("I" %in% names(cl_raw))) {
      cl_raw[, I := 1]
    }

    max_dp <- if (granularity == "output") {
      max_dp_o
    } else {
      max_dp_i
    }

    if (granularity == "output") {

      if (all(c("AP_o", "DP_o") %in% names(cl_raw))) {

        cl_dt <- cl_raw[
          ,
          .(I = sum(I, na.rm = TRUE)),
          by = .(
            AP = as.integer(AP_o),
            DP = as.integer(DP_o)
          )
        ]

      } else if (all(c("AP_o", "DP_rev_o") %in% names(cl_raw))) {

        cl_dt <- cl_raw[
          ,
          .(I = sum(I, na.rm = TRUE)),
          by = .(
            AP = as.integer(AP_o),
            DP_rev = as.integer(DP_rev_o)
          )
        ]

        cl_dt[
          ,
          DP := max_dp - DP_rev + 1L
        ]

        cl_dt[
          ,
          DP_rev := NULL
        ]

      } else if (all(c("AP_i", "DP_i") %in% names(cl_raw))) {

        cl_dt <- data.table::copy(cl_raw)

        cl_dt[
          ,
          `:=`(
            AP = as.integer(ceiling(AP_i * conversion_factor)),
            CP = as.integer(ceiling((AP_i + DP_i - 1L) * conversion_factor))
          )
        ]

        cl_dt[
          ,
          DP := CP - AP + 1L
        ]

        cl_dt <- cl_dt[
          ,
          .(I = sum(I, na.rm = TRUE)),
          by = .(AP, DP)
        ]

      } else {

        stop(
          "Could not construct the output-granularity chain-ladder benchmark. ",
          "Need either `AP_o`/`DP_o`, `AP_o`/`DP_rev_o`, or `AP_i`/`DP_i`.",
          call. = FALSE
        )
      }

    } else {

      if (all(c("AP_i", "DP_i") %in% names(cl_raw))) {

        cl_dt <- cl_raw[
          ,
          .(I = sum(I, na.rm = TRUE)),
          by = .(
            AP = as.integer(AP_i),
            DP = as.integer(DP_i)
          )
        ]

      } else if (all(c("AP_i", "DP_rev_i") %in% names(cl_raw))) {

        cl_dt <- cl_raw[
          ,
          .(I = sum(I, na.rm = TRUE)),
          by = .(
            AP = as.integer(AP_i),
            DP_rev = as.integer(DP_rev_i)
          )
        ]

        cl_dt[
          ,
          DP := max_dp - DP_rev + 1L
        ]

        cl_dt[
          ,
          DP_rev := NULL
        ]

      } else {

        stop(
          "Could not construct the input-granularity chain-ladder benchmark. ",
          "Need either `AP_i`/`DP_i` or `AP_i`/`DP_rev_i`.",
          call. = FALSE
        )
      }
    }

    cl_dt[
      ,
      CP := AP + DP - 1L
    ]

    cl_dt <- cl_dt[
      AP >= 1L &
        AP <= max_dp &
        DP >= 1L &
        DP <= max_dp &
        CP <= max_dp
    ]

    upper_grid <- data.table::CJ(
      AP = seq_len(max_dp),
      DP = seq_len(max_dp),
      sorted = FALSE
    )

    upper_grid[
      ,
      CP := AP + DP - 1L
    ]

    upper_grid <- upper_grid[
      CP <= max_dp
    ]

    upper_incremental <- cl_dt[
      ,
      .(I = sum(I, na.rm = TRUE)),
      by = .(AP, DP, CP)
    ][
      upper_grid,
      on = .(AP, DP, CP)
    ]

    upper_incremental[
      is.na(I),
      I := 0
    ]

    data.table::setorder(upper_incremental, AP, DP)

    upper_incremental[
      ,
      C := cumsum(I),
      by = AP
    ]

    dev_factors <- data.table::data.table(
      DP = seq_len(max_dp - 1L),
      f = rep(1, max_dp - 1L)
    )

    if (max_dp > 1L) {

      for (jj in seq_len(max_dp - 1L)) {

        den <- upper_incremental[
          DP == jj & AP + jj <= max_dp,
          sum(C, na.rm = TRUE)
        ]

        num <- upper_incremental[
          DP == jj + 1L & AP + jj <= max_dp,
          sum(C, na.rm = TRUE)
        ]

        if (is.finite(den) && den > 0) {
          dev_factors[
            DP == jj,
            f := num / den
          ]
        } else {
          dev_factors[
            DP == jj,
            f := 1
          ]
        }
      }
    }

    latest_observed <- upper_incremental[
      ,
      .SD[which.max(DP)],
      by = AP
    ][
      ,
      .(
        AP = AP,
        DP_max = DP,
        C_latest = C
      )
    ]

    cl_rows <- vector("list", nrow(latest_observed))

    for (ii in seq_len(nrow(latest_observed))) {

      aa <- latest_observed$AP[ii]
      latest_dp <- latest_observed$DP_max[ii]
      latest_c <- latest_observed$C_latest[ii]

      future_dp <- seq.int(latest_dp + 1L, max_dp)

      if (length(future_dp) == 0L) {
        next
      }

      tmp <- data.table::data.table(
        AP = aa,
        DP = future_dp
      )

      tmp[
        ,
        CP := AP + DP - 1L
      ]

      tmp[
        ,
        IBNR := 0
      ]

      prev_c <- latest_c

      for (rr in seq_len(nrow(tmp))) {

        dd <- tmp$DP[rr]

        ff <- dev_factors[
          DP == dd - 1L,
          f
        ]

        if (length(ff) == 0L || is.na(ff) || !is.finite(ff)) {
          ff <- 1
        }

        new_c <- prev_c * ff

        tmp$IBNR[rr] <- new_c - prev_c

        prev_c <- new_c
      }

      cl_rows[[ii]] <- tmp
    }

    cl_pred <- data.table::rbindlist(cl_rows, fill = TRUE)

    if (nrow(cl_pred) == 0L) {
      cl_pred <- data.table::data.table(
        AP = integer(),
        DP = integer(),
        CP = integer(),
        IBNR = numeric()
      )
    }

    cl_pred <- cl_pred[
      CP > max_dp
    ]

    cl_pred <- cl_pred[
      ,
      .(IBNR = sum(IBNR, na.rm = TRUE)),
      by = .(AP, DP, CP)
    ]

    prediction_list <- c(
      list(CL = cl_pred),
      prediction_list
    )
  }

  ## ------------------------------------------------------------------
  ## Deterministic metrics
  ## ------------------------------------------------------------------

  deterministic_metrics <- intersect(
    scoring_metrics,
    c("EI", "R-tot", "R-cell-wise", "R-cal-wise")
  )

  out <- list()

  for (metric in deterministic_metrics) {

    tmp_metric <- data.table::data.table(
      model = character(),
      metric = character(),
      score = numeric()
    )

    for (mm in names(prediction_list)) {

      pred_dt <- prediction_list[[mm]]

      score_dt <- merge(
        pred_dt,
        actual_dt,
        by = c("AP", "DP", "CP"),
        all = TRUE
      )

      score_dt[
        is.na(IBNR),
        IBNR := 0
      ]

      score_dt[
        is.na(actual),
        actual := 0
      ]

      actual_total <- score_dt[
        ,
        sum(actual, na.rm = TRUE)
      ]

      predicted_total <- score_dt[
        ,
        sum(IBNR, na.rm = TRUE)
      ]

      if (!is.finite(actual_total) || actual_total == 0) {

        score_value <- NA_real_

      } else {

        if (metric == "EI") {

          score_value <- predicted_total / actual_total

        } else if (metric == "R-tot") {

          score_value <- abs(predicted_total - actual_total) / actual_total

        } else if (metric == "R-cell-wise") {

          score_value <- score_dt[
            ,
            sum(abs(IBNR - actual), na.rm = TRUE)
          ] / actual_total

        } else if (metric == "R-cal-wise") {

          cal_dt <- score_dt[
            ,
            .(
              IBNR = sum(IBNR, na.rm = TRUE),
              actual = sum(actual, na.rm = TRUE)
            ),
            by = CP
          ]

          score_value <- cal_dt[
            ,
            sum(abs(IBNR - actual), na.rm = TRUE)
          ] / actual_total
        }
      }

      tmp_metric <- data.table::rbindlist(
        list(
          tmp_metric,
          data.table::data.table(
            model = mm,
            metric = metric,
            score = as.numeric(score_value)
          )
        ),
        use.names = TRUE
      )
    }

    out[[metric]] <- tmp_metric
  }

  ## ------------------------------------------------------------------
  ## CRPS placeholder guard
  ## ------------------------------------------------------------------

  ## ------------------------------------------------------------------
  ## CRPS for ReSurvFit models only
  ## ------------------------------------------------------------------

  if ("CRPS" %in% scoring_metrics) {

    crps_out <- data.table::data.table(
      model = character(),
      metric = character(),
      score = numeric()
    )

    for (mm in names(models)) {

      object <- models[[mm]]

      if (!inherits(object, "ReSurvFit")) {
        next
      }

      data_information <- object$data_information

      categorical_features <- data_information$categorical_features
      continuous_features <- data_information$continuous_features
      years <- data_information$years
      input_time_granularity <- data_information$input_time_granularity

      time_unit_string <- c("days", "months", "quarters", "semesters", "years")
      time_unit_numeric <- c(1 / 360, 1 / 12, 1 / 4, 1 / 2, 1)

      input_pos <- match(input_time_granularity, time_unit_string)

      if (is.na(input_pos)) {
        stop(
          "`input_time_granularity` must be one of: ",
          paste(time_unit_string, collapse = ", "),
          call. = FALSE
        )
      }

      max_dp_i_crps <- as.integer(years / time_unit_numeric[input_pos])

      hazard_frame <- data.table::copy(
        data.table::as.data.table(object$hazard_frame)
      )

      if (!("DP_rev_i" %in% names(hazard_frame))) {
        if ("DP_i" %in% names(hazard_frame)) {
          hazard_frame[
            ,
            DP_rev_i := max_dp_i_crps - as.integer(DP_i) + 1L
          ]
        } else {
          stop(
            "CRPS requires `object$hazard_frame` to contain `DP_rev_i` or `DP_i`.",
            call. = FALSE
          )
        }
      }

      if (!("S_i" %in% names(hazard_frame))) {

        if ("cum_f_i" %in% names(hazard_frame)) {
          hazard_frame[
            ,
            S_i := 1 / as.numeric(cum_f_i)
          ]
        } else if ("cum_dev_f_i" %in% names(hazard_frame)) {
          hazard_frame[
            ,
            S_i := 1 / as.numeric(cum_dev_f_i)
          ]
        } else {
          stop(
            "CRPS requires `S_i`, `cum_f_i`, or `cum_dev_f_i` in `object$hazard_frame`.",
            call. = FALSE
          )
        }
      }

      hazard_frame[
        ,
        `:=`(
          DP_rev_i = as.integer(DP_rev_i),
          S_i = as.numeric(S_i)
        )
      ]

      hazard_frame[
        S_i < 0,
        S_i := 0
      ]

      hazard_frame[
        S_i > 1,
        S_i := 1
      ]

      if (!is.null(categorical_features)) {
        for (cc in categorical_features) {
          if (cc %in% names(hazard_frame)) {
            hazard_frame[
              ,
              (cc) := as.character(get(cc))
            ]
          }
        }
      }

      test_for_crps <- data.table::copy(
        data.table::as.data.table(newdata)
      )

      if (nrow(test_for_crps) == 0L) {
        stop(
          "`newdata` contains no observations for CRPS.",
          call. = FALSE
        )
      }

      if (!("AP_i" %in% names(test_for_crps))) {
        if ("AP" %in% names(test_for_crps)) {
          test_for_crps[
            ,
            AP_i := as.integer(AP)
          ]
        } else {
          stop(
            "CRPS requires `newdata` to contain `AP_i` or raw `AP`.",
            call. = FALSE
          )
        }
      } else {
        test_for_crps[
          ,
          AP_i := as.integer(AP_i)
        ]
      }

      if (!("RP_i" %in% names(test_for_crps))) {
        if ("RP" %in% names(test_for_crps)) {
          test_for_crps[
            ,
            RP_i := as.integer(RP)
          ]
        }
      } else {
        test_for_crps[
          ,
          RP_i := as.integer(RP_i)
        ]
      }

      if (!("DP_i" %in% names(test_for_crps))) {

        if ("RP_i" %in% names(test_for_crps)) {

          test_for_crps[
            ,
            DP_i := RP_i - AP_i + 1L
          ]

        } else if ("DP_rev_i" %in% names(test_for_crps)) {

          test_for_crps[
            ,
            DP_i := max_dp_i_crps - as.integer(DP_rev_i) + 1L
          ]

        } else {

          stop(
            "CRPS requires input-scale development information in `newdata`: ",
            "`DP_i`, `DP_rev_i`, or `RP`/`RP_i`.",
            call. = FALSE
          )
        }

      } else {

        test_for_crps[
          ,
          DP_i := as.integer(DP_i)
        ]
      }

      if (!("RP_i" %in% names(test_for_crps))) {
        test_for_crps[
          ,
          RP_i := AP_i + DP_i - 1L
        ]
      }

      if (!("DP_rev_i" %in% names(test_for_crps))) {
        test_for_crps[
          ,
          DP_rev_i := max_dp_i_crps - DP_i + 1L
        ]
      } else {
        test_for_crps[
          ,
          DP_rev_i := as.integer(DP_rev_i)
        ]
      }

      if (!("TR_i" %in% names(test_for_crps))) {
        test_for_crps[
          ,
          TR_i := AP_i - 1L
        ]
      } else {
        test_for_crps[
          ,
          TR_i := as.integer(TR_i)
        ]
      }

      if (!("I" %in% names(test_for_crps))) {
        test_for_crps[
          ,
          I := 1
        ]
      } else {
        test_for_crps[
          ,
          I := as.numeric(I)
        ]
      }

      if (!is.null(categorical_features)) {
        for (cc in categorical_features) {
          if (cc %in% names(test_for_crps)) {
            test_for_crps[
              ,
              (cc) := as.character(get(cc))
            ]
          }
        }
      }

      test_for_crps <- test_for_crps[
        !is.na(DP_rev_i) &
          !is.na(TR_i) &
          DP_rev_i <= TR_i &
          DP_rev_i >= 1L &
          DP_rev_i <= max_dp_i_crps &
          !is.na(I) &
          I > 0
      ]

      if (nrow(test_for_crps) == 0L) {
        warning(
          "No lower-triangle observations are available for CRPS for model `",
          mm,
          "`.",
          call. = FALSE
        )

        crps_out <- data.table::rbindlist(
          list(
            crps_out,
            data.table::data.table(
              model = mm,
              metric = "CRPS",
              score = NA_real_
            )
          ),
          use.names = TRUE
        )

        next
      }

      group_cols <- unique(c(
        categorical_features,
        continuous_features
      ))

      group_cols <- group_cols[
        group_cols %in% names(hazard_frame)
      ]

      missing_group_cols <- setdiff(
        group_cols,
        names(test_for_crps)
      )

      if (length(missing_group_cols) > 0L) {
        stop(
          "CRPS could not match `newdata` to the ReSurv hazard groups. ",
          "Missing columns in `newdata`: ",
          paste(missing_group_cols, collapse = ", "),
          call. = FALSE
        )
      }

      if (length(group_cols) == 0L) {

        hazard_frame[
          ,
          crps_group := 1L
        ]

        test_for_crps[
          ,
          crps_group := 1L
        ]

      } else {

        hazard_frame[
          ,
          crps_group := .GRP,
          by = group_cols
        ]

        group_map <- unique(
          hazard_frame[
            ,
            c(group_cols, "crps_group"),
            with = FALSE
          ]
        )

        test_for_crps <- group_map[
          test_for_crps,
          on = group_cols
        ]

        if (anyNA(test_for_crps$crps_group)) {

          n_unmatched <- test_for_crps[
            is.na(crps_group),
            .N
          ]

          warning(
            n_unmatched,
            " observations in `newdata` could not be matched to a ReSurv hazard group ",
            "and are omitted from CRPS.",
            call. = FALSE
          )

          test_for_crps <- test_for_crps[
            !is.na(crps_group)
          ]
        }
      }

      if (nrow(test_for_crps) == 0L) {
        warning(
          "No observations could be matched to the ReSurv hazard groups for model `",
          mm,
          "`.",
          call. = FALSE
        )

        crps_out <- data.table::rbindlist(
          list(
            crps_out,
            data.table::data.table(
              model = mm,
              metric = "CRPS",
              score = NA_real_
            )
          ),
          use.names = TRUE
        )

        next
      }

      curve_dt <- hazard_frame[
        !is.na(crps_group) &
          !is.na(DP_rev_i) &
          !is.na(S_i),
        .(
          S_i = mean(S_i, na.rm = TRUE)
        ),
        by = .(crps_group, DP_rev_i)
      ]

      curve_dt <- curve_dt[
        DP_rev_i >= 1L &
          DP_rev_i <= max_dp_i_crps
      ]

      data.table::setorder(curve_dt, crps_group, DP_rev_i)

      curve_dt[
        ,
        delta := DP_rev_i - data.table::shift(
          DP_rev_i,
          fill = 0L
        ),
        by = crps_group
      ]

      curve_dt[
        !is.finite(delta) | is.na(delta) | delta <= 0,
        delta := 1
      ]

      test_for_crps[
        ,
        id_crps := .I
      ]

      test_for_crps_small <- test_for_crps[
        ,
        .(
          id_crps,
          crps_group,
          obs_DP_rev_i = as.integer(DP_rev_i),
          weight = as.numeric(I)
        )
      ]

      crps_long <- curve_dt[
        test_for_crps_small,
        on = "crps_group",
        allow.cartesian = TRUE
      ]

      crps_long[
        ,
        term := data.table::fifelse(
          DP_rev_i < obs_DP_rev_i,
          delta * (1 - S_i)^2,
          data.table::fifelse(
            DP_rev_i > obs_DP_rev_i,
            delta * S_i^2,
            0.5 * delta * ((1 - S_i)^2 + S_i^2)
          )
        )
      ]

      crps_by_obs <- crps_long[
        ,
        .(
          crps = sum(term, na.rm = TRUE),
          weight = unique(weight)[1L]
        ),
        by = id_crps
      ]

      crps_score <- crps_by_obs[
        ,
        stats::weighted.mean(
          x = crps,
          w = weight,
          na.rm = TRUE
        )
      ]

      crps_out <- data.table::rbindlist(
        list(
          crps_out,
          data.table::data.table(
            model = mm,
            metric = "CRPS",
            score = as.numeric(crps_score)
          )
        ),
        use.names = TRUE
      )
    }

    out[["CRPS"]] <- crps_out
  }

  attr(out, "granularity") <- granularity
  class(out) <- "Score_Reserving"

  out
}
