#' Predict IBNR frequency
#'
#' This function predicts the results from the ReSurv fits.
#'
#' @param object \code{ResurvFit} object specifying start time, end time and status.
#' @param newdata \code{IndividualDataPP} object that contains new data to predict.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#' @param lower_triangular_output \code{logical}, if set to \code{TRUE} we add the predicted lower triangle in input and output granularity to the \code{predict.ReSurvFit} output.
#' @param minimal_output \code{logical}, if set to \code{TRUE} return a reduced prediction object.
#' @param ... Additional arguments to pass to the predict function.
#'
#'
#' @return A \code{ReSurvPredict} object with fitted predictions, long triangle outputs, predicted counts, and optional lower-triangle outputs.
#' @importFrom dplyr bind_rows distinct relocate arrange
#' @export
#' @method predict ReSurvFit
predict.ReSurvFit <- function(object,
                              newdata = NULL,
                              lower_triangular_output = FALSE,
                              minimal_output = FALSE,
                              check_value = 1.85,
                              ...) {

  ## ------------------------------------------------------------------
  ## Metadata and validation
  ## ------------------------------------------------------------------

  if (!is.null(newdata)) {

    if (!identical(object$data_information$input_time_unit,
                   newdata$data_information$input_time_unit)) {
      stop("newdata must have the same input granularity as pastdata.")
    }

    if (!inherits(newdata, "IndividualDataPP")) {
      stop("newdata must be an IndividualDataPP object.")
    }

    newfeatures <- c(
      newdata$data_information$categorical_features,
      newdata$data_information$continuous_features
    )

    pastfeatures <- c(
      object$data_information$categorical_features,
      object$data_information$continuous_features
    )

    if (!identical(newfeatures, pastfeatures)) {
      stop("newdata must have the same features as pastdata.")
    }

    data_information <- newdata$data_information

  } else {

    data_information <- object$data_information
  }

  conversion_factor <- data_information$conversion_factor
  continuous_features <- data_information$continuous_features
  categorical_features <- data_information$categorical_features
  calendar_period_extrapolation <- data_information$calendar_period_extrapolation
  years <- data_information$years
  input_time_granularity <- data_information$input_time_granularity
  output_time_granularity <- data_information$output_time_granularity

  ## Inline pkg.env$maximum.time()

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

  max_dp_i <- as.integer(years / time_unit_numeric[input_pos])
  max_DP <- as.integer(years / time_unit_numeric[output_pos])

  ## ------------------------------------------------------------------
  ## Hazard frame preparation
  ## ------------------------------------------------------------------

  hazard_frame <- data.table::copy(
    data.table::as.data.table(object$hazard_frame)
  )

  if ("DP_i" %in% names(hazard_frame)) {
    hazard_frame[, DP_i := NULL]
  }

  old_nm <- intersect(c("f_i", "cum_f_i"), names(hazard_frame))

  if (length(old_nm) > 0L) {
    new_nm <- c(f_i = "dev_f_i", cum_f_i = "cum_dev_f_i")[old_nm]
    data.table::setnames(hazard_frame, old = old_nm, new = unname(new_nm))
  }

  ## ------------------------------------------------------------------
  ## Covariate mapping
  ## ------------------------------------------------------------------

  if (
    (length(continuous_features) == 1L && "AP_i" %in% continuous_features) ||
    (length(continuous_features) == 1L && "RP_i" %in% continuous_features) ||
    (length(continuous_features) == 2L &&
     sum(c("AP_i", "RP_i") %in% continuous_features) == 2L)
  ) {
    continuous_features_group <- NULL
  } else {
    continuous_features_group <- continuous_features[
      !(continuous_features %in% c("AP_i", "RP_i"))
    ]
  }

  data.table::setDT(hazard_frame)

  feature_cols <- unique(c(categorical_features, continuous_features_group))

  if (is.null(feature_cols) || length(feature_cols) == 0L) {
    hazard_frame[, covariate := "0"]
  } else {
    hazard_frame[
      ,
      covariate := do.call(paste, c(.SD, sep = "_")),
      .SDcols = feature_cols
    ]
  }

  time_features <- continuous_features[
    continuous_features %in% c("AP_i", "RP_i")
  ]

  if (length(time_features) > 0L) {

    group_cols <- c(time_features, "covariate")

    groups <- unique(
      hazard_frame[
        ,
        .SD,
        .SDcols = group_cols
      ]
    )

    groups[, group_i := .I]

    hazard_group <- merge(
      hazard_frame,
      groups,
      by = group_cols,
      all.x = TRUE
    )

  } else {

    groups <- unique(hazard_frame[, .(covariate)])
    groups[, group_i := .I]

    hazard_group <- merge(
      hazard_frame,
      groups,
      by = "covariate",
      all.x = TRUE
    )
  }

  groups[, group_o := group_i]

  if (
    conversion_factor != 1 &&
    sum(c("AP_i", "RP_i") %in% continuous_features) > 0L
  ) {

    output_time_cols <- paste0(substr(time_features, 1L, 2L), "_o")

    groups_o <- hazard_group[
      ,
      c(
        lapply(.SD, function(x) ceiling(x * conversion_factor)),
        list(covariate = covariate)
      ),
      .SDcols = time_features
    ]

    data.table::setnames(
      groups_o,
      old = time_features,
      new = output_time_cols
    )

    groups_o <- unique(groups_o)
    groups_o[, group_o := .I]

    groups[, group_o := NULL]

    groups[
      ,
      (output_time_cols) := lapply(
        .SD,
        function(x) ceiling(x * conversion_factor)
      ),
      .SDcols = time_features
    ]

    groups <- merge(
      groups,
      groups_o,
      by = c(output_time_cols, "covariate"),
      all.x = TRUE
    )
  }

  hazard_frame_grouped <- list(
    hazard_group = hazard_group,
    groups = groups
  )

  ## ------------------------------------------------------------------
  ## Latest observed values
  ## ------------------------------------------------------------------

  data_reserve <- data.table::copy(
    data.table::as.data.table(object$data_information$data_for_reserving)
  )

  max_dp_by_ap <- unique(
    data_reserve[
      ,
      .(DP_max_rev = min(max(DP_rev_i) - DP_i) + 1L),
      by = AP_i
    ]
  )

  cols_to_keep <- unique(c(
    "AP_i",
    "AP_o",
    "DP_rev_i",
    "DP_i",
    categorical_features,
    continuous_features,
    "I"
  ))

  cols_to_keep <- intersect(cols_to_keep, names(data_reserve))

  data_reserve2 <- data_reserve[, .SD, .SDcols = cols_to_keep]
  data_reserve2[, AP_i := as.numeric(AP_i)]

  data_reserve2 <- max_dp_by_ap[data_reserve2, on = "AP_i"]

  if (is.null(continuous_features)) {

    group_cols <- c(
      categorical_features,
      "AP_i",
      "AP_o",
      "DP_max_rev"
    )

    observed_so_far <- data_reserve2[
      ,
      .(latest_I = sum(I)),
      by = group_cols
    ]

    group_cols <- c(
      "AP_i",
      "AP_o",
      categorical_features,
      "DP_rev_i",
      "DP_i"
    )

    observed_dp_rev_i <- data_reserve2[
      ,
      .(I = sum(I)),
      by = group_cols
    ]

    feats <- c(categorical_features, continuous_features)

    if (is.null(feats) || length(feats) == 0L) {

      observed_so_far[, covariate := "0"]
      observed_dp_rev_i[, covariate := "0"]

    } else {

      observed_so_far[
        ,
        covariate := do.call(paste, c(.SD, sep = "_")),
        .SDcols = feats
      ]

      observed_dp_rev_i[
        ,
        covariate := do.call(paste, c(.SD, sep = "_")),
        .SDcols = feats
      ]
    }

    latest_cumulative <- hazard_frame_grouped$groups[
      observed_so_far,
      on = "covariate"
    ][
      ,
      .(AP_i, group_i, DP_max_rev, latest_I)
    ]

    observed_pr_dp <- hazard_frame_grouped$groups[
      observed_dp_rev_i,
      on = "covariate"
    ][
      ,
      .(AP_i, group_i, DP_rev_i, DP_i, I)
    ]

  } else {

    group_cols <- unique(c(
      "AP_i",
      "AP_o",
      categorical_features,
      continuous_features,
      "DP_max_rev"
    ))

    observed_so_far <- data_reserve2[
      ,
      .(latest_I = sum(I)),
      by = group_cols
    ]

    group_cols <- unique(c(
      "AP_i",
      "AP_o",
      categorical_features,
      continuous_features,
      "DP_rev_i",
      "DP_i"
    ))

    observed_dp_rev_i <- data_reserve2[
      ,
      .(I = sum(I)),
      by = group_cols
    ]

    feature_cols <- c(
      categorical_features,
      continuous_features[continuous_features != "AP_i"]
    )

    if (is.null(feature_cols) || length(feature_cols) == 0L) {

      observed_so_far[, covariate := "0"]
      observed_dp_rev_i[, covariate := "0"]

    } else {

      observed_so_far[
        ,
        covariate := do.call(paste, c(.SD, sep = "_")),
        .SDcols = feature_cols
      ]

      observed_dp_rev_i[
        ,
        covariate := do.call(paste, c(.SD, sep = "_")),
        .SDcols = feature_cols
      ]
    }

    time_features <- continuous_features[
      continuous_features %in% c("AP_i", "RP_i")
    ]

    latest_cumulative <- hazard_frame_grouped$groups[
      observed_so_far,
      on = c(time_features, "covariate")
    ]

    latest_cumulative <- latest_cumulative[
      ,
      c(unique(c("AP_i", time_features)), "group_i", "DP_max_rev", "latest_I"),
      with = FALSE
    ]

    observed_dp_rev_i <- hazard_frame_grouped$groups[
      observed_dp_rev_i,
      on = c(time_features, "covariate")
    ]

    observed_dp_rev_i <- observed_dp_rev_i[
      ,
      c(unique(c("AP_i", time_features)), "group_i", "DP_rev_i", "DP_i", "I"),
      with = FALSE
    ]

    observed_dp_rev_i <- observed_dp_rev_i[
      hazard_frame_grouped$groups[
        ,
        c(time_features, "group_i"),
        with = FALSE
      ],
      on = c(time_features, "group_i"),
      nomatch = 0
    ]

    observed_pr_dp <- observed_dp_rev_i
  }
  ## ------------------------------------------------------------------
  ## Expected IBNR on input granularity
  ## ------------------------------------------------------------------

  grouped_hazard_0 <- latest_cumulative[
    hazard_frame_grouped$hazard_group,
    on = c("group_i", "AP_i")
  ]

  expected_i <- grouped_hazard_0[
    ,
    .(
      DP_rev_i,
      AP_i,
      group_i,
      S_i,
      S_i_lag,
      DP_max_rev,
      latest_I,
      DP_max_rev_keep = DP_max_rev
    )
  ]

  hazard_tmp <- hazard_frame_grouped$hazard_group[
    ,
    .(
      DP_rev_i_key = DP_rev_i + 1L,
      AP_i,
      group_i,
      S_ultimate_i = S_i
    )
  ]

  expected_i <- hazard_tmp[
    expected_i,
    on = .(DP_rev_i_key = DP_max_rev, AP_i, group_i)
  ]

  expected_i[
    ,
    U := data.table::fcase(
      S_i_lag == 1,
      as.numeric(latest_I),

      DP_max_rev_keep == min(hazard_frame_grouped$hazard_group$DP_rev_i),
      as.numeric(latest_I),

      S_ultimate_i == 0,
      0.0,

      AP_i != 1,
      as.numeric(latest_I) / as.numeric(S_ultimate_i),

      default = as.numeric(latest_I)
    )
  ]

  expected_i_probability <- data.table::copy(expected_i)

  expected_i_probability[
    ,
    I_expected := 1 * (S_i_lag - S_i)
  ][
    ,
    IBNR := data.table::fifelse(DP_rev_i < DP_max_rev_keep, I_expected, NA_real_)
  ]

  expected_i[
    ,
    I_expected := U * (S_i_lag - S_i)
  ]

  expected_i[
    ,
    IBNR := data.table::fifelse(DP_rev_i < DP_max_rev_keep, I_expected, NA_real_)
  ]

  expected_i <- expected_i[
    ,
    .(AP_i, group_i, DP_rev_i, I_expected, IBNR)
  ]

  hazard_frame_input <- hazard_frame_grouped$hazard_group[
    expected_i,
    on = .(AP_i, group_i, DP_rev_i)
  ]

  hazard_frame_input[
    ,
    DP_i := max_dp_i - DP_rev_i + 1L
  ]

  ## ------------------------------------------------------------------
  ## Output granularity, if requested
  ## ------------------------------------------------------------------

  has_output_granularity <- conversion_factor != 1 && !isTRUE(minimal_output)

  if (has_output_granularity) {

    group_map <- unique(hazard_frame_grouped$groups[, .(group_i, group_o)])

    if (nrow(group_map) != data.table::uniqueN(group_map$group_i)) {
      stop(
        "Each `group_i` must map to exactly one `group_o` in output-granularity prediction.",
        call. = FALSE
      )
    }

    expected_o <- group_map[
      data.table::copy(expected_i),
      on = "group_i"
    ]

    if (anyNA(expected_o$group_o)) {
      stop(
        "Some `group_i` values in `expected_i` could not be mapped to `group_o`.",
        call. = FALSE
      )
    }

    expected_o[
      ,
      DP_i := max_dp_i - DP_rev_i + 1L
    ][
      ,
      `:=`(
        AP_o = ceiling(AP_i * conversion_factor),
        DP_rev_o = ceiling(max_dp_i * conversion_factor) -
          ceiling(
            (DP_i + (AP_i - 1L) %% (1 / conversion_factor)) *
              conversion_factor
          ) + 1L
      )
    ]

    expected_o <- expected_o[
      ,
      .(
        I_expected = sum(I_expected, na.rm = TRUE),
        IBNR = sum(IBNR, na.rm = TRUE)
      ),
      by = .(AP_o, DP_rev_o, group_o)
    ][
      ,
      .(AP_o, group_o, DP_rev_o, I_expected, IBNR)
    ]

    hazard_data_frame <- group_map[
      hazard_frame_grouped$hazard_group,
      on = "group_i"
    ]

    observed_pr_dp_o <- group_map[
      observed_pr_dp,
      on = "group_i"
    ]

    latest_cumulative_o <- group_map[
      latest_cumulative,
      on = "group_i"
    ][
      ,
      .(latest_I = sum(latest_I, na.rm = TRUE)),
      by = .(AP_i, group_o, DP_max_rev)
    ]

    expected_i_output_probability <- group_map[
      expected_i_probability,
      on = "group_i"
    ]

    grouped_hazard_0 <- data.table::copy(hazard_frame_grouped$hazard_group)[
      ,
      DP_i := max_dp_i - DP_rev_i + 1L
    ][
      ,
      DP_rev_o := floor(max_dp_i * conversion_factor) -
        ceiling(
          DP_i * conversion_factor +
            ((AP_i - 1L) %% (1 / conversion_factor)) *
            conversion_factor
        ) +
        1L
    ][
      DP_rev_o > 0
    ]

    grouped_hazard_0 <- object$data_information$dp_ranges[
      grouped_hazard_0,
      on = .(AP_i, DP_rev_o)
    ]

    grouped_hazard_0[
      group_map,
      on = "group_i",
      group_o := i.group_o
    ]

    grouped_hazard_0 <- latest_cumulative_o[
      ,
      .(group_o, AP_i, DP_max_rev, latest_I)
    ][
      grouped_hazard_0,
      on = .(group_o, AP_i),
      allow.cartesian = TRUE
    ]

    grouped_hazard_0 <- observed_pr_dp_o[
      ,
      .(group_o, AP_i, DP_rev_i, I)
    ][
      grouped_hazard_0,
      on = .(group_o, AP_i, DP_rev_i),
      allow.cartesian = TRUE
    ]

    cumulative_observed <- observed_pr_dp_o[
      order(DP_i),
      .(
        exposure = cumsum(data.table::fifelse(is.na(I), 0, I)),
        DP_rev_i = DP_rev_i - 1L,
        AP_i = AP_i,
        group_o = group_o
      ),
      by = .(AP_i, group_o)
    ][
      ,
      .(AP_i, group_o, DP_rev_i, exposure)
    ]

    exposures <- grouped_hazard_0[
      ,
      .SD[DP_rev_i == max(DP_rev_i)],
      by = .(AP_i, DP_rev_o, group_o)
    ]

    exposures[
      cumulative_observed,
      on = .(AP_i, group_o, max_dp = DP_rev_i),
      exposure := i.exposure
    ]

    no_exposure <- exposures[
      ,
      .(DP_rev_i, DP_rev_o, AP_i, group_o, S_i, DP_max_rev, latest_I)
    ]

    hazard_tmp <- hazard_data_frame[
      ,
      .(
        DP_max_rev = DP_rev_i + 1L,
        AP_i,
        group_o,
        S_ultimate_i = S_i
      )
    ]

    no_exposure <- hazard_tmp[
      no_exposure,
      on = .(DP_max_rev, AP_i, group_o)
    ]

    no_exposure[, U := 1]
    no_exposure[latest_I == 0, U := 0]
    no_exposure[, exposure_expected := U * S_i]

    no_exposure <- no_exposure[
      ,
      .(AP_i, group_o, DP_rev_o, DP_rev_i, exposure_expected)
    ]

    exposures_combined <- no_exposure[
      exposures,
      on = .(AP_i, DP_rev_o, DP_rev_i, group_o)
    ]

    exposures_combined[
      ,
      exposure_combined := data.table::fifelse(
        !is.na(exposure_expected),
        exposure_expected,
        0
      )
    ]

    expected_i_for_output <- expected_i_output_probability[
      ,
      .SD,
      .SDcols = setdiff(
        colnames(expected_i_output_probability),
        c("group_i", "S_i", "S_i_lag", "latest_I")
      )
    ]

    grouped_hazard_1 <- expected_i_for_output[
      grouped_hazard_0,
      on = .(AP_i, group_o, DP_rev_i),
      allow.cartesian = TRUE
    ]

    grouped_hazard_1[
      ,
      I_combined := data.table::fifelse(!is.na(I_expected), I_expected, 0)
    ]

    grouped_hazard_2 <- grouped_hazard_1[
      ,
      .(observed = sum(I_combined, na.rm = TRUE)),
      by = .(AP_i, DP_rev_o, group_o)
    ]

    grouped_hazard_2 <- exposures_combined[
      grouped_hazard_2,
      on = .(AP_i, group_o, DP_rev_o),
      allow.cartesian = TRUE
    ]

    grouped_hazard_2[
      latest_I == 0,
      observed := 0
    ]

    output_dev_factor <- grouped_hazard_2[
      ,
      .(
        dev_f_o = data.table::fifelse(
          sum(exposure_combined, na.rm = TRUE) == 0,
          1,
          (
            sum(observed, na.rm = TRUE) +
              sum(exposure_combined, na.rm = TRUE)
          ) /
            sum(exposure_combined, na.rm = TRUE)
        )
      ),
      by = .(DP_rev_o, group_o)
    ]

    output_dev_factor[
      ,
      DP_o := max_DP - DP_rev_o + 1L
    ]

    cols_to_keep <- unique(c(
      "AP_i",
      "covariate",
      "group_i",
      categorical_features,
      continuous_features
    ))

    cols_to_keep <- intersect(names(hazard_frame_input), cols_to_keep)

    dt1 <- data.table::copy(hazard_frame_input)[
      ,
      .SD,
      .SDcols = cols_to_keep
    ]

    dtm <- hazard_frame_grouped$groups[
      dt1,
      on = .(covariate, group_i)
    ]

    if (!("group_o" %in% names(dtm))) {
      dtm[, group_o := group_i]
    }

    if (!("AP_o" %in% names(dtm))) {
      dtm[, AP_o := ceiling(AP_i * conversion_factor)]
    }

    data.table::setorder(dtm, group_o, AP_i)

    output_feature_cols <- unique(c(
      "group_o",
      "AP_o",
      categorical_features,
      setdiff(continuous_features, c("AP_i"))
    ))

    output_feature_cols <- intersect(output_feature_cols, names(dtm))

    final_result <- unique(
      dtm[
        !is.na(group_o),
        .SD,
        .SDcols = output_feature_cols
      ]
    )

    if (!("AP_o" %in% names(final_result))) {

      ap_group_keys <- unique(expected_o[, .(AP_o, group_o)])

      final_result <- final_result[
        ap_group_keys,
        on = "group_o",
        allow.cartesian = TRUE
      ]

    } else {

      data.table::setorderv(
        final_result,
        intersect(c("AP_o", "group_o"), names(final_result))
      )

      final_result <- unique(
        final_result,
        by = intersect(c("AP_o", "group_o"), names(final_result))
      )
    }

    hazard_frame_output <- expected_o[
      final_result,
      on = .(AP_o, group_o),
      allow.cartesian = TRUE
    ]

    hazard_frame_output <- output_dev_factor[
      hazard_frame_output,
      on = .(DP_rev_o, group_o)
    ]

    hazard_frame_output[
      is.na(dev_f_o),
      dev_f_o := 1
    ]
  }

  ## ------------------------------------------------------------------
  ## Final output formatting
  ## ------------------------------------------------------------------

  old_nm <- intersect(c("dev_f_i", "I_expected"), names(hazard_frame_input))

  if (length(old_nm) > 0L) {
    new_nm <- c(dev_f_i = "f_i", I_expected = "expected_counts")[old_nm]
    data.table::setnames(hazard_frame_input, old = old_nm, new = unname(new_nm))
  }

  drop_cols_input <- intersect(
    c("expg", "baseline", "hazard", "DP_rev_i"),
    names(hazard_frame_input)
  )

  if (length(drop_cols_input) > 0L) {
    hazard_frame_input[, (drop_cols_input) := NULL]
  }

  long_triangle_format_out <- list(
    input_granularity = hazard_frame_input
  )

  if (has_output_granularity) {

    old_nm <- intersect(c("dev_f_o", "I_expected"), names(hazard_frame_output))

    if (length(old_nm) > 0L) {
      new_nm <- c(dev_f_o = "f_o", I_expected = "expected_counts")[old_nm]
      data.table::setnames(hazard_frame_output, old = old_nm, new = unname(new_nm))
    }

    drop_cols_output <- intersect("DP_rev_o", names(hazard_frame_output))

    if (length(drop_cols_output) > 0L) {
      hazard_frame_output[, (drop_cols_output) := NULL]
    }

    long_triangle_format_out$output_granularity <- hazard_frame_output
  }

  out <- list(
    ReSurvFit = object,
    long_triangle_format_out = long_triangle_format_out,
    predicted_counts = sum(hazard_frame_input$IBNR, na.rm = TRUE)
  )

  ## ------------------------------------------------------------------
  ## Lower triangle output
  ## ------------------------------------------------------------------

  if (lower_triangular_output) {

    ltr_input_long <- data.table::copy(hazard_frame_input)[
      ,
      .(value = sum(IBNR, na.rm = TRUE)),
      by = .(AP_i, DP_i)
    ]

    ltr_input <- matrix(
      0,
      nrow = max_dp_i,
      ncol = max_dp_i,
      dimnames = list(
        as.character(seq_len(max_dp_i)),
        as.character(seq_len(max_dp_i))
      )
    )

    idx_input <- !is.na(ltr_input_long$AP_i) &
      !is.na(ltr_input_long$DP_i) &
      ltr_input_long$AP_i >= 1L &
      ltr_input_long$AP_i <= max_dp_i &
      ltr_input_long$DP_i >= 1L &
      ltr_input_long$DP_i <= max_dp_i

    ltr_input[
      cbind(
        as.integer(ltr_input_long$AP_i[idx_input]),
        as.integer(ltr_input_long$DP_i[idx_input])
      )
    ] <- ltr_input_long$value[idx_input]

    input_mask <- outer(
      seq_len(max_dp_i),
      seq_len(max_dp_i),
      function(i, j) i + j - 1L <= max_dp_i
    )

    ltr_input[input_mask] <- NA_real_

    lower_triangle <- list(
      input_granularity = data.table::as.data.table(ltr_input)
    )

    if (has_output_granularity) {

      ltr_output_long <- data.table::copy(hazard_frame_output)[
        ,
        .(value = sum(IBNR, na.rm = TRUE)),
        by = .(AP_o, DP_o)
      ]

      ltr_output <- matrix(
        0,
        nrow = max_DP,
        ncol = max_DP,
        dimnames = list(
          as.character(seq_len(max_DP)),
          as.character(seq_len(max_DP))
        )
      )

      idx_output <- !is.na(ltr_output_long$AP_o) &
        !is.na(ltr_output_long$DP_o) &
        ltr_output_long$AP_o >= 1L &
        ltr_output_long$AP_o <= max_DP &
        ltr_output_long$DP_o >= 1L &
        ltr_output_long$DP_o <= max_DP

      ltr_output[
        cbind(
          as.integer(ltr_output_long$AP_o[idx_output]),
          as.integer(ltr_output_long$DP_o[idx_output])
        )
      ] <- ltr_output_long$value[idx_output]

      output_mask <- outer(
        seq_len(max_DP),
        seq_len(max_DP),
        function(i, j) i + j - 1L <= max_DP
      )

      ltr_output[
        output_mask & (is.na(ltr_output) | ltr_output == 0)
      ] <- NA_real_

      cut_point <- ceiling(max_dp_i * conversion_factor)

      ltr_output <- ltr_output[seq_len(cut_point), , drop = FALSE]

      lower_triangle$output_granularity <- data.table::as.data.table(ltr_output)
    }

    out[["lower_triangle"]] <- lower_triangle
  }

  class(out) <- c("ReSurvPredict")

  return(out)
}
