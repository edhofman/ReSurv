#' Individual Data Pre-Processing
#'
#' This function pre-processes the data for the application of a \code{ReSurv} model.
#'
#'
#'
#' Accident and reporting periods are indexed from one. Development time is
#' \code{DP_i = RP_i - AP_i + 1}; reverse development time is
#' \code{DP_rev_i = DP_max - DP_i + 1}, and truncation time is
#' \code{TR_i = AP_i - 1}. Training retains observed rows with
#' \code{DP_rev_i > TR_i}.
#'
#' The conversion factor is the input time unit divided by the output time
#' unit (for example, 1/3 for months to quarters). Accident and calendar
#' periods are grouped using ceiling; development periods also account for
#' the position of the accident period within each output period.
#' Days use a 360-day year for the development horizon.
#'
#' @param data \code{data.frame}, for the individual reserving. The number of development periods can be larger than the number of accident periods.
#' @param id \code{character}, \code{data} column that contains the policy identifier. If \code{NULL} (default), we assume that each row is an observation. We assume that each observation can only have one reporting time, if not null we take the reporting time of the first row for each \code{id}.
#' @param continuous_features \code{character}, continuous features columns to be scaled.
#' @param categorical_features \code{character}, categorical features columns to be one-hot encoded.
#' @param accident_period \code{character}, it contains the name of the column in data corresponding to the accident period.
#' @param calendar_period \code{character}, it contains the name of the column in data corresponding to the calendar period.
#' @param calendar_period_extrapolation \code{logical}, whether a spline for calendar extrapolation should be considered in the cox model fit.
#'                                       Default is `FALSE`.
#' @param input_time_granularity \code{character}, time unit of the input data. Granularity supported:
#' \itemize{
#' \item{\code{"days"}: the input data are daily.}
#' \item{\code{"months"}: the input data are monthly.}
#' \item{\code{"quarters"}: the input data are quarterly}
#' \item{\code{"semesters"}: six-month periods.}
#' \item{\code{"years"}: the input data are yearly.}
#' }
#' Default to \code{months}.
#'
#' @param output_time_granularity \code{character}, time unit of the output data. The granularity supported is the same as for the input data:
#'  \itemize{
#'  \item{\code{"days"}: the output data will be on a daily scale.}
#' \item{\code{"months"}: the output data will be on a monthly scale.}
#' \item{\code{"quarters"}: the output data will be on a quarterly scale.}
#' \item{\code{"semesters"}: six-month periods.}
#' \item{\code{"years"}: the output data will be on yearly scale.}
#' }
#' The output granularity must be equal to or coarser than the input granularity.
#' Also, the output granularity must be consistent with the input granularity, meaning that the time conversion must be possible.
#' E.g., it is possible to group quarters to years. Quarters can also be grouped to semesters.
#' Default to \code{quarters}.
#'
#' @param years \code{numeric}, number of development years in the study.
#' @param continuous_features_spline \code{character}, names of continuous features to model with splines; NULL uses linear terms. Use \code{"AP_i"} for a remapped accident-period feature.
#' @param degrees_cf \code{numeric}, degrees of the spline for smoothing continuous features.
#' @param degrees_of_freedom_cf \code{numeric}, degrees of freedom of the splines for smoothing continuous features.
#' @param degrees_cp \code{numeric}, degrees of the spline for smoothing the calendar period effect.
#' @param degrees_of_freedom_cp \code{numeric}, degrees of freedom of the splines for smoothing the calendar period effect.
#'
#'
#'
#'
#' @return An \code{IndividualDataPP} list containing \code{training.data}
#'   (observed rows), \code{full.data} (all encoded rows), and
#'   \code{data_information} (conversion factor, input/output formulas,
#'   feature names, time units, horizon, and original column names).
#'
#' After pre-processing, we provide a standard encoding for the time components. This regards the output in \code{training.data}.
#' In the \code{ReSurv} notation:
#'\itemize{
#'\item{\code{AP_i}: Input granularity accident period.}
#'\item{\code{AP_o}: Output granularity accident period.}
#'\item{\code{DP_i}: Input granularity development period in forward time.}
#'\item{\code{DP_rev_i}: Input granularity development period in reverse time.}
#'\item{\code{DP_rev_o}: Output granularity development period in reverse time.}
#'\item{\code{TR_i}: Input granularity truncation time.}
#'\item{\code{TR_o}: Output granularity truncation time.}
#'\item{\code{I}: event indicator, under this framework is equal to one for each entry. }
#'}
#'
#'
#'@examples
#'
#'input_data_0 <- data_generator(
#'random_seed = 1964,
#'scenario = "alpha",
#'time_unit = 1,
#'years = 2,
#'period_exposure = 100)
#'
#'individual_data <- IndividualDataPP(data = input_data_0,
#'categorical_features = "claim_type",
#'continuous_features = "AP",
#'accident_period = "AP",
#'calendar_period = "RP",
#'input_time_granularity = "years",
#'output_time_granularity = "years",
#'years = 2)
#'
#'
#'
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr slice_head
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr row_number
#' @importFrom dplyr summarize arrange relocate
#' @importFrom purrr map_df
#' @importFrom stats as.formula
#' @importFrom data.table data.table
#'
#'
#'
#' @references
#' Hiabu, M., Hofman, E., & Pittarello, G. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' @export
IndividualDataPP <- function(data,
                             id = NULL,
                             continuous_features = NULL,
                             categorical_features = NULL,
                             accident_period,
                             calendar_period,
                             input_time_granularity = "months",
                             output_time_granularity = "quarters",
                             years = NULL,
                             calendar_period_extrapolation = FALSE,
                             continuous_features_spline = NULL,
                             degrees_cf = 3,
                             degrees_of_freedom_cf = 4,
                             degrees_cp = 3,
                             degrees_of_freedom_cp = 4) {

  ## ------------------------------------------------------------------
  ## Basic validation
  ## ------------------------------------------------------------------

  tmp <- data.table::copy(data.table::as.data.table(data))

  if (!(accident_period %in% names(tmp))) {
    stop("`accident_period` is not a column of `data`.", call. = FALSE)
  }

  if (!(calendar_period %in% names(tmp))) {
    stop("`calendar_period` is not a column of `data`.", call. = FALSE)
  }

  if (!is.null(id) && !(id %in% names(tmp))) {
    stop("`id` is not a column of `data`.", call. = FALSE)
  }

  all_features <- unique(c(continuous_features, categorical_features))

  derived_time_features <- c(
    "AP_i", "DP_i", "RP_i", "CP_i", "DP_rev_i", "TR_i",
    "AP_o", "DP_o", "RP_o", "CP_o", "DP_rev_o", "TR_o"
  )

  features_to_check <- setdiff(
    all_features,
    derived_time_features
  )

  if (length(features_to_check) > 0L) {
    missing_features <- setdiff(features_to_check, names(tmp))

    if (length(missing_features) > 0L) {
      stop(
        "The following features are not columns of `data`: ",
        paste(missing_features, collapse = ", "),
        call. = FALSE
      )
    }
  }

  time_unit_string <- c(
    "days",
    "months",
    "quarters",
    "semesters",
    "years"
  )

  time_unit_numeric <- c(
    1 / 360,
    1 / 12,
    1 / 4,
    1 / 2,
    1
  )

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

  input_numeric <- time_unit_numeric[input_pos]
  output_numeric <- time_unit_numeric[output_pos]

  if (abs((1 / input_numeric) %% (1 / output_numeric)) > .Machine$double.eps^0.5) {
    stop(
      "The provided time granularities are not subsettable.",
      call. = FALSE
    )
  }

  conversion_factor <- input_numeric / output_numeric

  ## ------------------------------------------------------------------
  ## Inline date conversion and period encoding
  ## ------------------------------------------------------------------

  if (inherits(tmp[[accident_period]], "Date")) {
    ap1 <- lubridate::floor_date(
      min(tmp[[accident_period]], na.rm = TRUE),
      "year"
    )
  } else {
    ap1 <- min(tmp[[accident_period]], na.rm = TRUE)
  }

  if (inherits(ap1, "Date")) {

    if (input_time_granularity == "quarters") {
      ap1num <- floor(
        lubridate::time_length(ap1 - ap1, "months") / 3
      )
    } else if (input_time_granularity == "semesters") {
      ap1num <- floor(
        lubridate::time_length(ap1 - ap1, "months") / 6
      )
    } else {
      ap1num <- floor(
        lubridate::time_length(ap1 - ap1, input_time_granularity)
      )
    }

  } else {
    ap1num <- ap1
  }

  if (inherits(tmp[[accident_period]], "Date")) {

    if (input_time_granularity == "quarters") {
      x.ap <- floor(
        lubridate::time_length(tmp[[accident_period]] - ap1, "months") / 3
      )
    } else if (input_time_granularity == "semesters") {
      x.ap <- floor(
        lubridate::time_length(tmp[[accident_period]] - ap1, "months") / 6
      )
    } else {
      x.ap <- floor(
        lubridate::time_length(
          tmp[[accident_period]] - ap1,
          input_time_granularity
        )
      )
    }

  } else {
    x.ap <- tmp[[accident_period]]
  }

  if (inherits(tmp[[calendar_period]], "Date")) {

    if (input_time_granularity == "quarters") {
      x.cp <- floor(
        lubridate::time_length(tmp[[calendar_period]] - ap1, "months") / 3
      )
    } else if (input_time_granularity == "semesters") {
      x.cp <- floor(
        lubridate::time_length(tmp[[calendar_period]] - ap1, "months") / 6
      )
    } else {
      x.cp <- floor(
        lubridate::time_length(
          tmp[[calendar_period]] - ap1,
          input_time_granularity
        )
      )
    }

  } else {
    x.cp <- tmp[[calendar_period]]
  }

  seq_ap <- min(c(x.ap, ap1num), na.rm = TRUE):max(x.ap, na.rm = TRUE)
  tmp.ap <- seq_along(seq_ap)[match(x.ap, seq_ap)]

  seq_cp <- ap1num:max(x.cp, na.rm = TRUE)
  tmp.cp <- seq_along(seq_cp)[match(x.cp, seq_cp)]

  tmp.dp <- tmp.cp - tmp.ap + 1L

  ## ------------------------------------------------------------------
  ## Feature-name normalization
  ## ------------------------------------------------------------------

  if (!is.null(continuous_features)) {
    continuous_features[continuous_features == accident_period] <- "AP_i"
  }

  if (!is.null(categorical_features)) {
    categorical_features[categorical_features == accident_period] <- "AP_i"
  }

  continuous_features <- unique(continuous_features)
  categorical_features <- unique(categorical_features)

  ## ------------------------------------------------------------------
  ## Missing-period warnings
  ## ------------------------------------------------------------------

  ap_diff <- diff(as.integer(sort(unique(tmp.ap))))

  if (sum(ap_diff > 1L, na.rm = TRUE) > 0L) {
    warning(
      "Some accident periods are missing in the data",
      call. = FALSE
    )
  }

  cp_diff <- diff(as.integer(sort(unique(tmp.cp))))

  if (sum(cp_diff > 1L, na.rm = TRUE) > 0L) {
    warning(
      "Some calendar periods are missing in the data",
      call. = FALSE
    )
  }

  ## ------------------------------------------------------------------
  ## Time horizon
  ## ------------------------------------------------------------------

  if (is.null(years)) {
    years <- ceiling(max(tmp.dp, na.rm = TRUE) * input_numeric)
  }

  max_dp_i <- as.integer(round(years / input_numeric))
  max_dp_o <- as.integer(round(years / output_numeric))

  ## ------------------------------------------------------------------
  ## Main encoded variables
  ## ------------------------------------------------------------------

  tmp[
    ,
    `:=`(
      AP_i = as.integer(tmp.ap),
      DP_i = as.integer(tmp.dp),
      RP_i = as.integer(tmp.cp),
      DP_rev_i = as.integer(max_dp_i - tmp.dp + 1L),
      TR_i = as.integer(tmp.ap - 1L),
      I = 1L
    )
  ]

  ## In case there is an ID, keep only the first row per ID to avoid
  ## double counting. We assume one reporting time per claim.
  if (!is.null(id)) {
    tmp <- tmp[
      ,
      .SD[1L],
      by = id
    ]
  }

  tmp[
    ,
    `:=`(
      DP_rev_o = as.integer(
        floor(max_dp_i * conversion_factor) -
          ceiling(
            DP_i * conversion_factor +
              ((AP_i - 1L) %% (1 / conversion_factor)) *
              conversion_factor
          ) +
          1L
      ),
      AP_o = as.integer(ceiling(AP_i * conversion_factor))
    )
  ]

  tmp[
    ,
    `:=`(
      TR_o = as.integer(AP_o - 1L),
      DP_o = as.integer(max_dp_o - DP_rev_o + 1L)
    )
  ]

  if (isTRUE(calendar_period_extrapolation)) {
    tmp[
      ,
      RP_o := as.integer(ceiling(RP_i * conversion_factor))
    ]
  }

  if (!is.null(categorical_features)) {
    cat_cols_present <- intersect(categorical_features, names(tmp))

    if (length(cat_cols_present) > 0L) {
      tmp[
        ,
        (cat_cols_present) := lapply(.SD, as.factor),
        .SDcols = cat_cols_present
      ]
    }
  }

  ## ------------------------------------------------------------------
  ## Training data: observed upper triangle
  ## ------------------------------------------------------------------

  train <- data.table::copy(
    tmp[
      DP_rev_i > TR_i
    ]
  )

  train_cols <- unique(c(
    id,
    categorical_features,
    continuous_features,
    "AP_i",
    if (isTRUE(calendar_period_extrapolation)) "RP_i" else NULL,
    "AP_o",
    "DP_i",
    "DP_rev_i",
    "DP_rev_o",
    "TR_i",
    "TR_o",
    "I",
    if (isTRUE(calendar_period_extrapolation)) "RP_o" else NULL
  ))

  train_cols <- intersect(train_cols, names(train))

  train <- train[
    ,
    .SD,
    .SDcols = train_cols
  ]

  ## ------------------------------------------------------------------
  ## Inline formula editor
  ## ------------------------------------------------------------------

  continuous_features_i <- continuous_features
  continuous_features_o <- continuous_features

  if (!is.null(continuous_features_o)) {
    continuous_features_o[continuous_features_o == "AP_i"] <- "AP_o"
    continuous_features_o[continuous_features_o == "RP_i"] <- "RP_o"
  }

  make_formula <- function(continuous_features_local,
                           categorical_features_local,
                           continuous_features_spline_local,
                           calendar_period_local,
                           calendar_period_extrapolation_local,
                           input_output_local) {

    cat_terms <- NULL

    if (!is.null(categorical_features_local) &&
        length(categorical_features_local) > 0L) {
      cat_terms <- categorical_features_local
    }

    cont_terms <- NULL
    spline_terms <- NULL

    if (!is.null(continuous_features_local) &&
        length(continuous_features_local) > 0L) {

      spline_features <- intersect(
        continuous_features_local,
        continuous_features_spline_local
      )

      linear_features <- setdiff(
        continuous_features_local,
        spline_features
      )

      if (length(linear_features) > 0L) {
        cont_terms <- linear_features
      }

      if (length(spline_features) > 0L) {
        spline_terms <- paste0(
          "pspline(",
          spline_features,
          ",degree=",
          degrees_cf,
          ",df=",
          degrees_of_freedom_cf,
          ")"
        )
      }
    }

    calendar_term <- NULL

    if (isTRUE(calendar_period_extrapolation_local)) {
      calendar_term <- paste0(
        "pspline(",
        calendar_period_local,
        ",degree=",
        degrees_cp,
        ",df=",
        degrees_of_freedom_cp,
        ")"
      )
    }

    rhs_terms <- c(
      cat_terms,
      cont_terms,
      spline_terms,
      calendar_term
    )

    rhs_terms <- rhs_terms[
      !is.na(rhs_terms) &
        nzchar(rhs_terms)
    ]

    lhs <- paste0(
      "survival::Surv(TR_",
      input_output_local,
      ", DP_rev_",
      input_output_local,
      ", I) ~ "
    )

    if (length(rhs_terms) == 0L) {
      paste0(lhs, "1")
    } else {
      paste0(lhs, paste(rhs_terms, collapse = "+"))
    }
  }

  string_formula_i <- make_formula(
    continuous_features_local = continuous_features_i,
    categorical_features_local = categorical_features,
    continuous_features_spline_local = continuous_features_spline,
    calendar_period_local = "RP_i",
    calendar_period_extrapolation_local = calendar_period_extrapolation,
    input_output_local = "i"
  )

  string_formula_o <- make_formula(
    continuous_features_local = continuous_features_o,
    categorical_features_local = categorical_features,
    continuous_features_spline_local = continuous_features_spline,
    calendar_period_local = "RP_o",
    calendar_period_extrapolation_local = calendar_period_extrapolation,
    input_output_local = "o"
  )

  ## ------------------------------------------------------------------
  ## Output
  ## ------------------------------------------------------------------

  out <- list(
    training.data = train,
    full.data = tmp,
    data_information = list(
      conversion_factor = conversion_factor,
      string_formula_i = string_formula_i,
      string_formula_o = string_formula_o,
      continuous_features = continuous_features,
      categorical_features = categorical_features,
      calendar_period_extrapolation = calendar_period_extrapolation,
      years = years,
      accident_period = accident_period,
      calendar_period = calendar_period,
      input_time_granularity = input_time_granularity,
      output_time_granularity = output_time_granularity,
      input_time_unit = input_time_granularity,
      output_time_unit = output_time_granularity
    )
  )

  class(out) <- "IndividualDataPP"

  out
}





