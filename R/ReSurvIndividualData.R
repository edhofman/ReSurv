#' Fit \code{ReSurv} models on the individual data.
#'
#' This function fits and computes the reserves for the \code{ReSurv} models
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#' correspondence between hazard models and development factors:
#'
#' To be completed with final notation of the paper.
#'
#' The \code{ReSurv} package assumes proportional hazard models.
#' Given an i.i.d. sample \eqn{\left\{y_i,x_i\right\}_{i=1, \ldots, n}} the individual hazard at time \eqn{t} is:
#'
#' \eqn{\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}}
#'
#' Composed of a baseline \eqn{\lambda_0(t)} and a proportional effect \eqn{e^{y_i(x_i)}}.
#'
#' Currently, the implementation allows to optimize the partial likelihood (concerning the proportional effects) using one of the following statistical learning approaches:
#' \itemize{
#' \item{\href{https://github.com/therneau/survival}{COX}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' }
#'
#'
#' @param IndividualDataPP IndividualDataPP object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"COX"}: Standard Cox model for the hazard.}
#' \item{\code{"NN"}: Deep Survival Neural Network.}
#' \item{\code{"XGB"}: eXtreme Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed \code{integer}, random seed set for reproducibility
#' @param hparameters \code{list}, hyperparameters for the machine learning models. It will be disregarded for the cox approach.
#' @param percentage_data_training \code{numeric}, percentage of data used for training on the upper triangle.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#' @param eta \code{numeric}, Efron baseline and development-factor eta parameter.
#' @param simplifier \code{logical}, kept for compatibility. The simplified forecast frame is always used.
#'
#'
#' @return \code{ReSurv} fit. A list containing
#' \itemize{
#' \item{\code{model.out}: \code{list} containing the pre-processed covariates data for the fit (\code{data}) and the basic model output (\code{model.out};COX, XGB or NN).}
#' \item{\code{is_lkh}: \code{numeric} Training negative log likelihood.}
#' \item{\code{os_lkh}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{hazard_frame}: \code{data.frame} containing the fitted hazard model with the corresponding covariates. It contains:}
#'    \itemize{
#'    \item{\code{expg}: fitted risk score.}
#'    \item{\code{baseline}: fitted baseline.}
#'    \item{\code{hazard}: fitted hazard rate (\code{expg}*\code{baseline}).}
#'    \item{\code{f_i}: fitted development factors.}
#'    \item{\code{cum_f_i}: fitted cumulative development factors.}
#'    \item{\code{S_i}:fitted survival function.}
#'    \item{\code{S_i_lag}:fitted survival function (lag version, for further information see \code{?dplyr::lag}).}
#'    \item{\code{S_i_lead}:fitted survival function (lead version, for further information see \code{?dplyr::lead}).}
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
#'
#'
#'
#'@examples
#'
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 4,
#' period_exposure = 100)
#'
#' individual_data <- IndividualDataPP(data = input_data_0,
#' categorical_features = "claim_type",
#' continuous_features = "AP",
#' accident_period = "AP",
#' calendar_period = "RP",
#' input_time_granularity = "years",
#' output_time_granularity = "years",
#' years=4)
#'
#'
#' resurv_fit_cox <- ReSurv(individual_data,
#' hazard_model = "COX",
#' eta = 0)
#'
#'
#'
#'

#' @import xgboost
#' @import data.table
#' @importFrom dplyr reframe full_join
#' @importFrom tidyr replace_na
#'
#' @references
#' Munir, H., Emil, H., & Gabriele, P. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package â€˜survivalâ€™. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package â€˜xgboostâ€™. R version, 90, 1-66.
#'
#' @export
ReSurv <- function(IndividualDataPP,
                   hazard_model = "COX",
                   tie = "efron",
                   baseline = "spline",
                   continuous_features_scaling_method = "minmax",
                   random_seed = 1,
                   hparameters = list(),
                   percentage_data_training = .8,
                   grouping_method = "exposure",
                   check_value = 1.85,
                   eta=0.5,
                   simplifier=TRUE){

  UseMethod("ReSurv")

}
#' Fit \code{ReSurv} models on the individual data.
#'
#' This function fits and computes the reserves for the \code{ReSurv} models
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#' correspondence between hazard models and development factors:
#'
#' To be completed with final notation of the paper.
#'
#' The \code{ReSurv} package assumes proportional hazard models.
#' Given an i.i.d. sample \eqn{\left\{y_i,x_i\right\}_{i=1, \ldots, n}} the individual hazard at time \eqn{t} is:
#'
#' \eqn{\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}}
#'
#' Composed of a baseline \eqn{\lambda_0(t)} and a proportional effect \eqn{e^{y_i(x_i)}}.
#'
#' Currently, the implementation allows to optimize the partial likelihood (concerning the proportional effects) using one of the following statistical learning approaches:
#' \itemize{
#' \item{\href{https://github.com/therneau/survival}{COX}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' }
#'
#'
#' @param IndividualDataPP IndividualDataPP object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"COX"}: Standard Cox model for the hazard.}
#' \item{\code{"NN"}: Deep Survival Neural Network.}
#' \item{\code{"XGB"}: eXtreme Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed \code{integer}, random seed set for reproducibility
#' @param hparameters \code{list}, hyperparameters for the machine learning models. It will be disregarded for the cox approach.
#' @param percentage_data_training \code{numeric}, percentage of data used for training on the upper triangle.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#' @param eta \code{numeric}, Efron baseline and development-factor eta parameter.
#' @param simplifier \code{logical}, kept for compatibility. The simplified forecast frame is always used.
#'
#'
#' @return \code{ReSurv} fit. A list containing
#' \itemize{
#' \item{\code{model.out}: \code{list} containing the pre-processed covariates data for the fit (\code{data}) and the basic model output (\code{model.out};COX, XGB or NN).}
#' \item{\code{is_lkh}: \code{numeric} Training negative log likelihood.}
#' \item{\code{os_lkh}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{hazard_frame}: \code{data.frame} containing the fitted hazard model with the corresponding covariates. It contains:}
#'    \itemize{
#'    \item{\code{expg}: fitted risk score.}
#'    \item{\code{baseline}: fitted baseline.}
#'    \item{\code{hazard}: fitted hazard rate (\code{expg}*\code{baseline}).}
#'    \item{\code{f_i}: fitted development factors.}
#'    \item{\code{cum_f_i}: fitted cumulative development factors.}
#'    \item{\code{S_i}:fitted survival function.}
#'    \item{\code{S_i_lag}:fitted survival function (lag version, for further information see \code{?dplyr::lag}).}
#'    \item{\code{S_i_lead}:fitted survival function (lead version, for further information see \code{?dplyr::lead}).}
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
#'


#' @import xgboost

#'
#'
#'
#'
#'@examples
#'
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 4,
#' period_exposure = 100)
#'
#' individual_data <- IndividualDataPP(data = input_data_0,
#' categorical_features = "claim_type",
#' continuous_features = "AP",
#' accident_period = "AP",
#' calendar_period = "RP",
#' input_time_granularity = "years",
#' output_time_granularity = "years",
#' years=4)
#'
#'
#' resurv_fit_cox <- ReSurv(individual_data,
#' hazard_model = "COX",
#' eta = 0)
#'
#'
#'
#'
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package â€˜survivalâ€™. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package â€˜xgboostâ€™. R version, 90, 1-66.
#'
#' @export
ReSurv.default <- function(IndividualDataPP,
                           hazard_model = "COX",
                           tie = "efron",
                           baseline = "spline",
                           continuous_features_scaling_method = "minmax",
                           random_seed = 1,
                           hparameters = list(),
                           percentage_data_training = .8,
                           grouping_method = "exposure",
                           check_value = 1.85,
                           eta=0.5,
                           simplifier=TRUE){

  message('The object provided must be of class IndividualDataPP')

}



#' Fit \code{ReSurv} models on the individual data.
#'
#' This function fits and computes the reserves for the \code{ReSurv} models
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#' correspondence between hazard models and development factors:
#'
#' To be completed with final notation of the paper.
#'
#' The \code{ReSurv} package assumes proportional hazard models.
#' Given an i.i.d. sample \eqn{\left\{y_i,x_i\right\}_{i=1, \ldots, n}} the individual hazard at time \eqn{t} is:
#'
#' \eqn{\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}}
#'
#' Composed of a baseline \eqn{\lambda_0(t)} and a proportional effect \eqn{e^{y_i(x_i)}}.
#'
#' Currently, the implementation allows to optimize the partial likelihood (concerning the proportional effects) using one of the following statistical learning approaches:
#' \itemize{
#' \item{\href{https://github.com/therneau/survival}{COX}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' }
#'
#'
#' @param IndividualDataPP IndividualDataPP object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"COX"}: Standard Cox model for the hazard.}
#' \item{\code{"NN"}: Deep Survival Neural Network.}
#' \item{\code{"XGB"}: eXtreme Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed \code{integer}, random seed set for reproducibility
#' @param hparameters \code{list}, hyperparameters for the machine learning models. It will be disregarded for the cox approach.
#' @param percentage_data_training \code{numeric}, percentage of data used for training on the upper triangle.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#' @param eta \code{numeric}, Efron baseline and development-factor eta parameter.
#' @param simplifier \code{logical}, kept for compatibility. The simplified forecast frame is always used.
#'
#' @return \code{ReSurv} fit. A list containing
#' \itemize{
#' \item{\code{model.out}: \code{list} containing the pre-processed covariates data for the fit (\code{data}) and the basic model output (\code{model.out};COX, XGB or NN).}
#' \item{\code{is_lkh}: \code{numeric} Training negative log likelihood.}
#' \item{\code{os_lkh}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{hazard_frame}: \code{data.frame} containing the fitted hazard model with the corresponding covariates. It contains:}
#'    \itemize{
#'    \item{\code{expg}: fitted risk score.}
#'    \item{\code{baseline}: fitted baseline.}
#'    \item{\code{hazard}: fitted hazard rate (\code{expg}*\code{baseline}).}
#'    \item{\code{f_i}: fitted development factors.}
#'    \item{\code{cum_f_i}: fitted cumulative development factors.}
#'    \item{\code{S_i}:fitted survival function.}
#'    \item{\code{S_i_lag}:fitted survival function (lag version, for further information see \code{?dplyr::lag}).}
#'    \item{\code{S_i_lead}:fitted survival function (lead version, for further information see \code{?dplyr::lead}).}
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
#'


#' @import xgboost

#'
#'
#'
#'
#'@examples
#'
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 4,
#' period_exposure = 100)
#'
#' individual_data <- IndividualDataPP(data = input_data_0,
#' categorical_features = "claim_type",
#' continuous_features = "AP",
#' accident_period = "AP",
#' calendar_period = "RP",
#' input_time_granularity = "years",
#' output_time_granularity = "years",
#' years=4)
#'
#'
#' resurv_fit_cox <- ReSurv(individual_data,
#' hazard_model = "COX",
#' eta = 0)
#'
#'
#'
#'
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package â€˜survivalâ€™. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package â€˜xgboostâ€™. R version, 90, 1-66.
#'
#' @export
ReSurv.IndividualDataPP <- function(IndividualDataPP,
                                  hazard_model = "COX",
                                  tie = "efron",
                                  baseline = "spline",
                                  continuous_features_scaling_method = "minmax",
                                  random_seed = 1,
                                  hparameters = list(),
                                  percentage_data_training = .8,
                                  grouping_method = "exposure",
                                  check_value = 1.85,
                                  eta=0.5,
                                  simplifier=TRUE
){


  # validate eta
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta)) {
    stop("`eta` must be a single finite numeric value.", call. = FALSE)
  }

  if (eta < 0 || eta > 1) {
    stop("`eta` must lie in [0, 1].", call. = FALSE)
  }

  # covariates and formula extracted
  cont_f <- IndividualDataPP$data_information$continuous_features
  cat_f <- IndividualDataPP$data_information$categorical_features

  set.seed(random_seed)

  formula_ct <- as.formula(IndividualDataPP$data_information$string_formula_i)

  ## ------------------------------------------------------------------
  ## Inline simplified_df_2_fcst(), data.table-only
  ## ------------------------------------------------------------------

  train_dt <- data.table::as.data.table(IndividualDataPP$training.data)

  years <- IndividualDataPP$data_information$years
  input_time_granularity <- IndividualDataPP$data_information$input_time_granularity
  calendar_period_extrapolation <- IndividualDataPP$data_information$calendar_period_extrapolation

  time_features <- c("DP_i", "DP_rev_i", "RP_i")

  columns_for_grouping <- unique(c(
    setdiff(cont_f, time_features),
    setdiff(cat_f, time_features),
    "AP_i"
  ))

  columns_for_grouping <- intersect(columns_for_grouping, names(train_dt))

  out <- train_dt[
    ,
    .(.N),
    by = columns_for_grouping
  ][
    ,
    .SD,
    .SDcols = columns_for_grouping
  ]

  l4 <- data.table::CJ(
    DP_rev_i = min(train_dt[["DP_rev_i"]], na.rm = TRUE):
      max(train_dt[["DP_rev_i"]], na.rm = TRUE),
    sorted = FALSE
  )

  out <- data.table::setkey(
    out[, c(k = 1, .SD)],
    k
  )[
    l4[, c(k = 1, .SD)],
    allow.cartesian = TRUE
  ][
    ,
    k := NULL
  ]

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

  max_dp_i <- as.integer(years / time_unit_numeric[input_pos])

  out[
    ,
    DP_i := max_dp_i - DP_rev_i + 1L
  ]

  if (isTRUE(calendar_period_extrapolation) || "RP_i" %in% c(cont_f, cat_f)) {
    out[
      ,
      RP_i := AP_i + DP_i - 1L
    ]
  }

  if (!is.null(cat_f)) {
    time_cat_f <- intersect(cat_f, time_features)

    for (cft in time_cat_f) {
      if (cft %in% names(out) && cft %in% names(train_dt)) {
        out[
          ,
          (cft) := factor(
            get(cft),
            levels = levels(train_dt[[cft]])
          )
        ]
      }
    }
  }

  newdata <- as.data.frame(out)

  # logical: check if we work with a baseline model
  is_baseline_model = is.null(c(cont_f,
                                cat_f))


  # Proportional hazard model fitting -------

  data <- IndividualDataPP$training.data

  Y <- data[, .SD, .SDcols = c("DP_rev_i", "I", "TR_i")]

  if (!(hazard_model %in% c("COX", "NN", "XGB"))) {
    stop("`hazard_model` must be one of 'COX', 'NN', or 'XGB'.", call. = FALSE)
  }

  if (hazard_model %in% c("NN", "XGB")) {
    if (percentage_data_training > 1 || percentage_data_training < 0) {
      warning(
        paste0(
          "Traintestsplit has been put to ",
          percentage_data_training,
          ". The value needs to be between 0 and 1, defaulting to 0.8."
        )
      )
      training_test_split <- .8
    } else {
      training_test_split <- percentage_data_training
    }
  }

  X <- NULL
  X_base <- NULL
  datads_pp <- NULL
  pred <- NULL

  ## ------------------------------------------------------------
  ## 1. Model-specific fitting and prediction
  ## ------------------------------------------------------------

  if (hazard_model == "COX") {

    cox <- survival::coxph(formula_ct, data = data, ties = tie)

    cox_lp <- predict(
      cox,
      newdata = newdata,
      type = "lp",
      reference = "zero"
    )

    cox_training_lp <- predict(
      cox,
      newdata = data[order(DP_rev_i)],
      type = "lp",
      reference = "zero"
    )

    model.out <- list(
      cox = cox,
      cox_lp = cox_lp,
      expg = exp(cox_lp),
      train_expg = cox_training_lp
    )

    pred <- cox_lp
  }

  if (hazard_model %in% c("NN", "XGB")) {

    if (is_baseline_model) {

      X <- data.frame(intercept_1 = rep(1, nrow(Y)))

    } else {

      X_parts <- list()

      if (!is.null(cat_f)) {

        X_cat <- fastDummies::dummy_cols(
          data,
          select_columns = cat_f,
          remove_selected_columns = TRUE,
          remove_first_dummy = hazard_model == "XGB"
        )

        keep_cat_cols <- vapply(
          colnames(X_cat),
          function(z) {
            any(vapply(cat_f, function(p) grepl(pattern = p, x = z), logical(1)))
          },
          logical(1)
        )

        X_cat <- data.table::as.data.table(X_cat)
        X_parts[["categorical"]] <- X_cat[, .SD, .SDcols = colnames(X_cat)[keep_cat_cols]]
      }

      if (!is.null(cont_f)) {

        X_cont <- data.table::copy(data.table::as.data.table(data)[, .SD, .SDcols = cont_f])

        for (cft in cont_f) {
          if (continuous_features_scaling_method == "minmax") {
            X_cont[[cft]] <- 2 * (X_cont[[cft]] - min(X_cont[[cft]])) /
              (max(X_cont[[cft]]) - min(X_cont[[cft]])) - 1
          }

          if (continuous_features_scaling_method == "standard") {
            X_cont[[cft]] <- (X_cont[[cft]] - mean(X_cont[[cft]])) / sd(X_cont[[cft]])
          }
        }

        X_parts[["continuous"]] <- X_cont
      }

      X <- as.data.frame(do.call(cbind, X_parts))
    }
  }

  if (hazard_model == "NN") {

    ## Inline deep_surv_pp()

    tmp_order <- order(Y$DP_rev_i)

    X_ordered <- as.data.frame(X[tmp_order, , drop = FALSE])
    Y_ordered <- as.data.frame(Y[tmp_order, , drop = FALSE])

    tmp_ids <- data.frame(id = seq_len(nrow(X_ordered)))

    if (!is.numeric(training_test_split) ||
        length(training_test_split) != 1L ||
        !is.finite(training_test_split) ||
        training_test_split <= 0 ||
        training_test_split > 1) {
      stop("`training_test_split` must be a number in (0, 1].", call. = FALSE)
    }

    if (training_test_split == 1) {
      sampled_id <- tmp_ids$id
    } else {
      n_sample <- ceiling(nrow(tmp_ids) * training_test_split)
      n_sample <- max(1L, min(nrow(tmp_ids), n_sample))
      sampled_id <- sample(tmp_ids$id, size = n_sample, replace = FALSE)
    }

    id_train <- tmp_ids$id %in% sampled_id

    datads_pp <- list(
      x_train = as.matrix(X_ordered[id_train, , drop = FALSE]),
      y_train = as.matrix(Y_ordered[id_train, , drop = FALSE]),
      x_val = as.matrix(X_ordered[!id_train, , drop = FALSE]),
      y_val = as.matrix(Y_ordered[!id_train, , drop = FALSE]),
      lkh_eval_data = list(
        data_train = X_ordered[id_train, , drop = FALSE],
        data_val = X_ordered[!id_train, , drop = FALSE],
        y_train = Y_ordered[id_train, , drop = FALSE],
        y_val = Y_ordered[!id_train, , drop = FALSE]
      )
    )

    ## Inline nn_hparameter_nodes_grid()

    if ("num_layers" %in% names(hparameters)) {
      if (hparameters$num_layers == length(hparameters$num_nodes)) {
        for (i in seq_len(hparameters$num_layers)) {
          hparameters[[paste0("node_", i)]] <- hparameters$num_nodes[i]
        }
      } else if (length(hparameters$num_nodes) == 1L) {
        for (i in seq_len(hparameters$num_layers)) {
          hparameters[[paste0("node_", i)]] <- hparameters$num_nodes
        }
      } else {
        warning(
          "`num_nodes` was not supplied correctly. Using the first value for all layers.",
          call. = FALSE
        )

        for (i in seq_len(hparameters$num_layers)) {
          hparameters[[paste0("node_", i)]] <- hparameters$num_nodes[1]
        }
      }

      hparameters[["num_nodes"]] <- NULL
    }

    hparameters <- list(
      params = as.list.data.frame(hparameters),
      verbose = hparameters$verbose,
      epochs = hparameters$epochs,
      num_workers = hparameters$num_workers
    )

    ## Inline fit_deep_surv()

    input_dim <- ncol(datads_pp$x_train)

    activation_map <- list(
      "relu" = torch::nn_relu,
      "selu" = torch::nn_selu,
      "tanh" = torch::nn_tanh,
      "leakyrelu" = torch::nn_leaky_relu,
      "leaky_relu" = torch::nn_leaky_relu
    )

    activation <- tolower(as.character(hparameters$params$activation))
    act_fn <- activation_map[[activation]]

    if (is.null(act_fn)) {
      stop("Unknown activation: ", hparameters$params$activation, call. = FALSE)
    }

    layers <- list()
    in_features <- input_dim

    for (i in seq_len(hparameters$params$num_layers)) {
      node_i <- as.integer(hparameters$params[[paste0("node_", i)]])
      layers <- c(
        layers,
        list(
          torch::nn_linear(in_features, node_i),
          act_fn()
        )
      )
      in_features <- node_i
    }

    layers <- c(
      layers,
      list(torch::nn_linear(in_features, 1L, bias = FALSE))
    )

    net <- do.call(torch::nn_sequential, layers)

    torch::torch_manual_seed(random_seed)

    optimizer <- switch(
      hparameters$params$optim,
      "Adam" = torch::optim_adam(net$parameters, lr = hparameters$params$lr),
      "SGD" = torch::optim_sgd(net$parameters, lr = hparameters$params$lr),
      "AdamW" = if (exists("optim_adamw", envir = asNamespace("torch"))) {
        torch::optim_adamw(net$parameters, lr = hparameters$params$lr)
      } else {
        torch::optim_adam(net$parameters, lr = hparameters$params$lr)
      },
      stop("Unknown optimizer: ", hparameters$params$optim, call. = FALSE)
    )

    x_train_t <- torch::torch_tensor(
      as.matrix(datads_pp$x_train),
      dtype = torch::torch_float32()
    )

    x_val_t <- torch::torch_tensor(
      as.matrix(datads_pp$x_val),
      dtype = torch::torch_float32()
    )

    y_train_nn <- as.matrix(datads_pp$y_train)
    y_val_nn <- as.matrix(datads_pp$y_val)

    dur_train_t <- torch::torch_tensor(y_train_nn[, 1], dtype = torch::torch_float32())
    event_train_t <- torch::torch_tensor(y_train_nn[, 2], dtype = torch::torch_float32())
    trunc_train_t <- torch::torch_tensor(y_train_nn[, 3], dtype = torch::torch_float32())

    dur_val_t <- torch::torch_tensor(y_val_nn[, 1], dtype = torch::torch_float32())
    event_val_t <- torch::torch_tensor(y_val_nn[, 2], dtype = torch::torch_float32())
    trunc_val_t <- torch::torch_tensor(y_val_nn[, 3], dtype = torch::torch_float32())

    train_losses <- numeric(0)
    val_losses <- numeric(0)

    best_val <- Inf
    wait <- 0L

    for (epoch in seq_len(hparameters$epochs)) {

      net$train()
      optimizer$zero_grad()

      train_log_h <- net(x_train_t)$squeeze()

      cox_loss <- (function(log_h, durations, events, truncation) {
        log_h <- log_h$view(c(-1))
        durations <- durations$view(c(-1))
        events <- events$view(c(-1))
        truncation <- truncation$view(c(-1))

        durations_r <- as.numeric(durations)
        events_r <- as.numeric(events)

        event_times <- sort(unique(durations_r[events_r == 1]))

        if (length(event_times) == 0L) {
          stop("No events found in Cox torch loss.", call. = FALSE)
        }

        exp_log_h <- torch::torch_exp(log_h)
        loss <- log_h$sum() * 0
        n_events <- sum(events_r == 1)

        for (t_j in event_times) {
          event_ind <- ((durations == t_j) * (events == 1))$to(dtype = torch::torch_float32())
          risk_ind <- ((durations >= t_j) * (truncation < t_j))$to(dtype = torch::torch_float32())

          d_j <- sum(durations_r == t_j & events_r == 1)

          if (d_j <= 0L) {
            next
          }

          risk_sum <- torch::torch_sum(risk_ind * exp_log_h)
          event_sum_exp <- torch::torch_sum(event_ind * exp_log_h)
          event_sum_log <- torch::torch_sum(event_ind * log_h)

          efron_fraction <- torch::torch_tensor(
            seq(0, d_j - 1) / d_j,
            dtype = torch::torch_float32()
          )

          denominators <- risk_sum - efron_fraction * event_sum_exp

          if (any(as.numeric(denominators) <= 0)) {
            stop(
              "Non-positive denominator in Efron partial likelihood. ",
              "Check risk sets, truncation times, and event times.",
              call. = FALSE
            )
          }

          loss <- loss +
            torch::torch_sum(torch::torch_log(denominators)) -
            event_sum_log
        }

        loss / n_events
      })(train_log_h, dur_train_t, event_train_t, trunc_train_t)

      l2 <- torch::torch_tensor(0, dtype = torch::torch_float32())
      l1 <- torch::torch_tensor(0, dtype = torch::torch_float32())

      for (p in net$parameters) {
        if (length(p$shape) >= 2L) {
          l2 <- l2 + torch::torch_sum(p^2)
          l1 <- l1 + torch::torch_sum(torch::torch_abs(p))
        }
      }

      reg <- hparameters$params$eps *
        (hparameters$params$xi * l2 + (1 - hparameters$params$xi) * l1)

      total_loss <- cox_loss + reg

      total_loss$backward()
      optimizer$step()

      net$eval()

      torch::with_no_grad({

        train_log_h_eval <- net(x_train_t)$squeeze()

        t_loss <- (function(log_h, durations, events, truncation) {
          log_h <- log_h$view(c(-1))
          durations <- durations$view(c(-1))
          events <- events$view(c(-1))
          truncation <- truncation$view(c(-1))

          durations_r <- as.numeric(durations)
          events_r <- as.numeric(events)
          event_times <- sort(unique(durations_r[events_r == 1]))

          if (length(event_times) == 0L) {
            stop("No events found in Cox torch loss.", call. = FALSE)
          }

          exp_log_h <- torch::torch_exp(log_h)
          loss <- log_h$sum() * 0
          n_events <- sum(events_r == 1)

          for (t_j in event_times) {
            event_ind <- ((durations == t_j) * (events == 1))$to(dtype = torch::torch_float32())
            risk_ind <- ((durations >= t_j) * (truncation < t_j))$to(dtype = torch::torch_float32())

            d_j <- sum(durations_r == t_j & events_r == 1)

            risk_sum <- torch::torch_sum(risk_ind * exp_log_h)
            event_sum_exp <- torch::torch_sum(event_ind * exp_log_h)
            event_sum_log <- torch::torch_sum(event_ind * log_h)

            efron_fraction <- torch::torch_tensor(
              seq(0, d_j - 1) / d_j,
              dtype = torch::torch_float32()
            )

            denominators <- risk_sum - efron_fraction * event_sum_exp

            if (any(as.numeric(denominators) <= 0)) {
              stop(
                "Non-positive denominator in Efron partial likelihood.",
                call. = FALSE
              )
            }

            loss <- loss +
              torch::torch_sum(torch::torch_log(denominators)) -
              event_sum_log
          }

          loss / n_events
        })(train_log_h_eval, dur_train_t, event_train_t, trunc_train_t)

        val_log_h <- net(x_val_t)$squeeze()

        v_loss <- (function(log_h, durations, events, truncation) {
          log_h <- log_h$view(c(-1))
          durations <- durations$view(c(-1))
          events <- events$view(c(-1))
          truncation <- truncation$view(c(-1))

          durations_r <- as.numeric(durations)
          events_r <- as.numeric(events)
          event_times <- sort(unique(durations_r[events_r == 1]))

          if (length(event_times) == 0L) {
            stop("No events found in Cox torch loss.", call. = FALSE)
          }

          exp_log_h <- torch::torch_exp(log_h)
          loss <- log_h$sum() * 0
          n_events <- sum(events_r == 1)

          for (t_j in event_times) {
            event_ind <- ((durations == t_j) * (events == 1))$to(dtype = torch::torch_float32())
            risk_ind <- ((durations >= t_j) * (truncation < t_j))$to(dtype = torch::torch_float32())

            d_j <- sum(durations_r == t_j & events_r == 1)

            risk_sum <- torch::torch_sum(risk_ind * exp_log_h)
            event_sum_exp <- torch::torch_sum(event_ind * exp_log_h)
            event_sum_log <- torch::torch_sum(event_ind * log_h)

            efron_fraction <- torch::torch_tensor(
              seq(0, d_j - 1) / d_j,
              dtype = torch::torch_float32()
            )

            denominators <- risk_sum - efron_fraction * event_sum_exp

            if (any(as.numeric(denominators) <= 0)) {
              stop(
                "Non-positive denominator in Efron partial likelihood.",
                call. = FALSE
              )
            }

            loss <- loss +
              torch::torch_sum(torch::torch_log(denominators)) -
              event_sum_log
          }

          loss / n_events
        })(val_log_h, dur_val_t, event_val_t, trunc_val_t)
      })

      t_loss_val <- t_loss$item()
      v_loss_val <- v_loss$item()

      train_losses <- c(train_losses, t_loss_val)
      val_losses <- c(val_losses, v_loss_val)

      if (isTRUE(hparameters$params$early_stopping)) {
        if (v_loss_val < best_val) {
          best_val <- v_loss_val
          wait <- 0L
        } else {
          wait <- wait + 1L

          if (wait >= hparameters$params$patience) {
            if (hparameters$verbose) {
              message("Early stopping at epoch ", epoch)
            }
            break
          }
        }
      }

      if (hparameters$verbose) {
        message(sprintf(
          "Epoch %d | train_loss: %.6f | val_loss: %.6f",
          epoch,
          t_loss_val,
          v_loss_val
        ))
      }
    }

    model.out <- list(
      net = net,
      log = data.frame(
        train_loss = train_losses,
        val_loss = val_losses
      )
    )

    ## Forecast prediction

    if (is_baseline_model) {

      newdata.mx <- data.frame(intercept_1 = rep(1, nrow(newdata)))

    } else {

      X_parts_new <- list()

      if (!is.null(cat_f)) {
        X_cat_new <- fastDummies::dummy_cols(
          newdata,
          select_columns = cat_f,
          remove_selected_columns = TRUE,
          remove_first_dummy = FALSE
        )

        keep_cat_cols_new <- vapply(
          colnames(X_cat_new),
          function(z) {
            any(vapply(cat_f, function(p) grepl(pattern = p, x = z), logical(1)))
          },
          logical(1)
        )

        X_cat_new <- data.table::as.data.table(X_cat_new)
        X_parts_new[["categorical"]] <- X_cat_new[
          ,
          .SD,
          .SDcols = colnames(X_cat_new)[keep_cat_cols_new]
        ]
      }

      if (!is.null(cont_f)) {
        X_cont_new <- data.table::copy(
          data.table::as.data.table(newdata)[, .SD, .SDcols = cont_f]
        )

        for (cft in cont_f) {
          mnv <- min(data[[cft]])
          mxv <- max(data[[cft]])
          X_cont_new[[cft]] <- 2 * (X_cont_new[[cft]] - mnv) / (mxv - mnv) - 1
        }

        X_parts_new[["continuous"]] <- X_cont_new
      }

      newdata.mx <- do.call(cbind, X_parts_new)
    }

    net$eval()
    x_fc <- as.matrix(newdata.mx)

    x_tensor <- torch::torch_tensor(
      as.matrix(x_fc),
      dtype = torch::torch_float32()
    )

    pred <- as.numeric(torch::with_no_grad({
      net(x_tensor)
    }))
  }

  if (hazard_model == "XGB") {

    ## Inline xgboost_pp()

    xy <- data.frame(X, Y, check.names = FALSE)
    tmp <- xy[order(xy$DP_rev_i), , drop = FALSE]
    tmp$id <- seq_len(nrow(tmp))

    if (!is.numeric(training_test_split) ||
        length(training_test_split) != 1L ||
        !is.finite(training_test_split) ||
        training_test_split <= 0 ||
        training_test_split > 1) {
      stop("`training_test_split` must be a number in (0, 1].", call. = FALSE)
    }

    if (training_test_split == 1) {
      sampled_id <- tmp$id
    } else {
      n_sample <- ceiling(nrow(tmp) * training_test_split)
      n_sample <- max(1L, min(nrow(tmp), n_sample))
      sampled_id <- sample(tmp$id, size = n_sample, replace = FALSE)
    }

    samples_cn <- data.frame(id = sampled_id)

    tmp_train <- tmp[tmp$id %in% samples_cn$id, , drop = FALSE]
    tmp_train <- tmp_train[order(tmp_train$DP_rev_i), , drop = FALSE]

    tmp_train$efron_c <- ave(
      seq_along(tmp_train$DP_rev_i),
      tmp_train$DP_rev_i,
      FUN = function(ind) (seq_along(ind) - 1) / length(ind)
    )

    ds_train_m <- xgboost::xgb.DMatrix(
      data = as.matrix(tmp_train[, colnames(X), drop = FALSE]),
      label = tmp_train$I
    )

    event_times <- unique(tmp_train$DP_rev_i)

    risk_sets <- vector("list", length(event_times))
    event_sets <- vector("list", length(event_times))

    for (i in seq_along(event_times)) {
      stop_i <- event_times[i]
      risk_sets[[i]] <- which((tmp_train$TR_i < stop_i) & (tmp_train$DP_rev_i >= stop_i))
      event_sets[[i]] <- which(tmp_train$DP_rev_i == stop_i)
    }

    attr(ds_train_m, "truncation") <- tmp_train$TR_i
    attr(ds_train_m, "claim_arrival") <- tmp_train$DP_rev_i
    attr(ds_train_m, "risk_sets") <- risk_sets
    attr(ds_train_m, "event_sets") <- event_sets
    attr(ds_train_m, "efron_c") <- tmp_train$efron_c

    tie_table <- table(tmp_train$DP_rev_i)

    attr(ds_train_m, "tieid") <- unname(tie_table)
    attr(ds_train_m, "groups") <- rep(
      as.integer(names(tie_table)),
      unname(tie_table)
    )

    if (training_test_split < 1) {

      tmp_test <- tmp[!(tmp$id %in% samples_cn$id), , drop = FALSE]
      tmp_test <- tmp_test[order(tmp_test$DP_rev_i), , drop = FALSE]

      tmp_test$efron_c <- ave(
        seq_along(tmp_test$DP_rev_i),
        tmp_test$DP_rev_i,
        FUN = function(ind) (seq_along(ind) - 1) / length(ind)
      )

      ds_test_m <- xgboost::xgb.DMatrix(
        data = as.matrix(tmp_test[, colnames(X), drop = FALSE]),
        label = tmp_test$I
      )

      event_times_test <- unique(tmp_test$DP_rev_i)

      risk_sets_test <- vector("list", length(event_times_test))
      event_sets_test <- vector("list", length(event_times_test))

      for (i in seq_along(event_times_test)) {
        stop_i <- event_times_test[i]
        risk_sets_test[[i]] <- which((tmp_test$TR_i < stop_i) & (tmp_test$DP_rev_i >= stop_i))
        event_sets_test[[i]] <- which(tmp_test$DP_rev_i == stop_i)
      }

      attr(ds_test_m, "truncation") <- tmp_test$TR_i
      attr(ds_test_m, "claim_arrival") <- tmp_test$DP_rev_i
      attr(ds_test_m, "risk_sets") <- risk_sets_test
      attr(ds_test_m, "event_sets") <- event_sets_test
      attr(ds_test_m, "efron_c") <- tmp_test$efron_c

      tie_table_test <- table(tmp_test$DP_rev_i)

      attr(ds_test_m, "tieid") <- unname(tie_table_test)
      attr(ds_test_m, "groups") <- rep(
        as.integer(names(tie_table_test)),
        unname(tie_table_test)
      )

    } else {

      ds_test_m <- NULL
    }

    datads_pp <- list(
      ds_train_m = ds_train_m,
      ds_test_m = ds_test_m,
      samples_cn = samples_cn
    )

    ## Inline fit_xgboost()

    if (length(hparameters) == 0L) {
      hparameters <- list(
        params = list(
          booster = "gbtree",
          eta = .01,
          subsample = .5,
          alpha = 1,
          lambda = 1,
          min_child_weight = .2
        ),
        print_every_n = NULL,
        nrounds = 10,
        verbose = FALSE,
        early_stopping_rounds = 500
      )
    }

    evals <- list(train = datads_pp$ds_train_m)

    if (!is.null(datads_pp$ds_test_m)) {
      evals$eval <- datads_pp$ds_test_m
    }

    early_stopping_rounds <- hparameters$early_stopping_rounds

    if (is.null(datads_pp$ds_test_m)) {
      early_stopping_rounds <- NULL
    }

    model.out <- xgboost::xgb.train(
      params = hparameters$params,
      data = datads_pp$ds_train_m,
      obj = function(preds, dtrain) {

        Ti <- attr(dtrain, "truncation")
        Ei <- attr(dtrain, "claim_arrival")

        risk_sets <- attr(dtrain, "risk_sets")
        event_sets <- attr(dtrain, "event_sets")

        efron_c <- attr(dtrain, "efron_c")
        tieid <- attr(dtrain, "tieid")

        exp_p_sum <- sapply(
          risk_sets,
          FUN = function(x, ypred) sum(exp(ypred[x])),
          ypred = preds
        )

        exp_p_tie <- sapply(
          event_sets,
          FUN = function(x, ypred) sum(exp(ypred[x])),
          ypred = preds
        )

        tmp1 <- data.table::data.table(
          risks_s = rep(exp_p_sum, tieid),
          events_s = rep(exp_p_tie, tieid),
          efron_c = efron_c,
          ties = Ei
        )

        tmp_alpha_i <- tmp1[
          ,
          .(alpha_i = sum(1 / (risks_s - efron_c * events_s))),
          by = ties
        ]$alpha_i

        alpha_i <- vector("numeric", length = max(Ei))
        alpha_i[unique(Ei)] <- tmp_alpha_i

        tmp_beta_i <- tmp1[
          ,
          .(beta_i = sum(efron_c / (risks_s - efron_c * events_s))),
          by = ties
        ]$beta_i

        beta_i <- vector("numeric", length = max(Ei))
        beta_i[unique(Ei)] <- tmp_beta_i

        tmp_gamma_i <- tmp1[
          ,
          .(gamma_i = sum(1 / (risks_s - efron_c * events_s)^2)),
          by = ties
        ]$gamma_i

        gamma_i <- vector("numeric", length = max(Ei))
        gamma_i[unique(Ei)] <- tmp_gamma_i

        tmp_omega_i <- tmp1[
          ,
          .(omega_i = sum((1 - (1 - efron_c)^2) / (risks_s - efron_c * events_s)^2)),
          by = ties
        ]$omega_i

        omega_i <- vector("numeric", length = max(Ei))
        omega_i[unique(Ei)] <- tmp_omega_i

        exp_p <- exp(preds)

        cumsum_alpha <- c(0, cumsum(alpha_i))
        cumsum_gamma <- c(0, cumsum(gamma_i))

        alpha_i_lt <- cumsum_alpha[Ei + 1] - cumsum_alpha[Ti + 1]
        gamma_i_lt <- cumsum_gamma[Ei + 1] - cumsum_gamma[Ti + 1]
        beta_i_lt <- beta_i[Ei]
        omega_i_lt <- omega_i[Ei]

        grad <- exp_p * (alpha_i_lt - beta_i_lt) - 1
        hess <- grad - (exp_p^2) * (gamma_i_lt - omega_i_lt) + 1

        list(grad = grad, hess = hess)
      },
      nrounds = hparameters$nrounds,
      custom_metric = function(preds, dtrain) {
        risk_sets <- attr(dtrain, "risk_sets")
        event_sets <- attr(dtrain, "event_sets")
        efron_c <- attr(dtrain, "efron_c")
        tieid <- attr(dtrain, "tieid")

        exp_p_sum <- rep(
          sapply(risk_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds),
          tieid
        )

        exp_p_tie <- rep(
          sapply(event_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds),
          tieid
        )

        exp_p <- exp(preds)

        r_k <- exp_p_sum - efron_c * exp_p_tie
        lkh <- exp_p / r_k

        value <- -sum(log(lkh))

        list(metric = "log-partial likelihood", value = value / length(preds))
      },
      evals = evals,
      verbose = hparameters$verbose,
      print_every_n = hparameters$print_every_n,
      early_stopping_rounds = early_stopping_rounds,
      maximize = FALSE
    )

    ## Forecast prediction

    if (is_baseline_model) {

      newdata.mx <- xgboost::xgb.DMatrix(
        as.matrix(data.frame(intercept_1 = rep(1, nrow(newdata))))
      )

    } else {

      X_parts_new <- list()

      if (!is.null(cat_f)) {
        X_cat_new <- fastDummies::dummy_cols(
          newdata,
          select_columns = cat_f,
          remove_selected_columns = TRUE,
          remove_first_dummy = TRUE
        )

        keep_cat_cols_new <- vapply(
          colnames(X_cat_new),
          function(z) {
            any(vapply(cat_f, function(p) grepl(pattern = p, x = z), logical(1)))
          },
          logical(1)
        )

        X_cat_new <- data.table::as.data.table(X_cat_new)
        X_parts_new[["categorical"]] <- X_cat_new[
          ,
          .SD,
          .SDcols = colnames(X_cat_new)[keep_cat_cols_new]
        ]
      }

      if (!is.null(cont_f)) {
        X_cont_new <- data.table::copy(
          data.table::as.data.table(newdata)[, .SD, .SDcols = cont_f]
        )

        for (cft in cont_f) {
          mnv <- min(data[[cft]])
          mxv <- max(data[[cft]])
          X_cont_new[[cft]] <- 2 * (X_cont_new[[cft]] - mnv) / (mxv - mnv) - 1
        }

        X_parts_new[["continuous"]] <- X_cont_new
      }

      X_new_xgb <- as.matrix(do.call(cbind, X_parts_new))

      newdata.mx <- xgboost::xgb.DMatrix(
        X_new_xgb,
        label = rep(1, nrow(X_new_xgb))
      )
    }

    pred <- predict(model.out, newdata.mx)
  }

  ## ------------------------------------------------------------
  ## 2. Baseline design matrix
  ## ------------------------------------------------------------

  if (is_baseline_model) {

    X_base <- data.frame(intercept_1 = rep(1, nrow(Y)))

  } else if (hazard_model %in% c("NN", "XGB")) {

    X_base <- X

  } else if (hazard_model == "COX") {

    X_parts_base <- list()

    if (!is.null(cat_f)) {

      X_cat_base <- fastDummies::dummy_cols(
        data,
        select_columns = cat_f,
        remove_selected_columns = TRUE,
        remove_first_dummy = TRUE
      )

      keep_cat_cols_base <- vapply(
        colnames(X_cat_base),
        function(z) {
          any(vapply(cat_f, function(p) grepl(pattern = p, x = z), logical(1)))
        },
        logical(1)
      )

      X_cat_base <- data.table::as.data.table(X_cat_base)
      X_parts_base[["categorical"]] <- X_cat_base[
        ,
        .SD,
        .SDcols = colnames(X_cat_base)[keep_cat_cols_base]
      ]
    }

    if (!is.null(cont_f)) {

      X_cont_base <- data.table::copy(
        data.table::as.data.table(data)[, .SD, .SDcols = cont_f]
      )

      for (cft in cont_f) {
        if (continuous_features_scaling_method == "minmax") {
          X_cont_base[[cft]] <- 2 * (X_cont_base[[cft]] - min(X_cont_base[[cft]])) /
            (max(X_cont_base[[cft]]) - min(X_cont_base[[cft]])) - 1
        }

        if (continuous_features_scaling_method == "standard") {
          X_cont_base[[cft]] <- (X_cont_base[[cft]] - mean(X_cont_base[[cft]])) /
            sd(X_cont_base[[cft]])
        }
      }

      X_parts_base[["continuous"]] <- X_cont_base
    }

    X_base <- as.data.frame(do.call(cbind, X_parts_base))
  }

  ## ------------------------------------------------------------
  ## 3. Baseline calculation
  ## ------------------------------------------------------------

  xy_bsln <- data.frame(X_base, Y, check.names = FALSE)
  tmp_bsln <- xy_bsln[order(xy_bsln$DP_rev_i), , drop = FALSE]
  tmp_bsln$id <- seq_len(nrow(tmp_bsln))

  tmp_bsln$efron_c <- ave(
    seq_along(tmp_bsln$DP_rev_i),
    tmp_bsln$DP_rev_i,
    FUN = function(ind) (seq_along(ind) - 1) / length(ind)
  )

  ds_bsln <- xgboost::xgb.DMatrix(
    data = as.matrix(tmp_bsln[, colnames(X_base), drop = FALSE]),
    label = tmp_bsln$I
  )

  event_times_bsln <- unique(tmp_bsln$DP_rev_i)

  risk_sets_bsln <- vector("list", length(event_times_bsln))
  event_sets_bsln <- vector("list", length(event_times_bsln))

  for (i in seq_along(event_times_bsln)) {
    stop_i <- event_times_bsln[i]
    risk_sets_bsln[[i]] <- which((tmp_bsln$TR_i < stop_i) & (tmp_bsln$DP_rev_i >= stop_i))
    event_sets_bsln[[i]] <- which(tmp_bsln$DP_rev_i == stop_i)
  }

  attr(ds_bsln, "truncation") <- tmp_bsln$TR_i
  attr(ds_bsln, "claim_arrival") <- tmp_bsln$DP_rev_i
  attr(ds_bsln, "risk_sets") <- risk_sets_bsln
  attr(ds_bsln, "event_sets") <- event_sets_bsln
  attr(ds_bsln, "efron_c") <- tmp_bsln$efron_c

  tie_table_bsln <- table(tmp_bsln$DP_rev_i)

  attr(ds_bsln, "tieid") <- unname(tie_table_bsln)
  attr(ds_bsln, "groups") <- rep(
    as.integer(names(tie_table_bsln)),
    unname(tie_table_bsln)
  )

  if (hazard_model == "COX") {

    predict_bsln <- model.out$train_expg

  }

  if (hazard_model == "NN") {

    tmp_order_bsln <- order(Y$DP_rev_i)
    x_train_bsln <- as.matrix(as.data.frame(X_base[tmp_order_bsln, , drop = FALSE]))

    model.out$net$eval()

    x_train_bsln_t <- torch::torch_tensor(
      as.matrix(x_train_bsln),
      dtype = torch::torch_float32()
    )

    predict_bsln <- as.numeric(torch::with_no_grad({
      model.out$net(x_train_bsln_t)
    }))
  }

  if (hazard_model == "XGB") {

    predict_bsln <- predict(model.out, ds_bsln)

  }

  predict_bsln <- predict_bsln - predict_bsln[1]

  risk_sum <- vapply(
    attr(ds_bsln, "risk_sets"),
    FUN = function(x, ypred) sum(exp(ypred[x])),
    ypred = predict_bsln,
    FUN.VALUE = numeric(1)
  )

  event_sum <- vapply(
    attr(ds_bsln, "event_sets"),
    FUN = function(x, ypred) sum(exp(ypred[x])),
    ypred = predict_bsln,
    FUN.VALUE = numeric(1)
  )

  n_events <- lengths(attr(ds_bsln, "event_sets"))

  denom <- risk_sum - eta * event_sum

  if (any(!is.finite(denom)) || any(denom <= 0)) {
    stop(
      "Non-positive denominator in baseline hazard calculation. ",
      "Check `eta`, fitted risk scores, and event/risk sets.",
      call. = FALSE
    )
  }

  bsln <- n_events / denom

  bsln <- data.table::data.table(
    baseline = bsln,
    DP_rev_i = sort(as.integer(unique(data$DP_rev_i)))
  )

  ## ------------------------------------------------------------
  ## 4. Benchmark and relative predictions
  ## ------------------------------------------------------------

  if (is_baseline_model) {

    newdata.bs <- data.frame(intercept_1 = rep(1, nrow(newdata)))
    remove_first_dummy_benchmark <- FALSE

  } else {

    X_parts_benchmark <- list()

    if (!is.null(cat_f)) {
      X_cat_benchmark <- fastDummies::dummy_cols(
        newdata,
        select_columns = cat_f,
        remove_selected_columns = TRUE,
        remove_first_dummy = FALSE
      )

      keep_cat_cols_benchmark <- vapply(
        colnames(X_cat_benchmark),
        function(z) {
          any(vapply(cat_f, function(p) grepl(pattern = p, x = z), logical(1)))
        },
        logical(1)
      )

      X_cat_benchmark <- data.table::as.data.table(X_cat_benchmark)
      X_parts_benchmark[["categorical"]] <- X_cat_benchmark[
        ,
        .SD,
        .SDcols = colnames(X_cat_benchmark)[keep_cat_cols_benchmark]
      ]
    }

    if (!is.null(cont_f)) {
      X_cont_benchmark <- data.table::copy(
        data.table::as.data.table(newdata)[, .SD, .SDcols = cont_f]
      )

      for (cft in cont_f) {
        mnv <- min(data[[cft]])
        mxv <- max(data[[cft]])
        X_cont_benchmark[[cft]] <- 2 * (X_cont_benchmark[[cft]] - mnv) / (mxv - mnv) - 1
      }

      X_parts_benchmark[["continuous"]] <- X_cont_benchmark
    }

    newdata.bs <- do.call(cbind, X_parts_benchmark)
    remove_first_dummy_benchmark <- hazard_model %in% c("COX", "XGB")
  }

  DT_benchmark <- data.table::as.data.table(
    cbind(X_base, DP_rev_i = Y$DP_rev_i)
  )

  benchmark <- DT_benchmark[
    order(DP_rev_i)
  ][
    1,
    .SD,
    .SDcols = !"DP_rev_i"
  ]

  newdata_benchmark <- data.table::as.data.table(newdata.bs)

  if (isTRUE(remove_first_dummy_benchmark)) {
    newdata_benchmark <- newdata_benchmark[
      ,
      .SD,
      .SDcols = colnames(newdata_benchmark) %in% names(X_base)
    ]
  }

  benchmark_id <- newdata_benchmark[
    benchmark,
    on = names(newdata_benchmark),
    which = TRUE
  ][1]

  pred_relative <- pred - pred[benchmark_id]

  hazard_frame <- data.table::copy(newdata)
  data.table::setDT(hazard_frame)
  hazard_frame[, expg := exp(pred_relative)]

  ## ------------------------------------------------------------
  ## 5. Likelihood evaluation
  ## ------------------------------------------------------------

  if (hazard_model == "COX") {

    if (is_baseline_model) {
      X <- data.frame(intercept_1 = rep(1, nrow(data)))
    } else {
      X <- data[, .SD, .SDcols = c(cont_f, cat_f)]
    }

    xy_tr <- cbind(X, Y)

    tmp_train_lkh <- xy_tr |>
      dplyr::arrange(DP_rev_i) |>
      dplyr::group_by(DP_rev_i) |>
      dplyr::mutate(efron_c = (seq_along(DP_rev_i) - 1) / length(DP_rev_i)) |>
      as.data.frame()

    ds_lkh <- X

    event_times_lkh <- unique(tmp_train_lkh$DP_rev_i)

    risk_sets_lkh <- vector("list", length(event_times_lkh))
    event_sets_lkh <- vector("list", length(event_times_lkh))

    for (i in seq_along(event_times_lkh)) {
      stop_i <- event_times_lkh[i]
      risk_sets_lkh[[i]] <- which((tmp_train_lkh$TR_i < stop_i) & (tmp_train_lkh$DP_rev_i >= stop_i))
      event_sets_lkh[[i]] <- which(tmp_train_lkh$DP_rev_i == stop_i)
    }

    attr(ds_lkh, "truncation") <- tmp_train_lkh$TR_i
    attr(ds_lkh, "claim_arrival") <- tmp_train_lkh$DP_rev_i
    attr(ds_lkh, "risk_sets") <- risk_sets_lkh
    attr(ds_lkh, "event_sets") <- event_sets_lkh
    attr(ds_lkh, "efron_c") <- tmp_train_lkh$efron_c

    tie_table_lkh <- table(tmp_train_lkh$DP_rev_i)

    attr(ds_lkh, "tieid") <- unname(tie_table_lkh)
    attr(ds_lkh, "groups") <- rep(
      as.integer(names(tie_table_lkh)),
      unname(tie_table_lkh)
    )

    preds_tr <- predict(model.out$cox, ds_lkh)

    risk_sets <- attr(ds_lkh, "risk_sets")
    event_sets <- attr(ds_lkh, "event_sets")
    efron_c <- attr(ds_lkh, "efron_c")
    tieid <- attr(ds_lkh, "tieid")

    exp_p_sum <- rep(
      sapply(risk_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds_tr),
      tieid
    )

    exp_p_tie <- rep(
      sapply(event_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds_tr),
      tieid
    )

    exp_p <- exp(preds_tr)
    r_k <- exp_p_sum - efron_c * exp_p_tie
    lkh <- exp_p / r_k

    is_lkh <- list(
      metric = "log-partial likelihood",
      value = -sum(log(lkh)) / length(preds_tr)
    )

    os_lkh <- NULL
  }

  if (hazard_model == "NN") {

    for (lkh_set in c("is", "os")) {

      if (lkh_set == "is") {
        X_lkh <- datads_pp$lkh_eval_data$data_train
        Y_lkh <- datads_pp$lkh_eval_data$y_train
      } else {
        X_lkh <- datads_pp$lkh_eval_data$data_val
        Y_lkh <- datads_pp$lkh_eval_data$y_val
      }

      if (!inherits(X_lkh, "data.frame")) {
        X_lkh <- as.data.frame(X_lkh)
      }

      data_train_lkh <- cbind(X_lkh, DP_rev_i = Y_lkh$DP_rev_i) |>
        dplyr::arrange(DP_rev_i) |>
        dplyr::select(-DP_rev_i) |>
        as.matrix()

      model.out$net$eval()

      data_train_lkh_t <- torch::torch_tensor(
        as.matrix(data_train_lkh),
        dtype = torch::torch_float32()
      )

      preds_tr <- as.numeric(torch::with_no_grad({
        model.out$net(data_train_lkh_t)
      }))

      preds_tr <- preds_tr - preds_tr[1]

      xy_tr <- cbind(X_lkh, Y_lkh)

      tmp_train_lkh <- xy_tr |>
        dplyr::arrange(DP_rev_i) |>
        dplyr::group_by(DP_rev_i) |>
        dplyr::mutate(efron_c = (seq_along(DP_rev_i) - 1) / length(DP_rev_i)) |>
        as.data.frame()

      ds_lkh <- tmp_train_lkh |>
        dplyr::arrange(DP_rev_i) |>
        dplyr::group_by(DP_rev_i) |>
        dplyr::mutate(efron_c = (seq_along(DP_rev_i) - 1) / length(DP_rev_i)) |>
        as.data.frame()

      event_times_lkh <- unique(tmp_train_lkh$DP_rev_i)

      risk_sets_lkh <- vector("list", length(event_times_lkh))
      event_sets_lkh <- vector("list", length(event_times_lkh))

      for (i in seq_along(event_times_lkh)) {
        stop_i <- event_times_lkh[i]
        risk_sets_lkh[[i]] <- which((tmp_train_lkh$TR_i < stop_i) & (tmp_train_lkh$DP_rev_i >= stop_i))
        event_sets_lkh[[i]] <- which(tmp_train_lkh$DP_rev_i == stop_i)
      }

      attr(ds_lkh, "truncation") <- tmp_train_lkh$TR_i
      attr(ds_lkh, "claim_arrival") <- tmp_train_lkh$DP_rev_i
      attr(ds_lkh, "risk_sets") <- risk_sets_lkh
      attr(ds_lkh, "event_sets") <- event_sets_lkh
      attr(ds_lkh, "efron_c") <- tmp_train_lkh$efron_c

      tie_table_lkh <- table(tmp_train_lkh$DP_rev_i)

      attr(ds_lkh, "tieid") <- unname(tie_table_lkh)
      attr(ds_lkh, "groups") <- rep(
        as.integer(names(tie_table_lkh)),
        unname(tie_table_lkh)
      )

      risk_sets <- attr(ds_lkh, "risk_sets")
      event_sets <- attr(ds_lkh, "event_sets")
      efron_c <- attr(ds_lkh, "efron_c")
      tieid <- attr(ds_lkh, "tieid")

      exp_p_sum <- rep(
        sapply(risk_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds_tr),
        tieid
      )

      exp_p_tie <- rep(
        sapply(event_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds_tr),
        tieid
      )

      exp_p <- exp(preds_tr)
      r_k <- exp_p_sum - efron_c * exp_p_tie
      lkh <- exp_p / r_k

      tmp_lkh <- list(
        metric = "log-partial likelihood",
        value = -sum(log(lkh)) / length(preds_tr)
      )

      if (lkh_set == "is") {
        is_lkh <- tmp_lkh
      } else {
        os_lkh <- tmp_lkh
      }
    }
  }

  if (hazard_model == "XGB") {

    for (lkh_set in c("is", "os")) {

      xy_tr <- cbind(X, Y) |>
        dplyr::arrange(DP_rev_i) |>
        as.data.frame()

      id <- seq_len(nrow(X))
      cond <- id %in% datads_pp$samples_cn$id

      if (lkh_set == "os") {
        cond <- !cond
      }

      tmp_train_lkh <- xy_tr[cond, , drop = FALSE] |>
        dplyr::arrange(DP_rev_i) |>
        dplyr::group_by(DP_rev_i) |>
        dplyr::mutate(efron_c = (seq_along(DP_rev_i) - 1) / length(DP_rev_i)) |>
        as.data.frame()

      tmp_train_x <- data.matrix(tmp_train_lkh |> dplyr::select(colnames(X)))

      tmp_train_x <- matrix(
        as.numeric(tmp_train_x),
        nrow = nrow(tmp_train_x),
        ncol = ncol(tmp_train_x),
        dimnames = dimnames(tmp_train_x)
      )

      ds_lkh <- xgboost::xgb.DMatrix(
        tmp_train_x,
        label = tmp_train_lkh$I
      )

      event_times_lkh <- unique(tmp_train_lkh$DP_rev_i)

      risk_sets_lkh <- vector("list", length(event_times_lkh))
      event_sets_lkh <- vector("list", length(event_times_lkh))

      for (i in seq_along(event_times_lkh)) {
        stop_i <- event_times_lkh[i]
        risk_sets_lkh[[i]] <- which((tmp_train_lkh$TR_i < stop_i) & (tmp_train_lkh$DP_rev_i >= stop_i))
        event_sets_lkh[[i]] <- which(tmp_train_lkh$DP_rev_i == stop_i)
      }

      attr(ds_lkh, "truncation") <- tmp_train_lkh$TR_i
      attr(ds_lkh, "claim_arrival") <- tmp_train_lkh$DP_rev_i
      attr(ds_lkh, "risk_sets") <- risk_sets_lkh
      attr(ds_lkh, "event_sets") <- event_sets_lkh
      attr(ds_lkh, "efron_c") <- tmp_train_lkh$efron_c

      tie_table_lkh <- table(tmp_train_lkh$DP_rev_i)

      attr(ds_lkh, "tieid") <- unname(tie_table_lkh)
      attr(ds_lkh, "groups") <- rep(
        as.integer(names(tie_table_lkh)),
        unname(tie_table_lkh)
      )

      preds_tr <- predict(model.out, ds_lkh)
      preds_tr <- preds_tr - preds_tr[1]

      risk_sets <- attr(ds_lkh, "risk_sets")
      event_sets <- attr(ds_lkh, "event_sets")
      efron_c <- attr(ds_lkh, "efron_c")
      tieid <- attr(ds_lkh, "tieid")

      exp_p_sum <- rep(
        sapply(risk_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds_tr),
        tieid
      )

      exp_p_tie <- rep(
        sapply(event_sets, FUN = function(x, ypred) sum(exp(ypred[x])), ypred = preds_tr),
        tieid
      )

      exp_p <- exp(preds_tr)
      r_k <- exp_p_sum - efron_c * exp_p_tie
      lkh <- exp_p / r_k

      tmp_lkh <- list(
        metric = "log-partial likelihood",
        value = -sum(log(lkh)) / length(preds_tr)
      )

      if (lkh_set == "is") {
        is_lkh <- tmp_lkh
      } else {
        os_lkh <- tmp_lkh
      }
    }
  }

  # Development factors -----------


  hazard_frame <- data.table::as.data.table(hazard_frame)
  bsln <- data.table::as.data.table(bsln)

  hazard_frame <- merge(
    hazard_frame,
    bsln,
    by = "DP_rev_i",
    all = TRUE
  )

  # replace NA in column "baseline" with 0
  hazard_frame[is.na(baseline), baseline := 0]


  # hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']




  # From hazard to development factors

  hazard_frame[,hazard:=baseline*expg]

  hazard_frame[,dev_f_i := (1+(1-..eta)*hazard)/(1-..eta*hazard)]

  hazard_frame[is.na(dev_f_i), dev_f_i := 1]
  hazard_frame[dev_f_i < 0, dev_f_i := 1]

  hazard_frame <- hazard_frame[order(DP_rev_i)]

  # columns grouping
  group_cols <- c(unique(c("AP_i",IndividualDataPP$data_information$continuous_features)),IndividualDataPP$data_information$categorical_features)

  hazard_frame[, `:=`(
    cum_dev_f_i = cumprod(dev_f_i),
    S_i = fifelse(cumprod(dev_f_i) == 0, 0, 1 / cumprod(dev_f_i))
  ), by = group_cols]

  hazard_frame[, `:=`(
    S_i_lead = shift(S_i, type = "lead", fill = 0),
    S_i_lag  = shift(S_i, type = "lag",  fill = 1)
  ), by = group_cols]

  # hazard_frame[, c("expg", "baseline", "hazard") := NULL]


  cols_to_remove_na_from <- c("dev_f_i", "S_i", "S_i_lead", "S_i_lag", "cum_dev_f_i")

  hazard_frame[, (cols_to_remove_na_from) := lapply(.SD, function(x) fifelse(is.na(x), 1, x)), .SDcols = cols_to_remove_na_from]





  # Prepare software output -----

  data_information <- IndividualDataPP$data_information
  data <- data.table::as.data.table(IndividualDataPP$training.data)

  time_unit_string <- c("days", "months", "quarters", "semesters", "years")
  time_unit_numeric <- c(1 / 360, 1 / 12, 1 / 4, 1 / 2, 1)

  input_pos <- match(data_information$input_time_granularity, time_unit_string)

  if (is.na(input_pos)) {
    stop(
      "`input_time_granularity` must be one of: ",
      paste(time_unit_string, collapse = ", "),
      call. = FALSE
    )
  }

  max_dp_i <- as.integer(data_information$years / time_unit_numeric[input_pos])

  ## Hazard-frame post-processing

  hazard_frame[, DP_i := max_dp_i - DP_rev_i + 1L]

  cols_hazard <- names(hazard_frame)

  setcolorder(
    hazard_frame,
    append(
      cols_hazard[cols_hazard != "DP_i"],
      "DP_i",
      after = which(cols_hazard == "AP_i")
    )
  )

  setnames(
    hazard_frame,
    old = c("dev_f_i", "cum_dev_f_i"),
    new = c("f_i", "cum_f_i")
  )

  ## Data needed for reserving

  feature_cols <- unique(c(
    data_information$categorical_features,
    data_information$continuous_features,
    "AP_i"
  ))

  existing_cols <- unique(c(
    feature_cols,
    "DP_i"
  ))

  observed_condition <- (max_dp_i - data$DP_i + 1L) > (data$AP_i - 1L)

  tmp_grid <- unique(
    data[observed_condition, .SD, .SDcols = feature_cols]
  )[
    , .(DP_i = seq_len(max_dp_i)),
    by = feature_cols
  ]

  tmp_existing <- unique(
    data[observed_condition, .SD, .SDcols = existing_cols]
  )

  ## data.table equivalent of dplyr::setdiff()
  ## Requires same columns and same column order.

  tmp_existing <- tmp_existing[, .SD, .SDcols = names(tmp_grid)]

  tmp_missing <- data.table::fsetdiff(
    x = tmp_grid,
    y = tmp_existing,
    all = FALSE
  )

  if (nrow(tmp_missing) > 0L) {

    tmp_missing[
      ,
      `:=`(
        DP_rev_i = max_dp_i - DP_i + 1L,
        TR_i = AP_i - 1L,
        I = 0L
      )
    ]

    tmp_missing <- tmp_missing[DP_rev_i > TR_i]

    tmp_missing[
      ,
      `:=`(
        DP_rev_o =
          floor(max_dp_i * data_information$conversion_factor) -
          ceiling(
            DP_i * data_information$conversion_factor +
              ((AP_i - 1L) %% (1 / data_information$conversion_factor)) *
              data_information$conversion_factor
          ) +
          1L,
        AP_o = ceiling(AP_i * data_information$conversion_factor)
      )
    ][
      ,
      TR_o := AP_o - 1L
    ]

    if (!is.null(data_information$categorical_features)) {
      tmp_missing[
        ,
        (data_information$categorical_features) := lapply(.SD, as.factor),
        .SDcols = data_information$categorical_features
      ]
    }

    reserving_cols <- unique(c(
      data_information$categorical_features,
      data_information$continuous_features,
      "AP_i",
      "AP_o",
      "DP_i",
      "DP_rev_i",
      "DP_rev_o",
      "TR_i",
      "TR_o",
      "I"
    ))

    tmp_missing <- tmp_missing[
      ,
      .SD,
      .SDcols = intersect(reserving_cols, names(tmp_missing))
    ]

  } else {

    tmp_missing <- NULL
  }

  data_information$data_for_reserving <- data.table::rbindlist(
    list(data, tmp_missing),
    use.names = TRUE,
    fill = TRUE
  )

  ## Development-period ranges needed later

  development_periods <- unique(data[, .(AP_i, AP_o)])

  cf <- data_information$conversion_factor

  output_pos <- match(data_information$output_time_granularity, time_unit_string)

  if (is.na(output_pos)) {
    stop(
      "`output_time_granularity` must be one of: ",
      paste(time_unit_string, collapse = ", "),
      call. = FALSE
    )
  }

  max_dp_o <- as.integer(data_information$years / time_unit_numeric[output_pos])

  dp_ranges <- development_periods[
    ,
    .(DP_rev_o = seq_len(max_dp_o)),
    by = .(AP_i, AP_o)
  ][
    ,
    `:=`(
      min_dp = AP_i + (DP_rev_o - AP_o) / cf,
      max_dp = AP_i - 1 + (DP_rev_o - AP_o + 1) / cf
    )
  ]

  data_information$dp_ranges <- dp_ranges

  out <- list(
    model.out = list(
      data = X,
      model.out = model.out
    ),
    hazard_frame = hazard_frame,
    data_information = data_information,
    fit_information = list(
      hazard_model = hazard_model,
      is_lkh = is_lkh,
      os_lkh = os_lkh,
      eta = eta
    )
  )

  class(out) <- c("ReSurvFit")

  return(out)


}



