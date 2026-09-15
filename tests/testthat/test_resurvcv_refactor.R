test_that("cv_design_matrix uses data_information slots", {
  training_data <- data.frame(
    DP_rev_i = c(2, 2, 3, 3, 4, 4),
    I = 1,
    TR_i = c(0, 0, 1, 1, 2, 2),
    AP_i = c(1, 1, 2, 2, 3, 3),
    claim_type = factor(c("A", "B", "A", "B", "A", "B"))
  )

  x <- list(
    training.data = training_data,
    data_information = list(
      categorical_features = "claim_type",
      continuous_features = "AP_i"
    )
  )

  class(x) <- "IndividualDataPP"

  xy <- cv_design_matrix(
    IndividualDataPP = x,
    continuous_features_scaling_method = "minmax",
    remove_first_dummy = TRUE
  )

  testthat::expect_true(is.data.frame(xy$X))
  testthat::expect_true(is.data.frame(xy$Y))
  testthat::expect_equal(nrow(xy$X), nrow(training_data))
  testthat::expect_equal(nrow(xy$Y), nrow(training_data))
  testthat::expect_true(all(c("DP_rev_i", "I", "TR_i") %in% names(xy$Y)))
})


test_that("ReSurvCV validates inputs", {
  training_data <- data.frame(
    DP_rev_i = c(2, 2, 3, 3),
    I = 1,
    TR_i = c(0, 0, 1, 1),
    AP_i = c(1, 1, 2, 2)
  )

  x <- list(
    training.data = training_data,
    data_information = list(
      categorical_features = NULL,
      continuous_features = "AP_i"
    )
  )

  class(x) <- "IndividualDataPP"

  testthat::expect_error(
    ReSurvCV(
      x,
      model = "COX",
      hparameters_grid = list(eta = 0.1),
      folds = 2,
      random_seed = 1
    ),
    "`model`"
  )

  testthat::expect_error(
    ReSurvCV(
      x,
      model = "XGB",
      hparameters_grid = list(eta = 0.1),
      folds = 10,
      random_seed = 1
    ),
    "`folds` cannot exceed"
  )
})


test_that("ReSurvCV XGB runs on a tiny current-style IndividualDataPP object", {
  testthat::skip_if_not_installed("xgboost")

  training_data <- data.frame(
    DP_rev_i = c(2, 2, 3, 3, 4, 4, 5, 5),
    I = 1,
    TR_i = c(0, 0, 0, 1, 1, 2, 2, 3),
    AP_i = c(1, 1, 1, 2, 2, 3, 3, 4)
  )

  x <- list(
    training.data = training_data,
    data_information = list(
      categorical_features = NULL,
      continuous_features = "AP_i"
    )
  )

  class(x) <- "IndividualDataPP"

  out <- ReSurvCV(
    x,
    model = "XGB",
    hparameters_grid = list(
      booster = "gbtree",
      eta = 0.1,
      max_depth = 1,
      subsample = 1,
      alpha = 0,
      lambda = 1,
      min_child_weight = 0
    ),
    folds = 2,
    random_seed = 1,
    nrounds = 2,
    early_stopping_rounds = NULL,
    verbose = FALSE,
    verbose.cv = FALSE
  )

  testthat::expect_s3_class(out, "ReSurvCV")
  testthat::expect_true(all(c("out.cv", "out.cv.best.oos") %in% names(out)))
  testthat::expect_equal(nrow(out$out.cv), 1)
})
