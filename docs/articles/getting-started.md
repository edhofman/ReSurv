# Getting started with ReSurv

ReSurv estimates IBNR claim counts using reverse-time proportional
hazard models. The following small example runs during package checks.

## Prepare individual claims

``` r

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
head(individual$training.data)
#>    claim_type  AP_i  AP_o  DP_i DP_rev_i DP_rev_o  TR_i  TR_o     I
#>        <fctr> <int> <int> <int>    <int>    <int> <int> <int> <int>
#> 1:          0     2     2     3        2        2     1     1     1
#> 2:          0     2     2     1        4        4     1     1     1
#> 3:          0     2     2     3        2        2     1     1     1
#> 4:          1     1     1     3        2        2     0     0     1
#> 5:          1     1     1     4        1        1     0     0     1
#> 6:          1     1     1     4        1        1     0     0     1
```

Each row represents a claim. `accident_period` and `calendar_period`
name the accident and reporting columns. Supply `id` to retain only the
first row per identifier. Continuous and categorical features are
specified by column name. `continuous_features_spline` selects
continuous features to smooth in Cox models.

Supported time units are days, months, quarters, semesters, and years.
Output units must be equal to or coarser than input units and permit
exact grouping. The horizon is `years`; when omitted it is inferred from
development times. For a known valuation horizon, specify it explicitly.
`training.data` contains observed upper-triangle rows, `full.data`
contains all encoded rows, and `data_information` stores formulas,
features, and time-scale metadata.

## Fit and predict

``` r

fit <- ReSurv(individual, hazard_model = "COX", eta = 0)
prediction <- predict(fit)
summary(prediction)
#> 
#>  Hazard model:
#> "COX"
#> 
#> 
#> Categorical Features:
#> claim_type
#> Total IBNR level: 
#> [1]  1.16
head(predictReserve(fit, granularity = "output"))
#>       AP    DP    CP     IBNR
#>    <int> <int> <int>    <num>
#> 1:     2     4     5 1.159703
```

`eta` controls the baseline and hazard-to-development-factor conversion
and lies in \[0, 1\]; its default is 0.5. The example uses 0.
Predictions are counts, not monetary amounts.
[`predictReserve()`](https://edhofman.github.io/ReSurv/reference/predictReserve.md)
returns `AP`, `DP`, `CP`, and `IBNR` with `CP = AP + DP - 1`. Use
`granularity = "input"` for the original time scale.
[`predict()`](https://rdrr.io/r/stats/predict.html) provides fuller
results, including `predicted_counts` and `long_triangle_format_out`;
set `lower_triangular_output = TRUE` for matrices.

## Tune an XGBoost model

Cross-validation uses sequential folds and returns the best parameter
list; the final model is fitted separately. This short example uses only
two boosting rounds to demonstrate the API, rather than to select a
production model.

``` r

cv <- ReSurvCV(
  individual, model = "XGB",
  hparameters_grid = list(
    booster = "gbtree", eta = c(0.05, 0.1), max_depth = 1,
    subsample = 1, alpha = 0, lambda = 1, min_child_weight = 0,
    nthread = 1
  ),
  folds = 2, random_seed = 1, nrounds = 2
)
cv$out.cv.best.oos
#>   booster  eta max_depth subsample alpha lambda min_child_weight nthread
#> 1  gbtree 0.05         1         1     0      1                0       1
#>   train.lkh  test.lkh         time
#> 1 0.4583005 0.4719134 0.0002375166
xgb_fit <- ReSurv(
  individual, hazard_model = "XGB", eta = 0,
  hparameters = cv$hparameters.best
)
predict(xgb_fit)$predicted_counts
#> [1] 1.58113
```

The XGBoost grid’s `eta` is its learning rate; the `eta` argument to
[`ReSurv()`](https://edhofman.github.io/ReSurv/reference/ReSurv.md)
controls the development-factor conversion. `cv_data_subsample`
optionally reduces the rows used for cross-validation: 0.1 uses ten
percent, while 1 uses all rows. Refit on the full preprocessed data
after selection.

## Neural networks

NN models require the optional R package `torch` and its runtime.
Install these once with `install.packages("torch")` and
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html).
The example below is not run when building the vignette because runtime
installation is a separate user action.

``` r

nn_fit <- ReSurv(
  individual, hazard_model = "NN", eta = 0,
  hparameters = list(
    num_layers = 1, num_nodes = 8, activation = "relu",
    optim = "Adam", lr = 0.01, xi = 0.5, eps = 0,
    early_stopping = TRUE, patience = 5, epochs = 20,
    verbose = FALSE, num_workers = 0
  )
)
predictReserve(nn_fit)
```

## Evaluate realized counts

Use
`Score_Reserving(models = list(COX = fit), newdata = realized_claims)`
when the lower-triangle outcomes are available. `newdata` may contain
individual claims with the original accident/reporting columns, or an
aggregate table with `AP`, `DP`, `CP`, and `actual` counts on the
requested granularity scale. Metrics include `EI`, `R-tot`,
`R-cell-wise`, `R-cal-wise`, and optionally `CRPS`. The default includes
a chain-ladder benchmark; `clmplus_benchmark = c("ac", "apc")`
additionally requires the optional `clmplus` package.
