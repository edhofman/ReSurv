# K-fold cross-validation of a ReSurv model

K-fold cross-validation of a ReSurv model

## Usage

``` r
ReSurvCV(
  IndividualDataPP,
  model,
  hparameters_grid,
  folds,
  random_seed,
  continuous_features_scaling_method = "minmax",
  print_every_n = 1L,
  nrounds = NULL,
  early_stopping_rounds = NULL,
  epochs = 1,
  parallel = FALSE,
  ncores = 1,
  num_workers = 0,
  verbose = FALSE,
  verbose.cv = FALSE,
  cv_data_subsample = 1
)
```

## Arguments

- IndividualDataPP:

  An object of class `IndividualDataPP`.

- model:

  Character. Either `"NN"` or `"XGB"`.

- hparameters_grid:

  Named list defining the hyperparameter grid.

- folds:

  Integer. Number of folds.

- random_seed:

  Integer. Random seed.

- continuous_features_scaling_method:

  Character. Scaling method for continuous features.

- print_every_n:

  Integer. Passed to XGBoost.

- nrounds:

  Integer. Number of XGBoost boosting rounds.

- early_stopping_rounds:

  Integer. XGBoost early stopping.

- epochs:

  Integer. Number of NN epochs.

- parallel:

  Logical. Retained for compatibility; execution is sequential and
  setting this to TRUE produces a warning.

- ncores:

  Integer. Retained for compatibility; currently ignored.

- num_workers:

  Deprecated for the native torch backend. Ignored.

- verbose:

  Logical. Print model fitting output.

- verbose.cv:

  Logical. Print CV progress.

- cv_data_subsample:

  Fraction of training rows used for cross-validation, in (0, 1\].
  Values greater than 1 and at most 100 are interpreted as percentages.
  The default 1 uses all rows; use 0.01 for one percent.

## Value

An object of class `ReSurvCV` containing `out.cv` (all combinations and
mean training and validation losses), `out.cv.best.oos` (the row with
smallest validation loss), and `hparameters.best`, suitable for the
`hparameters` argument of
[`ReSurv()`](https://edhofman.github.io/ReSurv/reference/ReSurv.md).
Cross-validation does not refit the final model.
