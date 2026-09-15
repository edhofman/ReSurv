# Fit `ReSurv` models on the individual data.

This function fits and computes the reserves for the `ReSurv` models

## Usage

``` r
ReSurv(
  IndividualDataPP,
  hazard_model = "COX",
  tie = "efron",
  baseline = "spline",
  continuous_features_scaling_method = "minmax",
  random_seed = 1,
  hparameters = list(),
  percentage_data_training = 0.8,
  grouping_method = "exposure",
  check_value = 1.85,
  eta = 0.5,
  simplifier = TRUE
)

# Default S3 method
ReSurv(
  IndividualDataPP,
  hazard_model = "COX",
  tie = "efron",
  baseline = "spline",
  continuous_features_scaling_method = "minmax",
  random_seed = 1,
  hparameters = list(),
  percentage_data_training = 0.8,
  grouping_method = "exposure",
  check_value = 1.85,
  eta = 0.5,
  simplifier = TRUE
)

# S3 method for class 'IndividualDataPP'
ReSurv(
  IndividualDataPP,
  hazard_model = "COX",
  tie = "efron",
  baseline = "spline",
  continuous_features_scaling_method = "minmax",
  random_seed = 1,
  hparameters = list(),
  percentage_data_training = 0.8,
  grouping_method = "exposure",
  check_value = 1.85,
  eta = 0.5,
  simplifier = TRUE
)
```

## Arguments

- IndividualDataPP:

  IndividualDataPP object to use for the `ReSurv` fit.

- hazard_model:

  `character`, hazard model supported from our package, must be provided
  as a string. The model can be chosen from:

  - `"COX"`: Standard Cox model for the hazard.

  - `"NN"`: Deep Survival Neural Network.

  - `"XGB"`: eXtreme Gradient Boosting.

- tie:

  Handling of ties in the Cox fit, passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html):
  `"efron"` (default), `"breslow"`, or `"exact"`. NN and XGB fitting use
  Efron handling.

- baseline:

  Retained for compatibility. The baseline is estimated from the risk
  sets using the specified `eta`.

- continuous_features_scaling_method:

  Scaling method: `"minmax"` or `"standard"`.

- random_seed:

  `integer`, random seed set for reproducibility

- hparameters:

  `list`, hyperparameters for the machine learning models. It will be
  disregarded for the cox approach.

- percentage_data_training:

  `numeric`, fraction in (0, 1\] of observed data used for training; the
  remainder is used for validation in the NN and XGB models.

- grouping_method:

  Retained for compatibility; currently ignored.

- check_value:

  Retained for compatibility; currently ignored.

- eta:

  Numeric in \[0, 1\], controlling the baseline estimate and
  hazard-to-development-factor conversion. Default is 0.5.

- simplifier:

  `logical`, kept for compatibility. The simplified forecast frame is
  always used.

## Value

A `ReSurvFit` list containing `model.out` (design matrix and fitted
backend model), `hazard_frame` (risk scores, baseline, hazards,
development factors and survival probabilities), `data_information`
(preprocessing and reserving metadata), and `fit_information` (model
name, training loss `is_lkh`, validation loss `os_lkh`, and `eta`). Use
[`predict()`](https://rdrr.io/r/stats/predict.html) or
[`predictReserve()`](https://edhofman.github.io/ReSurv/reference/predictReserve.md)
to obtain claim-count predictions.

## Details

The model fit uses the theoretical framework of Hiabu et al. (2023),
that relies on the correspondence between hazard models and development
factors:

Neural networks use the native R torch backend. Install its runtime with
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html)
before fitting a neural network.

The `ReSurv` package assumes proportional hazard models. Given an i.i.d.
sample \\\left\\y_i,x_i\right\\\_{i=1, \ldots, n}\\ the individual
hazard at time \\t\\ is:

\\\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}\\

Composed of a baseline \\\lambda_0(t)\\ and a proportional effect
\\e^{y_i(x_i)}\\.

Currently, the implementation allows to optimize the partial likelihood
(concerning the proportional effects) using one of the following
statistical learning approaches:

- [COX](https://github.com/therneau/survival)

- [Neural
  Networks](https://link.springer.com/article/10.1186/s12874-018-0482-1)

- [eXtreme Gradient Boosting](https://xgboost.readthedocs.io/en/stable/)

## References

Hiabu, M., Hofman, E., & Pittarello, G. (2023). A machine learning
approach based on survival analysis for IBNR frequencies in non-life
reserving. arXiv preprint arXiv:2312.14549.

Therneau, T. M., & Lumley, T. (2015). Package â€˜survivalâ€™. R Top Doc,
128(10), 28-33.

Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., &
Kluger, Y. (2018). DeepSurv: personalized treatment recommender system
using a Cox proportional hazards deep neural network. BMC medical
research methodology, 18(1), 1-12.

Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package
â€˜xgboostâ€™. R version, 90, 1-66.

## Examples

``` r

input_data_0 <- data_generator(
random_seed = 1964,
scenario = "alpha",
time_unit = 1,
years = 4,
period_exposure = 100)

individual_data <- IndividualDataPP(data = input_data_0,
categorical_features = "claim_type",
continuous_features = "AP",
accident_period = "AP",
calendar_period = "RP",
input_time_granularity = "years",
output_time_granularity = "years",
years=4)


resurv_fit_cox <- ReSurv(individual_data,
hazard_model = "COX",
eta = 0)



```
