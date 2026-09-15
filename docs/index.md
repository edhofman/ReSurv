# ReSurv

ReSurv predicts incurred-but-not-reported (IBNR) claim **counts** from
individual claims data. It links reverse-time proportional hazards to
development factors, with Cox regression (`"COX"`), gradient boosting
(`"XGB"`), and neural networks (`"NN"`). It includes simulation,
preprocessing, cross-validation, prediction, and comparison with
realized claim counts.

## Installation

Install the development version with R 4.1 or later:

``` r

install.packages("remotes")
remotes::install_github("edhofman/ReSurv")
```

Neural networks use the native R **torch** backend. Before using
`hazard_model = "NN"`, install the optional package and its runtime:

``` r

install.packages("torch")
torch::install_torch()
```

Python and a virtual environment are no longer required. Cox and XGBoost
models do not require torch.

## Example

``` r

library(ReSurv)
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
fit <- ReSurv(individual, hazard_model = "COX", eta = 0)
prediction <- predict(fit)
summary(prediction)
head(predictReserve(fit, granularity = "output"))
```

The reserve table contains accident period (`AP`), development period
(`DP`), calendar period (`CP`), and predicted count (`IBNR`). Periods
start at one and `CP = AP + DP - 1`. Preprocessing retains the observed
upper triangle for fitting; the development horizon is controlled by
`years`.

[`ReSurvCV()`](https://edhofman.github.io/ReSurv/reference/ReSurvCV.md)
selects NN or XGB hyperparameters. Pass its `hparameters.best` component
to [`ReSurv()`](https://edhofman.github.io/ReSurv/reference/ReSurv.md)
to fit the selected model.
[`Score_Reserving()`](https://edhofman.github.io/ReSurv/reference/Score_Reserving.md)
compares models against realized claims and optionally adds chain-ladder
and `clmplus` benchmarks. See the getting-started article on the
[documentation website](https://edhofman.github.io/ReSurv/) and the R
help pages for parameters and return values. The vignette source is
`vignettes/getting-started.Rmd` in this repository.

## Development checks

GitHub Actions checks Windows, macOS, and Linux across R release,
development, and the previous release. A separate Linux release job
includes the PDF manual and vignettes and fails on notes. The audit
verifies generated documentation, runs tests, checks URLs, and exercises
all three model backends with torch installed. The website workflow
builds documentation on pull requests and publishes from `main` or
`master` to `gh-pages`.

Run the core checks locally with:

``` r

devtools::document()
testthat::test_local(stop_on_failure = TRUE)
devtools::check(document = FALSE, manual = TRUE)
```

The `articles/historical/` directory preserves earlier replication
sources. Its README explains their status. They are separate from the
current package vignettes.

## Reference

Hiabu, M., Hofman, E., and Pittarello, G. (2023). *A machine learning
approach based on survival analysis for IBNR frequencies in non-life
reserving.*
[doi:10.48550/arXiv.2312.14549](https://doi.org/10.48550/arXiv.2312.14549).
