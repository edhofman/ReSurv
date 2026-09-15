# Predict deterministic reserve table

Predict deterministic reserve table

Predict deterministic reserve table from a ReSurvFit object

## Usage

``` r
predictReserve(object, ...)

# S3 method for class 'ReSurvFit'
predictReserve(object, granularity = c("output", "input"), ...)
```

## Arguments

- object:

  A ReSurvFit object.

- ...:

  Additional arguments passed to predict.ReSurvFit().

- granularity:

  Character. Either "output" or "input".

## Value

A `data.table` with accident period `AP`, development period `DP`,
calendar period `CP`, and predicted claim count `IBNR`. These are claim
counts, not monetary reserves.

A data.table with columns AP, DP, CP, IBNR.

## See also

[`predict.ReSurvFit`](https://edhofman.github.io/ReSurv/reference/predict.ReSurvFit.md)
