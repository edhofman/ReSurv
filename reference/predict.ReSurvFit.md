# Predict IBNR frequency

This function predicts the results from the ReSurv fits.

## Usage

``` r
# S3 method for class 'ReSurvFit'
predict(
  object,
  newdata = NULL,
  lower_triangular_output = FALSE,
  minimal_output = FALSE,
  check_value = 1.85,
  ...
)
```

## Arguments

- object:

  A fitted `ReSurvFit` object.

- newdata:

  An optional `IndividualDataPP` object using the same features and time
  scale as the fitted data. NULL predicts for the fitted data.

- lower_triangular_output:

  `logical`, if set to `TRUE` we add the predicted lower triangle in
  input and output granularity to the `predict.ReSurvFit` output.

- minimal_output:

  `logical`, if set to `TRUE` return a reduced prediction object.

- check_value:

  Retained for compatibility; currently ignored.

- ...:

  Additional arguments to pass to the predict function.

## Value

A `ReSurvPredict` object with fitted predictions, long triangle outputs,
predicted counts, and optional lower-triangle outputs.
