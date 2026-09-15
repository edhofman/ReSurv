# Score reserving predictions

Compare predictions with realized claim counts in the lower triangle.

## Usage

``` r
Score_Reserving(
  models,
  newdata,
  scoring_metrics = c("EI", "R-tot", "R-cell-wise", "R-cal-wise"),
  granularity = c("output", "input"),
  chain_ladder = TRUE,
  clmplus_benchmark = NULL,
  ...
)
```

## Arguments

- models:

  A fitted `ReSurvFit` or a named list of models. At least one entry
  must be a `ReSurvFit`; other entries must support
  [`predictReserve()`](https://edhofman.github.io/ReSurv/reference/predictReserve.md).

- newdata:

  Realized claims, either individual data with the accident and
  reporting columns used for fitting, or an aggregate table with `AP`,
  `DP`, `CP`, and counts in `actual`, `I`, or `IBNR`. Aggregate periods
  must be on the requested `granularity` scale.

- scoring_metrics:

  Character vector selecting `"EI"`, `"R-tot"`, `"R-cell-wise"`,
  `"R-cal-wise"`, or `"CRPS"`. CRPS is available only for `ReSurvFit`
  models.

- granularity:

  Time scale for scoring: `"output"` (default) or `"input"`.

- chain_ladder:

  Logical; include the aggregate chain-ladder benchmark.

- clmplus_benchmark:

  Optional character vector containing any of \`"ac"\` or \`"apc"\`.
  Requested models are fitted with the \`clmplus\` package to the same
  aggregate triangle as the chain-ladder benchmark. \`NULL\` (the
  default) disables these benchmarks.

- ...:

  Additional arguments passed to
  [`predictReserve()`](https://edhofman.github.io/ReSurv/reference/predictReserve.md).

## Value

A list of score tables, one per requested metric, with class
`Score_Reserving` and a `granularity` attribute.
