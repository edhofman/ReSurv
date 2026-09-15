# Individual Data Pre-Processing

This function pre-processes the data for the application of a `ReSurv`
model.

## Usage

``` r
IndividualDataPP(
  data,
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
  degrees_of_freedom_cp = 4
)
```

## Arguments

- data:

  `data.frame`, for the individual reserving. The number of development
  periods can be larger than the number of accident periods.

- id:

  `character`, `data` column that contains the policy identifier. If
  `NULL` (default), we assume that each row is an observation. We assume
  that each observation can only have one reporting time, if not null we
  take the reporting time of the first row for each `id`.

- continuous_features:

  `character`, continuous features columns to be scaled.

- categorical_features:

  `character`, categorical features columns to be one-hot encoded.

- accident_period:

  `character`, it contains the name of the column in data corresponding
  to the accident period.

- calendar_period:

  `character`, it contains the name of the column in data corresponding
  to the calendar period.

- input_time_granularity:

  `character`, time unit of the input data. Granularity supported:

  - `"days"`: the input data are daily.

  - `"months"`: the input data are monthly.

  - `"quarters"`: the input data are quarterly

  - `"semesters"`: six-month periods.

  - `"years"`: the input data are yearly.

  Default to `months`.

- output_time_granularity:

  `character`, time unit of the output data. The granularity supported
  is the same as for the input data:

  - `"days"`: the output data will be on a daily scale.

  - `"months"`: the output data will be on a monthly scale.

  - `"quarters"`: the output data will be on a quarterly scale.

  - `"semesters"`: six-month periods.

  - `"years"`: the output data will be on yearly scale.

  The output granularity must be equal to or coarser than the input
  granularity. Also, the output granularity must be consistent with the
  input granularity, meaning that the time conversion must be possible.
  E.g., it is possible to group quarters to years. Quarters can also be
  grouped to semesters. Default to `quarters`.

- years:

  `numeric`, number of development years in the study.

- calendar_period_extrapolation:

  `logical`, whether a spline for calendar extrapolation should be
  considered in the cox model fit. Default is \`FALSE\`.

- continuous_features_spline:

  `character`, names of continuous features to model with splines; NULL
  uses linear terms. Use `"AP_i"` for a remapped accident-period
  feature.

- degrees_cf:

  `numeric`, degrees of the spline for smoothing continuous features.

- degrees_of_freedom_cf:

  `numeric`, degrees of freedom of the splines for smoothing continuous
  features.

- degrees_cp:

  `numeric`, degrees of the spline for smoothing the calendar period
  effect.

- degrees_of_freedom_cp:

  `numeric`, degrees of freedom of the splines for smoothing the
  calendar period effect.

## Value

An `IndividualDataPP` list containing `training.data` (observed rows),
`full.data` (all encoded rows), and `data_information` (conversion
factor, input/output formulas, feature names, time units, horizon, and
original column names).

After pre-processing, we provide a standard encoding for the time
components. This regards the output in `training.data`. In the `ReSurv`
notation:

- `AP_i`: Input granularity accident period.

- `AP_o`: Output granularity accident period.

- `DP_i`: Input granularity development period in forward time.

- `DP_rev_i`: Input granularity development period in reverse time.

- `DP_rev_o`: Output granularity development period in reverse time.

- `TR_i`: Input granularity truncation time.

- `TR_o`: Output granularity truncation time.

- `I`: event indicator, under this framework is equal to one for each
  entry.

## Details

Accident and reporting periods are indexed from one. Development time is
`DP_i = RP_i - AP_i + 1`; reverse development time is
`DP_rev_i = DP_max - DP_i + 1`, and truncation time is
`TR_i = AP_i - 1`. Training retains observed rows with
`DP_rev_i > TR_i`.

The conversion factor is the input time unit divided by the output time
unit (for example, 1/3 for months to quarters). Accident and calendar
periods are grouped using ceiling; development periods also account for
the position of the accident period within each output period. Days use
a 360-day year for the development horizon.

## References

Hiabu, M., Hofman, E., & Pittarello, G. (2023). A machine learning
approach based on survival analysis for IBNR frequencies in non-life
reserving. arXiv preprint arXiv:2312.14549.

## Examples

``` r

input_data_0 <- data_generator(
random_seed = 1964,
scenario = "alpha",
time_unit = 1,
years = 2,
period_exposure = 100)

individual_data <- IndividualDataPP(data = input_data_0,
categorical_features = "claim_type",
continuous_features = "AP",
accident_period = "AP",
calendar_period = "RP",
input_time_granularity = "years",
output_time_granularity = "years",
years = 2)



```
