# Individual data generator

This function generates individual claims data on the requested time
scale using the `SynthETIC` package. This simple function allows to
simulate from a sand-box to test out the `ReSurv` approach. Some
parameters of the simulation can be changed.

## Usage

``` r
data_generator(
  ref_claim = 2e+05,
  time_unit = 1/360,
  years = 4,
  random_seed = 1964,
  period_exposure = 200,
  period_frequency = 0.2,
  scenario = 1
)
```

## Arguments

- ref_claim:

  `integer`, reference claim size.

- time_unit:

  Numeric, fraction of a year per output period, for example 1/360 (the
  default) for days, 1/12 for months, or 1 for years.

- years:

  `integer`, number of years to be simulated.

- random_seed:

  `integer`, random seed for replicable code.

- period_exposure:

  `integer`, volume (number of policies) underwritten each period.

- period_frequency:

  `numeric`, expected frequency in each period.

- scenario:

  `character` or `numeric`, one of the scenarios described in the
  accompanying manuscript. Possible choices are 'alpha' (0), 'beta' (1),
  'gamma'(2), 'delta'(3),'epsilon'(4). Our simulated data are
  constituted of a mix of short tail claims (`claim_type 0`) and claims
  with longer resolution (`claim_type 1`). We chose the parameter of the
  simulator to resemble a mix of property damage (`claim_type 0`) and
  bodily injuries (`claim_type 1`). each scenario has distinctive
  characteristics. Scenario Alpha is a mix of `claim_type 0` and
  `claim_type 1` with same number of claims volume at each accident
  period. Differently from scenario Alpha, in scenario Beta the volumes
  of `claim_type 1` are decreasing in the most recent accident periods.
  In scenario Gamma we add an interaction between `claim_type 1` and
  accident period: in a real world setting this can be motivated by a
  change in consumer behavior or company policies resulted in different
  reporting patterns over time. In scenario Delta, we introduce a
  seasonality effect dependent on the accident period for `claim_type 0`
  and `claim_type 1`. In the real word, scenario Delta resembles
  seasonal changes in the workforce composition. Scenario Epsilon does
  not satisfy the proportionality assumption. Additional scenarios are
  `"zeta"` (5), with age, property value and business-use covariates,
  and `"eta"` (6), with business-use and accident-period effects.

## Value

Individual claims data. It contains the following columns:

- `claim_number`: Policy ID.

- `claim_type`: Type of claim (0 or 1) in scenarios alpha through
  epsilon. Zeta instead includes `age`, `property_value`, and
  `business_use`; eta includes `business_use`.

- `AP`: Accident period

- `RP`: Reporting period.

## References

Avanzi, B., Taylor, G., Wang, M., & Wong, B. (2021). SynthETIC: an
individual insurance claim simulator with feature control. Insurance:
Mathematics and Economics, 100, 296-308.

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



```
