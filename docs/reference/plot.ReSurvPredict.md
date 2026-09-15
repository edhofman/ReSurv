# Plot of the development factors

Plots the development factors by group code.

## Usage

``` r
# S3 method for class 'ReSurvPredict'
plot(
  x,
  granularity = "input",
  group_code = 1,
  color_par = "royalblue",
  linewidth_par = 2.5,
  ylim_par = NULL,
  ticks_by_par = NULL,
  base_size_par = NULL,
  title_par = NULL,
  x_text_par = NULL,
  plot.title.size_par = NULL,
  ...
)
```

## Arguments

- x:

  "ReSurvPredict" object specifying hazard and development factors.

- granularity:

  `character`, either `"input"` for `input_time_granularity` or
  `"output"` for `output_time_granularity`.

- group_code:

  `numeric`: Identifier for the group that will be plotted. Default
  is 1. The code identifiers can be find in the
  `ReSurvPredict$long_triangle_format_out` list. Depending on the
  granularity of interest, it will be either in
  `ReSurvPredict$long_triangle_format_out$input_granularity` for
  `input_time_granularity` or
  `ReSurvPredict$long_triangle_format_out$output_granularity` for
  `output_time_granularity`.

- color_par:

  `character`: `ggplot2` Colour of the line plot. Default is
  `'royalblue'`. Optional.

- linewidth_par:

  `numeric`: Line plot width. Optional.

- ylim_par:

  `numeric`: Highest intercept on the y-axis (development factors). The
  default is the highest predicted development factor. Optional.

- ticks_by_par:

  `numeric`: gap between each x-axis label (development period). Default
  is 2. Optional.

- base_size_par:

  `numeric`: base size of the plot. Default is 5. See `base_size` in the
  `?theme_bw` documentation. Optional.

- title_par:

  `character`: Title of the plot. Optional.

- x_text_par:

  `character`: Text on the x-axis. Optional.

- plot.title.size_par:

  `numeric`: size of the plot title. Default is 20. See `size` in the
  `?element_text` documentation. Optional.

- ...:

  Other arguments to be passed to Plot. Optional.

## Value

`ggplot2` of the development factors
