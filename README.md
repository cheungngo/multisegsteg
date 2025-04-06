# multisegsteg

<!-- badges: start -->
<!-- badges: end -->

The goal of multisegsteg is to implement methods for fitting and visualizing segmented regression models with multiple thresholds, including multi-segmented regression (continuous at thresholds) and multi-stegmented regression (allowing discontinuities).

## Installation

You can install the development version of multisegsteg from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("cheungngo/multisegsteg")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(multisegsteg)

# Generate example data
set.seed(123)
n <- 200
x <- runif(n, 0, 10)
z <- rnorm(n)
y <- 1 + 0.5 * z + 0 * x + 1 * pmax(x - 3, 0) + 2 * pmax(x - 7, 0) + rnorm(n)
data <- data.frame(y = y, x = x, z = z)

# Fit model with two thresholds
result <- multi_segmented(
  formula.1 = y ~ z,
  formula.2 = ~ x,
  family = "gaussian",
  data = data,
  n_thresholds = 2
)

# View the estimated thresholds
print(result$thresholds)
#> [1] 3.049755 6.954196

# Plot the fitted segmented regression line
plot_multi_segmented(result, data, formula.1 = y ~ z, formula.2 = ~ x, plot.type = "fit")
```

For models with discontinuities (stegmented regression):

``` r
# Generate data with both steps and slope changes
y_step <- 1 + 0.5 * z + 0 * x + 2 * (x > 3) + 1 * pmax(x - 3, 0) +
     3 * (x > 7) + 2 * pmax(x - 7, 0) + rnorm(n)
data_step <- data.frame(y = y_step, x = x, z = z)

# Fit stegmented model
steg_result <- multi_stegmented(
  formula.1 = y ~ z,
  formula.2 = ~ x,
  family = "gaussian",
  data = data_step,
  n_thresholds = 2
)

# Visualize the result
plot_multi_stegmented(steg_result, data_step, formula.1 = y ~ z, formula.2 = ~ x)
```
