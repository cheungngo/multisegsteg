# multisegsteg

<!-- badges: start -->
<!-- badges: end -->

The goal of multisegsteg is to provide tools for fitting, testing, and visualizing segmented regression models with multiple thresholds. It supports multi-segmented regression (continuous at thresholds) and multi-stegmented regression (allowing discontinuities at thresholds). Additionally, the package includes the segsteg.test function to perform statistical tests for detecting change points, offering both single and multiple threshold testing with robust p-value calculations.

## Installation

You can install the development version of multisegsteg from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("cheungngo/multisegsteg")
```

## Features

- **Model Fitting**: Fit segmented (`multi_segmented`) and stegmented (`multi_stegmented`) regression models with multiple thresholds.
- **Change Point Testing**: Test for the presence and significance of change points using `segsteg.test`, supporting both single and multiple thresholds.
- **Visualization**: Plot fitted models to visualize thresholds and regression lines with `plot_multi_segmented` and `plot_multi_stegmented`.

## Example: Fitting Segmented Models

Here's a basic example to demonstrate fitting segmented regression models with multiple thresholds:

```r
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

```r
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

## Example: Testing for Change Points with segsteg.test

The `segsteg.test` function allows you to test for the presence of change points in your regression model, providing statistical evidence for whether thresholds exist and whether additional thresholds improve the model fit. It supports both single and multiple threshold testing, with p-values calculated via Monte Carlo simulation (for single thresholds) or sequential likelihood ratio tests (for multiple thresholds).

### Single Threshold Test

To test for a single change point using the same dataset:

```r
# Test for a single threshold
single_test <- segsteg.test(
  formula.null = y ~ z,
  formula.chngpt = ~ x,
  family = "gaussian",
  data = data,
  type = "hinge",
  n_thresholds = 1,
  lb.quantile = 0.1,
  ub.quantile = 0.9,
  chngpts.cnt = 20,
  test.statistic = "lr",
  p.val.method = "MC",
  verbose = TRUE
)

# View results
print(single_test$chngpt)   # Estimated change point
print(single_test$p.value)  # P-value for presence of a change point
```

**What it does**: Tests whether there is at least one change point in the relationship between y and x, adjusting for z. The p-value indicates if the change point is statistically significant.

### Multiple Threshold Test

To test for multiple change points (e.g., two thresholds) and assess their significance:

```r
# Test for two thresholds
multi_test <- segsteg.test(
  formula.null = y ~ z,
  formula.chngpt = ~ x,
  family = "gaussian",
  data = data,
  type = "hinge",
  n_thresholds = 2,
  lb.quantile = 0.1,
  ub.quantile = 0.9,
  chngpts.cnt = 20,
  verbose = TRUE
)

# View results
print(multi_test$thresholds)  # Estimated thresholds
print(multi_test$p_values)    # Sequential p-values for each threshold
```

**What it does**: Fits models with 0, 1, and 2 thresholds, selecting the best thresholds for each. It then performs sequential likelihood ratio tests (LRTs) to calculate p-values:
- p_values[1]: Tests 1 vs. 0 thresholds (is there at least one change point?).
- p_values[2]: Tests 2 vs. 1 thresholds (does a second change point improve the fit?).

**P-Value Calculation**:
- For each step k, computes the LRT statistic: 2 * (loglik_k - loglik_{k-1}).
- Uses a chi-squared distribution with 1 degree of freedom (for "hinge") to compute the p-value: pchisq(lr_stat, df = 1, lower.tail = FALSE).
- A small p-value (< 0.05) suggests the additional threshold is significant.

### Stegmented Model Test

For the stegmented data (allowing discontinuities):

```r
# Test for two thresholds in stegmented model
multi_test_step <- segsteg.test(
  formula.null = y ~ z,
  formula.chngpt = ~ x,
  family = "gaussian",
  data = data_step,
  type = "stegmented",
  n_thresholds = 2,
  lb.quantile = 0.1,
  ub.quantile = 0.9,
  chngpts.cnt = 20,
  verbose = TRUE
)

# View results
print(multi_test_step$thresholds)
print(multi_test_step$p_values)
```

**Note**: For "stegmented" models, the LRT uses 2 degrees of freedom per threshold due to the step and segmented terms.

## Why Use segsteg.test?

The `segsteg.test` function complements the fitting functions (`multi_segmented`, `multi_stegmented`) by providing rigorous statistical tests to:

- Confirm the presence of change points.
- Determine if additional thresholds are justified.
- Support both continuous (hinge) and discontinuous (stegmented) models.

It's particularly useful when you want to validate the thresholds estimated by the fitting functions, ensuring your model is statistically sound.

## Notes

- **Data Requirements**: Ensure your data has sufficient variation in the change point variable (x in the examples) to detect thresholds reliably.
- **Computational Cost**: Testing multiple thresholds can be intensive, especially with many candidate thresholds (chngpts.cnt).
- **Interpretation**: For multiple thresholds, interpret p-values sequentially. A non-significant p-value at any step suggests stopping at the previous number of thresholds.
