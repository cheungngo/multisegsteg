#' Fit a multiple threshold step-segmented regression model
#'
#' @description
#' Fits a step-segmented (stegmented) regression model with multiple thresholds
#' that allows for both discontinuities (steps) and changes in slope at the threshold
#' points. The function employs a grid search approach to identify optimal threshold
#' values that maximize the model's log-likelihood.
#'
#' @details
#' The stegmented regression model implemented here follows the form:
#' \deqn{y = X\beta + \sum_{j=1}^{k} \alpha_j I(z > e_j) + \sum_{j=1}^{k} \gamma_j (z - e_j)_+ + \epsilon}
#' where \eqn{I(z > e_j)} is an indicator function, \eqn{(z - e_j)_+ = \max(0, z - e_j)},
#' \eqn{e_j} are the threshold values, \eqn{k} is the number of thresholds, and \eqn{X}
#' represents covariates from formula.1.
#'
#' This model extends the standard segmented regression by allowing jumps (discontinuities)
#' at the threshold points in addition to changes in slope. The \eqn{\alpha_j} parameters
#' control the size of the jumps, while the \eqn{\gamma_j} parameters control the changes in slope.
#'
#' The function evaluates all combinations of thresholds within the specified quantile
#' range of the threshold variable and selects the combination that maximizes the
#' log-likelihood. For models with clustered data, cluster-robust standard errors
#' are calculated using the sandwich package.
#'
#' @param formula.1 Formula. Specifies the relationship between the response variable
#'        and covariates (excluding the threshold terms).
#' @param formula.2 Formula. Specifies the threshold variable (right-hand side only).
#' @param family Character. Distribution family for the response variable.
#'        Currently supports "gaussian" for continuous outcomes and "binomial" for binary outcomes.
#' @param data Data frame. Contains all variables in the model.
#' @param n_thresholds Integer. Number of thresholds to estimate. Default is 1.
#' @param lb.quantile Numeric. Lower bound quantile for candidate thresholds. Default is 0.05.
#' @param ub.quantile Numeric. Upper bound quantile for candidate thresholds. Default is 0.95.
#' @param grid.search.max Numeric. Maximum number of candidate threshold values
#'        to consider. Default is Inf (all unique values in the specified range).
#' @param cluster_var Character. Name of the cluster variable for computing
#'        cluster-robust standard errors. Default is NULL.
#' @param verbose Logical. Whether to print progress information during the grid search. Default is FALSE.
#' @param ... Additional arguments passed to internal fitting functions.
#'
#' @return A list containing the following components:
#'   \item{thresholds}{Numeric vector of estimated threshold values}
#'   \item{fit}{Model fit object (lm for gaussian family or glm for binomial family)}
#'   \item{loglik}{Log-likelihood value for the optimal model}
#'   \item{all_logliks}{Numeric vector of log-likelihood values for all threshold combinations evaluated}
#'   \item{all_combinations}{Matrix of all threshold combinations evaluated}
#'   \item{vcov_cluster}{Variance-covariance matrix with cluster-robust standard errors (if cluster_var provided)}
#'   \item{se_cluster}{Cluster-robust standard errors (if cluster_var provided)}
#'
#' @examples
#'
#' # Generate example data with both steps and slope changes
#' set.seed(123)
#' n <- 200
#' x <- runif(n, 0, 10)
#' z <- rnorm(n)
#' # Data with steps at x=3 and x=7, and changes in slope at the same points
#' y <- 1 + 0.5 * z + 0 * x + 2 * (x > 3) + 1 * pmax(x - 3, 0) +
#'      3 * (x > 7) + 2 * pmax(x - 7, 0) + rnorm(n)
#' data <- data.frame(y = y, x = x, z = z, cluster = rep(1:20, each = 10))
#'
#' # Fit stegmented model with two thresholds
#' result <- multi_stegmented(
#'   formula.1 = y ~ z,
#'   formula.2 = ~ x,
#'   family = "gaussian",
#'   data = data,
#'   n_thresholds = 2
#' )
#'
#' # Print the estimated thresholds
#' print(result$thresholds)
#'
#' # Fit with cluster-robust standard errors
#' result_clustered <- multi_stegmented(
#'   formula.1 = y ~ z,
#'   formula.2 = ~ x,
#'   family = "gaussian",
#'   data = data,
#'   n_thresholds = 2,
#'   cluster_var = "cluster"
#' )
#'
#' # Access cluster-robust standard errors
#' print(result_clustered$se_cluster)
#'
#'
#' @seealso \code{\link{plot_multi_stegmented}} for visualizing the results
#' @seealso \code{\link{multi_segmented}} for fitting continuous segmented regression models
#'
#'
#' @importFrom sandwich vcovCL
#'
#' @export
multi_stegmented <- function(formula.1, formula.2, family, data, n_thresholds = 1,
                             lb.quantile = 0.05, ub.quantile = 0.95, grid.search.max = Inf,
                             cluster_var = NULL, verbose = FALSE, ...) {
  requireNamespace("sandwich", quietly = TRUE)

  # Extract outcome variable
  y <- model.frame(formula.1, data)[, 1]

  # Extract change point variable
  chngpt_var_name <- all.vars(formula.2)[1]
  chngpt_var <- data[[chngpt_var_name]]

  # Generate candidate thresholds based on quantiles
  chngpts <- unique(quantile(chngpt_var,
                             seq(lb.quantile, ub.quantile,
                                 length.out = min(grid.search.max, length(unique(chngpt_var))))))

  # Generate all ordered combinations of n_thresholds distinct thresholds
  comb <- t(combn(chngpts, n_thresholds))
  n_comb <- nrow(comb)

  # Initialize vector to store log-likelihoods
  logliks <- numeric(n_comb)

  # Base design matrix without threshold terms
  Z <- model.matrix(formula.1, data)

  # Evaluate each combination
  for (i in 1:n_comb) {
    thresholds <- comb[i, ]

    # Create step and segmented terms for each threshold
    X_step <- sapply(thresholds, function(e) as.numeric(chngpt_var > e))
    X_seg <- sapply(thresholds, function(e) pmax(chngpt_var - e, 0))

    # Combine base covariates with step and segmented terms
    X <- cbind(Z, X_step, X_seg)

    # Fit the model based on family
    if (family == "gaussian") {
      fit <- lm(y ~ X - 1)  # -1 because X includes intercept if in formula.1
      logliks[i] <- logLik(fit)
    } else if (family == "binomial") {
      fit <- glm(y ~ X - 1, family = binomial)
      logliks[i] <- logLik(fit)
    } else {
      stop("Only 'gaussian' and 'binomial' families are supported")
    }

    # Optional verbose output
    if (verbose) {
      cat(sprintf("Combination %d, thresholds: %s, log-likelihood: %.4f\n",
                  i, paste(thresholds, collapse = ", "), logliks[i]))
    }
  }

  # Find the best combination
  best_idx <- which.max(logliks)
  best_thresholds <- comb[best_idx, ]

  # Fit the final model with the best thresholds
  X_step <- sapply(best_thresholds, function(e) as.numeric(chngpt_var > e))
  X_seg <- sapply(best_thresholds, function(e) pmax(chngpt_var - e, 0))
  X <- cbind(Z, X_step, X_seg)
  if (family == "gaussian") {
    best_fit <- lm(y ~ X - 1)
  } else {
    best_fit <- glm(y ~ X - 1, family = binomial)
  }

  # Compute cluster-robust standard errors if cluster_var is provided
  if (!is.null(cluster_var)) {
    cluster <- data[[cluster_var]]
    vcov_cluster <- sandwich::vcovCL(best_fit, cluster = cluster)
    se_cluster <- sqrt(diag(vcov_cluster))
  } else {
    vcov_cluster <- NULL
    se_cluster <- NULL
  }

  # Return results as a list
  return(list(
    thresholds = best_thresholds,
    fit = best_fit,
    loglik = logliks[best_idx],
    all_logliks = logliks,
    all_combinations = comb,
    vcov_cluster = vcov_cluster,
    se_cluster = se_cluster
  ))
}

#' Plot the results of a multi-stegmented regression model
#'
#' @description
#' Creates visualizations for interpreting multi-stegmented regression models,
#' displaying both the fitted stegmented line (with discontinuities and slope changes)
#' and the log-likelihood surface or profile. For models with one threshold, the function
#' shows a profile of log-likelihood values. For models with two thresholds, it displays
#' a heatmap of the log-likelihood surface.
#'
#' @details
#' This function generates two types of plots:
#'
#' 1. **Fit Plot**: Shows the original data points and the fitted stegmented regression line,
#'    with vertical dashed lines indicating the estimated threshold locations. The plot
#'    illustrates both the discontinuities (steps) and changes in slope at each threshold.
#'
#' 2. **Log-likelihood Plot**:
#'    - For one threshold: Shows a profile of log-likelihood values across all candidate
#'      threshold values, with the optimal threshold marked.
#'    - For two thresholds: Displays a heatmap of the log-likelihood surface across
#'      all combinations of threshold values, with the optimal combination marked.
#'    - For more than two thresholds: Provides a text summary of the optimal thresholds
#'      and their log-likelihood value.
#'
#' The stegmented model differs from the standard segmented model by allowing
#' discontinuities at the threshold points in addition to changes in slope.
#'
#' @param result List. Output from the multi_stegmented function.
#' @param data Data frame. The original dataset used to fit the model.
#' @param formula.1 Formula. The same formula.1 used in the multi_stegmented function.
#' @param formula.2 Formula. The same formula.2 used in the multi_stegmented function.
#' @param plot.type Character. Type of plot to create:
#'        - "both" (default): Creates both the fit plot and log-likelihood plot
#'        - "fit": Creates only the fit plot with data points and fitted line
#'        - "loglik": Creates only the log-likelihood plot
#' @param n_points Integer. Number of points to use for generating the smooth fitted line.
#'        Higher values create a smoother curve. Default is 100.
#' @param ... Additional arguments passed to ggplot functions.
#'
#' @return If plot.type = "both", returns a grid of plots arranged using gridExtra.
#'         If plot.type = "fit" or "loglik", returns a single ggplot object.
#'
#' @examples
#'
#' # Generate example data with both steps and slope changes
#' set.seed(123)
#' n <- 200
#' x <- runif(n, 0, 10)
#' z <- rnorm(n)
#' # Data with steps at x=3 and x=7, and changes in slope at the same points
#' y <- 1 + 0.5 * z + 0 * x + 2 * (x > 3) + 1 * pmax(x - 3, 0) +
#'      3 * (x > 7) + 2 * pmax(x - 7, 0) + rnorm(n)
#' data <- data.frame(y = y, x = x, z = z)
#'
#' # Fit stegmented model with two thresholds
#' result <- multi_stegmented(
#'   formula.1 = y ~ z,
#'   formula.2 = ~ x,
#'   family = "gaussian",
#'   data = data,
#'   n_thresholds = 2
#' )
#'
#' library(ggplot2)
#' library(gridExtra)
#'
#' # Plot both the fitted line and log-likelihood surface
#' plot_multi_stegmented(result, data, formula.1 = y ~ z, formula.2 = ~ x)
#'
#' # Plot only the fitted line
#' fit_plot <- plot_multi_stegmented(
#'   result, data,
#'   formula.1 = y ~ z,
#'   formula.2 = ~ x,
#'   plot.type = "fit"
#' )
#' print(fit_plot)
#'
#' # Plot only the log-likelihood surface
#' loglik_plot <- plot_multi_stegmented(
#'   result, data,
#'   formula.1 = y ~ z,
#'   formula.2 = ~ x,
#'   plot.type = "loglik"
#' )
#' print(loglik_plot)
#'
#' @seealso \code{\link{multi_stegmented}} for fitting the stegmented regression model
#' @seealso \code{\link{plot_multi_segmented}} for plotting segmented regression models
#'
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_vline geom_tile geom_text
#'             scale_fill_gradient labs theme_minimal annotate theme_void
#' @importFrom gridExtra grid.arrange
#'
#' @export
plot_multi_stegmented <- function(result, data, formula.1, formula.2, plot.type = "both",
                                  n_points = 100, ...) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("gridExtra", quietly = TRUE)

  # Extract thresholds and fit
  thresholds <- result$thresholds
  n_thresholds <- length(thresholds)
  fit <- result$fit

  # Extract change point variable and response
  chngpt_var_name <- all.vars(formula.2)[1]
  x <- data[[chngpt_var_name]]
  y <- model.frame(formula.1, data)[, 1]

  # Generate smooth x values for plotting the fit
  x_range <- range(x)
  x_smooth <- seq(x_range[1], x_range[2], length.out = n_points)

  # Reconstruct the design matrix for prediction
  Z <- model.matrix(formula.1, data)
  Z_mean <- colMeans(Z)
  Z_pred <- matrix(rep(Z_mean, each = n_points), nrow = n_points,
                   dimnames = list(NULL, colnames(Z)))
  X_step_smooth <- sapply(thresholds, function(e) as.numeric(x_smooth > e))
  X_seg_smooth <- sapply(thresholds, function(e) pmax(x_smooth - e, 0))
  X_smooth <- cbind(Z_pred, X_step_smooth, X_seg_smooth)
  colnames(X_smooth) <- c(colnames(Z),
                          paste0("step_", seq_len(n_thresholds)),
                          paste0("seg_", seq_len(n_thresholds)))

  # Predict y values using the fitted model
  y_pred <- predict(fit, newdata = list(X = X_smooth))

  # Prepare data frames for plotting
  plot_data <- data.frame(x = x, y = y)
  fit_data <- data.frame(x = x_smooth, y = y_pred)

  # Plot 1: Data with Fitted Line and Thresholds
  p1 <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +  # Scatter plot of data
    geom_line(data = fit_data, aes(x = x, y = y), color = "blue", size = 1) +  # Fitted line
    geom_vline(xintercept = thresholds, linetype = "dashed", color = "red") +  # Thresholds
    labs(title = paste(n_thresholds, "Threshold(s) Stegmented Model Fit"),
         x = chngpt_var_name, y = "Response") +
    theme_minimal()

  # Plot 2: Log-likelihood Visualization
  if (n_thresholds == 2) {
    loglik_data <- data.frame(
      t1 = result$all_combinations[, 1],
      t2 = result$all_combinations[, 2],
      loglik = result$all_logliks
    )

    p2 <- ggplot(loglik_data, aes(x = t1, y = t2, fill = loglik)) +
      geom_tile() +  # Heatmap
      geom_point(data = data.frame(t1 = thresholds[1], t2 = thresholds[2], loglik = result$loglik),
                 aes(x = t1, y = t2), color = "red", size = 3, shape = 4) +  # Mark best thresholds
      scale_fill_gradient(low = "blue", high = "yellow") +  # Color gradient
      labs(title = "Log-Likelihood Surface",
           x = paste(chngpt_var_name, "Threshold 1"),
           y = paste(chngpt_var_name, "Threshold 2"),
           fill = "Log-Likelihood") +
      theme_minimal()
  } else if (n_thresholds == 1) {
    loglik_data <- data.frame(
      t = result$all_combinations,
      loglik = result$all_logliks
    )

    p2 <- ggplot(loglik_data, aes(x = t, y = loglik)) +
      geom_line() +
      geom_point(data = data.frame(t = thresholds, loglik = result$loglik),
                 aes(x = t, y = loglik), color = "red", size = 3) +
      labs(title = "Log-Likelihood vs Threshold",
           x = chngpt_var_name, y = "Log-Likelihood") +
      theme_minimal()
  } else {
    p2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Best thresholds:", paste(round(thresholds, 2), collapse = ", "),
                             "\nLog-Likelihood:", round(result$loglik, 2)),
               size = 5) +
      theme_void() +
      labs(title = "Log-Likelihood for Multiple Thresholds (>2)")
  }

  # Combine plots based on plot.type
  if (plot.type == "fit") {
    return(p1)
  } else if (plot.type == "loglik") {
    return(p2)
  } else if (plot.type == "both") {
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  } else {
    stop("plot.type must be 'fit', 'loglik', or 'both'")
  }
}
