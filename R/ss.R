#' Fit Segmented or Stegmented Models with Multiple Thresholds
#'
#' This function fits either a segmented or a "stegmented" model (a combination of step and segmented terms)
#' for multiple numbers of thresholds, selecting the best thresholds via grid search over quantiles of the
#' change point variable. It supports both Gaussian and binomial families and optionally computes
#' cluster-robust standard errors.
#'
#' @param formula.1 A formula specifying the outcome and base covariates (e.g., \code{y ~ x1 + x2}).
#' @param formula.2 A formula specifying the change point variable (e.g., \code{~ z}).
#' @param family A character string specifying the family of the model. Must be either \code{"gaussian"} or \code{"binomial"}.
#' @param data A data frame containing the variables specified in the formulas.
#' @param n_thresholds An integer or vector of integers specifying the number(s) of thresholds to fit.
#'   - If a single integer, fits models for 1 up to that number (e.g., \code{3} fits for 1, 2, and 3 thresholds).
#'   - If a vector, fits models for each specified number (e.g., \code{c(1, 3)} fits for 1 and 3 thresholds).
#' @param lb.quantile A numeric value between 0 and 1 specifying the lower bound quantile for threshold candidates. Default is 0.05.
#' @param ub.quantile A numeric value between 0 and 1 specifying the upper bound quantile for threshold candidates. Default is 0.95.
#' @param grid.search.max An integer specifying the maximum number of threshold candidates to consider. Default is \code{Inf}, meaning all unique quantiles are used.
#'   - For \code{n_thresholds = 1}, \code{grid.search.max} is set to 20.
#' @param cluster_var Optional. A character string specifying the name of the clustering variable for computing cluster-robust standard errors.
#' @param verbose Logical. If \code{TRUE}, prints detailed output for each threshold combination. Default is \code{FALSE}.
#' @param model_type A character string specifying the type of model to fit. Must be either \code{"segmented"} or \code{"stegmented"}.
#'   - \code{"segmented"}: Includes only segmented terms \code{(x - e_i)+}.
#'   - \code{"stegmented"}: Includes both step terms \code{I(x > e_i)} and segmented terms \code{(x - e_i)+}.
#' @param ... Additional arguments passed to the fitting functions (\code{lm} or \code{glm}).
#'
#' @return A list of model results, where each element corresponds to a specific number of thresholds.
#'   The list is named by the number of thresholds requested (e.g., "1", "2", etc.). Each element contains:
#'   \itemize{
#'     \item \code{n_thresholds_requested}: The number of thresholds requested for this model.
#'     \item \code{n_thresholds_used}: The number of thresholds actually used in the model (may be less than requested if data constraints apply).
#'     \item \code{thresholds}: The best thresholds selected based on log-likelihood.
#'     \item \code{fit}: The fitted model object (\code{lm} or \code{glm}).
#'     \item \code{loglik}: The log-likelihood of the best model.
#'     \item \code{all_logliks}: A vector of log-likelihoods for all threshold combinations evaluated.
#'     \item \code{all_combinations}: A matrix of all threshold combinations evaluated.
#'     \item \code{vcov_cluster}: The cluster-robust variance-covariance matrix, if \code{cluster_var} is specified.
#'     \item \code{se_cluster}: The cluster-robust standard errors, if \code{cluster_var} is specified.
#'   }
#'
#' @details
#' The function performs a grid search over combinations of thresholds within the specified quantile range of the change point variable.
#' For each number of thresholds, it selects the combination that maximizes the log-likelihood.
#' If the requested number of thresholds exceeds the number of unique thresholds available, the function reduces the number of thresholds
#' used and issues a warning. Only results where the number of thresholds used equals the number requested are included in the output.
#'
#' @note
#' The \code{sandwich} package is required for computing cluster-robust standard errors and must be installed.
#'
#' @examples
#' # Fit a segmented model for 1 to 3 thresholds
#' result <- multiss(formula.1 = y ~ x1 + x2, formula.2 = ~ z, family = "gaussian", data = mydata, n_thresholds = 3)
#'
#' # Fit a stegmented model for specific numbers of thresholds
#' result <- multiss(formula.1 = y ~ x1 + x2, formula.2 = ~ z, family = "binomial", data = mydata, n_thresholds = c(1, 3), model_type = "stegmented")
#'
#' @import sandwich
#' @export
multi_ss <- function(formula.1, formula.2, family, data, n_thresholds = 1,
                    lb.quantile = 0.05, ub.quantile = 0.95, grid.search.max = Inf,
                    cluster_var = NULL, verbose = FALSE, model_type = "segmented", ...) {
  # Validate and process n_thresholds
  if (!is.numeric(n_thresholds)) {
    stop("n_thresholds must be numeric")
  }

  # Convert to integer
  n_thresholds <- as.integer(n_thresholds)

  # Determine the sequence of thresholds to compute
  if (length(n_thresholds) == 1) {
    # If a single value, compute all thresholds from 1 to n_thresholds
    if (n_thresholds < 1) {
      stop("n_thresholds must be at least 1")
    }
    thresh_seq <- 1:n_thresholds
  } else {
    # If a vector, use the provided values directly
    if (any(n_thresholds < 1)) {
      stop("All values in n_thresholds must be at least 1")
    }
    thresh_seq <- n_thresholds
  }

  # Initialize results list
  results <- list()

  # Loop over each number of thresholds
  for (n in thresh_seq) {
    result <- fit_for_n_thresholds(n, formula.1, formula.2, family, data, lb.quantile, ub.quantile, grid.search.max, cluster_var, verbose, model_type, ...)
    results[[as.character(n)]] <- result
  }

  # Filter results to only include where n_thresholds_used equals n_thresholds_requested
  results <- results[sapply(results, function(x) x$n_thresholds_used == x$n_thresholds_requested)]

  return(results)
}

#' Plot Segmented or Stegmented Model Results
#'
#' @description
#' Visualizes the results of segmented or stegmented regression models, including
#' the data with the fitted model and/or the log-likelihood surface for the thresholds.
#'
#' @param result The output from a segmented or stegmented model fitting function,
#' containing thresholds, the fitted model, and log-likelihood information.
#' @param data The original dataset used for model fitting.
#' @param formula.1 The formula for the linear part of the model.
#' @param formula.2 The formula specifying the change point variable.
#' @param model.type A string indicating the type of model: "segmented" (default) or "stegmented".
#' @param plot.type A string specifying which plot(s) to produce: "fit", "loglik", or "both" (default).
#' @param n_points The number of points to use for generating the smooth fitted line (default is 100).
#' @param ... Additional arguments to be passed to plotting functions.
#'
#' @details
#' This function leverages the `ggplot2` package for creating plots and the `gridExtra`
#' package for arranging multiple plots. It adapts its visualization based on the number
#' of thresholds in the model:
#' - **One threshold**: Plots the log-likelihood against the threshold value.
#' - **Two thresholds**: Generates a heatmap of the log-likelihood surface.
#' - **More than two thresholds**: Displays a text summary of the best thresholds and
#'   their corresponding log-likelihood.
#'
#' @return
#' Depending on the `plot.type` argument, the function returns:
#' - `"fit"`: A plot of the data with the fitted model.
#' - `"loglik"`: A plot of the log-likelihood surface or summary.
#' - `"both"`: Both the fit and log-likelihood plots arranged side by side using `grid.arrange`.
#'
#' @examples
#' # Example with a segmented model
#' # Assuming 'result_seg' is the model output and 'data' is the dataset
#' # plot_ss(result_seg, data, formula.1 = y ~ x, formula.2 = ~ x,
#' #         model.type = "segmented", plot.type = "both")
#'
#' # Example with a stegmented model
#' # plot_ss(result_steg, data, formula.1 = y ~ x, formula.2 = ~ x,
#' #         model.type = "stegmented", plot.type = "fit")
#'
#' @author Your Name
#' @date 2023-10-01
#'
#' @import ggplot2
#' @import gridExtra
#'
#' @export
plot_ss <- function(result, data, formula.1, formula.2, model.type = "segmented",
                    plot.type = "both", n_points = 100, ...) {
  # Function body remains unchanged
}
plot_ss <- function(result, data, formula.1, formula.2, model.type = "segmented",
                    plot.type = "both", n_points = 100, ...) {
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

  # Construct design matrix based on model.type
  if (model.type == "segmented") {
    X_threshold_smooth <- sapply(thresholds, function(e) pmax(x_smooth - e, 0))
    X_smooth <- cbind(Z_pred, X_threshold_smooth)
    colnames(X_smooth) <- c(colnames(Z), paste0("X", seq_len(n_thresholds)))
  } else if (model.type == "stegmented") {
    X_step_smooth <- sapply(thresholds, function(e) as.numeric(x_smooth > e))
    X_seg_smooth <- sapply(thresholds, function(e) pmax(x_smooth - e, 0))
    X_smooth <- cbind(Z_pred, X_step_smooth, X_seg_smooth)
    colnames(X_smooth) <- c(colnames(Z),
                            paste0("step_", seq_len(n_thresholds)),
                            paste0("seg_", seq_len(n_thresholds)))
  } else {
    stop("model.type must be 'segmented' or 'stegmented'")
  }

  # Predict y values using the fitted model
  y_pred <- predict(fit, newdata = list(X = X_smooth))

  # Prepare data frames for plotting
  plot_data <- data.frame(x = x, y = y)
  fit_data <- data.frame(x = x_smooth, y = y_pred)

  # Plot 1: Data with Fitted Line and Thresholds
  p1 <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    geom_line(data = fit_data, aes(x = x, y = y), color = "blue", size = 1) +
    geom_vline(xintercept = thresholds, linetype = "dashed", color = "red") +
    labs(title = paste(n_thresholds, "Threshold(s)", model.type, "Model Fit"),
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
      geom_tile() +
      geom_point(data = data.frame(t1 = thresholds[1], t2 = thresholds[2], loglik = result$loglik),
                 aes(x = t1, y = t2), color = "red", size = 3, shape = 4) +
      scale_fill_gradient(low = "blue", high = "yellow") +
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
