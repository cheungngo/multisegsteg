#' multisegsteg: Regression Models with Multiple Break-Points Estimation
#'
#' @description
#' Implements methods for fitting and visualizing segmented regression models
#' with multiple thresholds, including multi-segmented regression (continuous
#' at thresholds) and multi-stegmented regression (allowing discontinuities).
#'
#' @details
#' The package provides tools for:
#' \itemize{
#'   \item Estimation of multiple thresholds in segmented and stegmented regression models
#'   \item Support for both Gaussian and binomial outcomes
#'   \item Grid search optimization for threshold locations
#'   \item Cluster-robust standard errors for correlated data
#'   \item Visualization of model fits and log-likelihood surfaces
#' }
#'
#' @docType package
#' @name multisegsteg-package
#' @aliases multisegsteg
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_vline geom_tile annotate scale_fill_gradient labs theme_minimal theme_void
#' @importFrom gridExtra grid.arrange
#' @importFrom sandwich vcovCL
#' @importFrom stats model.frame model.matrix logLik pmax predict quantile
NULL
