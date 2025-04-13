#' @title Change Point Test for Regression Models
#'
#' @description
#' Performs a test for change points in regression models, supporting single or multiple thresholds.
#' For multiple thresholds, it uses a sequential likelihood ratio test (LRT) to assess the significance
#' of adding each additional threshold. P-values are calculated using the chi-squared distribution
#' based on the LRT statistic.
#'
#' @param formula.null Formula for the null model (without change points).
#' @param formula.chngpt Formula specifying the change point variable.
#' @param family Family of the response variable ("binomial" or "gaussian").
#' @param data Data frame containing the variables.
#' @param type Type of change point model ("step", "hinge", "segmented", "stegmented").
#' @param test.statistic Test statistic for single threshold ("lr" or "score").
#' @param chngpts Optional vector of candidate thresholds.
#' @param lb.quantile Lower quantile for generating candidate thresholds.
#' @param ub.quantile Upper quantile for generating candidate thresholds.
#' @param chngpts.cnt Number of candidate thresholds to generate.
#' @param prec.weights Optional precision weights.
#' @param p.val.method Method for p-value calculation ("MC" or "param.boot").
#' @param mc.n Number of Monte Carlo simulations.
#' @param boot.B Number of bootstrap iterations.
#' @param robust Use robust standard errors.
#' @param keep.fits Keep model fits in the output.
#' @param verbose Print verbose output.
#' @param n_thresholds Number of thresholds to test.
#' @param cluster_var Optional cluster variable for robust SEs.
#'
#' @return
#' For \code{n_thresholds = 1}, a list of class "chngpt.test" containing:
#' \itemize{
#'   \item \code{chngpt}: Estimated change point.
#'   \item \code{statistic}: Test statistic.
#'   \item \code{p.value}: P-value.
#'   \item \code{loglik}: Log-likelihood of the alternative model.
#'   \item \code{fit.null}: Null model fit (if \code{keep.fits = TRUE}).
#'   \item \code{fit.alt}: Alternative model fit (if \code{keep.fits = TRUE}).
#' }
#' For \code{n_thresholds > 1}, a list of class "multi.segsteg.test" containing:
#' \itemize{
#'   \item \code{thresholds}: Estimated thresholds.
#'   \item \code{fit}: Final model fit with \code{n_thresholds} thresholds.
#'   \item \code{loglik}: Log-likelihood of the final model.
#'   \item \code{all_logliks}: Log-likelihoods for models with 0 to \code{n_thresholds} thresholds.
#'   \item \code{p_values}: Sequential LRT p-values for each additional threshold.
#' }
#'
#' @details
#' For single thresholds (\code{n_thresholds = 1}), the function tests for a change point using
#' likelihood ratio or score tests, with p-values computed via Monte Carlo or bootstrap methods.
#' For multiple thresholds (\code{n_thresholds > 1}), it supports "hinge" and "stegmented" types,
#' using a sequential LRT approach:
#' \enumerate{
#'   \item Fits a null model (0 thresholds).
#'   \item Fits models with 1 to \code{n_thresholds} thresholds, selecting the best combination.
#'   \item Computes LRT statistic: \code{2 * (loglik_k - loglik_{k-1})}.
#'   \item Uses chi-squared distribution for p-values (df = 1 for "hinge", 2 for "stegmented").
#' }
#'
#' @examples
#' # Generate sample data with two hinge thresholds
#' set.seed(123)
#' n <- 100
#' x <- runif(n, 0, 10)
#' z <- rnorm(n)
#' true_thresholds <- c(3, 7)
#' y <- 2 + 0.5 * z + 1 * (x > true_thresholds[1]) * (x - true_thresholds[1]) +
#'      1.5 * (x > true_thresholds[2]) * (x - true_thresholds[2]) + rnorm(n)
#' data <- data.frame(y = y, x = x, z = z)
#'
#' # Define formulas
#' formula.null <- y ~ z
#' formula.chngpt <- ~ x
#'
#' # Test for two thresholds
#' result_multi <- segsteg.test(
#'   formula.null = formula.null,
#'   formula.chngpt = formula.chngpt,
#'   family = "gaussian",
#'   data = data,
#'   type = "hinge",
#'   n_thresholds = 2,
#'   lb.quantile = 0.1,
#'   ub.quantile = 0.9,
#'   chngpts.cnt = 20,
#'   verbose = TRUE
#' )
#'
#' # Print results
#' print(result_multi$thresholds)  # Estimated thresholds
#' print(result_multi$p_values)    # P-values for each step
#'
#' @export
segsteg.test <- function(formula.null, formula.chngpt, family = c("binomial", "gaussian"),
                         data, type = c("step", "hinge", "segmented", "stegmented"),
                         test.statistic = c("lr", "score"), chngpts = NULL,
                         lb.quantile = 0.1, ub.quantile = 0.9, chngpts.cnt = 50,
                         prec.weights = NULL, p.val.method = c("MC", "param.boot"),
                         mc.n = 50000, boot.B = 10000, robust = FALSE, keep.fits = FALSE,
                         verbose = FALSE, n_thresholds = 1, cluster_var = NULL) {
  requireNamespace("sandwich", quietly = TRUE)
  requireNamespace("MASS", quietly = TRUE)  # For mvrnorm in MC p-value

  # Initial setup
  DNAME <- deparse(substitute(data))
  family <- match.arg(family)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  p.val.method <- match.arg(p.val.method)

  # Initialize .Random.seed if it doesn’t exist
  if (!exists(".Random.seed")) runif(1)
  save.seed <- .Random.seed
  set.seed(1)

  # Handle missing data
  subset.1 <- complete.cases(model.frame(formula.null, data, na.action = na.pass))
  subset.2 <- complete.cases(model.frame(formula.chngpt, data, na.action = na.pass))
  data <- data[subset.1 & subset.2, , drop = FALSE]

  # Extract outcome and change point variable
  y <- model.frame(formula.null, data)[, 1]
  chngpt_var_name <- all.vars(formula.chngpt)[1]
  chngpt_var <- data[[chngpt_var_name]]

  # Base design matrix
  Z <- model.matrix(formula.null, data)
  n <- nrow(Z)
  p.null <- ncol(Z)

  # Handle precision weights
  if (is.null(prec.weights)) prec.weights <- rep(1, n)
  data$prec.weights <- prec.weights
  prec.w.all.one <- all(prec.weights == 1)

  # Generate candidate thresholds
  if (is.null(chngpts)) {
    chngpts <- unique(quantile(chngpt_var, seq(lb.quantile, ub.quantile, length.out = chngpts.cnt)))
  }
  chngpts <- chngpts[chngpts > min(chngpt_var) & chngpts < max(chngpt_var)]
  M <- length(chngpts)

  if (n_thresholds == 1) {
    # --- Single-Threshold Logic ---

    # Define formulas
    tmp <- as.matrix(model.frame(formula.chngpt, data))
    col.Z <- colnames(Z)
    z.1.name <- intersect(colnames(tmp), col.Z)
    has.itxn <- length(z.1.name) > 0
    f.null <- if (type %in% c("segmented", "stegmented"))
      update(formula.null, as.formula(paste("~ . +", chngpt_var_name)))
    else formula.null
    f.alt <- update(f.null, get.f.alt(type, chngpt_var_name,
                                      modified.by = if (has.itxn) z.1.name else NULL))

    # Model dimensions
    p.alt <- switch(type, step = 1, hinge = 1, segmented = 1, stegmented = 2)
    p.alt.2 <- p.alt * ifelse(has.itxn, 2, 1)
    do.score <- p.alt.2 == 1 & test.statistic == "score"

    # Fit null model
    fit.null <- glm(formula = f.null, data = data, family = family, weights = prec.weights)
    glm.warn <- FALSE

    # Compute components needed for W.null and V.S.hat
    linear.predictors.null <- fit.null$linear.predictors
    mu.h <- fit.null$fitted.values
    D.h <- if (family == "binomial") {
      diag(mu.h * (1 - mu.h))
    } else {
      diag(n) * summary(fit.null)$dispersion
    }
    DW <- if (prec.w.all.one) D.h else diag(diag(D.h) * prec.weights)
    tryCatch({
      V.eta.h <- Z %*% solve(t(Z) %*% DW %*% Z) %*% t(Z)
    }, error = function(e) {
      stop("Failed to compute V.eta.h due to singular matrix. Check data for collinearity.")
    })
    R.h <- diag(n) - D.h %*% V.eta.h %*% diag(prec.weights)
    RDR <- if (robust) {
      V.y <- diag(resid(fit.null, type = "response")^2)
      R.h %*% V.y %*% t(R.h)
    } else {
      R.h %*% D.h %*% (if (prec.w.all.one) diag(n) else t(R.h))
    }

    if (do.score) {
      # Score test
      W.null <- matrix(0, nrow = n, ncol = p.alt * M)
      for (m in 1:M) {
        W.null[, 1:p.alt + p.alt * (m - 1)] <- make.chngpt.var(chngpt_var, chngpts[m], type, data)[, "chngpt.var", drop = FALSE]
      }
      B <- if (prec.w.all.one) W.null else diag(prec.weights) %*% W.null
      V.S.hat <- t(B) %*% RDR %*% B
      TT.0 <- c((t(W.null) %*% (y - mu.h)) / sqrt(diag(V.S.hat)))
      TT <- abs(TT.0)
      T.max <- max(TT, na.rm = TRUE)
      max.id <- which.max(TT)
      chngpt <- chngpts[max.id]
      Q.max <- T.max
      method <- "Maximum of Score Statistics"
      data <- make.chngpt.var(chngpt_var, chngpt, type, data)
      fit.alt <- suppressWarnings(glm(f.alt, data = data, family = family, weights = prec.weights))
    } else {
      # Likelihood ratio test
      W.null <- matrix(0, nrow = n, ncol = p.alt.2 * M)
      QQ <- numeric(M)
      valid.thresholds <- rep(TRUE, M)
      for (m in 1:M) {
        data_m <- make.chngpt.var(chngpt_var, chngpts[m], type, data)
        X <- model.matrix(f.alt, data_m)
        I.a <- t(X) %*% DW %*% X
        tryCatch({
          I.bb.a <- I.a[p.null + 1:p.alt.2, p.null + 1:p.alt.2] -
            I.a[p.null + 1:p.alt.2, 1:p.null] %*% solve(I.a[1:p.null, 1:p.null], I.a[1:p.null, p.null + 1:p.alt.2])
          I.bb.a.inv.sqrt <- if (p.alt.2 == 1) {
            if (I.bb.a <= 0) stop("Non-positive I.bb.a")
            I.bb.a^(-1/2)
          } else {
            eig <- eigen(solve(I.bb.a))
            if (any(eig$values <= 0)) stop("Non-positive eigenvalues")
            eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
          }
          W.null[, 1:p.alt.2 + p.alt.2 * (m - 1)] <- X[, p.null + 1:p.alt.2, drop = FALSE] %*% I.bb.a.inv.sqrt
        }, error = function(e) {
          if (verbose) message("Warning: Failed to compute W.null for threshold ", chngpts[m], ": ", e$message)
          glm.warn <<- TRUE
          valid.thresholds[m] <<- FALSE
        })
        fit.alt_m <- suppressWarnings(glm(f.alt, data = data_m, family = family, weights = prec.weights))
        if (any(is.na(coef(fit.alt_m)))) {
          glm.warn <- TRUE
          QQ[m] <- NA
        } else {
          QQ[m] <- fit.null$aic - fit.alt_m$aic + 2 * p.alt.2
        }
      }
      if (all(!valid.thresholds)) stop("W.null could not be computed for any threshold. Check data or model specification.")
      Q.max <- max(QQ, na.rm = TRUE)
      max.id <- which.max(QQ)
      chngpt <- chngpts[max.id]
      data <- make.chngpt.var(chngpt_var, chngpt, type, data)
      fit.alt <- glm(f.alt, data = data, family = family, weights = prec.weights)
      method <- "Maximum of Likelihood Ratio Statistics"
      B <- if (prec.w.all.one) W.null else diag(prec.weights) %*% W.null
      V.S.hat <- t(B) %*% RDR %*% B
    }

    # P-value computation
    if (p.val.method == "MC") {
      sam <- MASS::mvrnorm(mc.n, rep(0, ncol(W.null)), cov2cor(V.S.hat))
      if (do.score) {
        sam <- abs(sam)
        x.max <- apply(sam, 1, max)
        p.value <- mean(x.max > T.max, na.rm = TRUE)
      } else {
        sam <- sam^2
        x.max <- apply(sam, 1, function(aux) max(colSums(matrix(aux, nrow = p.alt.2))))
        p.value <- mean(x.max > Q.max, na.rm = TRUE)
      }
    } else if (p.val.method == "param.boot" && !do.score) {
      sd.null <- sqrt(summary(fit.null)$dispersion)
      linear.predictors.null.sorted <- linear.predictors.null[order(chngpt_var)]
      Q.max.boot <- sapply(1:boot.B, function(b) {
        y.b <- rnorm(n, linear.predictors.null.sorted, sd.null)
        llik.null.b <- logLik(lm(y.b ~ Z - 1))
        fit.alt.b <- suppressWarnings(segsteg.test(
          formula.null, formula.chngpt, family, data, type,
          test.statistic = "lr", chngpts = chngpts, prec.weights = prec.weights,
          p.val.method = "none", n_thresholds = 1, keep.fits = TRUE
        ))
        2 * (logLik(fit.alt.b$fit.alt) - llik.null.b)
      })
      p.value <- mean(Q.max.boot > Q.max, na.rm = TRUE)
    } else {
      p.value <- NA
    }
    .Random.seed <- save.seed

    # Cluster-robust SEs
    if (!is.null(cluster_var)) {
      cluster <- data[[cluster_var]]
      vcov_cluster <- sandwich::vcovCL(fit.alt, cluster = cluster)
      se_cluster <- sqrt(diag(vcov_cluster))
    } else {
      vcov_cluster <- NULL
      se_cluster <- NULL
    }

    # Results
    res <- list(
      chngpt = chngpt,
      statistic = Q.max,
      p.value = p.value,
      chngpts = chngpts,
      method = method,
      loglik = logLik(fit.alt),
      fit.null = if (keep.fits) fit.null else NULL,
      fit.alt = if (keep.fits) fit.alt else NULL,
      glm.warn = glm.warn,
      vcov_cluster = vcov_cluster,
      se_cluster = se_cluster,
      data.name = DNAME
    )
    class(res) <- c("chngpt.test", "htest")

  } else {
    # --- Multiple-Threshold Logic ---
    if (type != "hinge" && type != "stegmented") {
      stop("For n_thresholds > 1, only 'hinge' and 'stegmented' types are supported")
    }

    # Fit null model (0 thresholds)
    null_model <- glm(y ~ Z, family = switch(family, "gaussian" = gaussian(),
                                             "binomial" = binomial()),
                      weights = prec.weights)
    loglik_null <- logLik(null_model)

    # Initialize lists to store results
    fits <- list()
    logliks <- numeric(n_thresholds + 1)
    logliks[1] <- loglik_null
    thresholds_list <- list(NULL)
    all_combinations_list <- list(NULL)

    # Fit models with 1 to n_thresholds thresholds
    for (k in 1:n_thresholds) {
      comb <- t(combn(chngpts, k))
      n_comb <- nrow(comb)
      logliks_k <- numeric(n_comb)
      for (i in 1:n_comb) {
        thresholds <- comb[i, ]
        if (type == "hinge") {
          X_threshold <- sapply(thresholds, function(e) pmax(chngpt_var - e, 0))
        } else {
          X_step <- sapply(thresholds, function(e) as.numeric(chngpt_var > e))
          X_seg <- sapply(thresholds, function(e) pmax(chngpt_var - e, 0))
          X_threshold <- cbind(X_step, X_seg)
        }
        X <- cbind(Z, X_threshold)
        fit <- glm(y ~ 0 + X, family = switch(family, "gaussian" = gaussian(),
                                              "binomial" = binomial()),
                   weights = prec.weights)
        logliks_k[i] <- logLik(fit)
        if (verbose) {
          cat(sprintf("Combination %d, k=%d, thresholds: %s, log-likelihood: %.4f\n",
                      i, k, paste(thresholds, collapse = ", "), logliks_k[i]))
        }
      }
      best_idx <- which.max(logliks_k)
      best_thresholds <- comb[best_idx, ]
      logliks[k + 1] <- logliks_k[best_idx]
      thresholds_list[[k + 1]] <- best_thresholds
      all_combinations_list[[k + 1]] <- comb
      fits[[k]] <- fit
    }

    # Compute LRT p-values
    p_values <- numeric(n_thresholds)
    for (k in 1:n_thresholds) {
      lr_stat <- 2 * (logliks[k + 1] - logliks[k])
      df_diff <- if (type == "hinge") 1 else 2
      p_values[k] <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
    }

    # Fit final model
    k <- n_thresholds
    best_thresholds <- thresholds_list[[k + 1]]
    if (type == "hinge") {
      X_threshold <- sapply(best_thresholds, function(e) pmax(chngpt_var - e, 0))
    } else {
      X_step <- sapply(best_thresholds, function(e) as.numeric(chngpt_var > e))
      X_seg <- sapply(best_thresholds, function(e) pmax(chngpt_var - e, 0))
      X_threshold <- cbind(X_step, X_seg)
    }
    X <- cbind(Z, X_threshold)
    best_fit <- glm(y ~ 0 + X, family = switch(family, "gaussian" = gaussian(),
                                               "binomial" = binomial()),
                    weights = prec.weights)

    # Cluster-robust SEs
    if (!is.null(cluster_var)) {
      cluster <- data[[cluster_var]]
      vcov_cluster <- sandwich::vcovCL(best_fit, cluster = cluster)
      se_cluster <- sqrt(diag(vcov_cluster))
    } else {
      vcov_cluster <- NULL
      se_cluster <- NULL
    }

    # Results
    res <- list(
      thresholds = best_thresholds,
      fit = best_fit,
      loglik = logliks[n_thresholds + 1],
      all_logliks = logliks,
      all_combinations = all_combinations_list,
      p_values = p_values,
      vcov_cluster = vcov_cluster,
      se_cluster = se_cluster,
      data.name = DNAME
    )
    class(res) <- "multi.segsteg.test"
  }

  .Random.seed <- save.seed
  return(res)
}

#' @title Generate Alternative Model Formula for Change Point Analysis
#'
#' @description
#' Constructs the alternative model formula for change point analysis based on the specified type.
#' The formula is used in conjunction with the null model to test for change points.
#'
#' @param type Character string specifying the type of change point model. Supported types are:
#'   \itemize{
#'     \item \code{"step"}: Adds an indicator term \code{chngpt.var} to the model.
#'     \item \code{"hinge"}: Adds a hinge term \code{chngpt.var} to the model.
#'     \item \code{"segmented"}: Adds a segmented term \code{chngpt.var} to the model.
#'     \item \code{"stegmented"}: Adds both step and segmented terms \code{chngpt.var.step} and \code{chngpt.var.seg}.
#'   }
#' @param var Character string specifying the name of the change point variable.
#' @param modified.by Optional character string specifying a variable for interaction effects.
#'   If provided, the formula includes an interaction between \code{chngpt.var} and \code{modified.by}.
#'   Note: Interaction is not implemented for \code{type = "stegmented"}.
#'
#' @return
#' A formula object representing the alternative model for the specified change point type.
#'
#' @details
#' This function generates the formula for the alternative hypothesis in change point testing.
#' It dynamically constructs the formula based on the type of change point model and whether
#' an interaction term is specified. The resulting formula extends the null model by adding
#' change point terms, making it suitable for use in functions like \code{segsteg.test}.
#'
#' @examples
#' # Formula for a hinge model without interaction
#' get.f.alt(type = "hinge", var = "x")
#' # Returns: ~ . + chngpt.var
#'
#' # Formula for a step model with interaction
#' get.f.alt(type = "step", var = "x", modified.by = "z")
#' # Returns: ~ . + chngpt.var * z
#'
#' @export
get.f.alt <- function(type, var, modified.by = NULL) {
  if (is.null(modified.by)) {
    switch(type,
           step = as.formula("~ . + chngpt.var"),
           hinge = as.formula("~ . + chngpt.var"),
           segmented = as.formula("~ . + chngpt.var"),
           stegmented = as.formula("~ . + chngpt.var.step + chngpt.var.seg"))
  } else {
    if (type == "stegmented") {
      stop("Interaction with stegmented not implemented")
    } else {
      as.formula(paste("~ . + chngpt.var *", modified.by))
    }
  }
}

#' @title Create Change Point Variables in Data Frame
#'
#' @description
#' Modifies the input data frame by adding columns for change point variables based on the specified type.
#' These columns are used in the alternative model for change point analysis.
#'
#' @param chngpt_var Numeric vector representing the change point variable.
#' @param chngpt Numeric value specifying the threshold (change point) value.
#' @param type Character string specifying the type of change point model. Supported types are:
#'   \itemize{
#'     \item \code{"step"}: Creates a step indicator variable \code{chngpt.var}.
#'     \item \code{"hinge"}: Creates a hinge variable \code{chngpt.var}.
#'     \item \code{"segmented"}: Creates a segmented variable \code{chngpt.var}.
#'     \item \code{"stegmented"}: Creates both step and segmented variables \code{chngpt.var.step} and \code{chngpt.var.seg}.
#'   }
#' @param data Data frame to which the change point variables will be added.
#'
#' @return
#' The input data frame with additional columns for the change point variables based on the specified type.
#'
#' @details
#' This function computes the change point terms for a given threshold and adds them as new columns to the data frame.
#' The specific columns added depend on the type:
#' - For \code{"step"}, it adds \code{chngpt.var} as an indicator (1 if \code{chngpt_var > chngpt}, 0 otherwise).
#' - For \code{"hinge"} or \code{"segmented"}, it adds \code{chngpt.var} as the hinge term \code{pmax(chngpt_var - chngpt, 0)}.
#' - For \code{"stegmented"}, it adds two columns: \code{chngpt.var.step} (step indicator) and \code{chngpt.var.seg} (hinge term).
#' These variables align with the alternative model formula generated by \code{get.f.alt}.
#'
#' @examples
#' # Create a step change point variable
#' data <- data.frame(x = 1:10)
#' make.chngpt.var(chngpt_var = data$x, chngpt = 5, type = "step", data = data)
#' # Adds chngpt.var: c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#'
#' # Create hinge and step variables for stegmented model
#' make.chngpt.var(chngpt_var = data$x, chngpt = 5, type = "stegmented", data = data)
#' # Adds chngpt.var.step and chngpt.var.seg
#'
#' @export
make.chngpt.var <- function(chngpt_var, chngpt, type, data) {
  if (type == "step") {
    data$chngpt.var <- as.numeric(chngpt_var > chngpt)
  } else if (type %in% c("hinge", "segmented")) {
    data$chngpt.var <- pmax(chngpt_var - chngpt, 0)
  } else if (type == "stegmented") {
    data$chngpt.var.step <- as.numeric(chngpt_var > chngpt)
    data$chngpt.var.seg <- pmax(chngpt_var - chngpt, 0)
  }
  data
}
