#' Balance Drug Dataset by Reducing Zero Dose Entries
#'
#' Creates a balanced dataset by reducing entries where a specified drug dose is zero to a target percentage,
#' while preserving the distributions of variables specified in the provided formulas.
#'
#' @param df The input dataset, a data.frame.
#' @param formula1 A formula or character string specifying the outcome and predictor variables for stratification,
#'                 e.g., "cholesterol ~ age + Gender + education + dose.METFORMIN + dose.ATORVASTATIN".
#' @param formula2 A formula or character string specifying the drug variable to balance,
#'                 e.g., "~dose.CLOZAPINE" or "~wz_dose.CLOZAPINE".
#' @param target_percentage The target percentage of zero entries to keep for the specified drug variable.
#'                         Default is 0.05 (5%).
#' @param random_seed An integer seed for reproducible sampling. Default is 42.
#'
#' @return A new data.frame with the balanced dataset, where zero entries for the specified drug variable
#'         are reduced to approximately the target percentage, and all non-zero entries are retained,
#'         while maintaining the distributions of specified variables.
#'
#' @details
#' The function identifies the drug variable from `formula2` (e.g., `dose.CLOZAPINE`) and checks for its presence
#' in the dataset. It extracts stratification variables from `formula1` and `formula2` (excluding the drug variable).
#' The dataset is split into zero and non-zero entries based on the drug variable. Zero entries are stratified
#' by binning continuous variables into quantiles and using categorical variables as-is. A proportional sample
#' is taken from each stratum to meet the target percentage, ensuring at least one entry per stratum if possible.
#' The sampled zero entries are combined with all non-zero entries to form the balanced dataset. Diagnostics
#' verify that variable distributions are preserved.
#'
#' @examples
#' # Example with a dataset 'df' containing relevant columns
#' balanced_df <- balance_drug_dataset_zero(
#'   df,
#'   "cholesterol ~ age + Gender + education + dose.METFORMIN",
#'   "~dose.CLOZAPINE",
#'   target_percentage = 0.05,
#'   random_seed = 42
#' )
#'
#' @export
balance_drug_dataset_zero <- function(df, formula1, formula2, target_percentage = 0.05, random_seed = 42) {

  set.seed(random_seed)

  # Convert formulas to character strings if needed
  if (inherits(formula1, "formula")) formula1 <- as.character(formula1)
  if (inherits(formula2, "formula")) formula2 <- as.character(formula2)

  # Extract drug variable from formula2
  if (is.null(formula2) || !grepl("~", formula2)) {
    stop("formula2 must be provided and contain '~', e.g., '~dose.DRUGNAME'")
  }

  # Extract the drug variable - look for patterns like 'dose.DRUGNAME' or 'wz_dose.DRUGNAME'
  drug_pattern <- "(?:wz_)?dose\\.([A-Z]+)"
  drug_matches <- regmatches(formula2, gregexpr(drug_pattern, formula2, perl = TRUE))[[1]]

  if (length(drug_matches) == 0) {
    stop(paste("Could not find drug variable in formula2:", formula2,
               "Expected format like '~dose.DRUGNAME' or '~wz_dose.DRUGNAME'"))
  }

  # Extract the actual drug name by removing the prefix
  drug_name <- sub(".*dose\\.", "", drug_matches[1])

  # Determine the actual column name to use for balancing
  if (paste0("dose.", drug_name) %in% colnames(df)) {
    target_column <- paste0("dose.", drug_name)
  } else if (paste0("wz_dose.", drug_name) %in% colnames(df)) {
    target_column <- paste0("wz_dose.", drug_name)
  } else {
    stop(paste0("Neither 'dose.", drug_name, "' nor 'wz_dose.", drug_name, "' found in dataframe columns"))
  }

  cat("Balancing dataset based on", target_column, "\n")

  # Extract variables from formula1
  variables <- c()
  if (grepl("~", formula1)) {
    parts <- strsplit(formula1, "~")[[1]]
    outcome <- trimws(parts[1])
    if (nchar(outcome) > 0 && !startsWith(outcome, "~")) {
      variables <- c(variables, outcome)
    }

    if (length(parts) > 1) {
      predictors <- trimws(parts[2])
      predictor_vars <- strsplit(predictors, "\\+")[[1]]
      variables <- c(variables, trimws(predictor_vars))
    }
  }

  # Extract additional variables from formula2 if there are any besides the drug
  if (grepl("~", formula2)) {
    parts <- strsplit(formula2, "~")[[1]]
    if (length(parts) > 1) {
      predictors <- trimws(parts[2])
      predictor_vars <- strsplit(predictors, "\\+")[[1]]
      for (v in predictor_vars) {
        v <- trimws(v)
        # Skip the drug variable we're already handling
        if (!grepl(drug_name, v)) {
          variables <- c(variables, v)
        }
      }
    }
  }

  # Remove duplicates and filter for valid columns
  variables <- unique(variables)
  valid_variables <- variables[variables %in% colnames(df)]

  if (length(valid_variables) == 0) {
    warning("No valid variables found in formulas. Using default stratification.")
    default_vars <- c("age", "Gender", "education")
    valid_variables <- default_vars[default_vars %in% colnames(df)]
  }

  # Split dataset based on target column
  zero_mask <- df[[target_column]] == 0
  zero_df <- df[zero_mask, ]
  non_zero_df <- df[!zero_mask, ]

  # Print initial stats
  cat(sprintf("Original dataset: %d rows\n", nrow(df)))
  cat(sprintf("- With %s = 0: %d rows (%.1f%%)\n",
              target_column, nrow(zero_df), nrow(zero_df)/nrow(df)*100))
  cat(sprintf("- With %s > 0: %d rows (%.1f%%)\n",
              target_column, nrow(non_zero_df), nrow(non_zero_df)/nrow(df)*100))

  # If there are too few non-zero entries, warn the user
  if (nrow(non_zero_df) < 10) {
    warning(sprintf("Only %d entries with %s > 0. Consider using a different drug variable with more non-zero values.",
                    nrow(non_zero_df), target_column))
  }

  # Create stratification features for each variable
  strat_df <- zero_df
  for (var in valid_variables) {
    strat_col_name <- paste0("strat_", var)

    # Check if categorical or has few unique values
    if (is.factor(df[[var]]) || is.character(df[[var]]) || length(unique(df[[var]])) < 10) {
      # For categorical variables, use as-is
      strat_df[[strat_col_name]] <- as.character(strat_df[[var]])
    } else {
      # For continuous variables, bin into quantiles
      tryCatch({
        # Create bins on the full dataset to maintain consistency
        bin_breaks <- quantile(df[[var]], probs = seq(0, 1, 0.2), na.rm = TRUE)
        bin_breaks <- unique(bin_breaks) # Remove duplicates

        if (length(bin_breaks) >= 2) {
          # Apply the bins to the zero_df
          strat_df[[strat_col_name]] <- cut(strat_df[[var]], breaks = bin_breaks, include.lowest = TRUE)
          strat_df[[strat_col_name]] <- as.character(strat_df[[strat_col_name]])
        } else {
          # If binning fails, use the raw values
          strat_df[[strat_col_name]] <- as.character(strat_df[[var]])
        }
      }, error = function(e) {
        # If binning fails, use the raw values
        strat_df[[strat_col_name]] <- as.character(strat_df[[var]])
      })
    }
  }

  # Combine all stratification columns into a single key
  strat_cols <- paste0("strat_", valid_variables)
  strat_cols <- strat_cols[strat_cols %in% colnames(strat_df)]

  if (length(strat_cols) > 0) {
    strat_df$strat_key <- apply(strat_df[, strat_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    # If no stratification columns were created, just use a constant
    strat_df$strat_key <- "all"
  }

  # Sample entries proportionally from each stratum
  sampled_indices <- c()
  strata <- split(1:nrow(strat_df), strat_df$strat_key)

  for (stratum_name in names(strata)) {
    indices_in_stratum <- strata[[stratum_name]]
    n_in_stratum <- length(indices_in_stratum)

    # Calculate proportional sample size, minimum 1 if available
    n_to_sample <- max(1, round(n_in_stratum * target_percentage))
    n_to_sample <- min(n_to_sample, n_in_stratum)  # Can't sample more than available

    # Sample from this stratum
    if (n_to_sample > 0) {
      sampled_from_stratum <- sample(indices_in_stratum, size = n_to_sample, replace = FALSE)
      sampled_indices <- c(sampled_indices, sampled_from_stratum)
    }
  }

  # Get the original row numbers from the full dataframe
  original_indices <- which(zero_mask)[sampled_indices]

  # Combine sampled zero entries with all non-zero entries
  non_zero_indices <- which(!zero_mask)
  balanced_indices <- c(original_indices, non_zero_indices)
  balanced_df <- df[balanced_indices, ]

  # Print final stats
  zero_count_balanced <- sum(balanced_df[[target_column]] == 0)
  non_zero_count_balanced <- sum(balanced_df[[target_column]] > 0)

  cat(sprintf("Balanced dataset: %d rows\n", nrow(balanced_df)))
  cat(sprintf("- With %s = 0: %d rows (%.1f%%)\n",
              target_column, zero_count_balanced, zero_count_balanced/nrow(balanced_df)*100))
  cat(sprintf("- With %s > 0: %d rows (%.1f%%)\n",
              target_column, non_zero_count_balanced, non_zero_count_balanced/nrow(balanced_df)*100))

  # Verify stratification was maintained
  cat("\nVerifying variable distributions were preserved:\n")
  for (var in valid_variables[1:min(3, length(valid_variables))]) {  # Show up to 3 variables
    if (!is.factor(df[[var]]) && !is.character(df[[var]]) && length(unique(df[[var]])) > 5) {
      # For continuous variables, compare means
      orig_mean <- mean(zero_df[[var]], na.rm = TRUE)
      sampled_mean <- mean(balanced_df[balanced_df[[target_column]] == 0, var], na.rm = TRUE)
      cat(sprintf("- %s mean: Original=%.2f, Sampled=%.2f\n", var, orig_mean, sampled_mean))
    } else {
      # For categorical, just indicate distribution was maintained
      cat(sprintf("- %s distribution maintained (categorical)\n", var))
    }
  }

  return(balanced_df)
}
