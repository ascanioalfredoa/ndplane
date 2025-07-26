#' Run Niche Divergence Analysis
#'
#' Master function to run niche divergence analysis using specified method.
#'
#' @param method Character. Either "maxent" or "kde" to select analysis approach.
#' #@param params List of named parameters required by the selected method functions.
#' @return List containing results specific to the analysis method.
#' @export
run_niche_divergence_analysis <- function(method = c("maxent", "kde"), params = list()) {
  method <- match.arg(method)
  if (method == "maxent") {
    # Expecting params to contain all args needed by run_maxent_ndp_pipeline
    if (!requireNamespace("maxnet", quietly = TRUE)) {
      stop("maxnet package required for maxent analysis. Please install it.")
    }
    result <- do.call(run_maxent_ndp_pipeline, params)
  } else if (method == "kde") {
    # Expecting params to contain all args needed by run_sal_ndp_analysis
    if (!requireNamespace("pracma", quietly = TRUE)) {
      stop("pracma package required for kde analysis. Please install it.")
    }
    result <- do.call(run_sal_ndp_analysis, params)
  } else {
    stop("Unsupported method")
  }
  result
}

#' Plot Results
#'
#' Generic plot function dispatching based on method.
#' @param results List returned by run_niche_divergence_analysis
#' @param method Character. "maxent" or "kde" matching the analysis method.
#' @return Returns plot object if available.
#' @export
plot_niche_divergence_results <- function(results, method) {
  if (method == "kde") {
    results$plot
  } else if (method == "maxent") {
    # Plotting could be implemented to visualize maxent niche indices or response curves
    # Returning NULL for now
    message("Plotting for maxent not implemented yet.")
    return(NULL)
  } else {
    stop("Unsupported method")
  }
}


# Note: Remember to import and make visible run_maxent_ndp_pipeline and run_sal_ndp_analysis functions
# in your package namespace for this wrapper to work.
