#' Generate Response Curve Data from Model
#'
#' @param model A fitted model object (e.g., from glmnet_mx)
#' @param variable Character, name of the variable to generate the curve for
#' @param range Numeric vector of length 2, the range of the variable
#' @param n Numeric, number of points to generate
#' @export
get_response_curve <- function(model, variable, range = NULL, n = 100) {
  if (is.null(range)) {
    range <- c(model$varmin[variable], model$varmax[variable])
  }

  seq_x <- seq(range[1], range[2], length.out = n)

  # Create data frame with other variables at their sample means
  test_data <- as.data.frame(lapply(model$samplemeans, function(x) rep(x, n)))
  test_data[[variable]] <- seq_x

  pred <- predict.glmnet_mx(model, test_data, type = "exponent")
  data.frame(x = seq_x, y = pred / max(pred)) # Normalized to 1
}
