#' Trapezoidal Integration
#'
#' @param x Numerical vector of x values
#' @param y Numerical vector of y values
#'
#' @return Numeric area under the curve
#' @export
trapz <- function(x, y) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (length(x) < 2) return(0)
  idx <- 2:length(x)
  sum((x[idx] - x[idx - 1]) * (y[idx] + y[idx - 1])) / 2
}
