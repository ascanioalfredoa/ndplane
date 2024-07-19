#' Stores Mathematical Expression For 4-Parameter Beta Function
#'
#' @return An R expression
#' @export
#'
#' @examples
#' betafunction_max()
betafunction_max <- function(){
    expression(k*((x-a)^alpha)*((b-x)^gamma))
}

#' Creates X and Y values for a Beta Function Curve
#'
#' @param a Lower bound of beta curve
#' @param b Upper bound of beta curve
#' @param alpha \eqn{\alpha} shape parameter (right skew)
#' @param gamma \eqn{\gamma} shape parameter (left skew)
#' @param k Scale parameter, usually left at 1 by default
#' @param interval Interval length for steps between a and b. By default 0.001
#'
#' @return Data frame containing X and Y values for the desired intervals between a and b
#' @export
#'
#' @examples
#' a = 0; b = 1
#' alpha = gamma <- 2
#' betaPDF(a, b, alpha, gamma)

betaPDF <- function(a = 0, b = 1, alpha = 1, gamma = 1, k = 1, interval = 0.001) {
    x = round(seq(a, b, interval), 3)
    y = eval(betafunction_max())
    x = x[y >= 0]
    y = y[y >= 0]
    y = y/max(y)
    data.frame(x, y)
}

#' Performs betaPDF calculation on Parameter Space row
#'
#' @param psv (row) vector, or single row data frame from parameter space
#' @param k Scale parameter, usually left at 1 by default
#' @param interval Interval length for steps between a and b. By default 0.001
#'
#' @return Data frame containing X and Y values for the desired intervals between a and b
#' @export
#'
#' @examples
#' a <- 1; b <- 2
#' alpha = gamma <- 1:2
#' par_space <- set_parameter_space(a, b, alpha, gamma)
#' betaPDF_ps(psv = par_space[1, ])
betaPDF_ps <- function(psv, k = 1, interval = 0.001) {
    a <- psv[[1]]
    b <- psv[[2]]
    alpha <- psv[[3]]
    gamma <- psv[[4]]
    betaPDF(a, b, alpha, gamma, k = k, interval = interval)
}
