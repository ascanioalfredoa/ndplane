#' Create Parameter Space for Niche Divergence Plane Simulations
#'
#' @param a Numerical vector containing set of a values for beta function
#' @param b Numerical vector containing set of b values for beta function
#' @param alpha Numerical vector containing set of alpha values for beta function
#' @param gamma Numerical vector containing set of gamma values for beta function
#'
#' @return Data frame containing all combinations of parameters a, b, alpha, and gamma, where a < b
#' @export
#'
#' @examples
#' a = b <- round(seq(0, 1, 0.1), 3)
#' alpha = gamma <- round(c(0.001, 0.01, seq(0.1, 1, 0.2), seq(2, 10, 2), 100, 1000), 3)
#' set_parameter_space(a, b, alpha, gamma)
set_parameter_space <- function(a, b, alpha, gamma){
    par_space <- expand.grid(a = a, b = b, alpha = alpha, gamma = gamma)
    par_space <- par_space[par_space$a < par_space$b, ]
    par_space
}
