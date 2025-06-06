#' Wrapper for inverse cumulative hazard
#'
#' A wrapper around the Rcpp function `inverseScHaz`, used to find the inverse
#' by numeric methods in case of no simple analytical solution.
#'
#' @param p The random variable (typically `-log(U)`).
#' @param t Current time.
#' @param lower Lower bound for root finding.
#' @param upper Upper bound for root finding.
#' @param eta Numeric vector of shape parameters.
#' @param nu Numeric vector of scale parameters.
#' @param phi Numeric vector of multiplicative effects.
#' @param at_risk Numeric vector indicating at-risk indicators for each event type.
#' @param tol Numeric tolerance for root-finding. Default is 1e-9.
#' @param max_iter Maximum iterations. Default is 100.
#'
#' @return A numeric scalar, the root `u`.
#' @export
#'
#' @examples
#' eta <- c(0.1, 0.1)
#' nu <- c(1.1, 1.1)
#' phi <- c(1, 1)
#' at_risk <- c(1, 1)
#' inverseScHaz(0.5, t = 0, lower = 1e-15, upper = 200, eta, nu, phi, at_risk)
inverseScHaz <- function(p, t, lower, upper, eta, nu, phi, at_risk,
                              tol = 1e-9, max_iter = 100) {
  inverseScHazCpp(
    p = p,
    t = t,
    lower = lower,
    upper = upper,
    eta = eta,
    nu = nu,
    phi = phi,
    at_risk = at_risk,
    tol = tol,
    max_iter = max_iter
  )
}
