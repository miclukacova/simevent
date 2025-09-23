#' Wrapper for inverse cumulative hazard
#'
#' A wrapper around the Rcpp function `inverseScHazCppTV`, used to find the inverse
#' of the summed cumulative hazard.
#'
#' @param p The random variable (typically `-log(U)`).
#' @param t The time of the previous event
#' @param lower Lower bound for root finding.
#' @param upper Upper bound for root finding.
#' @param t_prime The time where the time varying effects change
#' @param eta Numeric vector of shape parameters.
#' @param nu Numeric vector of scale parameters.
#' @param phi Numeric vector of multiplicative effect bedfore time t_prime
#' @param phi_prime Numeric vector of multiplicative effects after time t_prime
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
#' phi_prime <- c(2, 2)
#' inverseScHazTV(p = 0.5, t= 1, t_prime = 2, eta = eta, nu = nu,
#'                        phi = phi, phi_prime = phi_prime, at_risk = at_risk)
#'
inverseScHazTV <- function(p, t, lower = 1e-15, upper = 200, t_prime, eta, nu,
                                   phi, phi_prime, at_risk, tol = 1e-9, max_iter = 100) {
  inverseScHazTVCpp(p = p,
                    t = t,
                    lower = lower,
                    upper = upper,
                    t_prime = t_prime,
                    eta = eta,
                    nu = nu,
                    phi = phi,
                    phi_prime = phi_prime,
                    at_risk = at_risk,
                    tol = tol,
                    max_iter = max_iter)
}
