#' Simulate Survival Data with Censoring and Event Times
#'
#' Simulates survival data for \eqn{N} individuals who are at risk for censoring (0) and an event (1).
#' The hazard functions for censoring and event times follow Weibull distributions parameterized by
#' shape parameters \eqn{\eta} and scale parameters \eqn{\nu}. Covariate effects on censoring and event
#' hazards are specified via a matrix \code{beta}.
#'
#' @title Simulate Survival Data
#'
#'@param N Numeric scalar. Number of individuals to simulate.
#' @param beta Numeric 2x2 matrix specifying effects of baseline covariates \code{L0} and \code{A0} on censoring and event hazards.
#'   - Rows correspond to covariates \code{L0} and \code{A0}.
#'   - Columns correspond to censoring (1st column) and event (2nd column).
#'   Defaults to zero matrix if \code{NULL}.
#' @param eta Numeric vector of length 2. Shape parameters for Weibull hazard with parameterization \eqn{\eta \nu t^{\nu - 1}}.
#'   Defaults to \code{rep(0.1, 2)}.
#' @param nu Numeric vector of length 2. Scale parameters for the Weibull hazard. Defaults to \code{rep(1.1, 2)}.
#' @param cens Numeric binary indicator (0 or 1) specifying if censoring is included (default 1).
#' @param ... Additional arguments passed to \code{simEventData}, including the argument \code{add_cov} to specify extra covariates.
#'
#' @return  Data frame containing the simulated survival data
#'
#' @examples
#' simSurvData(10)
#'
#' @export
#'
simSurvData <- function(N,
                        beta = NULL,
                        eta = rep(0.1,2),
                        nu = rep(1.1,2),
                        cens = 1,
                        ...
                        ){
    at_risk <- function(events) c(cens,1)

    # Set default beta matrix if NULL
    if(is.null(beta)){
        beta <- matrix(0, ncol = 2, nrow = 2)
    }

    # Pad beta matrix for simEventData function
    beta <- rbind(beta, matrix(0, nrow = 2, ncol = 2))

    # If additional covariates specified, add rows to beta accordingly
    dots <- list(...)
    has_add_cov <- "add_cov" %in% names(dots)
    if (has_add_cov) beta <- rbind(beta, matrix(c(rep(0, length(dots$add_cov)), rep(0.1, length(dots$add_cov))), nrow = length(dots$add_cov), ncol = ncol(beta)))

    # Simulate data using underlying simEventData function
    results <- simEventData(N, beta, eta = eta, nu = nu, at_risk = at_risk, ...)

    # Remove terminal event indicator columns
    results <- results[, !c("N0", "N1")]

    return(results)
}

