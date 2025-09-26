#' Simulate Competing Risks Data
#'
#' Simulates competing risks data for \eqn{N} individuals who are at risk of mutually exclusive event types.
#' Three event types are simulated, where one can be interpreted as censoring.
#'
#' The event intensities follow Weibull hazard models parameterized by shape and scale parameters \eqn{\eta} and \eqn{\nu}.
#' Covariate effects on the hazard are specified by the \code{beta} matrix, which models the effects of baseline covariates \code{L0} and \code{A0} on each event type.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param beta Numeric matrix of dimension 2x3. Covariate effects of \code{L0} and \code{A0} on the three competing processes (columns correspond to processes). Defaults to zero matrix if \code{NULL}.
#' @param eta Numeric vector of length 3. Shape parameters for Weibull hazards, parameterized as \eqn{\eta \nu t^{\nu - 1}}. Defaults to \code{rep(0.1, 3)}.
#' @param nu Numeric vector of length 3. Scale parameters for Weibull hazards. Defaults to \code{rep(1.1, 3)}.
#' @param cens Binary (0 or 1). Indicates if a censoring process is included. Default is 1.
#' @param ... Additional arguments passed to \code{\link{simEventData}}, including \code{add_cov} for extra covariates.
#'
#' @return A \code{data.frame} with simulated competing risk data including:
#' \itemize{
#'   \item \code{ID} - Individual identifier.
#'   \item \code{Time} - Event time.
#'   \item \code{Delta} - Event type (0, 1, or 2).
#'   \item \code{L0} - Baseline covariate.
#'   \item \code{A0} - Baseline treatment indicator.
#' }
#'
#' @export
#'
#' @examples
#' simCRdata(10)
simCRdata <- function(N,
                      beta = NULL,
                      eta = rep(0.1,3),
                      nu = rep(1.1,3),
                      cens = 1,
                      ...
){

    at_risk <- function(events) c(cens,1,1)

    if(is.null(beta)){
        beta <- matrix(0, ncol = 3, nrow = 2)
    }

    beta <- rbind(beta, matrix(0, ncol = 3, nrow = 3))

    dots <- list(...)
    has_add_cov <- "add_cov" %in% names(dots)
    if (has_add_cov) beta <- rbind(beta, matrix(c(rep(0, length(dots$add_cov)),
                                                  rep(0.1, length(dots$add_cov))),
                                                nrow = length(dots$add_cov),
                                                ncol = ncol(beta)))

    results <- simEventData(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk,
                            term_deltas = c(0,1,2),
                            ...)

    results <- results[, !c("N0", "N1", "N2")]

    return(results)
}
