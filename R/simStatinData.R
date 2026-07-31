#' @title Simulate Data in a Statin Setting
#'
#' @param N Numeric scalar. Number of individuals to simulate.
#' @param eta Numeric vector of length equal to number of processes. Shape parameters for Weibull intensities with parameterization
#' \eqn{\eta \nu t^{\nu - 1}}. Defaults to \code{rep(0.1, 8)}.
#' @param nu Numeric vector of length equal to number of processes. Scale parameters for the Weibull hazards. Defaults to \code{rep(1.1, 8)}.
#' @param followup Numeric scalar. Maximum follow-up (censoring) time. Defaults to \code{Inf}.
#' @param lower Numeric scalar. Lower bound for root-finding (inverse cumulative hazard) (default \code{1e-15}).
#' @param upper Numeric scalar. Upper bound for root-finding (default 200).
#' @param beta Numeric matrix. Of dimension p times 6. Regression coefficients matrix where columns correspond to event types (N0, ..., N5) and rows correspond to covariates (L0, A0, L1, L2, ...) followed by event counts (N0, ..., N5). Default is a zero matrix.
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0.Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#' @param gen_L0 Function. Function to generate the baseline covariate L0. Takes N as input. Default is a N(0,1) random variable.
#' @param add_cov Named list of functions. Functions generating additional baseline covariates. Each function takes integer N and returns a numeric vector of length N. Default is NULL.
#' @param at_risk Function. Function determining if an individual is at risk for each event type, given their current event counts. Takes a numeric vector of event counts and returns a binary vector. Default returns 1 for all events.
#' @param ... Additional arguments passed to \code{simEventData}
#'
#' @return A data frame containing the simulated data with columns:
#'  \item{ID}{Individual identifier}
#'  \item{Time}{Time of the event}
#'  \item{Delta}{Event type (0,...,5)}
#'  \item{L0}{Baseline covariate}
#'  \item{A0}{Baseline covariate}
#'  \item{L1,...Lp}{Additional baseline covariates}
#'
#' @examples
#' simDisease(10)
#'
#' @export
simStatinData <- function(N,
                          eta = NULL,
                          nu = NULL,
                          beta = NULL,
                          followup = 5,
                          lower = 10^(-15),
                          upper = 200,
                          gen_A0 = function(N, L0) pmin(stats::rexp(N, 0.3) + 70, 100),
                          gen_L0 = function(N) stats::rbinom(N, 1, 0.4),
                          add_cov = NULL,
                          at_risk = NULL,
                          ...){
  # Default values
  if(!is.null(beta)){
    n_proc <- ncol(beta)
  } else if(!is.null(eta)){
    n_proc <- length(eta)
  } else if(!is.null(nu)){
    n_proc <- length(nu)
  } else {
    n_proc <- 12
    beta <- matrix(0, nrow = n_proc + 2 + length(add_cov), ncol = n_proc)
  }

  if(is.null(eta)) eta <- rep(0.1, n_proc)
  if(is.null(nu)) nu <- rep(1.1, n_proc)

  if(is.null(at_risk)) {
    at_risk <- function(events) {
      return(rep(1, n_proc))
    }
  }

  data <- simEventData(N,
                       beta = beta,
                       eta = eta,
                       nu = nu,
                       at_risk = at_risk,
                       max_cens = followup,
                       lower = lower,
                       upper = upper,
                       term_deltas = c(0,1,2),
                       gen_A0 = gen_A0,
                       gen_L0 = gen_L0,
                       add_cov = add_cov,
                       ...)

  return(data)
}
