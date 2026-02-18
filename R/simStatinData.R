#' @title Simulate Data in a Statin Setting
#'
#' @param N Numeric scalar. Number of individuals to simulate.
#' @param eta Numeric vector of length 6. Shape parameters for Weibull intensities with parameterization
#' \eqn{\eta \nu t^{\nu - 1}}. Defaults to \code{rep(0.1, 6)}.
#' @param nu Numeric vector of length 6. Scale parameters for the Weibull hazards. Defaults to \code{rep(1.1, 6)}.
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{1}).
#' @param followup Numeric scalar. Maximum follow-up (censoring) time. Defaults to \code{Inf}.
#' @param lower Numeric scalar. Lower bound for root-finding (inverse cumulative hazard) (default \code{1e-15}).
#' @param upper Numeric scalar. Upper bound for root-finding (default 200).
#' @param beta Numeric matrix. Of dimension p times 6. Regression coefficients matrix where columns correspond to event types (N0, ..., N5) and rows correspond to covariates (L0, A0, L1, L2, ...) and event counts (N0, ..., N5). Default is a zero matrix.
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0.
#' Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
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
                          eta = rep(0.1,6),
                          nu = rep(1.1,6),
                          beta = beta,
                          cens = 1,
                          followup = 5,
                          lower = 10^(-15),
                          upper = 200,
                          gen_A0 = NULL,
                          ...){

  if(is.null(beta)) beta <- matrix(0, nrow = 30, ncol = 6)

  at_risk <- function(events) {
    return(c(cens,                          # If you have not yet been censored you are at risk (if there is a censoring process)
      1,                                    # If you have not you died yet are at risk
      1,                                    # If you have not had CVD you are at risk
      as.numeric(events[4] == 0),           # If you have not yet stopped statin use, you are at risk
      as.numeric(events[5] <= 3),           # Number of diseases
      as.numeric(events[6] <= 3)))          # Number of medicines
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
                       cens = cens,
                       ...)

  #colnames(data) <- c("Censoring", "Death", "CVD", "Off Stat", "Medicines", "Diseases")

  return(data)
}
