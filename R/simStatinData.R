#' @title Simulate Data in a Statin Setting
#'
#' @param N Numeric scalar. Number of individuals to simulate.
#' @param n_proc Number of processes.
#' @param eta Numeric vector of length equal to number of processes. Shape parameters for Weibull intensities with parameterization
#' \eqn{\eta \nu t^{\nu - 1}}. Defaults to \code{rep(0.1, 8)}.
#' @param nu Numeric vector of length equal to number of processes. Scale parameters for the Weibull hazards. Defaults to \code{rep(1.1, 8)}.
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{1}).
#' @param followup Numeric scalar. Maximum follow-up (censoring) time. Defaults to \code{Inf}.
#' @param lower Numeric scalar. Lower bound for root-finding (inverse cumulative hazard) (default \code{1e-15}).
#' @param upper Numeric scalar. Upper bound for root-finding (default 200).
#' @param beta Numeric matrix. Of dimension p times 6. Regression coefficients matrix where columns correspond to event types (N0, ..., N5) and rows correspond to covariates (L0, A0, L1, L2, ...) followed by event counts (N0, ..., N5). Default is a zero matrix.
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0.Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#' @param gen_L0 Function. Function to generate the baseline covariate L0. Takes N as input. Default is a N(0,1) random variable.
#' @param add_cov Named list of functions. Functions generating additional baseline covariates. Each function takes integer N and returns a numeric vector of length N. Default is NULL.
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
                          n_proc = 8,
                          eta = rep(0.1,8),
                          nu = rep(1.1,8),
                          beta = NULL,
                          cens = 1,
                          followup = 5,
                          lower = 10^(-15),
                          upper = 200,
                          gen_A0 = NULL,
                          gen_L0 = NULL,
                          add_cov = NULL,
                          ...){


  if(is.null(beta)) beta <- matrix(0, nrow = n_proc + 2 + length(add_cov), ncol = n_proc)

  at_risk <- function(events) {
    return(c(cens,                                 # If you have not yet been censored you are at risk (if there is a censoring process)
             1,                                    # If you have not you died yet are at risk
             1,                                    # If you have not had CVD you are at risk
             as.numeric(events[4] <= 10),          # Statin Stop
             as.numeric(events[5] <= 10),          # Increase in number of diseases
             as.numeric(events[6] <= 10),          # Increase in number of medicines
             as.numeric(events[7] <= 10),          # LDL increase
             as.numeric(events[8] <= 10)))         # LDL decrease

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

  colnames(data)[(6 + length(add_cov)):ncol(data)] <- c("C", "D", "CVD", "OS", "L", "A", "LDL1", "LDL2")

  return(data)
}
