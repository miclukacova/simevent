#' Simulate Data from a Disease Setting
#'
#' This function simulates event data representing three event types:
#' Censoring (0), Death (1), and Change in Covariate Process (2). Death and Censoring are terminal events,
#' while Change in Covariate Process can occur only once.
#'
#' Event intensities depend on previous events and predefined parameters \eqn{\nu} and \eqn{\eta}.
#'
#' The arguments \code{beta_X_Y} control how the X affects Y. A positive value means that a higher value of X
#' increases the intensity of Y, while a negative value decreases the intensity.
#'
#' The simulation uses an event history framework with terminal events (death, censoring) and a single recurrent covariate change.
#' The event intensities depend on covariates and previous events according to user-specified parameters.
#' Time-varying effects can be included via \code{beta_L_D_t_prime} and \code{t_prime}.
#'
#' @title Simulate Data in a Disease Setting
#'
#' @param N Numeric scalar. Number of individuals to simulate.
#' @param eta Numeric vector of length 3. Shape parameters for Weibull intensities with parameterization \eqn{\eta \nu t^{\nu - 1}}. Defaults to \code{rep(0.1, 3)}.
#' @param nu Numeric vector of length 3. Scale parameters for the Weibull hazards. Defaults to \code{rep(1.1, 3)}.
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{1}).
#' @param beta_L0_D Numeric scalar. Effect of baseline covariate L0 on death risk (default 1).
#' @param beta_L0_L Numeric scalar. Effect of baseline covariate L0 on covariate change risk (default 1).
#' @param beta_L_D Numeric scalar. Effect of covariate change (L = 1) on death risk (default 1).
#' @param beta_A0_D Numeric scalar. Effect of baseline treatment (A0 = 1) on death risk (default 0).
#' @param beta_A0_L Numeric scalar. Effect of baseline treatment (A0 = 1) on covariate change risk (default 0).
#' @param beta_L0_C Numeric scalar. Effect of baseline covariate L0 on censoring probability (default 0).
#' @param beta_A0_C Numeric scalar. Effect of baseline treatment A0 on censoring probability (default 0).
#' @param beta_L_C Numeric scalar. Effect of covariate change (L = 1) on censoring probability (default 0).
#' @param followup Numeric scalar. Maximum follow-up (censoring) time. Defaults to \code{Inf}.
#' @param lower Numeric scalar. Lower bound for root-finding (inverse cumulative hazard) (default \code{1e-15}).
#' @param upper Numeric scalar. Upper bound for root-finding (default 200).
#' @param beta_L_D_t_prime Numeric scalar or NULL. Additional effect of covariate change on death risk after time \code{t_prime} (optional).
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0. Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#' @param at_risk_cov Function. Function determining if an individual is at risk for each event type, given their covariates. Takes a numeric vector covariates and returns a binary vector. Default returns 1 for all events.
#'
#' @return A data frame containing the simulated data with columns:
#'  \item{ID}{Individual identifier}
#'  \item{Time}{Time of the event}
#'  \item{Delta}{Event type (0 = censoring, 1 = death, 2 = covariate change)}
#'  \item{L0}{Baseline covariate}
#'  \item{L}{Covariate indicating change in covariate process}
#'
#' @examples
#' simDisease(10)
#'
#' @export
simDisease <- function(N, eta = rep(0.1,3), nu = rep(1.1,3),  cens = 1,
                   beta_L0_D = 1, beta_L0_L = 1, beta_L_D = 1, beta_A0_D = 0,
                   beta_A0_L = 0, beta_L0_C = 0, beta_A0_C = 0, beta_L_C = 0,
                   followup = Inf, lower = 10^(-15), upper = 200,
                   beta_L_D_t_prime = NULL, t_prime = NULL, gen_A0 = NULL,
                   at_risk_cov = NULL){

  at_risk <- function(events) {
    return(c(
      cens,                                 # If you have not yet  been censored you are at risk
      1,                                    # If you have not died yet you are at risk
      as.numeric(events[3] == 0)))          # Only at risk for the covariate process if you have not experienced one yet
  }

  Time <- ID <- N0 <- N1 <- NULL
  beta <- matrix(ncol = 3, nrow = 5) # Three events and two baseline covariates

  # The effect of L0 on C, D, L
  beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_L)
  # The effect of A0 on C, D, L
  beta[2,] <- c(beta_A0_C, beta_A0_D, beta_A0_L)
  # Censoring process is a terminal process
  beta[3,] <- 0;
  # Death process is a terminal process
  beta[4,] <- 0;
  # The effect of L on C, D, L
  beta[5,] <- c(beta_L_C, beta_L_D, 0)

  if(!is.null(beta_L_D_t_prime) & !is.null(t_prime)){
    tv_eff <- matrix(0, ncol = 3, nrow = 5)
    tv_eff[5,2] <- beta_L_D_t_prime
    data <- simEventTV(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk,
                         max_cens = followup, lower = lower, upper = upper,
                         t_prime = t_prime, tv_eff = tv_eff, gen_A0 = gen_A0)
  }
  else{
    data <- simEventData(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk,
                         max_cens = followup, lower = lower, upper = upper,
                         gen_A0 = gen_A0, at_risk_cov = at_risk_cov)
  }

  # We don't need columns for terminal events
  data[, N0 := NULL]; data[, N1 := NULL]

  colnames(data)[6] <- c("L")

  return(data)
}
