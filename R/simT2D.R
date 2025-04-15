#' Function to simulate data from a T2D diabetes setting. 3 different types of
#' events are simulated. They can be interpreted as Censoring(0), Death (1) and
#' Change in Covariate Process (3). Death and Censoring are terminal events and Change in
#' Covariate Process can occur once. The intensities of the various events depend
#' upon previous events and the pre specified \eqn{\nu} and \eqn{\eta} parameters.
#' The dependence on previous events is controlled by parameters chosen so that a
#' large baseline covariate (L0) increases the probability of covariate change (L = 1).
#' Initial treatment (A0 = 1) reduces the risk of a change in the covariate process
#' (L= 1) and the risk of death. A large baseline covariate (L0) increase the risk of death.
#' Censoring does not depend on anything. The risk of death depends on change in
#' the covariate process (L = 1) through beta_L_D.
#'
#'
#' @title Simulate Data in a T2D Diabetes Setting
#'
#' @param N A double of the number of individuals
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 1.1 for all events.
#' @param followup A maximal censoring time. By default set to infinity.
#' @param beta_L_D Specifies how change in the covariate process affects risk of death. Is by default set to 1.
#' @param beta_L0_L Specifies how change in the covariate process affects risk of death. Is by default set to 1.
#' @param beta_A0_D Specifies how baseline treatment affects risk of death. Is by default set to 0.
#' @param beta_A0_L Specifies how baseline treatment affects risk of T2D. Is by default set to 0.
#' @param beta_L0_D Specifies how baseline covariate affects risk of Death. Is by default set to 1.
#' @param beta_L0_C The effect of covariate L0 on the probability of C = 1. Is by default set to 0.
#' @param beta_A0_C The effect of covariate A0 = 1 on the probability of C = 1. Is by default set to 0.
#' @param beta_L_C The effect of covariate L = 1 on the probability of C = 1. Is by default set to 0.
#' @param cens Specifies whether you are at risk of being censored
#'
#' @return  Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0) and additional covariate (L).
#' @export
#'
#' @examples
#' simT2D(10)
simT2D <- function(N, eta = rep(0.1,3), nu = rep(1.1,3),  cens = 1,
                   beta_L0_D = 1, beta_L0_L = 1, beta_L_D = 1, beta_A0_D = 0,
                   beta_A0_L = 0, beta_L0_C = 0, beta_A0_C = 0, beta_L_C = 0,
                   followup = Inf){

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
  # How A0 on C, D, L
  beta[2,] <- c(beta_A0_C, beta_A0_D, beta_A0_L)
  # Censoring process is a terminal process
  beta[3,] <- 0;
  # Death process is a terminal process
  beta[4,] <- 0;
  # The effect of L on C, D, L
  beta[5,] <- c(beta_L_C, beta_L_D, 0)

  data <- simEventData(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk,
                         max_cens = followup)

  # We don't need columns for terminal events
  data[, N0 := NULL]; data[, N1 := NULL]

  colnames(data)[6] <- c("L")

  return(data)
}
