#' Function to simulate data from a event history setting with a time dependent confounder.
#'
#' 4 different types of events are simulated, chosen to represent Censoring (0),Death (1),
#' Operation(2) and Change in Covariate Process(3). Death and Censoring are terminal events
#' and Operation and Change in Covariate Process can occur once.
#'
#' The intensities of the various events are controlled by the parameters beta parameters and
#' \eqn{\nu} and \eqn{\eta} parameters. For example `beta_L_A` determines the effect of L = 1
#' on the event A.
#'
#' @title Simulate Event History Data in a setting with a time-dependent confounder
#'
#' @param N A double of the number of individuals
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 0.1 for all events.
#' @param followup A maximal censoring time. By default set to infinity.
#' @param beta_L_A The effect of covariate L = 1 on the probability of A = 1
#' @param beta_L_D The effect of covariate L = 1 on the probability of D = 1
#' @param beta_A_D The effect of operation A = 1 on the probability of D = 1
#' @param beta_L0_A The effect of covariate L0 on the probability of A = 1
#' @param beta_A_L The effect of operation A = 1 on the probability of L = 1
#' @param beta_L0_L The effect of covariate L0 on the probability of L = 1
#' @param beta_L0_D The effect of covariate L0 on the probability of D = 1
#' @param beta_L0_C The effect of covariate L0 on the probability of C = 1
#' @param beta_L_C The effect of covariate L = 1 on the probability of C = 1
#' @param beta_A_C The effect of operation A = 1 on the probability of C = 1
#' @param cens Specifies whether you are at risk of being censored
#' @param op Specifies whether you are at risk of being operated
#'
#' @return  Data frame containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0), additional covariate (L) and  Treatment Process (A).
#' @export
#'
#' @examples
#' simConfounding(10)
#'
simConfounding <- function(N, beta_L_A = 1, beta_L_D = 1, beta_A_D = -1, beta_A_L = -0.5,
                           beta_L0_A = 1, beta_L0_L = 1, beta_L0_D = 1,
                           beta_L0_C = 0, beta_L_C = 0, beta_A_C = 0,
                           eta = rep(0.1,4), nu = rep(1.1,4),
                           followup = Inf, cens = 1, op = 1){

  Time <- A0 <- N0 <- N1 <- ID <- NULL

  if(op == 0){
    at_risk <- function(i, event_counts) {
      return(c(
        cens, # You might be at risk for censoring
        1, # If you have not died yet you are at risk for dying
        0, # You are never at risk for an operation
        as.numeric(event_counts[i,4] == 0))) # You are only at risk for a change in the covariate process if you have not experienced a change yet
    }
  }
  else{
    at_risk <- function(i, event_counts) {
      return(c(
        cens, # You might be at risk for censoring
        1, # If you have not died yet you are at risk for dying
        as.numeric(event_counts[i,3] == 0), # You are at risk for an operation if you have not had one yet
        as.numeric(event_counts[i,4] == 0))) # You are only at risk for a change in the covariate process if you have not experienced a change yet
    }
  }

  beta <- matrix(ncol = 4, nrow = 6)

  # The effect of L0 on the processes C, D, A, L
  beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_A, beta_L0_L)
  # A0 is 0
  beta[2,] <- 0
  # The processes C and D are terminal
  beta[c(3,4),] <- 0
  # The effect of A on the processes C, D, A, L
  beta[5,] <- c(beta_A_C, beta_A_D, 0, beta_A_L)
  # The effect of L on the processes C, D, A, L
  beta[6,] <- c(beta_L_C, beta_L_D, beta_L_A, 0)


  data <- simEventData(N, beta = beta, eta = eta, nu = nu, max_cens = followup,
                         at_risk = at_risk)
  data[, c("N0", "N1", "A0") := NULL]
  colnames(data)[c(5,6)] <- c("A", "L")

  return(data)
}
