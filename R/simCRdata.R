#' Function to simulate data from a Competing Risk setting. The function simulates data corresponding to $N$ individuals that
#' are at risk for mutually exclusive types of failure. 3 events can take place, one of which can be interpreted as censoring.
#'
#' @title Simulate Competing risk Data
#'
#' @param N A double of the number of individuals
#' @param beta A 2X3 matrix with the effect of L0 and A0 on the three processes.
#' The columns represent Process 1, Process 2 and Process 3, while the rows represent L0 and A0.
#' @param eta Vector of  length 3 of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of length 3 of scale parameters for the Weibull hazard.
#'
#' @return  Data frame containing the simulated competing risk data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0) and baseline Treatment (A0).
#' @export
#'
#' @examples
#' simCRdata(10)
simCRdata <- function(N,
                      beta = NULL,
                      eta = rep(0.1,3),
                      nu = rep(1.1,3)
){

  at_risk <- function(events) c(1,1,1)

  if(is.null(beta)){
    beta <- matrix(0, ncol = 3, nrow = 2)
  }

  beta <- rbind(beta, matrix(0, ncol = 3, nrow = 3))
  results <- simEventData(N, beta = beta, eta = eta, nu = nu, at_risk = at_risk,
                          term_deltas = c(0,1,2))
  results <- results[, !c("N0", "N1", "N2")]

  return(results)
}
