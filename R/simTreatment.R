#' Simulate Event History Data with Treatment and Time-Dependent Covariate
#'
#' Simulates event history data with four types of events representing censoring (0), death (1), treatment (2), and covariate change (3).
#' Death and censoring are terminal events; treatment and covariate events can occur only once.
#'
#' Event intensities are modeled using Weibull hazards with parameters \eqn{\nu} (scale) and \eqn{\eta} (shape),
#' and covariate effects controlled by specified \code{beta} parameters. For example, \code{beta_L_A} quantifies
#' the effect of covariate \code{L = 1} on the hazard of treatment.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param eta Numeric vector of length 4. Shape parameters for Weibull intensities, parameterized as
#'  \eqn{\eta \nu t^{\nu - 1}}. Default is \code{rep(0.1, 4)}.
#' @param nu Numeric vector of length 4. Scale parameters for Weibull hazards. Default is \code{rep(1.1, 4)}.
#' @param beta_L_A Numeric. Effect of covariate \code{L = 1} on treatment hazard. Default 1.
#' @param beta_L_D Numeric. Effect of covariate \code{L = 1} on death hazard. Default 1.
#' @param beta_A_D Numeric. Effect of treatment \code{A = 1} on death hazard. Default -1.
#' @param beta_L0_A Numeric. Effect of baseline covariate \code{L0} on treatment hazard. Default 1.
#' @param beta_A_L Numeric. Effect of treatment \code{A = 1} on covariate hazard. Default -0.5.
#' @param beta_L0_L Numeric. Effect of baseline covariate \code{L0} on covariate hazard. Default 1.
#' @param beta_L0_D Numeric. Effect of baseline covariate \code{L0} on death hazard. Default 1.
#' @param beta_L0_C Numeric. Effect of baseline covariate \code{L0} on censoring hazard. Default 0.
#' @param beta_L_C Numeric. Effect of covariate \code{L = 1} on censoring hazard. Default 0.
#' @param beta_A_C Numeric. Effect of treatment \code{A = 1} on censoring hazard. Default 0.
#' @param beta_L_A_prime Numeric. Additional effect of covariate \code{L = 1} on treatment hazard. Default 0.
#' @param beta_L_D_prime Numeric. Additionalffect of covariate \code{L = 1} on death hazard. Default 0.
#' @param beta_A_D_prime Numeric. Effect of treatment \code{A = 1} on death hazard. Default 0.
#' @param beta_L0_A_prime Numeric. Effect of baseline covariate \code{L0} on treatment hazard. Default 0.
#' @param beta_A_L_prime Numeric. Effect of treatment \code{A = 1} on covariate hazard. Default 0.
#' @param beta_L0_L_prime Numeric. Effect of baseline covariate \code{L0} on covariate hazard. Default 0.
#' @param beta_L0_D_prime Numeric. Effect of baseline covariate \code{L0} on death hazard. Default 0.
#' @param beta_L0_C_prime Numeric. Effect of baseline covariate \code{L0} on censoring hazard. Default 0.
#' @param beta_L_C_prime Numeric. Effect of covariate \code{L = 1} on censoring hazard. Default 0.
#' @param beta_A_C_prime Numeric. Effect of treatment \code{A = 1} on censoring hazard. Default 0.
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#' @param at_risk_cov Function. Function determining if an individual is at risk for each event type,
#' given their covariates. Takes a numeric vector covariates and returns a binary vector. Default returns 1 for all events.
#' @param cens Integer (0 or 1). Indicates if censoring is possible. Default 1.
#' @param op Integer (0 or 1). Indicates if treatment (operation) is possible. Default 1.
#' @param lower Numeric. Lower bound for root finding (inverse cumulative hazard). Default \code{1e-15}.
#' @param upper Numeric. Upper bound for root finding (inverse cumulative hazard). Default 200.
#' @param followup Numeric. Maximum censoring time. Defaults to \code{Inf} (no censoring).
#' @param ... Additional arguments passed to \code{simEventData} or \code{simEventTV}
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{ID} - Individual identifier.
#'   \item \code{Time} - Event time.
#'   \item \code{Delta} - Event type (0=censoring, 1=death, 2=treatment, 3=covariate change).
#'   \item \code{L0} - Baseline covariate.
#'   \item \code{L} - Time-dependent covariate.
#'   \item \code{A} - Treatment status.
#' }
#'
#' @export
#'
#' @examples
#' simTreatment(10)
simTreatment <- function(N,
                         eta = rep(0.1,4),
                         nu = rep(1.1,4),
                         beta_L_A = 1,
                         beta_L_D = 1,
                         beta_A_D = -1,
                         beta_A_L = -0.5,
                         beta_L0_A = 1,
                         beta_L0_L = 1,
                         beta_L0_D = 1,
                         beta_L0_C = 0,
                         beta_L_C = 0,
                         beta_A_C = 0,
                         beta_L_A_prime = 0,
                         beta_L_D_prime = 0,
                         beta_A_D_prime = 0,
                         beta_A_L_prime = 0,
                         beta_L0_A_prime = 0,
                         beta_L0_L_prime = 0,
                         beta_L0_D_prime = 0,
                         beta_L0_C_prime = 0,
                         beta_L_C_prime = 0,
                         beta_A_C_prime = 0,
                         t_prime = NULL,
                         at_risk_cov = NULL,
                         cens = 1,
                         op = 1,
                         lower = 10^(-15),
                         upper = 200,
                         followup = Inf,
                         ...){

  Time <- A0 <- N0 <- N1 <- ID <- NULL

  if(op == 0){
    at_risk <- function(events) {
      return(c(
        cens, # You might be at risk for censoring
        1, # If you have not died yet you are at risk for dying
        0, # You are never at risk for an treatment
        as.numeric(events[4] == 0))) # You are only at risk for a change in the covariate process if you have not experienced a change yet
    }
  }
  else{
    at_risk <- function(events) {
      return(c(
        cens, # You might be at risk for censoring
        1, # If you have not died yet you are at risk for dying
        as.numeric(events[3] == 0), # You are at risk for an treatment if you have not had one yet
        as.numeric(events[4] == 0))) # You are only at risk for a change in the covariate process if you have not experienced a change yet
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


  if(!is.null(t_prime)){
    tv_eff <- matrix(ncol = 4, nrow = 6)
    # The change in effect of L0 on the processes C, D, A, L
    tv_eff[1,] <- c(beta_L0_C_prime, beta_L0_D_prime, beta_L0_A_prime, beta_L0_L_prime)
    # A0 is 0
    tv_eff[2,] <- 0
    # The processes C and D are terminal
    tv_eff[c(3,4),] <- 0
    # The change in effect of A on the processes C, D, A, L
    tv_eff[5,] <- c(beta_A_C_prime, beta_A_D_prime, 0, beta_A_L_prime)
    # The change in effect of L on the processes C, D, A, L
    tv_eff[6,] <- c(beta_L_C_prime, beta_L_D_prime, beta_L_A_prime, 0)

    data <- simEventTV(N, beta = beta, eta = eta, nu = nu, max_cens = followup,
                       at_risk = at_risk, t_prime = t_prime, tv_eff = tv_eff,
                       at_risk_cov = at_risk_cov, upper = upper, lower = lower, ...)
  }
  else{
    data <- simEventData(N, beta = beta, eta = eta, nu = nu, max_cens = followup,
                         at_risk = at_risk, at_risk_cov = at_risk_cov, upper = upper,
                         lower = lower, ...)
  }


  data[, c("N0", "N1", "A0") := NULL]
  colnames(data)[c(5,6)] <- c("A", "L")

  return(data)
}
