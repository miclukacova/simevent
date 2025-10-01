#' Simulate Event Data from a "Drop In" Setting
#'
#' `simDropIn` is a function that simulates data corresponding to \code{N} individuals
#' that are at risk for 4 or 5 events. Censoring (C), Death (D), Drop In Initiation (Z),
#' Change in Covariate Process (L) and optionally Treatment (A).
#'
#' @title `simDropIn`
#'
#' @param N Integer. Number of individuals to simulate.
#' @param beta_L_A Numeric. Specifies how L affects A.
#' @param beta_L_Z Numeric. Specifies how L affects Z.
#' @param beta_L_D Numeric. Specifies how L affects D.
#' @param beta_L_C Numeric. Specifies how L affects C.
#' @param beta_A_L Numeric. Specifies how L affects A.
#' @param beta_A_Z Numeric. Specifies how L affects Z.
#' @param beta_A_D Numeric. Specifies how L affects D.
#' @param beta_A_C Numeric. Specifies how L affects C.
#' @param beta_Z_L Numeric. Specifies how L affects A.
#' @param beta_Z_A Numeric. Specifies how L affects Z.
#' @param beta_Z_D Numeric. Specifies how L affects D.
#' @param beta_Z_C Numeric. Specifies how L affects C.
#' @param beta_L0_L Numeric. Specifies how L affects A.
#' @param beta_L0_A Numeric. Specifies how L affects Z.
#' @param beta_L0_Z Numeric. Specifies how L affects Z.
#' @param beta_L0_D Numeric. Specifies how L affects D.
#' @param beta_L0_C Numeric. Specifies how L affects C.
#' @param beta_A0_L Numeric. Specifies how L affects A.
#' @param beta_A0_A Numeric. Specifies how L affects Z.
#' @param beta_A0_Z Numeric. Specifies how L affects Z.
#' @param beta_A0_D Numeric. Specifies how L affects D.
#' @param beta_A0_C Numeric. Specifies how L affects C.
#' @param eta Numeric vector of length 4 (or 5). Shape parameters of the Weibull baseline intensity for each event type.
#' \deqn{\eta \nu t^{\nu - 1}}.
#' @param nu Numeric vector of length 4 (or 5). Scale parameters for the Weibull hazard.
#' @param adherence Logical. Indicator of whether a Treatment process should be simulated.
#' @param followup Numeric. Maximum censoring time. Events occurring after this time are censored. Default is Inf (no censoring).
#' @param cens Logical. Indicator of whether there should be a censoring process.
#' @param generate.A0 Function. Function to generate the baseline treatment covariate A0. Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#' @param lower Numeric. Lower bound for root-finding in inverse cumulative hazard calculations. Default is \eqn{10^{-15}}.
#' @param upper Numeric. Upper bound for root-finding in inverse cumulative hazard calculations. Default is 200.
#' @param beta_L_A_prime Numeric. Specifies how L additionally affects A after time t_prime.
#' @param beta_L_Z_prime Numeric. Specifies how L additionally affects Z after time t_prime.
#' @param beta_L_D_prime Numeric. Specifies how L additionally affects D after time t_prime.
#' @param beta_L_C_prime Numeric. Specifies how L additionally affects C after time t_prime.
#' @param beta_A_L_prime Numeric. Specifies how L additionally affects A after time t_prime.
#' @param beta_A_Z_prime Numeric. Specifies how L additionally affects Z after time t_prime.
#' @param beta_A_D_prime Numeric. Specifies how L additionally affects D after time t_prime.
#' @param beta_A_C_prime Numeric. Specifies how L additionally affects C after time t_prime.
#' @param beta_Z_L_prime Numeric. Specifies how L additionally affects A after time t_prime.
#' @param beta_Z_A_prime Numeric. Specifies how L additionally affects Z after time t_prime.
#' @param beta_Z_D_prime Numeric. Specifies how L additionally affects D after time t_prime.
#' @param beta_Z_C_prime Numeric. Specifies how L additionally affects C after time t_prime.
#' @param beta_L0_L_prime Numeric. Specifies how L additionally affects after time A.
#' @param beta_L0_A_prime Numeric. Specifies how L additionally affects after time Z.
#' @param beta_L0_Z_prime Numeric. Specifies how L additionally affects after time Z.
#' @param beta_L0_D_prime Numeric. Specifies how L additionally affects after time D.
#' @param beta_L0_C_prime Numeric. Specifies how L additionally affects after time C.
#' @param beta_A0_L_prime Numeric. Specifies how L additionally affects after time A.
#' @param beta_A0_A_prime Numeric. Specifies how L additionally affects after time Z.
#' @param beta_A0_Z_prime Numeric. Specifies how L additionally affects after time Z.
#' @param beta_A0_D_prime Numeric. Specifies how L additionally affects after time D.
#' @param beta_A0_C_prime Numeric. Specifies how L additionally affects after time C.
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#'
#' @return  Data frame containing the simulated event history data
#' @export
#'
#' @examples
#' simDropIn(10)
simDropIn <- function(N,
                      beta_L_A = 1, beta_L_Z = 2, beta_L_D = 1.5, beta_L_C = 0,
                      beta_A_L = -0.5,  beta_A_Z = -0.5, beta_A_D = -1, beta_A_C = 0,
                      beta_Z_L = -1, beta_Z_A = 0, beta_Z_D = -1, beta_Z_C = 0,
                      beta_L0_L = 1, beta_L0_A = 1, beta_L0_Z = 1, beta_L0_D = 1, beta_L0_C = 0,
                      beta_A0_L = -1.5, beta_A0_A = 0, beta_A0_Z = 0, beta_A0_D = -2, beta_A0_C = 0,
                      eta = c(0.5, 0.5, 0.1, 0.25),
                      nu = c(1.1, 1.1, 1.1, 1.1),
                      adherence = FALSE,
                      followup = Inf,
                      cens = 1,
                      generate.A0 = function(N, L0) stats::rbinom(N, 1, 0.5),
                      lower = 1e-200, upper = 1e10,
                      beta_L_A_prime = 0, beta_L_Z_prime = 0, beta_L_D_prime = 0,
                      beta_L_C_prime = 0, beta_A_L_prime = 0, beta_A_Z_prime = 0,
                      beta_A_D_prime = 0, beta_A_C_prime = 0, beta_Z_L_prime = 0,
                      beta_Z_A_prime = 0, beta_Z_D_prime = 0, beta_Z_C_prime = 0,
                      beta_L0_L_prime = 0,beta_L0_A_prime = 0,beta_L0_Z_prime = 0,
                      beta_L0_D_prime = 0,beta_L0_C_prime = 0,beta_A0_L_prime = 0,
                      beta_A0_A_prime = 0,beta_A0_Z_prime = 0,beta_A0_D_prime = 0,
                      beta_A0_C_prime = 0, t_prime = NULL){

  N0 <- N1 <- ID <- NULL

  if (adherence) { #-- should add one more process
    if (length(eta) == 4) eta <- c(eta,eta[4])
    if (length(nu) == 4) nu <- c(nu,nu[4])
    at_risk <- function(events) {
      return(c(
        cens, # You might be at risk for censoring
        1, # If you have not died yet you are at risk for dying
        as.numeric(events[3] == 0), # You are at risk for drop-in initiation if you have not initiated yet
        as.numeric(events[4] == 0), # You are only at risk for a change in the covariate process if you have not experienced a change yet
        as.numeric(events[5] == 0))) # You are only at risk for a change in the treatment process if you have not experienced a change yet
    }
  } else {
    at_risk <- function(events) {
      return(c(
        cens, # You might be at risk for censoring
        1, # If you have not died yet you are at risk for dying
        as.numeric(events[3] == 0), # You are at risk for drop-in initiation if you have not initiated yet
        as.numeric(events[4] == 0))) # You are only at risk for a change in the covariate process if you have not experienced a change yet
    }
  }

  beta <- matrix(ncol = length(eta), nrow = 2+length(eta))
  # The effect of L0 on the processes C, D, Z, L
  beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_Z, beta_L0_L)
  # The effect of A0 on the processes C, D, Z, L
  beta[2,] <- c(beta_A0_C, beta_A0_D, beta_A0_Z, beta_A0_L)
  # The processes C and D are terminal
  beta[c(3,4),] <- 0
  # The effect of Z on the processes C, D, Z, L
  beta[5,] <- c(beta_Z_C, beta_Z_D, 0, beta_Z_L)
  # The effect of L on the processes C, D, Z, L
  beta[6,] <- c(beta_L_C, beta_L_D, beta_L_Z, 0)

  if (length(eta) == 5) {
    beta <- rbind(beta, c(beta_A_C, beta_A_D, beta_A_Z, beta_A_L))
    beta <- cbind(beta, c(beta_L0_A, beta_A0_A, 0, 0, beta_Z_A, beta_L_A, 0))
  }


  if(!is.null(t_prime)){
    beta_prime <- matrix(ncol = length(eta), nrow = 2 + length(eta))
    # The additional effect of L0 on the processes C, D, Z, L after t_prime
    beta_prime[1,] <- c(beta_L0_C_prime, beta_L0_D_prime, beta_L0_Z_prime, beta_L0_L_prime)
    # The additional effect of A0 on the processes C, D, Z, L after t_prime
    beta_prime[2,] <- c(beta_A0_C_prime, beta_A0_D_prime, beta_A0_Z_prime, beta_A0_L_prime)
    # The processes C and D are terminal
    beta_prime[c(3,4),] <- 0
    # The additional effect of Z on the processes C, D, Z, L after t_prime
    beta_prime[5,] <- c(beta_Z_C_prime, beta_Z_D_prime, 0, beta_Z_L_prime)
    # The additional effect of L on the processes C, D, Z, L after t_prime
    beta_prime[6,] <- c(beta_L_C_prime, beta_L_D_prime, beta_L_Z_prime, 0)
    if (length(eta) == 5) {
      beta_prime <- rbind(beta_prime, c(beta_A_C_prime, beta_A_D_prime, beta_A_Z_prime, beta_A_L_prime))
      beta_prime <- cbind(beta_prime, c(beta_L0_A_prime, beta_A0_A_prime, 0, 0, beta_Z_A_prime, beta_L_A_prime, 0))
    }

  data <- simEventTV(N,
                     beta = beta,
                     eta = eta,
                     nu = nu,
                     max_cens = followup,
                     at_risk = at_risk,
                     gen_A0 = generate.A0,
                     lower = lower,
                     upper = upper,
                     tv_eff = beta_prime,
                     t_prime = t_prime)
  }
  else{
    data <- simEventData(N,
                         beta = beta,
                         eta = eta,
                         nu = nu,
                         max_cens = followup,
                         at_risk = at_risk,
                         gen_A0 = generate.A0,
                         lower = lower,
                         upper = upper)
  }

  setnames(data, c("N2", "N3"), c("Z", "L"))
  if (length(eta)>4) setnames(data, c("N4"), c("A"))

  data[, N0 := NULL]
  data[, N1 := NULL]

  return(data)
}



