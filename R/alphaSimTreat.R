#' Simulation and Estimation in Treatment Setting with Modified Eta Parameter of Treatment Process
#'
#' This function simulates event history data from the Treatment scenario (see function simTreatment).
#' The \eqn{\eta} parameter of the treatment process is multiplied by a factor \code{alpha}.
#' The function computes the proportion or years lost of death and treatment by time \eqn{\tau}.
#'
#' @title Simulation and Estimation in Treatment Setting with Modified Eta Parameter of Treatment Process
#'
#' @param N Integer. Number of individuals to simulate.
#' @param alpha Numeric. Multiplicative factor applied to the \eqn{\eta} parameter of the treatment process.
#' @param tau Numeric. Time point at which event proportions or years lost are compared.
#' @param plot Logical. If TRUE, output plots of the first 250 events for respectively intervened and no invtervend data.
#' @param eta Numeric vector of length 4. Shape parameters for Weibull hazards.
#' @param nu Numeric vector of length 4. Scale parameters for Weibull hazards.
#' @param beta_L_A Numeric. Effect of covariate L on treatment process A.
#' @param beta_L_D Numeric. Effect of covariate L on death process D.
#' @param beta_A_D Numeric. Effect of treatment process A on death process D.
#' @param beta_A_L Numeric. Effect of treatment process A on covariate L.
#' @param beta_L0_A Numeric. Effect of L0 on treatment process A.
#' @param lower Numeric. Lower bound for root-finding algorithm to invert cumulative hazard.
#' @param upper Numeric. Upper bound for root-finding algorithm.
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{1}).
#' @param return_data Logical. If \code{TRUE} the simulated data is returned.
#' @param years_lost Logical. If \code{TRUE}, computes years lost instead of proportions.
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
#'
#' @return A list containing:
#' \describe{
#'   \item{effect_A}{Proportion or years lost of treated patients with and without intervention.}
#'   \item{effect_death}{Proportion or years lost of deaths with and without intervention.}
#' }
#' Or the simulated data.
#'
#' @export
#'
#' @examples
#' alphaSimTreat()
alphaSimTreat <- function(N = 1e4,alpha = 0.5, tau = 5, plot = FALSE,
                          nu = rep(1.1, 4), eta = rep(0.1, 4), beta_L_A = 1,
                          beta_L_D = 1, beta_A_D = -0.5, beta_A_L = -1,
                          beta_L0_A = 1, lower = 10^(-300), upper = 300,
                          cens = 0, return_data = FALSE, years_lost = FALSE,
                          beta_L_A_prime = 0, beta_L_D_prime = 0, beta_A_D_prime = 0,
                          beta_L0_A_prime = 0, beta_A_L_prime = 0, beta_L0_L_prime = 0,
                          beta_L0_D_prime = 0, beta_L0_C_prime = 0, beta_L_C_prime = 0,
                          beta_A_C_prime = 0, t_prime = 0){

  Delta <- Time <- A0 <- V1 <- tmp <- NULL
  # Generate large data
  data <- simTreatment(N = N, cens = cens, eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                       nu = nu, beta_L_A = beta_L_A, beta_L_D = beta_L_D,
                       beta_A_D = beta_A_D, beta_A_L = beta_A_L, beta_L0_A = beta_L0_A,
                       lower = lower, upper = upper, beta_L_A_prime = 0, beta_L_D_prime = 0,
                       beta_A_D_prime = 0, beta_L0_A_prime = 0, beta_A_L_prime = 0,
                       beta_L0_L_prime = 0, beta_L0_D_prime = 0, beta_L0_C_prime = 0,
                       beta_L_C_prime = 0, beta_A_C_prime = 0, t_prime = NULL)

  if (plot) plotEventData(data[1:250])
  if(return_data) return(data)

  #Proportion of subjects dying before some time $\tau$ in treatment group
  prop_D <- data[Delta == 1, mean(Delta == 1 & Time < tau)]
  #Proportion of subjects experiencing treatment before some time \tau in a0 group
  prop_A <- mean(data[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])

  if (years_lost) {
    data[, tmp := cumsum((Delta == 1)*(tau - pmin(tau, Time))), by = "ID"]
    prop_D <- data[, tmp[.N], by = "ID"][, mean(V1)]
    data[, tmp := cumsum((Delta == 2)*(tau - pmin(tau, Time))), by = "ID"]
    prop_A <- data[, tmp[.N], by = "ID"][, mean(V1)]
  }

  return(list(effect_A = prop_A,
              effect_D = prop_D))
}

