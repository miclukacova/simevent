#' Simulation and Estimation in Drop In Setting with Modified Eta Parameter of Drop In Process
#'
#' This function simulates event history data from the Drop In scenario. This is a
#' scenario with the events Censoring (C), Death (D), Drop In Initiation (Z),
#' Change in Covariate Process (L) and optionally Treatment (A). It simulates
#' under where \eqn{\eta} parameter of the Drop In process is multiplied by a factor \code{alpha}.
#' It evaluates the proportion of death and Drop In events by time \eqn{\tau} within
#' the subgroup defined by \code{A0 = a0}.
#'
#' @title Simulation and Estimation in Drop In Setting with Modified Eta Parameter of Drop In Process
#'
#' @param N Integer. Number of individuals to simulate.
#' @param alpha Numeric. Multiplicative factor applied to the \eqn{\eta} parameter of the Drop In process under intervention.
#' @param tau Numeric. Time point at which event proportions are compared.
#' @param a0 Binary (0 or 1). Group indicator to subset results.
#' @param plot Logical. If TRUE, plots of the first 250 events in each group are displayed.
#' @param eta Numeric vector of length 4. Shape parameters for the Weibull hazards (default length 4 for 4 processes).
#' @param nu Numeric vector of length 4. Scale parameters for the Weibull hazards.
#' @param adherence Logical. Whether to include a treatment adherence process (default FALSE).
#' @param lower Numeric. Lower bound for the root-finding algorithm to invert cumulative hazard.
#' @param upper Numeric. Upper bound for the root-finding algorithm to invert cumulative hazard.
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{1}).
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0. Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#' @param return_data Logical. If \code{TRUE} the simulated data is returned.
#' @param years_lost Logical. If \code{TRUE}, computes years lost instead of proportions.
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#' @param beta_L_A Numeric. Effect of process L on process A.
#' @param beta_L_Z Numeric. Effect of process L on process Z.
#' @param beta_L_D Numeric. Effect of process L on process D.
#' @param beta_L_C Numeric. Effect of process L on process C.
#' @param beta_A_L Numeric. Effect of process A on process L.
#' @param beta_A_Z Numeric. Effect of process A on process Z.
#' @param beta_A_D Numeric. Effect of process A on process D.
#' @param beta_A_C Numeric. Effect of process A on process C.
#' @param beta_Z_L Numeric. Effect of process Z on process L.
#' @param beta_Z_A Numeric. Effect of process Z on process A.
#' @param beta_Z_D Numeric. Effect of process Z on process D.
#' @param beta_Z_C Numeric. Effect of process Z on process C.
#' @param beta_L0_L Numeric. Effect of baseline covariate L0 on process L.
#' @param beta_L0_A Numeric. Effect of baseline covariate L0 on process A.
#' @param beta_L0_Z Numeric. Effect of baseline covariate L0 on process Z.
#' @param beta_L0_D Numeric. Effect of baseline covariate L0 on process D.
#' @param beta_L0_C Numeric. Effect of baseline covariate L0 on process C.
#' @param beta_A0_L Numeric. Effect of baseline covariate A0 on process L.
#' @param beta_A0_A Numeric. Effect of baseline covariate A0 on process A.
#' @param beta_A0_Z Numeric. Effect of baseline covariate A0 on process Z.
#' @param beta_A0_D Numeric. Effect of baseline covariate A0 on process D.
#' @param beta_A0_C Numeric. Effect of baseline covariate A0 on process C.
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
#' @param beta_L0_L_prime Numeric. Specifies how L0 additionally affects L after time t_prime.
#' @param beta_L0_A_prime Numeric. Specifies how L0 additionally affects A after time t_prime.
#' @param beta_L0_Z_prime Numeric. Specifies how L0 additionally affects Z after time t_prime.
#' @param beta_L0_D_prime Numeric. Specifies how L0 additionally affects D after time t_prime.
#' @param beta_L0_C_prime Numeric. Specifies how L0 additionally affects C after time t_prime.
#' @param beta_A0_L_prime Numeric. Specifies how A0 additionally affects L after time t_prime.
#' @param beta_A0_A_prime Numeric. Specifies how A0 additionally affects A after time t_prime.
#' @param beta_A0_Z_prime Numeric. Specifies how A0 additionally affects Z after time t_prime.
#' @param beta_A0_D_prime Numeric. Specifies how A0 additionally affects D after time t_prime.
#' @param beta_A0_C_prime Numeric. Specifies how A0 additionally affects C after time t_prime.
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#'
#' @return A list containing:
#' \describe{
#'   \item{effect_Z}{Proportion of subjects experiencing Drop In by time \eqn{\tau} with and without intervention.}
#'   \item{effect_death}{Proportion of subjects dying by time \eqn{\tau} with and without intervention.}
#' }
#' Or the simulated data.
#'
#' @export
#'
#' @examples
#' alphaSimDropIn()
alphaSimDropIn <- function(N = 1e4, alpha = 0.5, tau = 5, a0 = 1, plot = FALSE,
                           eta = rep(0.1, 4), nu = rep(1.1, 4), adherence = FALSE,
                           lower = 10^(-30), upper = 200, cens = 0, gen_A0 = NULL,
                           return_data = FALSE, years_lost = FALSE,
                           beta_L_A = 0, beta_L_Z = 1, beta_L_D = 0.5, beta_L_C = 0,
                           beta_A_L = -0.5,  beta_A_Z = -0.5, beta_A_D = -1, beta_A_C = 0,
                           beta_Z_L = -1, beta_Z_A = 0, beta_Z_D = -1, beta_Z_C = 0,
                           beta_L0_L = 1, beta_L0_A = 1, beta_L0_Z = 0.1, beta_L0_D = 1, beta_L0_C = 0,
                           beta_A0_L = -1, beta_A0_A = 0, beta_A0_Z = 0, beta_A0_D = -0.5, beta_A0_C = 0,
                           beta_L_A_prime = 0, beta_L_Z_prime = 0, beta_L_D_prime = 0,
                           beta_L_C_prime = 0, beta_A_L_prime = 0, beta_A_Z_prime = 0,
                           beta_A_D_prime = 0, beta_A_C_prime = 0, beta_Z_L_prime = 0,
                           beta_Z_A_prime = 0, beta_Z_D_prime = 0, beta_Z_C_prime = 0,
                           beta_L0_L_prime = 0,beta_L0_A_prime = 0,beta_L0_Z_prime = 0,
                           beta_L0_D_prime = 0,beta_L0_C_prime = 0,beta_A0_L_prime = 0,
                           beta_A0_A_prime = 0,beta_A0_Z_prime = 0,beta_A0_D_prime = 0,
                           beta_A0_C_prime = 0, t_prime = NULL){

  Delta <- Time <- A0 <- V1 <- tmp <- NULL

  # Generate large data set under the intervened intensity
  data <- simDropIn(N = N,
                    cens = cens,
                    generate.A0 = gen_A0,
                    eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                    nu = nu,
                    beta_L_A = beta_L_A, beta_L_Z = beta_L_Z, beta_L_D = beta_L_D, beta_L_C = beta_L_C,
                    beta_A_L = beta_A_L, beta_A_Z = beta_A_Z, beta_A_D = beta_A_D, beta_A_C = beta_A_C,
                    beta_Z_L = beta_Z_L, beta_Z_A = beta_Z_A, beta_Z_D = beta_Z_D, beta_Z_C = beta_Z_C,
                    beta_L0_L = beta_L0_L, beta_L0_A = beta_L0_A, beta_L0_Z = beta_L0_Z, beta_L0_D = beta_L0_D, beta_L0_C = beta_L0_C,
                    beta_A0_L = beta_A0_L, beta_A0_A = beta_A0_A, beta_A0_Z = beta_A0_Z, beta_A0_D = beta_A0_D, beta_A0_C = beta_A0_C,
                    beta_L_A_prime = beta_L_A_prime, beta_L_Z_prime = beta_L_Z_prime, beta_L_D_prime = beta_L_D_prime,
                    beta_L_C_prime = beta_L_C_prime, beta_A_L_prime = beta_A_L_prime, beta_A_Z_prime = beta_A_Z_prime,
                    beta_A_D_prime = beta_A_D_prime, beta_A_C_prime = beta_A_C_prime, beta_Z_L_prime = beta_Z_L_prime,
                    beta_Z_A_prime = beta_Z_A_prime, beta_Z_D_prime = beta_Z_D_prime, beta_Z_C_prime = beta_Z_C_prime,
                    lower = lower, upper = upper, t_prime = t_prime)

  if (plot) plotEventData(data[1:250])
  if(return_data) return(data)

  #Proportion of subjects dying before some time $\tau$ in a0 group
  prop_D <- data[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)]
  #Proportion of subjects experiencing Drop In before some time \tau in a0 group
  prop_Z <- mean(data[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])

  if (years_lost) {
    data[, tmp := cumsum((Delta == 1)*(tau - pmin(tau, Time))), by = "ID"]
    prop_D <- data[, tmp[.N], by = "ID"][, mean(V1)]
    data[, tmp := cumsum((Delta == 2)*(tau - pmin(tau, Time))), by = "ID"]
    prop_Z <- data[, tmp[.N], by = "ID"][, mean(V1)]
  }


  return(list(effect_Z = prop_Z, effect_D = prop_D))

}
