#' Calculate the effect of an intervention modifying the Drop In process intensity.
#'
#' This function simulates event history data from the Drop In scenario. This is a
#' scenario with the events Censoring (C), Death (D), Drop In Initiation (Z),
#' Change in Covariate Process (L) and optionally Treatment (A). It simulates
#' under the intervnetion where \eqn{\eta} parameter of the Drop In process is multiplied
#' by a factor \code{alpha}. It evaluates the proportion of death and Drop In events
#' by time \eqn{\tau} within the subgroup defined by \code{A0 = a0}, comparing intervened
#' and non-intervened scenarios.
#'
#' @title Perform intervention on Drop In process intensity
#'
#' @param N Integer. Number of individuals to simulate.
#' @param alpha Numeric. Multiplicative factor applied to the \eqn{\eta} parameter of the Drop In process under intervention.
#' @param tau Numeric. Time point at which event proportions are compared.
#' @param a0 Binary (0 or 1). Group indicator to subset results.
#' @param plot Logical. If TRUE, plots of the first 250 events in each group are displayed.
#' @param eta Numeric vector of length 4. Shape parameters for the Weibull hazards (default length 4 for 4 processes).
#' @param nu Numeric vector of length 4. Scale parameters for the Weibull hazards.
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
#' @param adherence Logical. Whether to include a treatment adherence process (default FALSE).
#' @param lower Numeric. Lower bound for the root-finding algorithm to invert cumulative hazard.
#' @param upper Numeric. Upper bound for the root-finding algorithm to invert cumulative hazard.
#'
#' @return A list containing:
#' \describe{
#'   \item{effect_Z}{Proportion of subjects experiencing Drop In by time \eqn{\tau} with and without intervention.}
#'   \item{effect_death}{Proportion of subjects dying by time \eqn{\tau} with and without intervention.}
#' }
#'
#' @export
#'
#' @examples
#' intEffectAlphaDropIn()
intEffectAlphaDropIn <- function(N = 1e4,
                                 alpha = 0.5,
                                 tau = 5,
                                 a0 = 1,
                                 plot = FALSE,
                                 eta = rep(0.1, 4),
                                 nu = rep(1.1, 4),
                                 beta_L_A = 0, beta_L_Z = 1, beta_L_D = 0.5, beta_L_C = 0,
                                 beta_A_L = -0.5,  beta_A_Z = -0.5, beta_A_D = -1, beta_A_C = 0,
                                 beta_Z_L = -1, beta_Z_A = 0, beta_Z_D = -1, beta_Z_C = 0,
                                 beta_L0_L = 1, beta_L0_A = 1, beta_L0_Z = 0.1, beta_L0_D = 1, beta_L0_C = 0,
                                 beta_A0_L = -1, beta_A0_A = 0, beta_A0_Z = 0, beta_A0_D = -0.5, beta_A0_C = 0,
                                 adherence = FALSE,
                                 lower = 10^(-30),
                                 upper = 200){

  Delta <- Time <- A0 <- NULL

  # Generate large data set under the intervened intensity
  data_G1 <- simDropIn(N = N,
                       cens = 0,
                       generate.A0 = function(N,L0) rep(1,N),
                       eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                       nu = nu,
                       #followup = tau+1,
                       beta_L_A = beta_L_A, beta_L_Z = beta_L_Z, beta_L_D = beta_L_D, beta_L_C = beta_L_C,
                       beta_A_L = beta_A_L, beta_A_Z = beta_A_Z, beta_A_D = beta_A_D, beta_A_C = beta_A_C,
                       beta_Z_L = beta_Z_L, beta_Z_A = beta_Z_A, beta_Z_D = beta_Z_D, beta_Z_C = beta_Z_C,
                       beta_L0_L = beta_L0_L, beta_L0_A = beta_L0_A, beta_L0_Z = beta_L0_Z, beta_L0_D = beta_L0_D, beta_L0_C = beta_L0_C,
                       beta_A0_L = beta_A0_L, beta_A0_A = beta_A0_A, beta_A0_Z = beta_A0_Z, beta_A0_D = beta_A0_D, beta_A0_C = beta_A0_C,
                       lower = lower, upper = upper)

  # Generate large data set without intervened intensity
  data_G2 <- simDropIn(N = N,
                       cens = 0,
                       generate.A0 = function(N,L0) rep(0,N),
                       eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                       nu = nu,
                       #followup = tau+1,
                       beta_L_A = beta_L_A, beta_L_Z = beta_L_Z, beta_L_D = beta_L_D, beta_L_C = beta_L_C,
                       beta_A_L = beta_A_L, beta_A_Z = beta_A_Z, beta_A_D = beta_A_D, beta_A_C = beta_A_C,
                       beta_Z_L = beta_Z_L, beta_Z_A = beta_Z_A, beta_Z_D = beta_Z_D, beta_Z_C = beta_Z_C,
                       beta_L0_L = beta_L0_L, beta_L0_A = beta_L0_A, beta_L0_Z = beta_L0_Z, beta_L0_D = beta_L0_D, beta_L0_C = beta_L0_C,
                       beta_A0_L = beta_A0_L, beta_A0_A = beta_A0_A, beta_A0_Z = beta_A0_Z, beta_A0_D = beta_A0_D, beta_A0_C = beta_A0_C,
                       lower = lower, upper = upper)

  if (plot) gridExtra::grid.arrange(plotEventData(data_G1[1:250]),
                         plotEventData(data_G2[1:250]), nrow = 1)

  #Proportion of subjects dying before some time $\tau$ in a0 group
  prop_G1 <- data_G1[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)] # with intervention
  prop_G2 <- data_G2[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)] # without intervention
  #Proportion of subjects experiencing Drop In before some time \tau in a0 group
  prop_G1_Z <- mean(data_G1[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention
  prop_G2_Z <- mean(data_G2[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # without intervention

  return(list(effect_Z = c(G1 = prop_G1_Z, G2 = prop_G2_Z),
              effect_death = c(G1 = prop_G1, G2 = prop_G2)))

}
