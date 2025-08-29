#' Function to calculate the effect of an intervention where the eta parameter of
#' the Drop In process is multiplied by the factor alpha. We consider the outcome
#' variable death and Drop In at time \eqn{\tau} - that is the output of the function is
#' the proportion of dead and T2D diagnosed patients by time \eqn{\tau} in A0 = a0
#' group in respectively intervened data and non intervened data.
#'
#' @title Perform intervention
#'
#' @param N A double of the number of individuals to simulate
#' @param alpha Double corresponding to the factor to be multiplied onto eta
#' parameter of the T2D process
#' @param tau A double deciding the time we compare proportions.
#' @param a0 A binary variable deciding whether we are comparing treatment group
#' or placebo group.
#' @param plot Logical indicating whether a plot should be outputed
#' @param eta Vector of  length 3 of shape parameters for the Weibull hazard
#' with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}.
#' @param nu Vector of length 3 of scale parameters for the Weibull hazard.
#' @param beta_L_A Specifies how L affects A.
#' @param beta_L_Z Specifies how L affects Z.
#' @param beta_L_D Specifies how L affects D.
#' @param beta_L_C Specifies how L affects C.
#' @param beta_A_L Specifies how L affects A.
#' @param beta_A_Z Specifies how L affects Z.
#' @param beta_A_D Specifies how L affects D.
#' @param beta_A_C Specifies how L affects C.
#' @param beta_Z_L Specifies how L affects A.
#' @param beta_Z_A Specifies how L affects Z.
#' @param beta_Z_D Specifies how L affects D.
#' @param beta_Z_C Specifies how L affects C.
#' @param beta_L0_L Specifies how L affects A.
#' @param beta_L0_A Specifies how L affects Z.
#' @param beta_L0_Z Specifies how L affects Z.
#' @param beta_L0_D Specifies how L affects D.
#' @param beta_L0_C Specifies how L affects C.
#' @param beta_A0_L Specifies how L affects A.
#' @param beta_A0_A Specifies how L affects Z.
#' @param beta_A0_Z Specifies how L affects Z.
#' @param beta_A0_D Specifies how L affects D.
#' @param beta_A0_C Specifies how L affects C.
#' @param adherence Logical indication whether a Treatment process should be added.
#' @param lower Lower bound for the uniroot function used to find the inverse
#' cumulative hazard.
#' @param upper Upper bound for the uniroot function used to find the inverse
#' cumulative hazard.
#'
#' @return   A list of proportions of dead and T2D patients by time \eqn{\tau}
#' in the A0 = a0 group in respectively intervened data and non intervened data.
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
                                 beta_L_A = 0, beta_L_Z = 3, beta_L_D = 0.5, beta_L_C = 0,
                                 beta_A_L = -0.5,  beta_A_Z = -0.5, beta_A_D = -1, beta_A_C = 0,
                                 beta_Z_L = -2, beta_Z_A = 0, beta_Z_D = -3, beta_Z_C = 0,
                                 beta_L0_L = 1, beta_L0_A = 1, beta_L0_Z = 0.1, beta_L0_D = 1, beta_L0_C = 0,
                                 beta_A0_L = -2.5, beta_A0_A = 0, beta_A0_Z = 0, beta_A0_D = -0.5, beta_A0_C = 0,
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
