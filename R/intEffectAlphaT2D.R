#' Function to calculate the effect of an intervention where the eta parameter of
#' the T2D process is multiplied by the factor alpha. We consider the outcome
#' variable death and T2D at time \eqn{\tau} - that is the output of the function is
#' the proportion of dead and T2D diagnosed patients by time \eqn{\tau} in A0 = a0
#' group in respectively intervened data and non intervened data.
#'
#' @title Perform intervention
#'
#' @param N A double of the number of individuals to simulate
#' @param alpha Double corresponding to the factor to be multiplied onto eta
#' parameter of the T2D process
#' @param tau A double deciding the time we compare proportions.
#' @param rmst Logical
#' @param a0 A binary variable deciding whether we are comparing treatment group
#' or placebo group.
#' @param plot Logical indicating whether a plot should be outputed
#' @param eta Vector of  length 3 of shape parameters for the Weibull hazard
#' with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}.
#' @param nu Vector of length 3 of scale parameters for the Weibull hazard.
#' @param beta_L_D Specifies how change in the covariate process affects risk of death. Is by default set to 1.
#' @param beta_L0_L Specifies how change in the covariate process affects risk of death. Is by default set to 1.
#' @param beta_A0_D Specifies how baseline treatment affects risk of death. Is by default set to 0.
#' @param beta_A0_L Specifies how baseline treatment affects risk of T2D. Is by default set to 0.
#' @param beta_L0_D Specifies how baseline covariate affects risk of Death. Is by default set to 1.
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
#' intEffectAlphaT2D()

intEffectAlphaT2D <- function(N = 1e4,
                      alpha = 0.5,
                      tau = 5,
                      rmst = FALSE,
                      a0 = 1,
                      plot = FALSE,
                      eta = rep(0.1, 3),
                      nu = rep(1.1, 3),
                      beta_L0_D = 0.5,
                      beta_L0_L = 2,
                      beta_L_D = 1,
                      beta_A0_D = -0.1,
                      beta_A0_L = -1,
                      lower = 10^(-30),
                      upper = 200){

  Delta <- Time <- A0 <- NULL

  # Generate large data set under the intervened intensity
  data_G1 <- simT2D(N = N,
                    cens = 0,
                    eta = c(eta[1:2],eta[3]*alpha),
                    nu = nu,
                    beta_A0_D = beta_A0_D,
                    beta_L_D = beta_L_D,
                    beta_L0_D = beta_L0_D,
                    beta_L0_L = beta_L0_L,
                    beta_A0_L = beta_A0_L,
                    upper = upper,
                    lower = lower)

  # Generate large data set without intervened intensity
  data_G2 <- simT2D(N = N,
                    cens = 0,
                    eta = eta,
                    nu = nu,
                    beta_A0_D = beta_A0_D,
                    beta_L0_L = beta_L0_L,
                    beta_A0_L = beta_A0_L,
                    beta_L_D = beta_L_D,
                    beta_L0_D = beta_L0_D,
                    upper = upper,
                    lower = lower)

  if (plot) gridExtra::grid.arrange(plotEventData(data_G1[1:250]),
                         plotEventData(data_G2[1:250]), nrow = 1)

  #Proportion of subjects dying before some time $\tau$ in a0 group
  prop_G1 <- data_G1[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)] # with intervention
  prop_G2 <- data_G2[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)] # without intervention
  #Proportion of subjects experiencing T2D before some time \tau in a0 group
  prop_G1_L <- mean(data_G1[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention
  prop_G2_L <- mean(data_G2[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # without intervention

  if (rmst) {
    tgrid <- seq(0, tau, length = 100)
    prop_G1 <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      data_G1[A0 == a0 & Delta == 1, mean(Delta == 1 & Time <= t)])) # with intervention
    prop_G2 <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      data_G2[A0 == a0 & Delta == 1, mean(Delta == 1 & Time <= t)])) # without intervention

    prop_G1_L <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      mean(data_G1[A0 == a0, any(Delta == 2 & Time <= t)[1], by = "ID"][[2]]))) # with intervention
    prop_G2_L <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      mean(data_G2[A0 == a0, any(Delta == 2 & Time <= t)[1], by = "ID"][[2]]))) # without intervention

  }

  return(list(effect_L = c(G1 = prop_G1_L, G2 = prop_G2_L),
              effect_death = c(G1 = prop_G1, G2 = prop_G2)))
}
