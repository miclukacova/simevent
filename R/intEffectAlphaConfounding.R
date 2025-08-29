#' Function to calculate the effect of an intervention where the eta parameter of
#' the operation process is multiplied by the factor alpha. We consider the outcome
#' variable death and operation at time \eqn{\tau}.
#'
#' @title Perform intervention
#'
#' @param N A double of the number of individuals to simulate
#' @param alpha Double corresponding to the factor to be multiplied onto eta
#' parameter of the operation process
#' @param tau A double deciding the time we compare proportions.
#' @param rmst ???
#' @param plot Logical indicating whether a plot should be outputed
#' @param eta Vector of  length 3 of shape parameters for the Weibull hazard
#' with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}.
#' @param nu Vector of length 3 of scale parameters for the Weibull hazard.
#' @param beta_L_A Specifies how L affects A.
#' @param beta_L_D Specifies how L affects D.
#' @param beta_A_D Specifies how A affects D.
#' @param beta_A_L Specifies how A affects L.
#' @param beta_L0_A Specifies how L0 affects A.
#' @param lower Lower bound for the uniroot function used to find the inverse
#' cumulative hazard.
#' @param upper Upper bound for the uniroot function used to find the inverse
#' cumulative hazard.
#'
#' @return  A list of proportions of dead and operated patients by time \eqn{\tau}
#' in the intervened data and non intervened data.
#' @export
#'
#' @examples
#' intEffectAlphaConf()

intEffectAlphaConf <- function(N = 1e4,
                               alpha = 0.5,
                               tau = 5,
                               rmst = FALSE,
                               plot = FALSE,
                               nu = rep(1.1, 4),
                               eta = rep(0.1, 4),
                               beta_L_A = 1,
                               beta_L_D = 1,
                               beta_A_D = -0.5,
                               beta_A_L = -1,
                               beta_L0_A = 1,
                               lower = 10^(-300),
                               upper = 300){

  Delta <- Time <- A0 <- NULL
  # Generate large data
  data_G1 <- simConfounding(N = N,
                            cens = 0,
                            eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                            nu = nu,
                            beta_L_A = beta_L_A, beta_L_D = beta_L_D,
                            beta_A_D = beta_A_D, beta_A_L = beta_A_L,
                            beta_L0_A = beta_L0_A,
                            lower = lower, upper = upper)

  data_G2 <- simConfounding(N = N,
                            cens = 0,
                            eta = eta,
                            nu = nu,
                            beta_L_A = beta_L_A, beta_L_D = beta_L_D,
                            beta_A_D = beta_A_D, beta_A_L = beta_A_L,
                            beta_L0_A = beta_L0_A,
                            lower = lower, upper = upper)

  if (plot) gridExtra::grid.arrange(plotEventData(data_G1[1:250]),
                       plotEventData(data_G2[1:250]), nrow = 1)

  #Proportion of subjects dying before some time $\tau$ in treatment group
  prop_G1 <- data_G1[Delta == 1, mean(Delta == 1 & Time < tau)] # with intervention
  prop_G2 <- data_G2[Delta == 1, mean(Delta == 1 & Time < tau)] # without intervention
  #Proportion of subjects experiencing T2D before some time \tau in a0 group
  prop_G1_A <- mean(data_G1[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention
  prop_G2_A <- mean(data_G2[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # without intervention

  if (rmst) {
    tgrid <- seq(0, tau, length = 100)
    prop_G1 <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      data_G1[Delta == 1, mean(Delta == 1 & Time <= t)])) # with intervention
    prop_G2 <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      data_G2[Delta == 1, mean(Delta == 1 & Time <= t)])) # without intervention

    prop_G1_A <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      mean(data_G1[, any(Delta == 2 & Time <= t)[1], by = "ID"][[2]]))) # with intervention
    prop_G2_A <- sum(diff(tgrid)[1]*sapply(tgrid[-1], function(t)
      mean(data_G2[, any(Delta == 2 & Time <= t)[1], by = "ID"][[2]]))) # without intervention
    }

  return(list(effect_A = c(G1 = prop_G1_A, G2 = prop_G2_A),
              effect_death = c(G1 = prop_G1, G2 = prop_G2)))
}

