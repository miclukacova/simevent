#' Calculate the effect of an intervention modifying the treatment process intensity.
#'
#' This function simulates event history data from the Treatment scenario (see function simTreatment)
#' under intervention and without intervention. We consider the intervention where
#' the \eqn{\eta} parameter of the treatment process is multiplied by a factor \code{alpha}.
#' It computes the proportion or years lost of death and treatment by time \eqn{\tau}
#' comparing intervened and non-intervened data.
#'
#' @title Perform intervention on treatment process intensity
#'
#' @param N Integer. Number of individuals to simulate.
#' @param alpha Numeric. Multiplicative factor applied to the \eqn{\eta} parameter of the treatment process.
#' @param tau Numeric. Time point at which event proportions or years lost are compared.
#' @param years_lost Logical. If TRUE, compute years lost estimand instead of simple proportions.
#' @param plot Logical. If TRUE, output plots of the first 250 events for respectively intervened and no invtervened data.
#' @param eta Numeric vector of length 4. Shape parameters for Weibull hazards.
#' @param nu Numeric vector of length 4. Scale parameters for Weibull hazards.
#' @param beta_L_A Numeric. Effect of covariate L on treatment process A.
#' @param beta_L_D Numeric. Effect of covariate L on death process D.
#' @param beta_A_D Numeric. Effect of treatment process A on death process D.
#' @param beta_A_L Numeric. Effect of treatment process A on covariate L.
#' @param beta_L0_A Numeric. Effect of L0 on treatment process A.
#' @param lower Numeric. Lower bound for root-finding algorithm to invert cumulative hazard.
#' @param upper Numeric. Upper bound for root-finding algorithm.
#'
#' @return A list containing:
#' \describe{
#'   \item{effect_A}{Proportion or years lost of treated patients with and without intervention.}
#'   \item{effect_death}{Proportion or years lost of deaths with and without intervention.}
#' }
#'
#' @export
#'
#' @examples
#' intEffectAlphaTreat()

intEffectAlphaTreat <- function(N = 1e4,
                               alpha = 0.5,
                               tau = 5,
                               years_lost = FALSE,
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
  data_G1 <- simTreatment(N = N,
                            cens = 0,
                            eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                            nu = nu,
                            beta_L_A = beta_L_A, beta_L_D = beta_L_D,
                            beta_A_D = beta_A_D, beta_A_L = beta_A_L,
                            beta_L0_A = beta_L0_A,
                            lower = lower, upper = upper)

  data_G2 <- simTreatment(N = N,
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
  #Proportion of subjects experiencing treatment before some time \tau in a0 group
  prop_G1_A <- mean(data_G1[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention
  prop_G2_A <- mean(data_G2[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # without intervention

  if (years_lost) {
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

