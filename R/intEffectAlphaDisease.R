#' Estimate Effect of Intervention: Modifying Eta Parameter of Disease Process
#'
#' This function simulates data from the disease setting in two scenarios. Under
#' intervention on the shape parameter \eqn{\eta} of the disease process is multiplied by
#' \code{alpha}, and a baseline (non-intervened) scenario. It computes the proportion
#' of individuals who experience death or disease by a specified time \eqn{\tau}
#' in the group \code{A0 = a0}, optionally returning years_lost.
#' The function can also plot a sample of the event data for each scenario for comparison.
#'
#' @param N Integer. Number of individuals to simulate. Default is 10,000.
#' @param alpha Numeric scalar. Multiplicative factor applied to the disease process shape parameter \eqn{\eta}.
#' @param tau Numeric scalar. Time horizon at which proportions are computed.
#' @param years_lost Logical. If \code{TRUE}, computes years lost instead of proportions.
#' @param a0 Binary (0/1). Specifies the group for comparison.
#' @param plot Logical. If \code{TRUE}, plots timelines for sample of intervention and non intervention data.
#' @param eta Numeric vector of length 3. Shape parameters for Weibull hazards (default \code{rep(0.1, 3)}).
#' @param nu Numeric vector of length 3. Scale parameters for Weibull hazards (default \code{rep(1.1, 3)}).
#' @param beta_L0_D Numeric scalar. Effect of baseline covariate on death risk (default 0.5).
#' @param beta_L0_L Numeric scalar. Effect of baseline covariate on covariate risk (default 2).
#' @param beta_L_D Numeric scalar. Effect of covariate process on death risk (default 1).
#' @param beta_A0_D Numeric scalar. Effect of baseline treatment on death risk (default -0.1).
#' @param beta_A0_L Numeric scalar. Effect of baseline treatment on covariate risk (default -1).
#' @param lower Numeric scalar. Lower bound for root-finding in hazard inversion (default 1e-30).
#' @param upper Numeric scalar. Upper bound for root-finding in hazard inversion (default 200).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{effect_L}}{Proportion (or years lost) of individuals diagnosed with disease by time \eqn{\tau} in group \code{A0 = a0}, under intervention (\code{G1}) and without intervention (\code{G2}).}
#'   \item{\code{effect_death}}{Proportion (or years lost) of individuals who died by time \eqn{\tau} in group \code{A0 = a0}, under intervention (\code{G1}) and without intervention (\code{G2}).}
#' }
#' @export
#'
#' @examples
#' intEffectAlphaDisease(N = 1000, alpha = 0.7, tau = 5, years_lost = FALSE, a0 = 1)
intEffectAlphaDisease <- function(N = 1e4,
                      alpha = 0.5,
                      tau = 5,
                      years_lost = FALSE,
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
  data_G1 <- simDisease(N = N,
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
  data_G2 <- simDisease(N = N,
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
  #Proportion of subjects experiencing Disease before some time \tau in a0 group
  prop_G1_L <- mean(data_G1[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention
  prop_G2_L <- mean(data_G2[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # without intervention

  if (years_lost) {
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
