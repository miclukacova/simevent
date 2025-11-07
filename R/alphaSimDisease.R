#' Simulation and Estimation in Setting Disease with Modified Shape Parameter of Disease Process
#'
#' This function simulates event history data from the Disease setting (see simDisease)
#' where the shape parameter \eqn{\eta} of the disease process is multiplied by
#' \code{alpha}. The function either
#' * returns number of years lost before \eqn{\tau} of death and disease
#' * returns simulated data.
#'
#' @param N Integer. Number of individuals to simulate. Default is 10,000.
#' @param alpha Numeric scalar. Multiplicative factor applied to the disease process shape parameter \eqn{\eta}.
#' @param tau Numeric scalar. Time horizon at which proportions are computed.
#' @param years_lost Logical. If \code{TRUE}, computes years lost instead of proportions.
#' @param a0 Binary (0/1). Specifies the group for comparison.
#' @param eta Numeric vector of length 3. Shape parameters for Weibull hazards (default \code{rep(0.1, 3)}).
#' @param nu Numeric vector of length 3. Scale parameters for Weibull hazards (default \code{rep(1.1, 3)}).
#' @param beta_L0_D Numeric scalar. Effect of baseline covariate on death risk (default 0.5).
#' @param beta_L0_L Numeric scalar. Effect of baseline covariate on covariate risk (default 2).
#' @param beta_L_D Numeric scalar. Effect of covariate process on death risk (default 1).
#' @param beta_A0_D Numeric scalar. Effect of baseline treatment on death risk (default -0.1).
#' @param beta_A0_L Numeric scalar. Effect of baseline treatment on covariate risk (default -1).
#' @param lower Numeric scalar. Lower bound for root-finding in hazard inversion (default 1e-30).
#' @param upper Numeric scalar. Upper bound for root-finding in hazard inversion (default 200).
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{1}).
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0. Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#' @param return_data Logical. If \code{TRUE} the simulated data is returned.
#' @param beta_L_D_t_prime Numeric scalar or NULL. Additional effect of covariate change on death risk after time \code{t_prime} (optional).
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{effect_L}}{Proportion (or years lost) of individuals diagnosed with disease by time \eqn{\tau} in group \code{A0 = a0}, under intervention (\code{G1}) and without intervention (\code{G2}).}
#'   \item{\code{effect_death}}{Proportion (or years lost) of individuals who died by time \eqn{\tau} in group \code{A0 = a0}, under intervention (\code{G1}) and without intervention (\code{G2}).}
#' }
#' Or the simulated data.
#' @export
#'
#' @examples
#' alphaSimDisease(N = 1000, alpha = 0.7, tau = 5, years_lost = FALSE, a0 = 1)
alphaSimDisease <- function(N = 1e4,
                            alpha = 0.5,
                            tau = 5,
                            years_lost = FALSE,
                            a0 = 1,
                            eta = rep(0.1, 3),
                            nu = rep(1.1, 3),
                            beta_L0_D = 0.5,
                            beta_L0_L = 2,
                            beta_L_D = 1,
                            beta_A0_D = -0.1,
                            beta_A0_L = -1,
                            lower = 10^(-30),
                            upper = 200,
                            cens = 1,
                            gen_A0 = NULL,
                            return_data = FALSE,
                            beta_L_D_t_prime = 0,
                            t_prime = NULL){

  Delta <- Time <- A0 <- tmp <- V1 <- NULL

  # Generate large data set under the intervened intensity
  data <- simDisease(N = N,
                     cens = 0,
                     eta = c(eta[1:2],eta[3]*alpha),
                     nu = nu,
                     beta_A0_D = beta_A0_D,
                     beta_L_D = beta_L_D,
                     beta_L0_D = beta_L0_D,
                     beta_L0_L = beta_L0_L,
                     beta_A0_L = beta_A0_L,
                     upper = upper,
                     lower = lower,
                     gen_A0 = gen_A0,
                     beta_L_D_t_prime = beta_L_D_t_prime,
                     t_prime = t_prime)

  if(return_data) return(data)

  #Proportion of subjects dying before some time $\tau$ in a0 group
  prop_D <- data[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)]
  #Proportion of subjects experiencing Disease before some time \tau in a0 group
  prop_L <- mean(data[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])

  if (years_lost) {
    data[, tmp := cumsum((Delta == 1)*(tau - pmin(tau, Time))), by = "ID"]
    prop_D <- data[, tmp[.N], by = "ID"][, mean(V1)]
    data[, tmp := cumsum((Delta == 2)*(tau - pmin(tau, Time))), by = "ID"]
    prop_L <- data[, tmp[.N], by = "ID"][, mean(V1)]
  }

  return(list(effect_D = prop_D,
              effect_L = prop_L))
}
