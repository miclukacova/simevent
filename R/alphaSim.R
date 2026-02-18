#' Simulation and Estimation with Modified Shape Parameter
#'
#' This function simulates event history data from either the Disease, Treatment or Drop In setting.
#' See simDisease, simTreatment and simDropIn. The shape parameter \eqn{\eta} of respectively the disease process,
#' the Drop In process and the Treatment process is multiplied by \code{alpha}. The function either
#' * returns the proportion of individuals who experience death and the proportion of individuals who experience disease/drop in/treatment
#' by a specified time \eqn{\tau} (in group \code{A0 = a0} for drop in and disease).
#' * returns number of years lost before \eqn{\tau} of death and disease/drop in/treatment
#' * returns simulated data.
#' One can specify all the same parameters as in the functions \code{simDisease}, \code{simTreatment} and \code{simDropIn}.
#'
#' @param N Integer. Number of individuals to simulate. Default is 10,000.
#' @param alpha Numeric scalar. Multiplicative factor applied to the disease process shape parameter \eqn{\eta}.
#' @param tau Numeric scalar. Time horizon at which proportions are computed.
#' @param years_lost Logical. If \code{TRUE}, computes years lost instead of proportions.
#' @param a0 Binary (0/1). Specifies the group for comparison in setting Drop In and Disease.
#' @param eta Numeric vector. Shape parameters for Weibull hazards. Length of the vector should
#' match number of events. For the Disease and Drop In setting this is 4. For the Treatment setting,
#' this is 3. (default \code{rep(0.1, 4)}).
#' @param nu Numeric vector. Scale parameters for Weibull hazards. Length of the vector should
#' match number of events. For the Disease and Drop In setting this is 4.  For the Treatment setting,
#' this is 3. (default \code{rep(1.1, 4)}).
#' @param cens Binary scalar. Indicates whether individuals are at risk of censoring (default \code{0}).
#' @param setting Character string. Must be either "Disease", "Drop In" or "Treatment". Depending on the simulation setting.
#' @param return_data Logical. If \code{TRUE} the simulated data is returned.
#' @param ... Additional arguments passed to respectively simDisease, simTreatment and simDropIn.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{effect_L}}{Proportion (or years lost) of individuals diagnosed with disease by time \eqn{\tau}, under intervention.}
#'   \item{\code{effect_death}}{Proportion (or years lost) of individuals who died by time \eqn{\tau} under intervention.}
#' }
#' Or the simulated data.
#' @export
#'
#' @examples
#' alphaSim(N = 100, eta = rep(0.1,3), nu = rep(1.1,3), alpha = 0.5, setting = "Disease")
#' alphaSim(N = 100, setting = "Drop In", beta_A0_Z = 1)
alphaSim <- function(N = 1e4,
                     eta = rep(0.1,4),
                     nu = rep(1.1,4),
                     alpha = 0.5,
                     tau = 5,
                     a0 = 1,
                     years_lost = FALSE,
                     setting = "Disease",
                     return_data = FALSE,
                     cens = 0,
                     ...){

  Delta <- Time <- A0 <- tmp <- V1 <- NULL

  # Generate large data set under the intervened intensity
  if (setting == "Disease"){
    if(length(eta) != 3 | length(nu) != 3) stop ("eta and nu must be of length 3 in the Disease setting")
    data <- simDisease(N = N,
                       eta = c(eta[1:2],eta[3]*alpha),
                       cens = cens,
                       ...)
  } else if (setting == "Drop In"){
    if(length(eta) != 4 | length(nu) != 4) stop ("eta and nu must be of length 3 in the Drop In setting")
    data <- simDropIn(N = N,
                      eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                      nu = nu,
                      cens = cens,
                      ...)
  } else if (setting == "Treatment"){
    if(length(eta) != 4 | length(nu) != 4) stop ("eta and nu must be of length 4 in the Treatment setting")
    data <- simTreatment(N = N,
                         eta = c(eta[1:2],eta[3]*alpha, eta[4]),
                         cens = cens,
                         nu = nu,
                         ...)
  } else {
    stop ("Setting must be either Disease, Drop In or Treatment")
  }

  if(return_data) return(data)

  # Proportion of subjects dying before some time $\tau$
  if(setting == "Treatment") prop_D <- data[Delta == 1, mean(Delta == 1 & Time < tau)]
  else prop_D <- data[A0 == a0 & Delta == 1, mean(Delta == 1 & Time < tau)]

  # Proportion of subjects experiencing Treatment/Drop In/Disease
  if(setting == "Treatment") prop2 <- mean(data[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])
  if(setting == "Drop In") prop2 <- mean(data[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])
  if(setting == "Disease") prop2 <- mean(data[A0 == a0, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])

  if (years_lost) {
    data[, tmp := cumsum((Delta == 1)*(tau - pmin(tau, Time))), by = "ID"]
    prop_D <- data[, tmp[.N], by = "ID"][, mean(V1)]
    data[, tmp := cumsum((Delta == 2)*(tau - pmin(tau, Time))), by = "ID"]
    prop2 <- data[, tmp[.N], by = "ID"][, mean(V1)]
  }

  return(list(effectDeath = prop_D,
              effectSetting = prop2))
}
