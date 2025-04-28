#' Function to simulate data from a Survival setting. The function simulates data corresponding to $N$ individuals that
#' are at risk for censoring (0) and an event (1).
#'
#' @title Simulate Survival Data
#'
#' @param N A double of the number of individuals
#' @param beta A 2X2 matrix with the effect of L0 and A0 on Censoring and Event.
#' The columns represent Censoring and Death, while the rows represent L0 and A0.
#' @param eta Vector of  length 2 of shape parameters for the Weibull hazard with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of length 2 of scale parameters for the Weibull hazard.
#' @param cens Binary variable indicating whether there should be a censoring rpocess
#' @param ... Lets you specify additional covariates. See the arfument `add_cov`
#' in the function `simEventData`.
#'
#' @return  Data frame containing the simulated survival data
#' @export
#'
#' @examples
#' simSurvData(10)
simSurvData <- function(N,
                        beta = NULL,
                        eta = rep(0.1,2),
                        nu = rep(1.1,2),
                        cens = 1,
                        ...
                        ){


    at_risk <- function(events) c(cens,1)

    if(is.null(beta)){
        beta <- matrix(0, ncol = 2, nrow = 2)
    }

    beta <- rbind(beta, matrix(0, nrow = 2, ncol = 2))

    dots <- list(...)
    has_add_cov <- "add_cov" %in% names(dots)
    if (has_add_cov) beta <- rbind(beta, matrix(c(rep(0, length(dots$add_cov)), rep(0.1, length(dots$add_cov))), nrow = length(dots$add_cov), ncol = ncol(beta)))

    results <- simEventData(N, beta, eta = eta, nu = nu, at_risk = at_risk, ...)

    # We don't need columns for terminal events
    results <- results[, !c("N0", "N1")]

    return(results)
}

