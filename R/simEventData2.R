#' `simEventData2` is a function to simulate event data, e.g. observational healthcare data. The number of events
#' simulated corresponds to the length of the \eqn{\eta} and \eqn{\nu} vector, and the number of columns in the
#' \eqn{\beta} matrix. By default 4 different types of events are simulated, which can be chosen to represent Censoring (0),
#' Death (1), Operation(2) and Change in Covariate Process(3). Death and Censoring are terminal events and Operation
#' and Change in Covariate Process can occur once. The intensities of the various events depend upon previous events,
#' baseline covariates and the pre specified \eqn{\beta}, \eqn{\nu} and \eqn{\eta} parameters.
#'
#' Different settings can be specified: Survival setting, competing risk setting and operation setting.
#' The different settings can be specified by the at_risk function. Default is the operation setting.
#'
#' @title Simulate Event Data
#'
#' @param N A double for the number of simulated individuals
#' @param beta A matrix of doubles for the effects on the intensities. The columns
#' represent the events. In the default case Censoring, Death, Operation, and Covariate Change.
#' The rows represent the baseline covariate \eqn{L0} (which can be interpreted
#' as age of participant), baseline treatment \eqn{A0}, the indicator for change
#' in the covariate process  \eqn{L}, and the indicator for operation \eqn{A}. If additional
#' covariates are specified, each of these corresponds to a row in the beta matrix. The
#' \eqn{\beta} matrix is by default set to 0.
#' @param eta Vector of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}} Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull intensity. Default is set
#' to 1.1 for all events.
#' @param at_risk At risk function. Default is set to the operation setting. A
#' survival or competing risk setting can be specified as well. The \code{at_risk}
#' function is an indicator of whether an individual is at risk for a specific event.
#' The function takes as input i (the index belonging to a particular individual),
#' L (indicator for a change in covariate) and A (indicator for change in treatment process).
#' The function returns a vector of 0's  and 1's corresponding to whether individual
#' i is at risk for a particular event.
#' @param term_deltas Terminal events. Default is set so that event 1 and 2 are terminal events.
#' @param max_cens A maximum censoring time. By default set to infinity.
#' @param add_cov List of random generator functions for the distributions of
#' additional baseline covariates. The functions should take the number of observations
#' as input. By default set to NULL.
#'
#' @return data.table containing the simulated data. There is a column for ID, time of event (Time),
#' event type (Delta), baseline covariate (L0), indicator for change in covariate process (L), Baseline Treatment (A0),
#' indicator for Operation (A), and potentially an additional baseline covariate (L1) representing sex.
#' @export
#'
#' @examples
#' simEventData2(N = 10)

simEventData2 <- function(N,                     # Number of individuals
                         beta = NULL,            # Effects
                         eta = rep(0.1,4),       # Shape parameters
                         nu = rep(1.1,4),        # Scale parameters
                         at_risk = NULL,         # Function defining the setting
                         term_deltas = c(0,1),   # Terminal events
                         max_cens = Inf,         # Followup time
                         add_cov = NULL          # Additional baseline covariates
){
  ID <- NULL

  # Checks
  if(!(is.null(add_cov) | is.list(add_cov))){
    stop("add_cov needs to be list of random functions")
  }

  # Number of additional baseline covariates
  num_add_cov <- length(add_cov)

  if(is.null(beta)){
    beta <- matrix(0, nrow = 4 + num_add_cov, ncol = length(eta))
  }

  # Events
  x <- 1:ncol(beta)

  if(is.null(at_risk)){
    at_risk <- function(i, L, A) {
      return(c(
        1,1, # If you have not died or been censored yet, you are at risk for dying or being censored
        as.numeric(A[i] == 0), # You are only at risk for an operation if you have not had an operation yet
        as.numeric(L[i] == 0))) # You are only at risk for a change in the covariate process if you have not had a change yet
      }}

  # Intensities
  phi <- function(i) {
    exp(L0[i] / 50 * beta[1,] +
        A0[i] * beta[2,] +
        L[i] * beta[3,] +
        A[i] * beta[4,] +
        if(num_add_cov > 0) L1[i,] %*% beta[5:(4 + num_add_cov),] else 0)
  }

  lambda <- function(t, i) {
    at_risk(i, L, A) * eta * nu * t ^ (nu - 1) * phi(i)
  }

  # If all events have same parameter, the inverse simplifies
  if(all(nu[1] == nu)){
    inverse_sc_haz <- function(p, t, i) {
      denom <- sum(at_risk(i, L, A) * eta * phi(i))
      (p / denom + t^nu[1])^(1 / nu[1]) - t
    }
  } else{
    # Summed cumulative hazard
    sum_cum_haz <- function(u, t, i) {
      sum(at_risk(i, L, A) * eta * phi(i) * ((t + u) ^ nu - t ^ nu))
    }

    # Inverse summed cumulative hazard function
    inverse_sc_haz <- function(p, t, i) {
      root_function <- function(u) sum_cum_haz(u, t, i) - p
      stats::uniroot(root_function, lower = 10^-15, upper = 200)$root
    }
  }

  # Event probabilities
  probs <- function(t, i){
    if(t > max_cens){
      probs <- c(1,0,0,0)
    }
    else{
      probs <- lambda(t, i)
      summ <- sum(probs)
      probs / summ
    }
  }

  # Draw baseline covariates
  L0 <- stats::runif(N, 30, 70)
  A0 <- stats::rbinom(N, 1, 0.5)

  # Generate additional covariates if distributions are specified
  if (num_add_cov != 0) {
    L1 <- sapply(add_cov, function(f) f(N))
    colnames(L1) <- paste0("L", seq_len(ncol(L1)))
  } else {
    L1 <- NULL
  }

  # Initialize
  T_k <- rep(0,N)
  Delta <- -1
  L <- rep(0, N)
  A <- rep(0, N)
  alive <- 1:N

  res <- data.table()

  while(length(alive) != 0){
    # Simulate time
    V <- stats::runif(N)
    W <- sapply(alive, function(i) inverse_sc_haz(-log(V)[i], T_k[i], i))
    T_k[alive] <- T_k[alive] + W

    # Simulate event
    Deltas <- sapply(alive, function(i) sample(x, size = 1, prob = probs(T_k[i], i)) - 1)

    kth_event <- data.table(ID = alive,
                            Time = T_k[alive],
                            Delta = Deltas,
                            L0 = L0[alive],
                            L = L[alive],
                            A0 = A0[alive],
                            A = A[alive])
    kth_event <- cbind(kth_event, L1[alive,])

    res <- rbind(res, kth_event)

    # Update treatment process and covariate change indicators
    A[alive][Deltas == 2] <- 1
    L[alive][Deltas == 3] <- 1

    # Who is still alive and uncensored?
    alive <- alive[! Deltas %in% term_deltas]
  }

  setkey(res, ID)

  return(res)
}
