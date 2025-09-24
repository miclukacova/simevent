#' Simulate Continuous Time-to-Event Data with Multiple Event Types
#'
#' `simEventData` simulates event times and types for a cohort of individuals in a
#' counting process framework. It supports multiple event types (by default 4),
#' including terminal events, with intensities influenced by baseline covariates
#' and previous event history.
#'
#' The event intensities for event type \( x \) at time \( t \) are given by
#' \deqn{
#' \lambda^x(t) = \lambda_0^x(t) \exp(\beta_x^T L),
#' }
#' where the baseline intensity follows a Weibull hazard function:
#' \deqn{
#' \lambda_0^x(t) = \eta^x \nu^x t^{\nu^x - 1}.
#' }
#' Here, \( L \) is the vector of covariates and event counts, and \( \beta \) is the matrix
#' of regression coefficients representing the effect of covariates and previous events
#' on the intensity.
#'
#' @title Simulate Event Data with Multiple Event Types and Covariates
#'
#' @param N Integer. Number of individuals to simulate.
#' @param beta Numeric matrix. Regression coefficients matrix where columns correspond to event types (N0, N1, ...) and rows correspond to covariates (L0, A0, L1, L2, ...) and event counts (N0, N1, ...). Default is a zero matrix.
#' @param eta Numeric vector. Shape parameters of the Weibull baseline intensity for each event type. Default is 0.1 for all events.
#' @param nu Numeric vector. Scale parameters of the Weibull baseline intensity for each event type. Default is 1.1 for all events.
#' @param at_risk Function. Function determining if an individual is at risk for each event type, given their current event counts. Takes a numeric vector and returns a binary vector. Default returns 1 for all events.
#' @param term_deltas Integer vector. Event types considered terminal (after which no further events occur). Default is c(0, 1).
#' @param max_cens Numeric. Maximum censoring time. Events occurring after this time are censored. Default is Inf (no censoring).
#' @param add_cov Named list of functions. Functions generating additional baseline covariates. Each function takes integer N and returns a numeric vector of length N. Default is NULL.
#' @param override_beta Named list. Used to specify entries of the \code{beta} matrix to override defaults. For example, \code{list("L0" = c("N1" = 2))} sets the effect of L0 on N1 to 2.
#' @param max_events Integer. Maximum number of events to simulate per individual. Default is 10.
#' @param lower Numeric. Lower bound for root-finding in inverse cumulative hazard calculations. Default is \(10^{-15}\).
#' @param upper Numeric. Upper bound for root-finding in inverse cumulative hazard calculations. Default is 200.
#' @param gen_A0 Function. Function to generate the baseline treatment covariate A0. Takes N and L0 as inputs. Default is a Bernoulli(0.5) random variable.
#'
#' @return A \code{data.table} with columns:
#' \item{ID}{Individual identifier}
#' \item{Time}{Time of event}
#' \item{Delta}{Event type at time}
#' \item{L0}{Baseline covariate}
#' \item{A0}{Baseline treatment}
#' \item{L1, L2, ...}{Additional baseline covariates if specified}
#' \item{N0, N1, ...}{Event counts up to the current event}
#'
#' @examples
#' # Simulate data for 10 individuals with default settings
#' sim_data <- simEventData(N = 10)
#' head(sim_data)
#'
#' @export

simEventData <- function(N,                      # Number of individuals
                         beta = NULL,            # Effects
                         eta = NULL,             # Shape parameters
                         nu = NULL,              # Scale parameters
                         at_risk = NULL,         # Function defining the setting
                         term_deltas = c(0,1),   # Terminal events
                         max_cens = Inf,         # Followup time
                         add_cov = NULL,         # Additional baseline covariates
                         override_beta = NULL,   # Override beta
                         max_events = 10,        # Maximal events per individual
                         lower = 10^(-15),       # Lower bound for ICH
                         upper = 200,            # Upper bound for ICH
                         gen_A0 = NULL           # Generation of A0
){
  ID <- NULL

  ############################ Check and useful quantities #####################
  # Check of add_cov
  if(!(is.null(add_cov) | is.list(add_cov))){
    stop("add_cov needs to be list of random functions")
  }

  # Number of additional baseline covariates
  num_add_cov <- length(add_cov)

  # Determine number of events
  num_events <- if (!is.null(eta)) length(eta) else
    if (!is.null(nu)) length(nu) else
      if (!is.null(beta)) ncol(beta) else 4

  # Useful indices
  N_start <- 3 + num_add_cov
  N_stop <- 2 + num_add_cov + num_events

  ############################ Default values ##################################

  # Set default values for beta, eta, and nu
  beta <- if (!is.null(beta)) beta else matrix(0, nrow = N_stop, ncol = num_events)
  colnames(beta) <- paste0("N", seq(0, num_events -1))

  if((N_stop) != nrow(beta)){
    stop("Number of rows in beta should equal the sum of number of event and
         number of additional covariates + 2")
  }

  eta  <- if (!is.null(eta)) eta   else rep(0.1, num_events)
  nu   <- if (!is.null(nu)) nu     else rep(1.1, num_events)

  # Check of dimensions
  if(num_events != length(nu) || num_events != ncol(beta)){
    stop("Length of eta should be equal to nu and number of columns of beta")
  }

  # Default at_risk
  if(is.null(at_risk)){
    riskss <- rep(1, num_events)
    at_risk <- function(events) return(riskss)
  }
  # Default A0 generation
  if(is.null(gen_A0)){
    gen_A0 <- function(N, L0) stats::rbinom(N, 1, 0.5)
  }

  # Matrix for storing values
  simmatrix <- matrix(0, nrow = N, ncol = (2 + num_events + num_add_cov))

  # Generate additional covariates if distributions are specified
  if (num_add_cov != 0) {
    simmatrix[,3:(2+length(add_cov))] <- sapply(add_cov, function(f) f(N))
  }

  # Naming of matrices
  if (is.null(names(add_cov)) && num_add_cov != 0) {
      colnames(simmatrix) <- c("L0", "A0", paste0("L", seq_len(num_add_cov)), colnames(beta))
  } else {
      colnames(simmatrix) <- c("L0", "A0", names(add_cov), colnames(beta))
  }

  rownames(beta) <- colnames(simmatrix)

  # Filling out beta matrix
  if(!is.null(override_beta)){
      for (bb in 1:length(override_beta)) {
          if (names(override_beta)[bb] %in% rownames(beta)) {
              beta[names(override_beta)[bb], names(override_beta[[bb]])] <- override_beta[[bb]]
          } else {
              beta <- rbind(beta, matrix(0, nrow = 1, ncol = ncol(beta)))
              beta[nrow(beta), names(override_beta[[bb]])] <- override_beta[[bb]]
              rownames(beta)[nrow(beta)] <- names(override_beta)[bb]
          }
      }
  }

  ############################ Functions #######################################

  # Proportional hazard
  calculate_phi <- function(simmatrix) {
    if(nrow(beta) == N_stop) return(exp(simmatrix %*% beta)) else {
      obj <- as.data.frame(simmatrix)
      X <- sapply(rownames(beta), function(expr) {
        eval(parse(text = expr), envir = obj)
      })
      effects <- as.matrix(X) %*% beta
      return(exp(effects))
    }
  }

  # Intensities
  lambda <- function(t, i) {
    risk_vec <- at_risk(simmatrix[i, N_start:N_stop])
    risk_vec * eta * nu * t^(nu - 1) * phi[i,]
  }

  # If all events have the same parameter, the inverse of the cumulative hazard simplifies
  if(all(nu[1] == nu) && all(eta[1] == eta)){
    inverse_sc_haz <- function(p, t, i) {
      denom <- sum(at_risk(simmatrix[i, N_start:N_stop]) * eta * phi[i,])
      (p / denom + t^nu[1])^(1 / nu[1]) - t
    }
  # Otherwise we use a numerical inverse coded in rcpp
  } else{
    inverse_sc_haz <- function(p, t, i) {
      inverseScHaz(p, t, lower = lower, upper = upper, eta = eta, nu = nu,
                     phi = phi[i,], at_risk = at_risk(simmatrix[i, N_start:N_stop]))

    }
  }

  # Event probabilities
  probs <- function(t, i){
    if(t == max_cens) return(c(1, rep(0, (num_events - 1))))
    probs <- lambda(t, i)
    summ <- sum(probs)
    probs / summ
  }

  ############################ Initializing Simulations ########################

  # Draw baseline covariates
  simmatrix[,1] <- stats::runif(N)                   # L0
  simmatrix[,2] <- gen_A0(N, simmatrix[,1])          # A0

  # Initialize
  T_k <- rep(0,N)                                    # Time 0
  alive <- 1:N                                       # Keeping track of who is alive
  res_list <- vector("list", max_events)             # For results
  idx <- 1                                           # Index


  ############################ Simulations #####################################

  while(length(alive) != 0){
    # Simulate time
    V <- -log(stats::runif(N))
    phi <- calculate_phi(simmatrix)
    W <- sapply(alive, function(i) inverse_sc_haz(V[i], T_k[i], i))
    T_k[alive] <- T_k[alive] + W

    # Maximal censoring time
    T_k[T_k > max_cens] <- max_cens

    # Simulate event
    probs_mat <- sapply(alive, function(i) probs(T_k[i], i), simplify = "array")
    Deltas <- sampleEvents(probs_mat)

    # Update event counts
    simmatrix[cbind(alive, 2 + num_add_cov + Deltas + 1)] <-
      simmatrix[cbind(alive, 2 + num_add_cov + Deltas + 1)] + 1

    # Store data
    kth_event <- data.table(ID = alive,
                            Time = T_k[alive],
                            Delta = Deltas)

    res_list[[idx]] <- cbind(kth_event, data.table::as.data.table(simmatrix[alive, , drop = FALSE]))
    idx <- idx + 1

    # Who is still alive and uncensored?
    alive <- alive[!Deltas %in% term_deltas]
  }

  res <- data.table::rbindlist(res_list)
  setkey(res, ID)

  return(res)
}




