#' Simulate Event Data with Time-Varying Effects
#'
#' `simEventTV` is a function that simulates event data, with the option of
#' adding time varying effects. The function is build up in the same way as `simEventData`,
#' with the additional arguments `tv_eff` and `t_prime`, which specify the change
#' of the beta matrix at time `t_prime`.
#'
#' @title simEventTV
#'
#' @param N Integer. Number of individuals to simulate.
#' @param beta Matrix. Coefficients for covariates and processes. Columns correspond to events (`N0`, `N1`, ...),
#'             rows correspond to covariates (`L0`, `A0`, ..., and past event counts).
#' @param tv_eff Matrix. Time-varying changes to `beta`, applied at time `t_prime`.
#'               Must have same dimensions as `beta`.
#' @param t_prime Numeric. Time at which `tv_eff` is applied to `beta`.
#' @param eta Numeric vector. Shape parameters of Weibull intensities for each event.
#' @param nu Numeric vector. Scale parameters of Weibull intensities for each event.
#' @param at_risk Function. Determines which events an individual is at risk for, based on event history.
#' @param term_deltas Integer vector. Event types considered terminal (e.g., death).
#' @param max_cens Numeric. Maximum censoring time. Defaults to `Inf`.
#' @param add_cov Named list of functions for generating additional baseline covariates.
#'                Each function takes one argument `N` and returns a vector of length `N`.
#' @param override_beta Named list to override elements of `beta`. Format:
#'                      `list("covariate" = c("event" = value))`.
#' @param max_events Integer. Maximum number of events allowed per individual.
#' @param lower Numeric. Lower bound for the root-finding algorithm used in inverse cumulative hazard computation.
#' @param upper Numeric. Upper bound for the root-finding algorithm used in inverse cumulative hazard computation.
#' @param gen_A0 Function. Generates baseline treatment assignment. Takes arguments `N` and `L0`.
#' @param tv_eff A matrix of the same dimensions as beta, specifying the
#' change of the effects at time t_prime.the matrix is in the same format as beta.
#' @param t_prime The time where the effects change.
#'
#' @return A `data.table` with columns:
#'   \item{ID:}{Individual identifier}
#'   \item{Time:}{Time of event}
#'   \item{Delta:}{Type of event}
#'   \item{L0:}{Baseline covariate}
#'   \item{A0:}{Baseline treatment}
#'   \item{N0, N1, ...:}{Cumulative event counts}
#'   \item{L1, L2, ...:}{Additional covariates (if specified)}
#'
#' @examples
#' eta <- rep(0.1, 2)
#' simEventTV(N = 100, t_prime = 1, eta = eta, term_deltas = c(0, 1))
#'
#' @export
#'
simEventTV <- function(N,                      # Number of individuals
                       beta = NULL,            # Effects
                       tv_eff = NULL,          # Time varying effects
                       t_prime = Inf,          # Time of change in effects
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

  # Check of dimensions of tv_eff
  if (!is.null(tv_eff) && any(dim(tv_eff) != c(N_stop, num_events))) {
    stop(sprintf("Dimensions of tv_eff must be (%d, %d)", N_stop, num_events))
  }

  ############################ Default values ##################################

  # Set default values for beta, eta, and nu
  beta <- if(!is.null(beta)) beta else matrix(0, nrow = N_stop, ncol = num_events)
  tv_eff <- if(!is.null(tv_eff)) tv_eff else matrix(0, nrow = N_stop, ncol = num_events)
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

  # In case of time varying effects
  beta_prime <- beta
  beta_prime[1:N_stop,] <- beta[1:N_stop,] + tv_eff

  ############################ Functions #######################################

  # Proportional hazard
  calculate_phi <- function(simmatrix, beta = beta){
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
  lambda <- function(t, i, phi) {
    risk_vec <- at_risk(simmatrix[i, N_start:N_stop])
    risk_vec * eta * nu * t^(nu - 1) * phi[i,]
  }

  # We calculate the inverse numerically
  inverse_sc_haz <- function(p, t, i) {
    inverseScHazTV(p, t, lower = lower, upper = upper, t_prime = t_prime,
                   eta = eta, nu = nu, phi = phi[i,], phi_prime = phi_prime[i,],
                   at_risk = at_risk(simmatrix[i, N_start:N_stop]))
  }

  # Event probabilities
  probs <- function(t, i){
    if(t == max_cens) return(c(1, rep(0, (num_events - 1))))
    probs <- if(t > t_prime) lambda(t, i, phi_prime) else lambda(t, i, phi)
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
    phi <- calculate_phi(simmatrix, beta)
    phi_prime <- calculate_phi(simmatrix, beta_prime)
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




