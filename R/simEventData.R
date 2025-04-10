#' `simEventData` is a function to simulate continuous time to event data, e.g. observational
#' healthcare data. The number of events simulated corresponds to the length of the \eqn{\eta}
#' and \eqn{\nu} vector, as well as the number of columns in the \eqn{\beta} matrix. By
#' default 4 different types of events are simulated. The first two being terminal
#' processes. The simulations build upon a counting process framework, where the events
#' have intensities given by
#' \deqn{\lambda^x(t) = \lambda_0^x \exp(\beta^T_x L)}
#' L here represents covariates and event counts. The baseline intensity \eqn{\lambda_0^x}
#' is given by
#' \deqn{\lambda_0^x(t)=\eta^x \nu^x t^{\nu^x - 1}}
#' The intensities of the various events depend upon previous events and
#' baseline covariates through the \eqn{\beta} matrix. The pre specified \eqn{\nu}
#' and \eqn{\eta} parameters also influence the intensities of the various events
#' through the baseline intensity.
#'
#' @title Simulate Event Data
#'
#' @param N A double for the number of simulated individuals
#' @param beta A matrix of doubles for the effects of covariates and events on the
#' intensities. The columns represent the events N0, N1, .... In the default case 4 events.
#' The rows represent covariates and processes: the first two rows determine the
#' effects of the baseline covariate \eqn{L0} and \eqn{A0} on the processes. The
#' next rows determine the effects of the additional baseline covariates, by default
#'  named L1, L2,... on the processes, followed by the effects of the processes N0,
#'  N1,.... The \eqn{\beta} matrix is by  default set to 0.
#' @param eta Vector of shape parameters for the baseline intensity. Default is set
#' to 0.1 for all events.
#' @param nu Vector of scale parameters for the baseline intensity. Default is set
#' to 1.1 for all events.
#' @param at_risk At risk function. The \code{at_risk} function determines whether
#' an individual is at risk for a specific event. The function takes as input a vector
#' `events` of event counts. The function returns a vector of 0's and 1's indicating
#' which events the subject is at risk for. Default is set to a setting where you
#' are always at risk for all events.
#' @param term_deltas Terminal events. Default is set so that event 0 and 1 are
#' terminal events.
#' @param max_cens A maximum censoring time. By default set to infinity. If the event
#' time is larger than this maximimal censoring time the event becomes event 0 with
#' probability 1.
#' @param add_cov Named list of random generator functions for the distributions of
#' additional baseline covariates. The functions should take the number of observations
#' as input. By default set to NULL.
#' @param override_beta This argument is a named list. The argument has two applications.
#' One possibility is instead of specifying the whole beta matrix, the user can
#' specify the relevant entries, and the rest will by default be 0. Imagine you want
#' to specify the effect of L0 on N1 to be equal to 2, this could be done by
#' override_beta = list("L0" = c("N1" = 2))
#' In general if you want the effect \eqn{\beta_{x,y}= z}, you can specify the
#' override_beta argument as
#' override_beta = list("x" = c("y" = z))
#' Where x and y are names of processes or covariates.
#' @param max_events Number of maximal events per individual
#' @param lower Number of maximal events per individual
#' @param upper Number of maximal events per individual
#' @param gen_A0 Function for generation of A0 covariate. Function of N (number
#' of individuals) and L0 (baseline covariate).
#'
#'
#' @return data.table containing the simulated data. There is a column for ID, time
#' of event (Time), event (Delta), baseline covariate (L0), Baseline Treatment (A0),
#' the count of the various events: N1, N2, .... In case of additional covariates
#' these are included in the data as well, named L1, L2, ....
#' @export
#'
#' @examples
#' simEventData(N = 10)

simEventData <- function(N,                      # Number of individuals
                         beta = NULL,            # Effects
                         eta = NULL,             # Shape parameters
                         nu = NULL,              # Scale parameters
                         at_risk = NULL,         # Function defining the setting
                         term_deltas = c(0,1),   # Terminal events
                         max_cens = Inf,         # Followup time
                         add_cov = NULL,         # Additional baseline covariates
                         override_beta = NULL,   # Argument to easily override entries in beta
                         max_events = 10,        # Number of maximal events per individual
                         lower = 10^(-15),       # Lower bound for inverse cumulative hazard
                         upper = 200,            # Upper bound for inverse cumulative hazard
                         gen_A0 = NULL           # Function for generation of A0
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

  # Filing out beta matrix
  if(!is.null(override_beta)){
    for (bb in 1:length(override_beta)) {
      if (names(override_beta)[bb] %in% rownames(beta)) {
        beta[bb, names(override_beta[[bb]])] <- override_beta[[bb]]
      } else {
        beta <- rbind(beta, matrix(0, nrow = 1, ncol = ncol(beta)))
        beta[nrow(beta), names(override_beta[[bb]])] <- override_beta[[bb]]
        rownames(beta)[nrow(beta)] <- names(override_beta)[bb]
      }
    }
  }

  ############################ Functions #######################################

  # Proportional hazard
  if(nrow(beta) == N_stop){
    calculate_phi <- function(simmatrix) {
      exp(simmatrix %*% beta)
    }
  } else {
    calculate_phi <- function(simmatrix) {
      obj <- data.frame(simmatrix)
      effects <- matrix(nrow = N, ncol = num_events)
      for(i in 1:N){
        terms <- sapply(1:nrow(beta), function(bb) {
          with(obj[i,], eval(parse(text = rownames(beta)[bb]))) * beta[bb,]
        })
        effects[i,] <- rowSums(terms)
      }
      return(exp(effects))
    }
  }

  # Intensities
  lambda <- function(t, i) {
    risk_vec <- at_risk(simmatrix[i, N_start:N_stop])
    risk_vec * eta * nu * t^(nu - 1) * phi[i,]
  }

  # If all events have same parameter, the inverse simplifies
  if(all(nu[1] == nu) && all(eta[1] == eta)){
    inverse_sc_haz <- function(p, t, i) {
      denom <- sum(at_risk(simmatrix[i, N_start:N_stop]) * eta * phi[i,])
      (p / denom + t^nu[1])^(1 / nu[1]) - t
    }
  } else{
    # Summed cumulative hazard
    sum_cum_haz <- function(u, t, i) {
      sum(at_risk(simmatrix[i, N_start:N_stop]) * eta * phi[i,] * ((t + u) ^ nu - t ^ nu))
    }

    # Inverse summed cumulative hazard function
    inverse_sc_haz <- function(p, t, i) {
      root_function <- function(u) sum_cum_haz(u, t, i) - p
      stats::uniroot(root_function, lower = lower, upper = upper)$root
    }
  }

  # Event probabilities
  probs <- function(t, i){
    probs <- lambda(t, i)
    summ <- sum(probs)
    probs / summ
  }

  ############################ Initializing Simulations ########################

  # Draw baseline covariates
  simmatrix[,1] <- stats::runif(N)           # L0
  simmatrix[,2] <- gen_A0(N, L0)             # A0

  # Initialize
  T_k <- rep(0,N)                         # Time 0
  alive <- 1:N                            # Keeping track of who is alive
  res_list <- vector("list", max_events)  # For results
  idx <- 1                                # Index


  ############################ Simulations #####################################

  while(length(alive) != 0){
    # Simulate time
    V <- -log(stats::runif(N))
    phi <- calculate_phi(simmatrix)
    W <- sapply(alive, function(i) inverse_sc_haz(V[i], T_k[i], i))
    T_k[alive] <- T_k[alive] + W

    # Maximal censoring time
    alive[T_k[alive] > max_cens] <- 0
    T_k[T_k[alive] > max_cens] <- max_cens

    # Simulate event
    probs_mat <- sapply(alive, function(i) probs(T_k[i], i), simplify = "array")
    Deltas <- apply(probs_mat, 2, function(p) sample(seq_along(p), 1, prob = p)) - 1

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
    alive <- alive[! Deltas %in% term_deltas]
  }

  res <- data.table::rbindlist(res_list)
  setkey(res, ID)

  return(res)
}

#N = 10                      # Number of individuals
#beta = NULL            # Effects
#eta = NULL             # Shape parameters
#nu = NULL              # Scale parameters
#at_risk = NULL        # Function defining the setting
#term_deltas = c(0,1)   # Terminal events
#max_cens = Inf        # Followup time
#add_cov = NULL        # Additional baseline covariates
#override_beta = NULL   # Argument to easily override entries in beta
#max_events = 10      # Number of maximal events per individual
#lower = 10^(-10)      # Lower bound for inverse cumulative hazard
#upper = 200          # Upper bound for inverse cumulative hazard
#gen_A0 = NULL           # Function for generation of A0
