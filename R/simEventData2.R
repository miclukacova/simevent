#' `simEventData2` is a function to simulate continuous time to event data, e.g. observational
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
#' next rows determine the effects of the processes N0, N1,..., followed by additional
#' baseline covariates, by default named L1, L2,.... The \eqn{\beta} matrix is by
#' default set to 0.
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
#' @param override_beta This argument is a named list containing named vectors.
#' override_beta has two applications. One possibility is instead of specifying
#' the whole beta matrix, the user can specify the relevant entries, and the rest
#' will by default be 0. Imagine you want to specify the effect of L0 on N1 to be
#' equal to 2, this could be done by override_beta = list("L0" = c("N1" = 2))
#' In general if you want the effect \eqn{\beta_{x,y}= z}, you can specify the
#' override_beta argument as
#' override_beta = list("x" = c("y" = z))
#' Where x and y are names of processes or covariates.
#' The second application is for specifying non linearities in case of recurrent
#' events.
#' @param max_events Number of maximal events per individual
#'
#' @return data.table containing the simulated data. There is a column for ID, time
#' of event (Time), event (Delta), baseline covariate (L0), Baseline Treatment (A0),
#' the count of the various events: N1, N2, .... In case of additional covariates
#' these are included in the data as well, named L1, L2, ....
#' @export
#'
#' @examples
#' #simEventData2(10)

simEventData2 <- function(N,                     # Number of individuals
                         beta = NULL,            # Effects
                         eta = NULL,             # Shape parameters
                         nu = NULL,              # Scale parameters
                         at_risk = NULL,         # Function defining the setting
                         term_deltas = c(0,1),   # Terminal events
                         max_cens = Inf,         # Followup time
                         add_cov = NULL,         # Additional baseline covariates
                         override_beta = NULL,   # For overriding beta entries
                         max_events = 10         # Maximal number of events
){
  ID <- NULL

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

  # Set default values for beta, eta, and nu
  beta <- if (!is.null(beta)) beta else matrix(0, nrow = 2 + num_events + num_add_cov, ncol = num_events)
  eta  <- if (!is.null(eta)) eta   else rep(0.1, num_events)
  nu   <- if (!is.null(nu)) nu     else rep(1.1, num_events)

  # Check of dimensions
  if(num_events != length(nu) || num_events != ncol(beta)){
    stop("Length of eta should be equal to nu and number of columns of beta")
  }
  if((num_events + num_add_cov + 2) != nrow(beta)){
    stop("Number of rows in beta should equal the sum of number of event and
         number of additional covariates + 2")
  }

  # Default at_risk
  if(is.null(at_risk)){
    at_risk <- function(events) return(rep(1,num_events))
  }

  # Intensities
  phi <- function(i) {
    browser()
    effects <- numeric(num_events)

    for (bb in seq_len(nrow(beta))) {
      rowname_expr <- rownames(beta)[bb]

      # Retrieve values efficiently
      if (rowname_expr %in% names(L1)) {
        value <- L1[i, rowname_expr]
      } else if (exists(rowname_expr)) {
        value <- get(rowname_expr)[i]
      } else {
        modified_expr <- gsub("N([0-9]+)", "event_counts[i, 'N\\1']", rowname_expr, perl = TRUE)
        value <- eval(parse(text = modified_expr))
      }
      effects <- effects + value * beta[bb, ]
    }
    exp(effects)
  }

  lambda <- function(t, i) {
    risk_vec <- at_risk(event_counts[i, ])
    risk_vec * eta * nu * t^(nu - 1) * phi(i)
  }

  # If all events have same parameter, the inverse simplifies
  if(all(nu[1] == nu)){
    inverse_sc_haz <- function(p, t, i) {
      denom <- sum(at_risk(event_counts[i,]) * eta * phi(i))
      (p / denom + t^nu[1])^(1 / nu[1]) - t
    }
  } else{
    # Summed cumulative hazard
    sum_cum_haz <- function(u, t, i) {
      sum(at_risk(event_counts[i,]) * eta * phi(i) * ((t + u) ^ nu - t ^ nu))
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
      probs <- c(1, rep(0, num_events - 1))
    }
    else{
      probs <- lambda(t, i)
      summ <- sum(probs)
      probs / summ
    }
  }

  # Draw baseline covariates
  L0 <- stats::runif(N)
  A0 <- stats::rbinom(N, 1, 0.5)

  # Generate additional covariates if distributions are specified
  if (num_add_cov != 0) {
    L1 <- as.data.frame(lapply(add_cov, function(f) f(N)))
    colnames(L1) <- paste0("L", seq_len(ncol(L1)))
  } else {
    L1 <- NULL
  }

  # Naming of beta matrix
  colnames(beta) <- paste0("N", seq(0,num_events -1))
  rownames(beta) <- c("L0", "A0", colnames(beta), colnames(L1))

  # Filing out beta matrix
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

  # Initialize
  T_k <- rep(0,N)                         # Time 0
  alive <- 1:N                            # Keeping track of who is alive
  res_list <- vector("list", max_events)  # For results
  idx <- 1                                # Index

  # Create a matrix for the event counts
  event_counts <- matrix(0, nrow = N, ncol = num_events)
  colnames(event_counts) <- colnames(beta)

  while(length(alive) != 0){
    # Simulate time
    V <- -log(stats::runif(N))
    W <- sapply(alive, function(i) inverse_sc_haz(V[i], T_k[i], i))
    T_k[alive] <- T_k[alive] + W

    # Simulate event
    Deltas <- sapply(alive, function(i) sample(seq_len(num_events), size = 1, prob = probs(T_k[i], i)) - 1)

    # Update event counts
    for (j in seq_len(num_events)) {
      event_counts[alive, j] <- event_counts[alive, j] + (Deltas == (j - 1))
    }

    # Store data
    kth_event <- data.table(ID = alive,
                            Time = T_k[alive],
                            Delta = Deltas,
                            L0 = L0[alive],
                            A0 = A0[alive])

    res_list[[idx]] <- cbind(kth_event, L1[alive, ], event_counts[alive, ])
    idx <- idx + 1

    # Who is still alive and uncensored?
    alive <- alive[! Deltas %in% term_deltas]
  }

  res <- rbindlist(res_list)
  setkey(res, ID)

  return(res)
}
