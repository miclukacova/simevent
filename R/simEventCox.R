#' Simulate Event History Data Based on Cox Models
#'
#' Simulates recurrent and terminal event data for a cohort of individuals based on a list
#' of fitted Cox proportional hazards models. Each event type is governed by its own model,
#' and simulation proceeds by iteratively sampling event times until a terminal event occurs.
#'
#' @param N Integer. The number of individuals to simulate.
#' @param cox_fits A named list of fitted Cox proportional hazards models (`coxph` objects),
#'   one for each event type. The names are used as event type labels.
#' @param L0_old A vector of previously observed baseline covariate values for L0,
#'   used for resampling baseline covariates.
#' @param A0_old A vector of previously observed baseline covariate values for A0,
#'   used for resampling baseline covariates.
#' @param max_events Integer. The maximum number of events to simulate per individual.
#' @param n_event_max Integer vector. Maximum number of times each event type can occur
#'   per individual.
#' @param term_events Integer or integer vector. Indices of event types that are terminal,
#'   i.e., events that stop further simulation for an individual.
#' @param intervention1 Optional function. Takes arguments `(j, sim_matrix)` and returns
#'   an updated simulation matrix. Used to modify covariates dynamically at each event iteration.
#' @param intervention2 Optional function. Takes arguments `(j, H_j)` and returns a modified
#'   baseline cumulative hazard vector for event type `j`. Allows dynamic hazard modification.
#'
#' @details
#' The function simulates individual event histories by:
#' \enumerate{
#'   \item Sampling initial baseline covariates (`L0`, `A0`) by resampling observed values.
#'   \item Extracting baseline cumulative hazard functions from the Cox models.
#'   \item Iteratively sampling event times using inverse transform sampling based on hazards.
#'   \item Updating covariate histories and event counts.
#'   \item Stopping simulation per individual after a terminal event or maximum events reached.
#' }
#'
#' @return A `data.table` with one row per event per individual containing:
#' \itemize{
#'   \item `ID` — Individual identifier.
#'   \item `Time` — Event time.
#'   \item `Delta` — Event type indicator.
#'   \item Baseline covariates `L0`, `A0`.
#'   \item Columns for each event type indicating cumulative event counts.
#' }
#'
#' @import survival
#' @import data.table
#'
#' @examples
#' # The observed data
#' data_obs <- simDisease(N = 1000)
#' data_obs <- IntFormatData(data_obs, N_cols = 6)
#'
#' # Fit Cox models
#' cox_death <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1)
#' ~ L0 + A0 + L, data = data_obs)
#' cox_Disease <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2)
#' ~ L0 + A0, data = data_obs[L == 0])
#'
#' # Then simulate new data:
#' cox_fits <- list("D" = cox_death, "L" = cox_Disease)
#' new_data <- simEventCox(100, cox_fits = cox_fits, L0_old = data_obs$L0, A0_old = data_obs$A0)
#'
#' @export
simEventCox <- function(N,
                        cox_fits,
                        L0_old,
                        A0_old,
                        max_events = 5,
                        n_event_max = c(1,1),
                        term_events = 1,
                        intervention1 = NULL,
                        intervention2 = NULL) {

  ID <- NULL

  # Initialize
  num_events <- length(cox_fits)                          # Number of events
  alive <- 1:N                                            # Vector for keeping track of who is alive
  num_alive <- N                                          # Number of alive individuals
  T_k <- rep(0, N)                                        # Last event time

  # Data frame for storing data
  sim_data <- data.frame(L0 = sample(L0_old, N, TRUE),
                         A0 = sample(A0_old, N, TRUE))
  for (name in names(cox_fits)) sim_data[[name]] <- 0

  # List for results
  res_list <- vector("list", max_events)                  # For results
  idx <- 1                                                # Index

  # Base hazard
  basehazz_list <- lapply(cox_fits, function(model) basehaz(model, centered = FALSE))

  # The cumulative hazard and inverse cumulative hazard
  cumhaz_fn <- vector("list", num_events)
  invhaz_fn <- vector("list", num_events)
  for(j in seq_len(num_events)) {
    H_j <- c(0, basehazz_list[[j]][["hazard"]])
    t_j <- c(0, basehazz_list[[j]][["time"]])
    if(!is.null(intervention2)) H_j <- intervention2(j, H_j)
    cumhaz_fn[[j]] <- stats::approxfun(t_j,       H_j,
                                       method="linear", yright = Inf)
    # We choose ties = max to ensure that event times are strictly increasing
    invhaz_fn[[j]] <- stats::approxfun(H_j,       t_j,
                                       method="linear", rule=2, ties = max)
  }

  # Loop
  while(num_alive != 0){
    # Intervention Cox term
    if(!is.null(intervention1)){
      cox_term <- list()
      for(j in seq_len(num_events)){
        sim_data_cox <- intervention1(j, sim_data)
        cox_term[[j]] <- exp(stats::predict(cox_fits[[j]], newdata = sim_data_cox, type="lp", reference = "zero"))
        }
      # Calculate the non intervention Cox term
      } else{
      cox_term <- lapply(cox_fits, function(model)
        exp(stats::predict(model, newdata = sim_data, type="lp", reference = "zero")))
    }

    # Calculate the cumulative intensity per individual per event
    cum_int_Tk <- sapply(seq_len(num_events), function(j) {
      cumhaz_fn[[j]](T_k) * cox_term[[j]]
    })

    # Simulate the uniform random variable
    U <- matrix(-log(stats::runif(num_alive * num_events)), ncol = num_events)  # matrix for the random draws
    V <- U + cum_int_Tk

    # Find the event times
    event_times <- sapply(seq_len(num_events), function(j) {
      invhaz_fn[[j]](V[,j] / cox_term[[j]])
    })

    # How many times can you experience the various events?
    for(j in seq_len(num_events)){
      event_times[sim_data[, (2+j)] == n_event_max[j], j] <- Inf
    }

    # The next event is the minimum of these events
    T_k <- apply(event_times, 1, min)
    Deltas <- apply(event_times, 1, which.min)

    # Update event counts
    sim_data[cbind(seq_len(num_alive), Deltas + 2)] <- sim_data[cbind(seq_len(num_alive), Deltas + 2)] + 1

    # Store data
    kth_event <- data.table(ID = alive,
                            Time = T_k,
                            Delta = Deltas)

    res_list[[idx]] <- cbind(kth_event, data.table::as.data.table(sim_data))
    idx <- idx + 1

    # Who is still alive?
    alive <- alive[!(Deltas %in% term_events)]
    num_alive <- length(alive)
    # For the next iteration we only keep data from individuals alive
    T_k <- T_k[!(Deltas %in% term_events)]
    sim_data <- sim_data[!(Deltas %in% term_events), , drop = FALSE]
  }

  res <- data.table::rbindlist(res_list)
  setkey(res, ID)
  return(res)
}
