#' Simulate Event History Data Based on Cox Models
#'
#' Simulates recurrent and terminal event data for a cohort of individuals based on a list
#' of fitted Cox proportional hazards models. Each event type is governed by its own model,
#' and simulation proceeds by iteratively sampling event times until a terminal event occurs.
#'
#' @param N Integer. The number of individuals to simulate.
#' @param cox_fits A list of fitted Cox models (`coxph` objects), one for each type of event. Should be named, so that the
#' names can be used as column names.
#' @param L0_old A vector of previously observed baseline covariate values for L0. Used for resampling.
#' @param A0_old A vector of previously observed baseline covariate values for A0. Used for resampling.
#' @param max_events Integer. The maximum number of events to simulate (i.e., simulation rounds).
#' @param n_event_max Integer vector. Maximum number of times each event type can occur per individual.
#' @param term_events Integer or vector of integers. Indices of event types that are terminal (i.e., stop further simulation for an individual).
#' @param intervention description
#'
#' @details
#' The function simulates individual event histories by:
#' \enumerate{
#'   \item Sampling initial covariate values (L0, A0) from provided vectors.
#'   \item Iteratively computing cumulative hazard functions based on provided Cox models.
#'   \item Drawing event times using inverse transform sampling.
#'   \item Updating covariate histories and censoring individuals after terminal events.
#' }
#' Simulation stops for an individual once a terminal event (as defined by `term_events`) occurs.
#'
#' The function supports multiple event types, and different maximum occurrences for each via `n_event_max`.
#'
#' @return A `data.table` with one row per event per individual, including:
#' \itemize{
#'   \item `ID` — Individual identifier.
#'   \item `Time` — Event time.
#'   \item `Delta` — Type of event that occurred.
#'   \item `L0`, `A0` — Baseline covariates.
#'   \item One column for each event type, indicating cumulative counts.
#' }
#'
#' @import survival
#' @import data.table
#' @export
#'
#' @examples
#' # The observed data
#' data_obs <- simT2D(N = 1000)
#' data_obs <- IntFormatData(data_obs, N_cols = 6)
#'
#' # Fit some Cox models
#' cox_death <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1)
#' ~ L0 + A0 + L, data = data_obs)
#' cox_t2d <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2)
#' ~ L0 + A0, data = data_obs[L == 0])
#'
#' # Then simulate new data:
#' cox_fits <- list("D" = cox_death, "L" = cox_t2d)
#' new_data <- simEventCox(100, cox_fits = cox_fits, L0_old = data_obs$L0, A0_old = data_obs$A0)
#'
simEventCox <- function(N,
                        cox_fits,
                        L0_old,
                        A0_old,
                        max_events = 5,
                        n_event_max = c(1,1),
                        term_events = 1,
                        intervention = list("No Int" = NULL)) {

  ID <- NULL

  # Initialize
  num_events <- length(cox_fits)                          # Number of events
  alive <- 1:N                                            # Vector for keeping track of who is alive
  num_alive <- N                                          # Number of alive individuals
  T_k <- rep(0, N)                                        # Last event time

  # Matrix for storing data
  sim_data <- data.frame(L0 = sample(L0_old, N, TRUE),
                         A0 = sample(A0_old, N, TRUE))

  for (name in names(cox_fits)) sim_data[[name]] <- 0


  # List for results
  res_list <- vector("list", max_events)                  # For results
  idx <- 1                                                # Index

  # Function to find the cumulative intensity at time t
  cum_int <- function(int_vals, time_vals, t) {
    idx <- findInterval(t, time_vals)
    ifelse(idx == 0, 0, int_vals[idx])
  }

  # Function to find the inverse cumulative intensity
  inv_cum_int <- function(int_vals, time_vals, p) {
    idx <- findInterval(p, int_vals)
    ifelse(idx == length(time_vals), Inf,  time_vals[idx])
  }

  # Base hazard
  basehazz_list <- lapply(cox_fits, function(model) basehaz(model, centered = FALSE))
  # Jump times of the base hazard
  jump_times <- lapply(basehazz_list, function(df) df[["time"]])
  basehazz_list <- lapply(basehazz_list, function(df) df[["hazard"]])
  if(names(intervention) == "basehaz"){
    basehazz_list <- lapply(basehazz_list, intervention[[1]])
  }

  while(length(alive) != 0){

    # Intervention Cox term
    if(names(intervention) == "cox"){
      cox_term <- list()
      for(j in seq_len(num_events)){
        sim_data_cox <- intervention[[1]](j, sim_data)
        cox_term[[j]] <- exp(stats::predict(cox_fits[[j]], newdata = sim_data_cox, type="lp", reference = "zero"))
      }
      # Calculate the non intervention Cox term
      } else{
      cox_term <- lapply(cox_fits, function(model)
        exp(stats::predict(model, newdata = sim_data, type="lp", reference = "zero")))
    }

    # Calculate the cumulative intensity
    cumInt_list <- mapply(function(base, cox) {
      tcrossprod(base, cox)
    }, basehazz_list, cox_term, SIMPLIFY = FALSE)

    # We find the cumulative intensity of event j and the ith individual at the event time T_k
    cum_int_Tk <- matrix(ncol = num_events, nrow = num_alive)
    for(j in seq_len(num_events)){
      cum_int_Tk[,j] <- sapply(1:num_alive, function(i) cum_int(cumInt_list[[j]][,i], jump_times[[j]], T_k[i]))
    }

    # Simulate the uniform random variable
    U <- matrix(-log(stats::runif(num_alive * num_events)), ncol = num_events)  # matrix for the random draws
    V <- U + cum_int_Tk

    # Find the event times
    event_times <- matrix(ncol = num_events, nrow = num_alive)                  # matrix for event times
    for(j in seq_len(num_events)){
      event_times[,j] <- sapply(1:num_alive, function(i) inv_cum_int(cumInt_list[[j]][,i], jump_times[[j]], V[i,j]))
    }

    # How many times can you experience the various events?
    for(j in seq_len(num_events)){
      event_times[sim_data[, (2+j)] == n_event_max[j], j] <- Inf
    }

    # The next event is the minimum of these events
    T_k <- apply(event_times, 1, min)
    Deltas <- apply(event_times, 1, which.min)

    # Update event counts
    for(i in 1:num_alive){
      sim_data[i, (Deltas[i] + 2)] <- sim_data[i, (Deltas[i] + 2)] + 1
    }

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
