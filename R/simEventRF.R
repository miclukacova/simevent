#' Simulate Survival and Competing Risk Data Based on a fitted RF
#'
#' Simulates survival or competing risk data for a cohort of individuals based
#' on a random forest fitted using the `randomForestSRC` package. Simulation proceeds
#' by iteratively sampling event times until a terminal event occurs.
#'
#' @param N Integer. The number of individuals to simulate.
#' @param RF_fit A rfsrc object. The object contains estimates of the cumulative hazard of each of the processes.
#' @param event_names A character vector. Containing the names of the various processes. The argument is optional.
#' @param L0_old A vector of previously observed baseline covariate values for L0,
#'   used for resampling baseline covariates.
#' @param A0_old A vector of previously observed baseline covariate values for A0,
#'   used for resampling baseline covariates.
#' @param n_event_max Integer vector. Maximum number of times each event type can occur
#'   per individual.
#' @param term_events Integer or integer vector. Indices of event types that are terminal,
#'   i.e., events that stop further simulation for an individual.
#' @param intervention1 Optional function. Not implemented.
#' @param intervention2 Optional function. Not implemented.
#'
#' @details
#' The function simulates individual event histories by:
#' \enumerate{
#'   \item Sampling initial baseline covariates (`L0`, `A0`) by resampling observed values.
#'   \item Extracting cumulative hazard functions from the RF model.
#'   \item Iteratively sampling event times.
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
#' @import randomForestSRC
#' @import data.table
#'
#' @examples
#' # The observed data
#' beta = matrix(c(0.5,-1,-0.5,0.5,0,0.5), ncol = 3, nrow = 2)
#' data <- simCRdata(N = 100, beta = beta)
#' # Random forest model
#' RF_fit <- randomForestSRC::rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data)
#' # New data
#' new_data <- simEventRF(100, RF_fit, L0_old = data$L0, A0_old = data$A0, term_events = c(1,2))
#'
#' @export
simEventRF <- function(N,
                       RF_fit,
                       event_names = NULL,
                       L0_old,
                       A0_old,
                       n_event_max = c(1,1),
                       term_events = 1,
                       intervention1 = NULL,
                       intervention2 = NULL) {

  ID <- NULL

  # Initialize
  alive <- 1:N                                            # Vector for keeping track of who is alive
  num_alive <- N                                          # Number of alive individuals
  T_k <- rep(0, N)                                        # Last event time

  # Data frame for storing data
  sim_data <- data.frame(L0 = sample(L0_old, N, TRUE),
                         A0 = sample(A0_old, N, TRUE))
  if(!is.null(event_names)) for (name in event_names) sim_data[[name]] <- 0

  # List for results
  res_list <- vector("list", sum(n_event_max))            # For results
  idx <- 1                                                # Index

  # The cumulative hazard and inverse cumulative hazard
  y.pred <- predict(RF_fit, sim_data)                     # Cumulative hazard for simulated data
  times <-  c(0,y.pred$time.interest)                     # Time points

  num_events <- if(is.na(dim(y.pred$chf)[3])) 1 else dim(y.pred$chf)[3]         # Number of events
  if(!is.null(event_names)) for (name in event_names) sim_data[[name]] <- 0 else
    for (name in paste0("N", 1:num_events)) sim_data[[name]] <- 0

  # Defining the cumulativ hazard and the inverse cumulative hazard
  cumhaz_fn <- array(vector("list", N * num_events), dim = c(N, num_events))
  invhaz_fn <- array(vector("list", N * num_events), dim = c(N, num_events))

  for(j in seq_len(num_events)) {
    cumhazz <-  if(num_events == 1) y.pred$chf else y.pred$chf[,,j]
    if(num_events)
    for(i in alive){
      cumhaz_fn[[i,j]] <- stats::approxfun(times, c(0, cumhazz[i,]),
                                         method="linear", yright = Inf)
      # We choose ties = max to ensure that event times are strictly increasing
      invhaz_fn[[i,j]] <- stats::approxfun(c(0,cumhazz[i,]), times,
                                         method="linear", rule=2, ties = max)
    }
  }

  # Loop
  while(num_alive != 0){
    # Calculate the cumulative intensity per individual per event
    cum_int_Tk <- matrix(nrow = num_alive, ncol = num_events)
    for(j in seq_len(num_events)) {
      for(i in seq_len(num_alive)){
      cum_int_Tk[i,j] <- cumhaz_fn[[alive[i],j]](T_k[i])
      }
    }

    # Simulate the uniform random variable
    U <- matrix(-log(stats::runif(num_alive * num_events)), ncol = num_events)  # matrix for the random draws
    V <- U + cum_int_Tk

    # Find the event times
    event_times <- matrix(nrow = num_alive, ncol = num_events)
    for(j in seq_len(num_events)) {
      for(i in seq_len(num_alive)){
        event_times[i,j] <- invhaz_fn[[alive[i],j]](V[i,j])
      }
    }

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
