#' Simulate Survival and Competing Risk Data Based on a General Model
#'
#' The `simEventObj` function simulates survival or competing risk data for a cohort
#' of individuals based on a general model with a `predict2` method. The function is useful
#' for simulating additional data under the same distribution as an original data set.
#' The procedure consists of fitting a model, such as a random forest or Cox Proportional
#' Hazards model on an original data set. Next the model is equipped with a `predict2` method,
#' and passed as an argument to the `simEventObj` function, which simulates new data using
#' the `predict2` method. The method should output the cumulative hazard array and the
#' jump times of the cumulative hazard. Simulation proceeds by sampling from the uniform
#' distribution and obtaining event times using the inverse of the cumulative hazard function(s).
#'
#' @param N Integer. The number of individuals to simulate.
#' @param obj An object of class `simevent`. The object should have a predict2 method.
#' The method should return a list containing `chf` and `time`. `chf` should be an array
#' of dimension Individuals x Times x Events. The array should contain the cumulative
#' hazard values for each individual, at each time for each event. `time` should be
#' a vector of times where the cumulative hazard function jumps.
#' @param event_names A character vector. Containing the names of the various processes.
#' The argument is optional. By default events will be named `N1`, `N2`, ....
#' @param list_old_vars A named list containing the old covariates. New covariates will
#' be simulated by drawing from the old covariates with replacement.
#'
#' @details
#' The function simulates individual event histories by:
#' \enumerate{
#'   \item Sampling initial baseline covariates by resampling observed values.
#'   \item Extracting cumulative hazard functions from the object.
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
#'   \item Baseline covariates.
#'   \item Columns for each event type indicating cumulative event counts.
#' }
#'
#' @import data.table
#'
#' @export
simEventObj <- function(N,
                        obj,
                        event_names = NULL,
                        list_old_vars = NULL) {

  ID <- predict2 <- NULL

  # Initialize last event time
  T0 <- rep(0, N)

  # Sampling new covariates
  if(is.null(names(list_old_vars))) warning("list_old_vars must be named list")
  num_cov <- length(list_old_vars)
  sim_data <- data.frame(matrix(ncol = num_cov, nrow = N))
  for(j in 1:num_cov){
    sim_data[,j] <- sample(list_old_vars[[j]], N, TRUE)
  }
  colnames(sim_data) <- names(list_old_vars)

  # The cumulative hazard and inverse cumulative hazard
  y.pred <- predict2(obj, sim_data)                                             # Predictions
  times <-  c(0,y.pred$time)                                                    # Time points
  chf_mat <- y.pred$chf                                                         # Cumulative hazard for simulated data
  num_events <- dim(chf_mat)[3]                                                 # Number of events

  # Creating columns for event counts
  if(!is.null(event_names)) for (name in event_names) sim_data[[name]] <- 0 else
    for (name in paste0("N", 1:num_events)) sim_data[[name]] <- 0

  # Defining the cumulativ hazard and the inverse cumulative hazard
  cumhaz_fn <- function(t, i, j){
    idx <- findInterval(t, times)
    t1 <- times[idx]; t2 <- times[idx + 1]
    y1 <- chf_mat[cbind(i,idx, j)]; y2 <- chf_mat[cbind(i,(idx + 1), j)]
    y1 + (t - t1) * (y2 - y1) / (t2 - t1)
  }

  invcumhaz_fn <- function(p, i, j){
    cumhazz <-  chf_mat[,,j]
    idx <- sapply(1:length(i), FUN = function(k) findInterval(p[k], cumhazz[i[k],]))
    idx[idx == ncol(cumhazz)] <- ncol(cumhazz) - 1
    p1 <- cumhazz[cbind(i, idx)]; p2 <- cumhazz[cbind(i,idx + 1)]
    y1 <- times[idx]; y2 <- times[idx + 1]
    # If the cumulative hazard flattens, we choose the smallest time
    y1 + ifelse(p1 == p2, 0, (p - p1) * (y2 - y1) / (p2 - p1))
  }

  # Calculate the cumulative intensity per individual per event
  cum_int_Tk <- matrix(nrow = N, ncol = num_events)
  for(j in seq_len(num_events)) {
    cum_int_Tk[,j] <- cumhaz_fn(T0, 1:N, j)
  }

  # Simulate the uniform random variable
  U <- matrix(-log(stats::runif(N * num_events)), ncol = num_events)  # matrix for the random draws
  V <- U + cum_int_Tk

  # Find the event times
  event_times <- matrix(nrow = N, ncol = num_events)
  for(j in seq_len(num_events)) {
    event_times[,j] <- invcumhaz_fn(V[,j], 1:N, j)
  }

  # The next event is the minimum of these events
  T_k <- apply(event_times, 1, min)
  Deltas <- apply(event_times, 1, which.min)

  # Update event counts
  for(i in 1:num_events){
    sim_data[seq_len(N), num_cov + i] <- sim_data[seq_len(N), num_cov + i] + ifelse(Deltas == i, 1, 0)
  }

  # Store data
  kth_event <- data.table(ID = 1:N,
                          Time = T_k,
                          Delta = Deltas)

  res <- cbind(kth_event, data.table::as.data.table(sim_data))
  data.table::setkey(res, ID)

  return(res)
}
