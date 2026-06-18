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
#' @param old_vars A named matrix containing the old covariates. New covariates will
#' be simulated by drawing rows from the old covariates with replacement.
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
                        old_vars = NULL,
                        useOldVars = FALSE) {

  ID <- predict2 <- NULL

  # Naming
  if(is.null(colnames(old_vars))) colnames(old_vars) <- paste0("L", 1:ncol(old_vars))
  # Number of covariates
  num_cov <- ncol(old_vars)

  # Simulation matrix
  if(useOldVars){
    sim_data <- old_vars
    N <- nrow(sim_data)
  } else{
    # Sampling new covariates
    sim_data <- data.frame(matrix(ncol = num_cov, nrow = N))
    sim_data[,1:num_cov] <- old_vars[sample(1:nrow(old_vars), N, TRUE),]
  }
  colnames(sim_data) <- colnames(old_vars)

  # Initialize last event time
  T0 <- rep(0, N)

  # The cumulative hazard and inverse cumulative hazard
  y.pred <- predict2(obj, sim_data)

  # Original output
  times <- c(0, y.pred$time)
  chf_raw <- y.pred$chf

  num_events <- dim(chf_raw)[3]

  # Add H(0)=0
  chf_mat <- array(
    0,
    dim = c(dim(chf_raw)[1], dim(chf_raw)[2] + 1, dim(chf_raw)[3])
  )

  chf_mat[, -1, ] <- chf_raw

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
    # We select the relevant chf
    H <-  chf_mat[,,j]

    # If the simulated value is larger than any observed value
    H_max <- H[cbind(seq_along(i), rep(ncol(H), length(i)))]
    t_max <- tail(times, 1)
    too_large <- p >= H_max[i]

    # We find the hazard intervals into which the simulated times fall
    idx <- sapply(1:length(i), FUN = function(k) findInterval(p[k], H[i[k],]))
    idx[idx == ncol(H)] <- ncol(H) - 1

    i_idx <- cbind(i, idx)
    i_idx2 <- cbind(i, idx + 1)

    p1 <- H[i_idx]
    p2 <- H[i_idx2]

    t1 <- times[idx]
    t2 <- times[idx + 1]

    out <- t1 + ifelse(
      p2 == p1,
      0,
      (p - p1) * (t2 - t1) / (p2 - p1)
    )

    out[too_large] <- t_max
    out

  }

  # Calculate the cumulative intensity per individual per event
  # this is kind of not necessary since we currently do not have recurrent events
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
