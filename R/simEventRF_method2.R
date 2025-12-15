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
#' @param intervention1 Optional function. Not implemented.
#' @param intervention2 Optional function. Not implemented.
#'
#' @details
#' The function simulates individual event histories by:
#' \enumerate{
#'   \item Sampling initial baseline covariates (`L0`, `A0`) by resampling observed values.
#'   \item Extracting cumulative hazard functions from the RF model.
#'   \item Sampling event times.
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
#' new_data <- simEventRF2(100, RF_fit, L0_old = data$L0, A0_old = data$A0)
#'
#' @export
simEventRF2 <- function(N,
                       RF_fit,
                       event_names = NULL,
                       L0_old,
                       A0_old,
                       intervention1 = NULL,
                       intervention2 = NULL) {

  ID <- NULL

  # Data frame for storing data
  sim_data <- data.frame(ID = 1:N,
                         L0 = sample(L0_old, N, TRUE),
                         A0 = sample(A0_old, N, TRUE))

  # The cumulative hazard and inverse cumulative hazard
  y.pred <- predict(RF_fit, sim_data)                     # Cumulative hazard for simulated data
  times <-  c(0,y.pred$time.interest)                     # Time points
  tp <- length(times) - 1

  # Number of events
  num_events <- if(is.na(dim(y.pred$chf)[3])) 1 else dim(y.pred$chf)[3]
  if(!is.null(event_names)) for (name in event_names) sim_data[[name]] <- 0 else
    for (name in paste0("N", 1:num_events)) sim_data[[name]] <- 0

  # ------------------------- ICF -----------------------------------------

  # We need to be able to tell the rows from one another
  distinction <- matrix(rep(seq(1, N)*10, tp), nrow = N, ncol = tp)
  # For the adherent timepoint
  mat_add <- matrix(c(rep(0,N), rep(1,N)), ncol = 2)

  vec2mat_idx <- function(idx) {
    j <- idx %% tp
    j[j == 0] <- tp
    i <- ceiling(idx / tp)
    return(cbind(i,j))
  }

  invcumhaz_fn <- function(p, j){
    cumhazz <-  if(num_events == 1) y.pred$chf else y.pred$chf[,,j]
    cumhazz_idx <- c(t(cumhazz + distinction))
    p_idx <- p + seq(1,N)*10
    idx <- findInterval(p_idx, cumhazz_idx) |> vec2mat_idx()
    idx[idx[,2] == ncol(cumhazz)] <- ncol(cumhazz) - 1
    idx2 <- (idx + mat_add)
    p1 <- cumhazz[idx]; p2 <- cumhazz[idx2]
    y1 <- times[idx[,2]]; y2 <- times[idx2[,2]]
    # If the cumulative hazard flattens, we choose the smallest time
    y1 + ifelse(p1 == p2, 0, (p - p1) * (y2 - y1) / (p2 - p1))
  }

  # ------------------------- SAMPLING -----------------------------------------

  # Simulate the uniform random variable
  U <- matrix(-log(stats::runif(N * num_events)), ncol = num_events)  # matrix for the random draws

  # Find the event times
  event_times <- matrix(nrow = N, ncol = num_events)

  for(j in seq_len(num_events)) {
    event_times[,j] <- invcumhaz_fn(U[,j], j)
  }

  # The next event is the minimum of these events
  T_k <- apply(event_times, 1, min)
  Deltas <- apply(event_times, 1, which.min)

  # Update event counts
  sim_data[cbind(1:N, Deltas + 2)] <- sim_data[cbind(1:N, Deltas + 2)] + 1

  # Store data
  res <- data.table(ID = 1:N,
                    Time = T_k,
                    Delta = Deltas,
                    sim_data)

  setkey(res, ID)
  return(res)
}
