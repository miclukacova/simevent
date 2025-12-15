#' Simulate Survival and Competing Risk Data Based on a fitted RF
#'
#' Simulates survival or competing risk data for a cohort of individuals based
#' on a random forest fitted using the `randomForestSRC` package. Simulation proceeds
#' by iteratively sampling event times until a terminal event occurs.
#'
#' @param N Integer. The number of individuals to simulate.
#' @param RF_fit A rfsrc object. The object contains estimates of the cumulative hazard of each of the processes.
#' @param event_names A character vector. Containing the names of the various processes. The argument is optional.
#' @param list_old_vars A named list containing the old covariates
#' @param n_event_max Integer vector. Maximum number of times each event type can occur
#'   per individual.
#' @param term_events Integer or integer vector. Indices of event types that are terminal,
#'   i.e., events that stop further simulation for an individual.
#' @param intervention1 Optional function. Not implemented.
#' @param intervention2 Optional function. Not implemented.
#' @param parallel Whether to parallelize
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
#' list_old_vars <- list(L0 = data$L0, A0 = data$A0)
#' new_data <- simEventRF(100, RF_fit, list_old_vars = list_old_vars, term_events = c(1,2))
#'
#' @export
simEventRF <- function(N,
                       RF_fit,
                       event_names = NULL,
                       list_old_vars,
                       n_event_max = c(1,1),
                       term_events = 1,
                       intervention1 = NULL,
                       intervention2 = NULL,
                       parallel = FALSE) {

  ID <- chf_mat <- NULL

  # Initialize
  alive <- 1:N                                            # Vector for keeping track of who is alive
  num_alive <- N                                          # Number of alive individuals
  T_k <- rep(0, N)                                        # Last event time

  # Sampling new covariates
  if(is.null(names(list_old_vars))) warning("list_old_vars must be named list")
  num_cov <- length(list_old_vars)
  sim_data <- matrix(ncol = num_cov, nrow = N)
  for(j in 1:num_cov){
    sim_data[,j] <- sample(list_old_vars[[j]], N, TRUE)
  }
  colnames(sim_data) <- names(list_old_vars)
  sim_data <- data.frame(sim_data)

  if(!is.null(event_names)) for (name in event_names) sim_data[[name]] <- 0

  # List for results
  res_list <- vector("list", sum(n_event_max))            # For results
  idx <- 1                                                # Index

  if(parallel == TRUE){
    # Number of cores to use
    n_cores <- parallel::detectCores() - 1   # Leave one core free

    # Split sim_data into roughly equal chunks
    split_data <- split(sim_data, cut(1:nrow(sim_data), n_cores, labels = FALSE))

    # Cluster setup
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist = c("RF_fit"))  # Export the model to the cluster
    parallel::clusterEvalQ(cl, library(randomForestSRC))  # Load required packages on each worker

    # Predict in parallel
    pred_list <- parallel::parLapply(cl, split_data, function(chunk) {
      predict(RF_fit, chunk)
    })

    parallel::stopCluster(cl)

    # Combine results
    chf_mat <- do.call(rbind, lapply(pred_list, function(x) x$chf))
    times <-  c(0,pred_list[[1]]$time.interest)                                   # Time points

  } else {
    # The cumulative hazard and inverse cumulative hazard
    y.pred <- predict(RF_fit, sim_data)                     # Cumulative hazard for simulated data
    times <-  c(0,y.pred$time.interest)                     # Time points
    chf_mat <- y.pred$chf                                   # Cumulative hazard for simulated data
  }

  num_events <- if(is.na(dim(chf_mat)[3])) 1 else dim(chf_mat)[3]               # Number of events
  if(!is.null(event_names)) for (name in event_names) sim_data[[name]] <- 0 else
    for (name in paste0("N", 1:num_events)) sim_data[[name]] <- 0

  # Defining the cumulativ hazard and the inverse cumulative hazard
  cumhaz_fn <- function(t, i, j){
    cumhazz <-  if(num_events == 1) chf_mat else chf_mat[,,j]
    idx <- findInterval(t, times, checkSorted = F)
    t1 <- times[idx]; t2 <- times[idx + 1]
    y1 <- cumhazz[cbind(i,idx)]; y2 <- cumhazz[cbind(i,(idx + 1))]
    y1 + (t - t1) * (y2 - y1) / (t2 - t1)
  }

  invcumhaz_fn <- function(p, i, j){
    cumhazz <-  if(num_events == 1) chf_mat else chf_mat[,,j]
    idx <- sapply(1:length(i), FUN = function(k) findInterval(p[k], cumhazz[i[k],], checkSorted = F))
    idx[idx == ncol(cumhazz)] <- ncol(cumhazz) - 1
    p1 <- cumhazz[cbind(i, idx)]; p2 <- cumhazz[cbind(i,idx + 1)]
    y1 <- times[idx]; y2 <- times[idx + 1]
    # If the cumulative hazard flattens, we choose the smallest time
    y1 + ifelse(p1 == p2, 0, (p - p1) * (y2 - y1) / (p2 - p1))
  }

  # Loop
  while(num_alive != 0){
    # Calculate the cumulative intensity per individual per event
    cum_int_Tk <- matrix(nrow = num_alive, ncol = num_events)
    for(j in seq_len(num_events)) {
      cum_int_Tk[,j] <- cumhaz_fn(T_k, alive, j)
    }

    # Simulate the uniform random variable
    U <- matrix(-log(stats::runif(num_alive * num_events)), ncol = num_events)  # matrix for the random draws
    V <- U + cum_int_Tk

    # Find the event times
    event_times <- matrix(nrow = num_alive, ncol = num_events)
    for(j in seq_len(num_events)) {
      event_times[,j] <- invcumhaz_fn(V[,j], alive, j)
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
  data.table::setkey(res, ID)
  return(res)
}
