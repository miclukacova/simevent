#' Transform Event Data into Interval Format for Classical Inference
#'
#' Converts simulated event history data with time-dependent covariates into an interval (start-stop) format,
#' suitable for classical survival analysis functions like \code{coxph}.
#' Adds interval start and stop times (\code{tstart}, \code{tstop}) and a counting variable \code{k} indexing events.
#' Optionally, the function can split intervals at a specified time point to accomadate estimation of time-varying effects.
#'
#' @param data A \code{data.frame} or \code{data.table} containing event data with columns \code{ID}, \code{Time}, \code{Delta},
#'   and counting process columns specified by \code{N_cols}.
#' @param N_cols Integer vector. Column indices of \code{data} that correspond to counting process variables.
#'   Defaults to \code{6:9}.
#' @param timeVar Logical. If \code{TRUE}, the intervals are split at \code{t_prime} to allow time-varying covariate effects. Default is \code{FALSE}.
#' @param t_prime Numeric scalar. Time point at which to split intervals if \code{timeVar = TRUE}.
#'
#' @return A \code{data.table} with columns \code{tstart}, \code{tstop}, \code{k}, and other original variables,
#'   formatted for survival analysis.
#' @export
#'
#' @examples
#' data <- simEventData(10)
#' IntFormatData(data)

IntFormatData <- function(data, N_cols = 6:9, timeVar = FALSE, t_prime = NULL) {

  k <- ID <- tstart <- tstop <- Time <- t_group <- NULL
  data <- copy(data)

  # Shifting values of counting processes
  data[, (N_cols) := lapply(.SD, data.table::shift, fill = 0), by = ID, .SDcols = N_cols]

  # Creating k variable
  data[, k:= stats::ave(ID, ID, FUN = seq_along)]
  max_k <- max(data$k)

  # Starting the new data format
  data_k <- list()
  data_k[[1]] <- data[data$k == 1,]

  data_k[[1]][, tstart := 0]
  data_k[[1]][, tstop := Time]

  # Going through all events
  for(i in 2:max_k){
    data_k[[i]] <- data[data$k == i,]

    data_k[[i]][, tstart := data_k[[i-1]][data_k[[i-1]]$ID %in% data_k[[i]]$ID,]$tstop]
    data_k[[i]][, tstop := Time]
  }

  res <- do.call(rbind, data_k)

  if(timeVar == TRUE){
    # We select the rows we need to split
    rows_to_split <- res[tstart <= t_prime & tstop > t_prime]

    # We create a new row marking the time point change
    data1 <- copy(rows_to_split)[, `:=` (tstop = t_prime, Time = t_prime, Delta = -1)]
    data2 <- copy(rows_to_split)[, `:=` (tstart = t_prime)]

    # We take all the remaining rows
    data_new <- res[!(tstart <= t_prime & tstop > t_prime)]
    res <- rbind(data_new, data1, data2)

    # We create a new variable indicating the time group
    res[, t_group := (1+(Time > t_prime))]
  }

  setorder(res, ID, tstart)
  return(res)
}
