#' The function transforms the simulated event data containing time dependent covariates
#' into a format that can be used for classical inference. Three variables are added,
#' called `tstart`, `tstop` and `k`. The two former are used in for example the coxph
#' function, and the latter counts the number of events.
#'
#' @title Transformation of Event Data into Interval Format
#'
#' @param data Event data: a data frame containing an ID, Time, Delta, Covariates and Processes columns.
#' @param N_cols The indices of the columns in data that belong to counting processes
#' @param timeVar Logical indicating whether the effects are time varying
#' @param t_prime The time where the effects change
#'
#' @return Event data in a tstart tstop format
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
