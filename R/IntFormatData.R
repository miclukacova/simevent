#' The function transforms the simulated event data containing time dependent covariates
#' into the tstart tstop format used by the coxph function.
#'
#' @title Transformation of Event Data into Interval Format
#'
#' @param data Event data: a data frame containing an ID, Time, Delta, Covariates and Processes columns.
#' @param N_cols The indices of the columns in data that belong to counting processes
#'
#' @return Event data in a tstart tstop format
#' @export
#'
#' @examples
#' data <- simEventData(10)
#' IntFormatData(data)

IntFormatData <- function(data, N_cols = 6:9) {

  k <- ID <- tstart <- tstop <- Time <- NULL
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

  setkey(res, ID)

  return(res)
}
