#' Function to visualize simulated event data.
#'
#' @title Plot Event Data
#'
#' @param data Event data: a data frame containing an ID, Time, Delta, L0 and L column.
#' @param title Title: a string that will be the title of the plot.
#'
#' @return A plot of the event data
#' @export
#'
#' @examples
#' data <- simEventData(10)
#' plotEventData(data)

plotEventData <- function(data, title = "Event Data") {
  max_time <- Time <- ID <- Delta <- NULL

  data <- data.table::copy(data)[, c("ID", "Time", "Delta")]
  n <- length(unique(data$ID))

  # We order according to Time
  ordering <- data[, list(max_time = max(Time)), by = ID]
  setkey(ordering, max_time)

  # We add the start time to the data set
  data[, ID := factor(ID, levels = ordering$ID)]
  plotdata <- rbind(data, data.table(ID = unique(data$ID),
                                     Time = rep(0, n),
                                     Delta = rep("start", n)))

  # Shapes and color for the plot
  diff_events <- length(unique(plotdata$Delta))
  cols <- c("green4", "blue1", "orange", "red2", "lightgreen", "purple1", "yellow")
  shapess <- rep(20, diff_events)

  ggplot2::ggplot(plotdata) +
    ggplot2::geom_line(ggplot2::aes(x = Time, y = ID, group = ID), color = "grey60", size = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = Time, y = ID, shape = factor(Delta), color = factor(Delta)),
                        size = 2.5, data = data, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::scale_shape_manual(values = shapess[1:diff_events]) +
    ggplot2::scale_color_manual(values = cols[1:diff_events]) +
    ggplot2::labs(
      title = title,
      x = "Time",
      y = "Patient ID",
      shape = "Event Type",
      color = "Event Type") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.text.y = ggplot2::element_blank(),
      legend.position = "top"
    )
}
