#' Plot Simulated Event History Data
#'
#' Visualizes event history data by plotting individual event times colored and shaped by event type.
#' Each individual's timeline is displayed horizontally with events marked along it.
#'
#' @param data A \code{data.frame} or \code{data.table} containing at least the columns \code{ID}, \code{Time}, and \code{Delta}.
#' @param title Character string specifying the plot title. Defaults to \code{"Event Data"}.
#'
#' @return A \code{ggplot} object representing the event data visualization.
#' @export
#'
#' @examples
#' data <- simEventData(10)
#' plotEventData(data)

plotEventData <- function(data, title = "Event Data") {
  max_time <- Time <- ID <- Delta <- NULL

  data <- data.table::copy(data)[, c("ID", "Time", "Delta")]

  # Extract unique patient IDs and number of patients
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
  shapes <- c(20, 17, 15, 3, 7, 8, 13)  # varied shapes for clarity

  # Make sure color and shape vectors are long enough
  if (diff_events > length(cols)) {
    warning("More event types than available colors; colors will recycle.")
    cols <- rep(cols, length.out = diff_events)
    shapes <- rep(shapes, length.out = diff_events)
  } else {
    cols <- cols[1:diff_events]
    shapes <- shapes[1:diff_events]
  }

  p <- ggplot2::ggplot(plotdata) +
    ggplot2::geom_line(ggplot2::aes(x = Time, y = ID, group = ID), color = "grey60", size = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = Time, y = ID, shape = factor(Delta), color = factor(Delta)),
                        size = 2.5, data = data, alpha = 0.8) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::scale_shape_manual(values = shapes[1:diff_events]) +
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

  return(p)
}
