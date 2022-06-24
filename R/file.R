#' Write ggplot object to a file
#'
#' @param plot ggplot object.
#' @param filename a file name.
#' @param filepath a directory.
#' @param width image width.
#' @param height image height.
#' @param units Plot size in units ("in", "cm", "mm", or "px"). If not supplied, uses the size of current graphics device.
#' @importFrom ggplot2 ggsave
#' @export
save2pdf <- function(plot, filename, filepath, width = 4, height = 4, units = "in") {
    ggsave(filename = filename, plot = plot, path = filepath, width = width, height = height, units = units)
}
