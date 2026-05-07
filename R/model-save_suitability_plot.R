#' Save a ggplot suitability plot to file
#'
#' Save a ggplot object to disk with configurable dimensions, resolution, and
#' typography. File format is inferred from the filename extension or defaults
#' to PNG. Relative paths are resolved against the project directory.
#'
#' @details
#' Width and height are specified in pixels and converted to inches using
#' \code{res} (DPI) for \code{ggplot2::ggsave()}. The plot theme is augmented
#' with the specified font sizes before saving. Supported extensions: \code{png},
#' \code{pdf}, \code{tiff}, \code{jpg}, \code{jpeg}, \code{svg}, \code{eps}.
#'
#' @param p ggplot. The ggplot object to save.
#' @param fn Character. Output path with or without extension. If no extension
#'   is provided, \code{.png} is appended. Relative paths are resolved against
#'   the project directory.
#' @param width Numeric. Image width in pixels. Default 1000.
#' @param height Numeric. Image height in pixels. Default 750.
#' @param res Numeric. Resolution in DPI. Default 72.
#' @param pointsize Numeric. Base point size passed to \code{ggsave()}.
#'   Default 15.
#' @param title_size Numeric. Font size for the plot title. Default 28.
#' @param subtitle_size Numeric. Font size for the plot subtitle. Default 20.
#' @param legend_size Numeric. Font size for legend text. Default 15.
#' @param legend_title_size Numeric. Font size for the legend title. Default 15.
#' @param axis_title_size Numeric. Font size for axis titles. Default 20.
#' @param axis_text_size Numeric. Font size for axis tick labels. Default 18.
#'
#' @return Invisibly returns the file path written. Side effect: writes an
#'   image file to disk.
#'
#' @seealso \code{\link{plot_suitability}}, \code{\link[ggplot2]{ggsave}}
#' @importFrom ggplot2 ggsave theme element_text
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#'   p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
#'   save_suitability_plot(p, "plots/mtcars_scatter.png")
#'
#'   # No extension -> defaults to PNG:
#'   save_suitability_plot(p, "plots/mtcars_scatter")
#' }
#'
#' @export
save_suitability_plot <- function(p,
                                  fn,
                                  width            = 1000,
                                  height           = 750,
                                  res              = 72,
                                  pointsize        = 15,
                                  title_size       = 28,
                                  subtitle_size    = 20,
                                  legend_size      = 15,
                                  legend_title_size = 15,
                                  axis_title_size  = 20,
                                  axis_text_size   = 18) {

  if (!inherits(p, "ggplot")) {
    stop("`p` must be a ggplot object.", call. = FALSE)
  }
  if (!is.character(fn) || length(fn) != 1L || !nzchar(fn)) {
    stop("`fn` must be a non-empty single character string.", call. = FALSE)
  }

  project_dir <- rENM_project_dir()

  # Resolve relative paths against project directory
  if (!grepl("^(/|[A-Za-z]:)", fn)) {
    fn <- file.path(project_dir, fn)
  }

  ext <- tools::file_ext(fn)
  if (identical(ext, "")) {
    fn  <- paste0(fn, ".png")
    ext <- "png"
  }
  ext <- tolower(ext)

  supported <- c("png", "pdf", "tiff", "jpg", "jpeg", "svg", "eps")
  if (!ext %in% supported) {
    stop(
      "Unsupported file extension: '", ext, "'.\n",
      "Supported: ", paste(supported, collapse = ", "),
      call. = FALSE
    )
  }

  dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)

  p <- p +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = title_size,  face = "bold"),
      plot.subtitle = ggplot2::element_text(size = subtitle_size),
      legend.title  = ggplot2::element_text(size = legend_title_size),
      legend.text   = ggplot2::element_text(size = legend_size),
      axis.title    = ggplot2::element_text(size = axis_title_size),
      axis.text     = ggplot2::element_text(size = axis_text_size)
    )

  ggplot2::ggsave(
    filename  = fn,
    plot      = p,
    width     = width  / res,
    height    = height / res,
    dpi       = res,
    units     = "in",
    pointsize = pointsize
  )

  message("Plot saved to: ", fn)
  invisible(fn)
}
