#' Save a ggplot suitability plot to file
#'
#' Save a ggplot object to disk with configurable size, resolution, and
#' typography. File format is inferred from the extension or defaults to PNG.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' This function saves ggplot outputs generated during the rENM workflow,
#' resolving relative paths against the project directory defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item A ggplot object.
#'   \item An output file path, with or without extension.
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item A saved image file written to disk.
#'   \item If no extension is provided, a PNG file is written.
#' }
#'
#' \strong{Methods}
#' Width and height are specified in pixels and converted to inches using
#' the resolution (DPI) for ggplot2::ggsave().
#'
#' Supported extensions:
#' \itemize{
#'   \item "png", "pdf", "tiff", "jpg", "jpeg", "svg", "eps"
#' }
#'
#' Before saving, the plot theme is enhanced to apply specified font sizes
#' to the title, subtitle, legend, and axis text elements.
#'
#' \strong{Data requirements}
#' The output path is resolved relative to the project directory if a
#' non-absolute path is provided.
#'
#' @param p ggplot. A ggplot object to save.
#' @param fn Character. Full output path and file name, with or without an
#' extension. If no extension is provided, ".png" is appended.
#' @param width Numeric. Image width in pixels. Default 1000.
#' @param height Numeric. Image height in pixels. Default 750.
#' @param res Numeric. Resolution (DPI). Default 72.
#' @param pointsize Numeric. Base point size for text passed to ggsave().
#' Default 15.
#' @param title_size Numeric. Font size for the plot title. Default 28.
#' @param subtitle_size Numeric. Font size for the plot subtitle. Default 20.
#' @param legend_size Numeric. Font size for legend text. Default 15.
#' @param legend_title_size Numeric. Font size for the legend title. Default 15.
#' @param axis_title_size Numeric. Font size for axis titles. Default 20.
#' @param axis_text_size Numeric. Font size for axis tick labels. Default 18.
#'
#' @return
#' Invisibly returns the file path written.
#' \itemize{
#'   \item Character. File path to the saved plot.
#'   \item Side effect: Writes an image file to disk.
#' }
#'
#' @importFrom ggplot2 ggsave theme element_text
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#'   p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
#'          ggplot2::geom_point()
#'   save_suitability_plot(p, "plots/mtcars_scatter.png")
#'
#'   # No extension provided -> defaults to PNG:
#'   save_suitability_plot(p, "plots/mtcars_scatter")
#' }
#'
#' @seealso [ggplot2::ggsave()], [ggplot2::theme()],
#'   [ggplot2::element_text()]
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
  # ---- Argument checks ----
  if (!inherits(p, "ggplot")) {
    stop("`p` must be a ggplot object.", call. = FALSE)
  }
  if (!is.character(fn) || length(fn) != 1L || !nzchar(fn)) {
    stop("`fn` must be a non-empty single character string.", call. = FALSE)
  }

  if (!exists("rENM_project_dir", mode = "function")) {
    stop("Function `rENM_project_dir()` must be available (e.g., exported by this package).")
  }

  project_dir <- rENM_project_dir()

  # Resolve relative paths against project directory
  if (!grepl("^(/|[A-Za-z]:)", fn)) {
    fn <- file.path(project_dir, fn)
  }

  # ---- Determine extension ----
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

  # ---- Ensure output directory exists ----
  dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)

  # ---- Apply font sizes via theme tweaks ----
  p <- p +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = title_size,  face = "bold"),
      plot.subtitle = ggplot2::element_text(size = subtitle_size),
      legend.title  = ggplot2::element_text(size = legend_title_size),
      legend.text   = ggplot2::element_text(size = legend_size),
      axis.title    = ggplot2::element_text(size = axis_title_size),
      axis.text     = ggplot2::element_text(size = axis_text_size)
    )

  # ---- Save (ggsave width/height expect inches) ----
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
