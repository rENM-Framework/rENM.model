#' Create and save a binary presence-absence range map
#'
#' Build a ggplot2 map for a binary (0/1) presence-absence raster, overlaid
#' with state boundaries and the species GAP range polygon. Writes PNG, ASC,
#' and GeoTIFF outputs to the per-year model directory.
#'
#' @details
#' Output files are written to:
#' \code{<project_dir>/runs/<alpha_code>/TimeSeries/<year>/model/}
#'
#' The species range shapefile is resolved via \code{get_species_info(alpha_code)},
#' which looks up \code{GAP.RANGE} in the species table and reads the shapefile
#' from \code{<project_dir>/data/shapefiles/<gap_range>/<gap_range>.shp}.
#' State and range vectors are reprojected to the raster CRS if necessary.
#' Values in the raster are coerced to binary by \code{round()}.
#' Sidecar files (\code{.prj}, \code{.aux.xml}) are removed after writing.
#'
#' @param pa SpatRaster, Character, or Raster. A presence-absence raster.
#' @param alpha_code Character. Species code (e.g., \code{"CASP"}), used to
#'   resolve output paths and to look up \code{GAP.RANGE} via
#'   \code{get_species_info()}.
#' @param year Integer or Character. Year inserted into the output path and name.
#' @param state_shp Character. Path to a states/administrative-boundaries
#'   shapefile. Default resolves to
#'   \code{<project_dir>/data/shapefiles/tl_2012_us_state/tl_2012_us_state.shp}.
#' @param title Character. Plot title. Default \code{"Range Map"}.
#' @param layer Integer. Layer index when \code{pa} has multiple layers.
#'   Default 1L.
#' @param minor_div Integer. Minor grid divisions between major axis breaks.
#'   Default 2L.
#' @param width_px,height_px Numeric. PNG dimensions in pixels. Defaults
#'   1000 and 750.
#' @param dpi Numeric. Plot DPI for \code{ggsave()}. Default 72.
#' @param title_size,legend_size,legend_title_size,axis_title_size,axis_text_size
#'   Numeric. Font sizes for plot elements.
#'
#' @return Invisible list with \code{plot} (ggplot2 object), \code{img_path},
#'   \code{asc_path}, \code{tif_path}. Side effects are raster and image files
#'   written to disk and sidecar files removed.
#'
#' @seealso \code{\link{plot_suitability}}, \code{\link{rENM_project_dir}}
#' @importFrom terra rast vect crs same.crs project ext nlyr writeRaster as.data.frame
#' @importFrom sf st_as_sf st_crs
#' @importFrom ggplot2 ggplot geom_raster aes_string scale_fill_manual geom_sf labs
#'   coord_sf scale_x_continuous scale_y_continuous theme_bw theme element_text ggsave
#'
#' @examples
#' \dontrun{
#'   create_range_map("rasters/casp_1980_pa.tif", "casp", 1980)
#'
#'   r <- terra::rast("rasters/casp_stack.tif")
#'   create_range_map(r, "casp", 1990, layer = 2L, title = "CASP 1990 Range")
#' }
#'
#' @export
create_range_map <- function(pa,
                             alpha_code,
                             year,
                             state_shp  = NULL,
                             title      = "Range Map",
                             layer      = 1L,
                             minor_div  = 2L,
                             width_px   = 1000,
                             height_px  = 750,
                             dpi        = 72,
                             title_size        = 25,
                             legend_size       = 15,
                             legend_title_size = 15,
                             axis_title_size   = 20,
                             axis_text_size    = 18) {

  if (!requireNamespace("terra",   quietly = TRUE)) stop("Package 'terra' is required.", call. = FALSE)
  if (!requireNamespace("sf",      quietly = TRUE)) stop("Package 'sf' is required.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)

  if (missing(alpha_code) || !is.character(alpha_code) || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character string (e.g., 'casp').")
  }

  project_dir <- rENM_project_dir()

  # ---- Resolve default state shapefile ----
  if (is.null(state_shp)) {
    state_shp <- file.path(project_dir, "data", "shapefiles",
                           "tl_2012_us_state", "tl_2012_us_state.shp")
  }

  fn <- file.path(
    project_dir, "runs", alpha_code, "TimeSeries", year,
    "model", paste0(alpha_code, "-", year, "-Range")
  )
  dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)

  img_path <- paste0(fn, ".png")
  asc_path <- paste0(fn, ".asc")
  tif_path <- paste0(fn, ".tif")

  code <- toupper(alpha_code)

  si <- get_species_info(code)
  if (is.null(si)) stop("get_species_info() returned NULL for alpha_code: ", code)
  if (!("GAP.RANGE" %in% names(si))) stop("`GAP.RANGE` column not found for alpha_code: ", code)

  gap_range <- as.character(si[["GAP.RANGE"]])[1]
  if (!nzchar(gap_range) || is.na(gap_range)) stop("Missing GAP.RANGE for alpha_code: ", code)

  range_shp <- file.path(project_dir, "data", "shapefiles",
                          gap_range, paste0(gap_range, ".shp"))

  # ---- Load PA raster ----
  if (inherits(pa, "SpatRaster")) {
    r <- pa
  } else if (is.character(pa) && length(pa) == 1L) {
    if (!file.exists(pa)) stop("Raster file not found: ", pa)
    r <- terra::rast(pa)
  } else if (inherits(pa, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    r <- terra::rast(pa)
  } else {
    stop("`pa` must be a terra::SpatRaster, a file path, or a raster::Raster*.")
  }

  if (terra::nlyr(r) < layer) stop("Requested layer ", layer, " but raster has ", terra::nlyr(r), " layer(s).")
  if (terra::nlyr(r) > 1L) r <- r[[layer]]

  if (!file.exists(state_shp)) stop("State shapefile not found: ", state_shp)
  if (!file.exists(range_shp)) {
    stop("Species range shapefile not found at: ", range_shp,
         "\n(derived from GAP.RANGE = '", gap_range, "')")
  }

  stl <- terra::vect(state_shp)
  spv <- terra::vect(range_shp)

  if (!is.na(terra::crs(r))) {
    if (!terra::same.crs(stl, r)) stl <- terra::project(stl, r)
    if (!terra::same.crs(spv, r)) spv <- terra::project(spv, r)
  }

  stl_sf <- sf::st_as_sf(stl)
  spv_sf <- sf::st_as_sf(spv)

  r_df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(r_df)[3] <- "value"
  r_df <- r_df[!is.na(r_df$value), , drop = FALSE]
  r_df$value <- factor(
    as.integer(round(r_df$value)),
    levels = c(1, 0),
    labels = c("Presence", "Absence")
  )

  e    <- terra::ext(r)
  xlim <- c(e$xmin, e$xmax)
  ylim <- c(e$ymin, e$ymax)

  pretty_with_minors <- function(lim, n_major = 6, minor_div = 2L) {
    majors <- pretty(lim, n = n_major)
    step   <- unique(round(diff(majors), 10))[1]
    if (is.na(step) || step == 0) return(list(majors = majors, minors = NULL))
    minor_step <- step / max(1L, minor_div + 1L)
    minors <- c(vapply(
      majors[-length(majors)],
      function(m0) m0 + seq_len(minor_div) * minor_step,
      numeric(minor_div)
    ))
    list(majors = majors, minors = minors)
  }

  bx <- pretty_with_minors(xlim, 6, minor_div)
  by <- pretty_with_minors(ylim, 6, minor_div)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data    = r_df,
      mapping = ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_manual(
      values = c("Presence" = "#F2C655", "Absence" = "gray95"),
      breaks = c("Presence", "Absence"),
      name   = "Range"
    ) +
    ggplot2::geom_sf(data = stl_sf, fill = NA, color = "gray40", linewidth = 0.3) +
    ggplot2::geom_sf(data = spv_sf, fill = NA, color = "black", linewidth = 0.7) +
    ggplot2::labs(title = title, x = "Longitude", y = "Latitude") +
    ggplot2::coord_sf(
      xlim = xlim, ylim = ylim, expand = FALSE,
      crs  = if (!is.na(sf::st_crs(stl_sf))) sf::st_crs(stl_sf) else NA
    ) +
    ggplot2::scale_x_continuous(
      breaks       = bx$majors,
      minor_breaks = bx$minors,
      expand       = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks       = by$majors,
      minor_breaks = by$minors,
      expand       = c(0, 0)
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = title_size),
      legend.title     = ggplot2::element_text(size = legend_title_size),
      legend.text      = ggplot2::element_text(size = legend_size),
      axis.title       = ggplot2::element_text(size = axis_title_size),
      axis.text        = ggplot2::element_text(size = axis_text_size),
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = "grey85"),
      panel.grid.minor = ggplot2::element_line(linewidth = 0.2, colour = "grey92"),
      legend.position  = "right"
    )

  ggplot2::ggsave(
    filename = img_path,
    plot     = p,
    width    = width_px  / dpi,
    height   = height_px / dpi,
    dpi      = dpi,
    units    = "in"
  )

  r_out <- round(r)
  terra::writeRaster(r_out, asc_path, overwrite = TRUE, NAflag = -9999)
  terra::writeRaster(r_out, tif_path, overwrite = TRUE, NAflag = -9999)

  sidecars <- c(
    sub("\\.asc$", ".prj",         asc_path),
    sub("\\.asc$", ".asc.aux.xml", asc_path),
    sub("\\.tif$", ".prj",         tif_path),
    sub("\\.tif$", ".tif.aux.xml", tif_path)
  )
  for (f in sidecars) {
    if (file.exists(f)) try(file.remove(f), silent = TRUE)
  }

  message(
    "Saved:",
    "\n  Plot:   ", img_path,
    "\n  ASCII:  ", asc_path,
    "\n  GeoTIFF:", tif_path
  )

  invisible(list(
    plot     = p,
    img_path = img_path,
    asc_path = asc_path,
    tif_path = tif_path
  ))
}
