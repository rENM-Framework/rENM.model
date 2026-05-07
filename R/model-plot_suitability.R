#' Plot a climatic suitability raster
#'
#' Build a ggplot2 map of continuous climatic suitability values from a raster,
#' overlaid with state boundaries and the species GAP range polygon.
#'
#' @details
#' The species range shapefile is resolved via \code{get_species_info(alpha_code)},
#' which looks up \code{GAP.RANGE} and reads the shapefile from
#' \code{<project_dir>/data/shapefiles/<gap_range>/<gap_range>.shp}.
#' State and range vectors are reprojected to the raster CRS if necessary.
#' NA raster cells render transparent.
#'
#' @param en SpatRaster, Character, or Raster. A suitability raster.
#' @param alpha_code Character. Species code used to look up \code{GAP.RANGE}
#'   via \code{get_species_info()}.
#' @param state_shp Character. Path to a states/administrative-boundaries
#'   shapefile. Default resolves to
#'   \code{<project_dir>/data/shapefiles/tl_2012_us_state/tl_2012_us_state.shp}.
#' @param title Character. Plot title. Default \code{"Climatic Suitability Map"}.
#' @param limits Numeric. Length-2 vector for the color scale limits; default
#'   \code{c(0, 1)}.
#' @param layer Integer. Band index when \code{en} has multiple layers;
#'   default 1L.
#' @param minor_div Integer. Minor grid divisions between major axis breaks;
#'   default 2L.
#'
#' @return A ggplot2 object. No files are written.
#'
#' @seealso \code{\link{save_suitability_plot}}, \code{\link{rENM_project_dir}}
#' @importFrom terra rast vect crs same.crs project ext nlyr as.data.frame
#' @importFrom sf st_as_sf st_crs
#' @importFrom ggplot2 ggplot geom_raster aes scale_fill_gradientn geom_sf labs
#'   coord_sf scale_x_continuous scale_y_continuous theme_bw theme element_text
#'   element_line
#' @importFrom scales squish
#'
#' @examples
#' \dontrun{
#'   p <- plot_suitability(en, "CASP")
#'   p <- plot_suitability(en, "CASP", title = "CASP 2000 Suitability")
#' }
#'
#' @export
plot_suitability <- function(en,
                             alpha_code,
                             state_shp  = NULL,
                             title      = "Climatic Suitability Map",
                             limits     = c(0, 1),
                             layer      = 1L,
                             minor_div  = 2L) {

  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required but not installed.", call. = FALSE)
  }

  project_dir <- rENM_project_dir()

  if (missing(alpha_code) || !is.character(alpha_code) || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character string (e.g., 'casp').")
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits))) {
    stop("`limits` must be a numeric length-2 vector of finite values, e.g., c(0, 1).")
  }
  layer     <- as.integer(layer)
  minor_div <- as.integer(minor_div)

  # ---- Resolve default state shapefile ----
  if (is.null(state_shp)) {
    state_shp <- file.path(project_dir, "data", "shapefiles",
                           "tl_2012_us_state", "tl_2012_us_state.shp")
  }

  # ---- Look up species info ----
  code <- toupper(alpha_code)
  si <- get_species_info(code)
  if (is.null(si)) stop("get_species_info() returned NULL for alpha_code: ", code)
  if (!("GAP.RANGE" %in% names(si))) stop("`GAP.RANGE` not found for: ", code)

  gap_range <- as.character(si[["GAP.RANGE"]])[1]
  if (!nzchar(gap_range) || is.na(gap_range)) stop("Missing GAP.RANGE for alpha_code: ", code)

  range_shp <- file.path(project_dir, "data", "shapefiles",
                          gap_range, paste0(gap_range, ".shp"))

  # ---- Coerce / load the raster ----
  if (inherits(en, "SpatRaster")) {
    r <- en
  } else if (is.character(en) && length(en) == 1L) {
    if (!file.exists(en)) stop("Raster file not found: ", en)
    r <- terra::rast(en)
  } else if (inherits(en, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    r <- terra::rast(en)
  } else {
    stop("`en` must be a terra::SpatRaster, a file path, or a raster::Raster*.")
  }

  if (terra::nlyr(r) < layer) {
    stop("Requested layer index ", layer, " but raster has only ", terra::nlyr(r), " layer(s).")
  }
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

  pal <- grDevices::colorRampPalette(c(
    "white", "lightblue", "blue", "green",
    "yellow", "orange", "red", "darkred", "black"
  ))(400)

  r_df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(r_df)[3] <- "value"

  e    <- terra::ext(r)
  xlim <- c(e$xmin, e$xmax)
  ylim <- c(e$ymin, e$ymax)

  pretty_with_minors <- function(lim, n_major = 6, minor_div = 2L) {
    majors <- pretty(lim, n = n_major)
    step   <- unique(round(diff(majors), 10))[1]
    if (is.na(step) || step == 0) {
      return(list(majors = majors, minors = NULL))
    }
    minor_step <- step / max(1L, minor_div + 1L)
    minors <- c(vapply(
      majors[-length(majors)],
      function(m0) m0 + seq_len(minor_div) * minor_step,
      numeric(minor_div)
    ))
    list(majors = majors, minors = minors)
  }

  bx <- pretty_with_minors(xlim, n_major = 6, minor_div = minor_div)
  by <- pretty_with_minors(ylim, n_major = 6, minor_div = minor_div)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data    = r_df,
      mapping = ggplot2::aes(x = x, y = y, fill = value)
    ) +
    ggplot2::scale_fill_gradientn(
      colours  = pal,
      limits   = limits,
      oob      = scales::squish,
      na.value = NA,
      name     = "Suitability"
    ) +
    ggplot2::geom_sf(data = stl_sf, fill = NA, color = "gray40", linewidth = 0.3) +
    ggplot2::geom_sf(data = spv_sf, fill = NA, color = "black",  linewidth = 0.7) +
    ggplot2::labs(title = title, x = "Longitude", y = "Latitude") +
    ggplot2::coord_sf(
      xlim   = xlim,
      ylim   = ylim,
      expand = FALSE,
      crs    = if (!is.na(sf::st_crs(stl_sf))) sf::st_crs(stl_sf) else NA
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
      panel.grid.major = ggplot2::element_line(),
      panel.grid.minor = ggplot2::element_line(),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold")
    )

  p
}
