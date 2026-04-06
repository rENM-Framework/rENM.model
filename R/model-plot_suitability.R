#' Plot a climatic suitability raster with state boundaries and species range
#'
#' Build a ggplot2 map of continuous climatic suitability values from a
#' raster and overlay state or administrative boundaries plus a species
#' range polygon. If the raster has multiple layers, layer selects which
#' to plot.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' Build a ggplot2 map of continuous climatic suitability values from a
#' raster (e.g., 0-1) and overlay state or administrative boundaries
#' plus a species range polygon. If the raster has multiple layers,
#' \code{layer} selects which to plot.
#'
#' The species range shapefile is looked up via an internal helper
#' \code{get_species_info(alpha_code)}, which must return a component
#' \code{GAP.RANGE}. The shapefile is then read from:
#' \code{<project_dir>/data/shapefiles/<gap_range>/}
#' \code{<gap_range>.shp}
#'
#' \strong{Inputs}
#' \itemize{
#'   \item A suitability raster provided as a SpatRaster, file path, or
#'   raster object.
#'   \item Species code used to look up GAP.RANGE via
#'   \code{get_species_info(alpha_code)}.
#'   \item Path to a states or administrative boundaries shapefile.
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item A ggplot object representing the climatic suitability surface
#'   with state boundaries and species range overlays.
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item The state and species-range vectors are reprojected to the
#'   raster CRS (if set) to ensure alignment.
#'   \item The color scale uses a smooth multi-hue palette via
#'   \code{grDevices::colorRampPalette()}.
#'   \item NA cell values render transparent.
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Species range shapefile located at:
#'   \code{<project_dir>/data/shapefiles/<gap_range>/}
#'   \code{<gap_range>.shp}
#'   \item Default state shapefile located at:
#'   \code{<project_dir>/data/shapefiles/tl_2012_us_state/}
#'   \code{tl_2012_us_state.shp}
#' }
#'
#' @param en SpatRaster, Character, Raster. A suitability raster. One of:
#'   a terra SpatRaster, a single-file path (character) to a raster
#'   readable by terra rast(), or a raster::Raster* object converted
#'   internally to SpatRaster.
#' @param alpha_code Character. Species code used to look up GAP.RANGE via
#'   get_species_info(alpha_code) and, from that, the species range
#'   shapefile under
#'   \code{<project_dir>/data/shapefiles/<gap_range>/}
#'   \code{<gap_range>.shp}
#' @param state_shp Character. Path to a states or administrative
#'   boundaries shapefile (.shp). Defaults to
#'   \code{<project_dir>/data/shapefiles/tl_2012_us_state/}
#'   \code{tl_2012_us_state.shp}
#' @param title Character. Plot title. Default "Climatic Suitability Map".
#' @param limits Numeric. Length-2 vector giving the lower and upper
#'   limits of the color scale (default c(0, 1)).
#' @param layer Integer. Band index to plot when en has multiple layers
#'   (default 1L).
#' @param minor_div Integer. Number of minor grid divisions between major
#'   axis breaks (default 2L).
#'
#' @return
#' A ggplot object.
#' Side effects:
#' \itemize{
#'   \item None. The function returns a plot object without writing files.
#' }
#'
#' @importFrom terra rast vect crs same.crs project ext nlyr as.data.frame
#' @importFrom sf st_as_sf st_crs
#' @importFrom ggplot2 ggplot geom_raster aes scale_fill_gradientn geom_sf labs
#' @importFrom ggplot2 coord_sf scale_x_continuous scale_y_continuous theme_bw
#' @importFrom ggplot2 theme element_text element_line
#' @importFrom scales squish
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' plot_suitability(en, "CASP")
#' }
#'
#' @export
plot_suitability <- function(en,
                             alpha_code,
                             state_shp  = "~/rENM/data/shapefiles/tl_2012_us_state/tl_2012_us_state.shp",
                             title      = "Climatic Suitability Map",
                             limits     = c(0, 1),
                             layer      = 1L,
                             minor_div  = 2L) {
  # ---- Dependency checks ----
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required but not installed.", call. = FALSE)
  }
  
  if (!exists("rENM_project_dir", mode = "function")) {
    stop("Function `rENM_project_dir()` must be available (e.g., exported by this package).")
  }
  project_dir <- rENM_project_dir()
  
  # ---- Argument validation ----
  if (missing(alpha_code) || !is.character(alpha_code) || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character string (e.g., 'casp').")
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits))) {
    stop("`limits` must be a numeric length-2 vector of finite values, e.g., c(0, 1).")
  }
  layer     <- as.integer(layer)
  minor_div <- as.integer(minor_div)
  
  # ---- Look up species info via internal helper ----
  code <- toupper(alpha_code)
  
  if (!exists("get_species_info", mode = "function")) {
    stop(
      "Internal helper `get_species_info()` not found in the package namespace.\n",
      "Ensure it is defined in rENM.core and that the package is loaded."
    )
  }
  
  si <- get_species_info(code)
  if (is.null(si)) {
    stop("get_species_info() returned NULL for alpha_code: ", code)
  }
  if (!("GAP.RANGE" %in% names(si))) {
    stop("`GAP.RANGE` not found in result of get_species_info() for: ", code)
  }
  
  gap_range <- as.character(si[["GAP.RANGE"]])[1]
  if (!nzchar(gap_range) || is.na(gap_range)) {
    stop("Missing GAP.RANGE for alpha_code: ", code)
  }
  
  # Construct expected shapefile path from GAP.RANGE
  shp_root_ex <- file.path(project_dir, "data", "shapefiles")
  range_dir   <- file.path(shp_root_ex, gap_range)
  range_shp   <- file.path(range_dir, paste0(gap_range, ".shp"))
  
  # ---- Coerce / load the raster ----
  if (inherits(en, "SpatRaster")) {
    r <- en
  } else if (is.character(en) && length(en) == 1L) {
    if (!file.exists(en)) stop("Raster file not found: ", en)
    r <- terra::rast(en)
  } else if (inherits(en, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    r <- terra::rast(en)
  } else {
    stop("`en` must be a terra::SpatRaster, a file path to a raster, or a raster::Raster*.")
  }
  
  # ---- Select layer if multi-band ----
  if (terra::nlyr(r) < layer) {
    stop("Requested layer index ", layer, " but raster has only ", terra::nlyr(r), " layer(s).")
  }
  if (terra::nlyr(r) > 1L) r <- r[[layer]]
  
  # ---- Read vectors & project to raster CRS if needed ----
  if (identical(state_shp, "~/rENM/data/shapefiles/tl_2012_us_state/tl_2012_us_state.shp")) {
    state_shp_ex <- file.path(project_dir, "data", "shapefiles", "tl_2012_us_state", "tl_2012_us_state.shp")
  } else {
    state_shp_ex <- state_shp
  }
  
  if (!file.exists(state_shp_ex)) {
    stop("State shapefile not found: ", state_shp_ex)
  }
  if (!file.exists(range_shp)) {
    stop(
      "Species range shapefile not found at: ", range_shp,
      "\n(derived from GAP.RANGE = '", gap_range, "')"
    )
  }
  
  stl <- terra::vect(state_shp_ex)
  spv <- terra::vect(range_shp)
  if (!is.na(terra::crs(r))) {
    if (!terra::same.crs(stl, r)) stl <- terra::project(stl, r)
    if (!terra::same.crs(spv, r)) spv <- terra::project(spv, r)
  }
  
  stl_sf <- sf::st_as_sf(stl)
  spv_sf <- sf::st_as_sf(spv)
  
  # ---- Palette ----
  pal <- grDevices::colorRampPalette(c(
    "white", "lightblue", "blue", "green",
    "yellow", "orange", "red", "darkred", "black"
  ))(400)
  
  # ---- Raster -> data.frame for ggplot ----
  r_df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(r_df)[3] <- "value"
  
  # ---- Axis breaks with minors ----
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
  
  # ---- Plot ----
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