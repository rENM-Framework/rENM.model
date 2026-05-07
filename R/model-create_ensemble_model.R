#' Build an ensemble ENM for a species and time slice
#'
#' Constructs an ensemble ecological niche model (ENM) for a given species and
#' 5-year time bin using the \pkg{sdm} framework, producing ensemble prediction
#' rasters, a binary range map, model statistics, and variable importance files.
#'
#' @details
#' Input paths follow the rENM run structure:
#' \preformatted{
#' <project_dir>/runs/<alpha_code>/TimeSeries/<year>/occs/of-<year>.csv
#' <project_dir>/runs/<alpha_code>/TimeSeries/<year>/vars/*.asc
#' <project_dir>/runs/<alpha_code>/TimeSeries/<year>/model/  (outputs)
#' }
#'
#' Models are trained via \code{sdm::sdm()} with replicated subsampling
#' (\code{reps}) and a held-out test partition (\code{tp} percent). Background
#' points are drawn with \code{method = "gRandom"}.
#'
#' The default ensemble uses \code{"maxnet"} (native R, no Java dependency).
#' A binary range map is derived by thresholding at the mean evaluation
#' threshold returned by \code{sdm::getEvaluation(opt = op)}.
#'
#' @param alpha_code Character. Four-letter banding code, case-insensitive
#'   (must match \code{^[A-Z]\{4\}$}; stored upper-case).
#' @param year Integer. A 5-year bin start in \{1980, ..., 2020\}
#'   divisible by 5.
#' @param methods Character. SDM algorithms to fit. Default
#'   \code{c("maxnet", "rf", "brt", "glm", "mars")}.
#' @param reps Integer. Number of replicated subsampling runs. Default 3.
#' @param bg Integer. Number of background points. Default 2500.
#' @param tp Numeric. Test partition percent (0--100). Default 40.
#' @param op Integer. Threshold optimization option for
#'   \code{sdm::getEvaluation}; \code{2} = max(sensitivity + specificity).
#'   Default 2.
#' @param ensemble_method Character. Ensemble combination method passed to
#'   \code{sdm::ensemble()}. Default \code{"unweighted"}.
#' @param io Character. Raster I/O backend for writing outputs; one of
#'   \code{"raster"} or \code{"terra"}. Default \code{"raster"}.
#' @param overwrite Logical. Overwrite existing output files. Default TRUE.
#' @param verbose Logical. Emit timestamped console progress. Default TRUE.
#'
#' @return Invisible list with elements:
#' \itemize{
#'   \item \code{data}: sdmData object
#'   \item \code{model}: sdmModel object
#'   \item \code{ensemble}: RasterLayer of predicted suitability
#'   \item \code{range}: RasterLayer of binary presence/absence
#'   \item \code{paths}: Named list of written file paths
#'   \item \code{stats}: List with \code{ms} (model stats), \code{pi_raw},
#'         \code{pi_ranked}
#' }
#' Side effects include writing raster outputs, statistics and variable
#' importance files, plots, and a log entry appended to the run log.
#'
#' @seealso \code{\link[sdm]{sdm}}, \code{\link[sdm]{sdmData}},
#'   \code{\link[sdm]{ensemble}}, \code{\link[sdm]{getEvaluation}}
#'
#' @importFrom utils read.csv write.table capture.output
#' @importFrom raster stack writeRaster nlayers values raster
#' @importFrom sdm sdm sdmData ensemble getVarImp getEvaluation
#' @importFrom sp coordinates proj4string CRS
#'
#' @examples
#' \dontrun{
#' result <- create_ensemble_model("CASP", 2000)
#'
#' result <- create_ensemble_model(
#'   alpha_code = "CASP", year = 2000,
#'   methods = c("maxnet", "glm", "brt", "rf", "mars"),
#'   reps = 5, bg = 5000, tp = 30, io = "terra"
#' )
#' }
#'
#' @export
create_ensemble_model <- function(
    alpha_code, year,
    methods = c("maxnet", "rf", "brt", "glm", "mars"),
    reps = 3,
    bg = 2500,
    tp = 40,
    op = 2,
    ensemble_method = "unweighted",
    io = c("raster", "terra"),
    overwrite = TRUE,
    verbose = TRUE
) {
  io <- match.arg(io)

  stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  say <- function(...) if (isTRUE(verbose)) {
    message(paste0("[", stamp(), "] ", paste(..., collapse = "")))
  }
  ensure_dir <- function(path) if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  # Validate required packages without library() calls
  req_pkgs <- c("sdm", "raster", "terra", "sp", "Hmisc", "dismo", "maxnet")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Required package(s) not installed: ",
         paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }

  if (!is.character(alpha_code) || length(alpha_code) != 1) {
    stop("alpha_code must be a single character value")
  }
  alpha_code <- toupper(trimws(alpha_code))
  if (!grepl("^[A-Z]{4}$", alpha_code)) {
    stop("alpha_code must be exactly four letters (A-Z)")
  }

  if (length(year) != 1 || is.na(suppressWarnings(as.integer(year)))) {
    stop("year must be a single integer")
  }
  year <- as.integer(year)
  if (year < 1980 || year > 2020 || (year %% 5) != 0) {
    stop("year must be in [1980, 2020] and divisible by 5 (e.g., 1980, 1985, ..., 2020)")
  }

  # ---- Resolve paths ----
  project_dir <- rENM_project_dir()
  base_dir  <- file.path(project_dir, "runs", alpha_code)
  ts_dir    <- file.path(base_dir, "TimeSeries", year)
  occs_fn   <- file.path(ts_dir, "occs", sprintf("of-%d.csv", year))
  vars_dir  <- file.path(ts_dir, "vars")
  model_dir <- file.path(ts_dir, "model")
  ensure_dir(model_dir)
  log_fn    <- file.path(base_dir, "_log.txt")

  say("Starting create_ensemble_model for ", alpha_code, " / ", year, " ...")

  if (!file.exists(occs_fn)) {
    stop(sprintf("Occurrences file not found: %s", occs_fn))
  }
  say("Loading occurrences: ", occs_fn)
  sp_df <- utils::read.csv(occs_fn, stringsAsFactors = FALSE)

  required_cols <- c("longitude", "latitude")
  missing_cols <- setdiff(required_cols, names(sp_df))
  if (length(missing_cols)) {
    stop(sprintf("Occurrence file missing columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  sp_df$species <- 1
  sp::coordinates(sp_df) <- ~ longitude + latitude
  sp::proj4string(sp_df) <- sp::CRS(NA_character_)

  if (!dir.exists(vars_dir)) {
    stop(sprintf("Predictor directory not found: %s", vars_dir))
  }
  asc_files <- list.files(vars_dir, pattern = "\\.asc$", full.names = TRUE)
  if (!length(asc_files)) {
    stop(sprintf("No .asc predictors found in %s", vars_dir))
  }
  say("Stacking predictors (", length(asc_files), ") from ", vars_dir)
  pr <- raster::stack(asc_files)

  say("Preparing sdmData() with bg=", bg, ", test.p=", tp)
  d <- sdm::sdmData(
    species ~ .,
    train      = sp_df,
    predictors = pr,
    bg         = list(n = bg, method = "gRandom")
  )
  utils::capture.output(print(d))

  say("Training sdm() with methods=", paste(methods, collapse = ","), ", reps=", reps)
  m <- sdm::sdm(
    species ~ .,
    d,
    methods     = methods,
    replication = "sub",
    test.p      = tp,
    n           = reps
  )
  utils::capture.output(print(m))

  pred_tif <- file.path(model_dir, sprintf("%s-%d-Prediction.tif", alpha_code, year))
  pred_asc <- file.path(model_dir, sprintf("%s-%d-Prediction.asc", alpha_code, year))

  say("Creating ensemble prediction (method=", ensemble_method, ", io=", io, ") ...")
  en <- sdm::ensemble(
    m, pr,
    filename  = if (io == "raster") pred_tif else "",
    setting   = list(method = ensemble_method),
    overwrite = overwrite
  )
  utils::capture.output(print(en))

  if (io == "raster") {
    raster::writeRaster(en, filename = pred_asc, format = "ascii", overwrite = overwrite)
  } else {
    en_tr <- terra::rast(en)
    terra::writeRaster(en_tr, pred_tif, overwrite = overwrite)
    terra::writeRaster(en_tr, pred_asc, overwrite = overwrite)
  }
  say("Wrote ensemble rasters:\n  ", pred_tif, "\n  ", pred_asc)

  say("Computing model statistics and variable importance ...")
  ms  <- utils::capture.output(print(m))
  pi0 <- utils::capture.output(print(sdm::getVarImp(m)))

  pi1 <- rank_variable_importance(pi0, combine_metrics = TRUE)

  stats_txt   <- file.path(model_dir, sprintf("%s-%d-Stats.txt",      alpha_code, year))
  pi_graph_fn <- file.path(model_dir, sprintf("%s-%d-PI-Graphed.txt", alpha_code, year))
  pi_rank_fn  <- file.path(model_dir, sprintf("%s-%d-PI-Ranked.txt",  alpha_code, year))

  writeLines(ms,  con = stats_txt)
  writeLines(pi0, con = pi_graph_fn)
  utils::write.table(pi1, file = pi_rank_fn, row.names = FALSE, quote = FALSE, sep = "\t")
  say("Wrote stats:\n  ", stats_txt, "\n  ", pi_graph_fn, "\n  ", pi_rank_fn)

  say("Plotting climatic suitability map ...")
  title1 <- sprintf("%s Climatic Suitability (%d)", alpha_code, year)
  p1     <- plot_suitability(en, alpha_code, title = title1)
  utils::capture.output(print(p1))

  save_prefix <- file.path(model_dir, sprintf("%s-%d-Prediction", alpha_code, year))
  save_suitability_plot(p1, save_prefix)

  say("Deriving binary range map using threshold op=", op, " ...")
  ev <- sdm::getEvaluation(m, stat = c("threshold"), opt = op)
  pa <- raster::raster(en)
  raster::values(pa) <- ifelse(
    raster::values(en) >= mean(ev$threshold, na.rm = TRUE),
    1, 0
  )
  title2 <- sprintf("%s Range (%d)", alpha_code, year)
  p2     <- create_range_map(pa, alpha_code, year = year, title = title2)
  utils::capture.output(print(p2))

  n_pres <- tryCatch(nrow(as.data.frame(sp_df)), error = function(e) NA_integer_)
  n_pred <- tryCatch(raster::nlayers(pr),        error = function(e) NA_integer_)

  sep_line <- strrep("-", 72)
  kv <- function(label, value, width = 20) sprintf("%-*s %s", width, paste0(label, ":"), value)
  per_year_line <- sprintf(
    "  %d  occs=%s   preds=%s   methods=%s   reps=%d   ensemble=tif+asc",
    year,
    ifelse(is.na(n_pres), "-", format(n_pres, big.mark = ",")),
    ifelse(is.na(n_pred), "-", format(n_pred, big.mark = ",")),
    paste(methods, collapse = "/"),
    reps
  )
  legacy_lines <- c(
    sep_line,
    "Processing summary (create_ensemble_model)",
    kv("Timestamp",      stamp()),
    kv("Alpha code",     alpha_code),
    kv("Year processed", year),
    kv("Occurrences",    ifelse(is.na(n_pres), "-", format(n_pres, big.mark = ","))),
    kv("Predictors",     ifelse(is.na(n_pred), "-", format(n_pred, big.mark = ","))),
    kv("Methods",        paste(methods, collapse = ", ")),
    kv("Replicates",     reps),
    kv("BG points",      bg),
    kv("Test percent",   tp),
    kv("Threshold op",   sprintf("%d (max(se+sp))", op)),
    "Outputs:",
    paste0("  ", pred_tif),
    paste0("  ", pred_asc),
    paste0("  ", stats_txt),
    paste0("  ", pi_graph_fn),
    paste0("  ", pi_rank_fn),
    "Per-year status:",
    per_year_line,
    ""
  )
  cat(paste0(paste(legacy_lines, collapse = "\n"), "\n"), file = log_fn, append = TRUE)
  say("Run complete. Summary appended to ", log_fn)

  invisible(list(
    data     = d,
    model    = m,
    ensemble = en,
    range    = pa,
    paths    = list(
      base_dir    = base_dir,
      ts_dir      = ts_dir,
      model_dir   = model_dir,
      pred_tif    = pred_tif,
      pred_asc    = pred_asc,
      stats_txt   = stats_txt,
      pi_graph_fn = pi_graph_fn,
      pi_rank_fn  = pi_rank_fn,
      log         = log_fn
    ),
    stats    = list(ms = ms, pi_raw = pi0, pi_ranked = pi1)
  ))
}
