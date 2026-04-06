#' Build an ensemble ENM and derived products for a species
#'
#' Constructs an ensemble ecological niche model (ENM) for a given
#' species (four-letter banding code) and a 5-year time slice using
#' the sdm framework, producing ensemble predictions and derived
#' products including range maps and model diagnostics.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}:
#' Inputs/paths must follow the rENM run structure:
#' \preformatted{
#' <rENM_project_dir()>/runs/<ALPHA>/TimeSeries/<YEAR>/occs/of-<YEAR>.csv
#' <rENM_project_dir()>/runs/<ALPHA>/TimeSeries/<YEAR>/vars/*.asc
#' <rENM_project_dir()>/runs/<ALPHA>/TimeSeries/<YEAR>/model/   (outputs written here)
#' }
#'
#' \strong{Inputs}:
#' The occurrence CSV must include \code{longitude}, \code{latitude}.
#'
#' \strong{Methods}:
#' Models are trained via \code{sdm::sdm()} with replicated subsampling
#' (\code{reps}) and a test partition (\code{tp}, percent). Background
#' points (\code{bg}) are drawn using \code{method = "gRandom"}.
#'
#' The default modeling ensemble uses \code{"maxnet"} (a native R
#' implementation of maximum entropy modeling) instead of
#' \code{"maxent"} to avoid Java dependencies and improve
#' reproducibility in CRAN-compliant environments.
#'
#' \strong{Outputs}:
#' \itemize{
#'   \item Ensemble prediction rasters (GeoTIFF and ASCII)
#'   \item Model statistics and variable importance summaries
#'   \item Binary range map derived via threshold optimization
#' }
#'
#' \strong{Logging}:
#' A human-readable, append-only summary is written to:
#' \code{<rENM_project_dir()>/runs/<alpha_code>/}
#' \code{_log.txt}
#' with key parameters, file outputs, and per-year status. Console
#' messages include timestamps for traceability.
#'
#' @param alpha_code Character. Four-letter banding code,
#' case-insensitive (validated as exactly four A-Z letters;
#' stored upper-case).
#' @param year Integer. A 5-year bin start in (1980, 2020)
#' divisible by 5.
#' @param methods Character. SDM algorithms to fit. Defaults to
#' \code{c("maxnet","rf","brt","glm","mars")}.
#' @param reps Integer. Number of replicated subsampling runs.
#' Default \code{3}.
#' @param bg Integer. Number of background points. Default
#' \code{2500}.
#' @param tp Numeric. Test partition percent (0-100). Default
#' \code{40}.
#' @param op Integer. Threshold optimization option passed to
#' \code{sdm::getEvaluation} for binarization; \code{2} =
#' max(se+sp). Default \code{2}.
#' @param ensemble_method Character. Ensemble combination method
#' passed to \code{sdm::ensemble} in
#' \code{setting = list(method = ...)}. Default
#' \code{"unweighted"}.
#' @param io Character. Raster I/O backend for writing outputs;
#' one of \code{"raster"} or \code{"terra"}. Predictors are read
#' via raster either way to maintain sdm compatibility. Default
#' \code{"raster"}.
#' @param overwrite Logical. Overwrite existing output files.
#' Default \code{TRUE}.
#' @param verbose Logical. Emit timestamped console progress.
#' Default \code{TRUE}.
#'
#' @return Invisibly returns a List with structured elements:
#' \itemize{
#'   \item \code{data}: sdmData object containing training data
#'   \item \code{model}: sdmModel object representing fitted models
#'   \item \code{ensemble}: RasterLayer of predicted suitability
#'   \item \code{range}: RasterLayer of binary presence/absence
#'   \item \code{paths}: List of important written file paths
#'   \item \code{stats}: List of captured model statistics and
#'   variable importance outputs
#' }
#' Side effects include writing raster outputs, statistics files,
#' variable importance summaries, plots, and appending a log entry
#' to the run log file.
#'
#' @importFrom utils read.csv write.table capture.output
#' @importFrom raster stack writeRaster nlayers values raster
#' @importFrom sdm sdm sdmData ensemble getVarImp getEvaluation
#' @importFrom sp coordinates proj4string CRS
#'
#' @examples
#' \dontrun{
#' # Default knobs; raster I/O
#' result <- create_ensemble_model("CASP", 2000)
#'
#' # Heavier replication, custom methods, terra I/O
#' result <- create_ensemble_model(
#'   alpha_code = "CASP", year = 2000,
#'   methods = c("maxnet","glm","brt","rf","mars"),
#'   reps = 5, bg = 5000, tp = 30, op = 2,
#'   ensemble_method = "unweighted",
#'   io = "terra", overwrite = TRUE, verbose = TRUE
#' )
#' }
#'
#' @seealso \code{\link[sdm]{sdm}}, \code{\link[sdm]{sdmData}},
#' \code{\link[sdm]{ensemble}}, \code{\link[sdm]{getEvaluation}}
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
  must_have_pkg <- function(p) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Required package '%s' is not installed.", p))
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
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

  invisible(lapply(c("sdm", "raster", "terra", "sp", "Hmisc", "dismo", "maxnet"), must_have_pkg))

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

  rvi_fun <- get("rank_variable_importance", mode = "function")
  pi1     <- rvi_fun(pi0, combine_metrics = TRUE)

  stats_txt   <- file.path(model_dir, sprintf("%s-%d-Stats.txt",      alpha_code, year))
  pi_graph_fn <- file.path(model_dir, sprintf("%s-%d-PI-Graphed.txt", alpha_code, year))
  pi_rank_fn  <- file.path(model_dir, sprintf("%s-%d-PI-Ranked.txt",  alpha_code, year))

  writeLines(ms,  con = stats_txt)
  writeLines(pi0, con = pi_graph_fn)
  utils::write.table(pi1, file = pi_rank_fn, row.names = FALSE, quote = FALSE, sep = "\t")
  say("Wrote stats:\n  ", stats_txt, "\n  ", pi_graph_fn, "\n  ", pi_rank_fn)

  say("Plotting climatic suitability map ...")
  title1 <- sprintf("%s Climatic Suitability (%d)", alpha_code, year)

  ps_fun <- get("plot_suitability", mode = "function")
  p1     <- ps_fun(en, alpha_code, title = title1)
  utils::capture.output(print(p1))

  ssp_fun     <- get("save_suitability_plot", mode = "function")
  save_prefix <- file.path(model_dir, sprintf("%s-%d-Prediction", alpha_code, year))
  ssp_fun(p1, save_prefix)

  say("Deriving binary range map using threshold op=", op, " ...")
  ev <- sdm::getEvaluation(m, stat = c("threshold"), opt = op)
  pa <- raster::raster(en)
  raster::values(pa) <- ifelse(
    raster::values(en) >= mean(ev$threshold, na.rm = TRUE),
    1, 0
  )
  title2 <- sprintf("%s Range (%d)", alpha_code, year)

  crm_fun <- get("create_range_map", mode = "function")
  p2      <- crm_fun(pa, alpha_code, year = year, title = title2)
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
