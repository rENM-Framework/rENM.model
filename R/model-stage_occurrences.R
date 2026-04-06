#' Stage occurrence CSVs into time series bins
#'
#' Copies per-bin occurrence files into the TimeSeries structure for
#' downstream modeling while ensuring consistency with preprocessing
#' outputs and recording a detailed audit log.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' Copies per-bin occurrence files (\code{of-<year>.csv}) from:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/of-<year>.csv}
#'
#' \strong{Outputs}:
#' Files are staged into:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{TimeSeries/<year>/occs/of-<year>.csv}
#'
#' Existing files are overwritten to ensure consistency with the latest
#' preprocessing outputs.
#'
#' \strong{Methods}:
#' \itemize{
#'   \item Validates requested years (must be (1980, 2020), divisible by 5)
#'   \item Creates destination directories as needed
#'   \item Copies files using \code{file.copy(..., overwrite = TRUE)}
#'   \item Records per-year status: Copied, Missing, or Error
#'   \item Aggregates results into a structured log block
#' }
#'
#' \strong{Logging}:
#' A standardized summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' including:
#' \itemize{
#'   \item Timestamp and alpha code
#'   \item Per-year staging table
#'   \item Counts (total, copied, missing)
#'   \item Output paths (multi-line)
#'   \item Total elapsed time
#' }
#'
#' Each run is separated by a 72-character divider for readability.
#'
#' @param alpha_code Character. Species alpha code (e.g., \code{"CASP"}).
#' @param year Integer. Vector of 5-year bins. Default:
#'   \code{seq(1980, 2020, 5)}.
#' @param project_dir Character, NULL. Optional project root directory.
#'   If NULL, resolved using \code{rENM_project_dir()}.
#'
#' @return List (invisibly returned) with elements:
#' \itemize{
#'   \item \code{alpha_code}: Character species alpha code
#'   \item \code{years}: Integer vector of processed years
#'   \item \code{copied_paths}: Character vector of copied file paths
#'   \item \code{counts}: List with elements total_years, valid_files,
#'         copied_files, missing_files
#'   \item \code{elapsed_seconds}: Numeric elapsed time in seconds
#'   \item \code{log_file}: Character path to the run log file
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Copies occurrence CSV files into TimeSeries directories
#'   \item Overwrites existing files in destination directories
#'   \item Appends a processing summary to the run log
#'   \item Writes audit output to the console
#' }
#'
#' @examples
#' \dontrun{
#'   stage_occurrences("CASP")
#'   stage_occurrences("CASP", year = 1990)
#'   stage_occurrences("CASP", year = c(1980, 2000, 2015))
#' }
#'
#' @export
stage_occurrences <- function(alpha_code,
                              year = seq(1980, 2020, 5),
                              project_dir = NULL) {

  ## ---- timing ---------------------------------------------------------------
  t_start <- Sys.time()

  ## ---- validate inputs ------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }
  alpha_code <- toupper(trimws(alpha_code))

  valid_steps <- seq(1980L, 2020L, by = 5L)
  if (!is.numeric(year) || any(!year %in% valid_steps)) {
    stop("`year` must be values in \\{1980, 1985, ..., 2020\\}.", call. = FALSE)
  }
  years <- sort(unique(as.integer(year)))

  ## ---- resolve paths --------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  base_dir <- file.path(project_root, "runs", alpha_code)
  occs_dir <- file.path(base_dir, "_occs")
  log_file <- file.path(base_dir, "_log.txt")

  ## ---- helpers --------------------------------------------------------------
  sep72 <- paste(rep("-", 72), collapse = "")
  now   <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  ## ---- counters -------------------------------------------------------------
  total_years   <- length(years)
  valid_files   <- 0L
  copied_files  <- 0L
  missing_files <- 0L
  copied_paths  <- character(0)
  per_year_lines <- character(0)

  ## ---- processing -----------------------------------------------------------
  for (yr in years) {

    src  <- file.path(occs_dir, sprintf("of-%d.csv", yr))
    dest_dir <- file.path(base_dir, "TimeSeries", yr, "occs")
    dest <- file.path(dest_dir, sprintf("of-%d.csv", yr))

    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (!file.exists(src)) {
      warning(sprintf("Missing source file: %s", src))
      missing_files <- missing_files + 1L
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6d %-10s %s", yr, "Missing", "-"))
      next
    }

    ok <- tryCatch(file.copy(src, dest, overwrite = TRUE),
                   error = function(e) FALSE)

    if (isTRUE(ok)) {
      valid_files  <- valid_files + 1L
      copied_files <- copied_files + 1L
      copied_paths <- c(copied_paths, dest)

      per_year_lines <- c(per_year_lines,
                          sprintf("%-6d %-10s %s", yr, "Copied", dest))
    } else {
      warning(sprintf("Copy failed: %s", src))
      missing_files <- missing_files + 1L
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6d %-10s %s", yr, "Error", "-"))
    }
  }

  ## ---- timing ---------------------------------------------------------------
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  timestamp <- now()

  ## ---- log formatting -------------------------------------------------------
  yearly_block <- c(
    "",
    sep72,
    "Yearly staging details (stage_occurrences)",
    sprintf("%-6s %-10s %s", "Year", "Status", "Destination"),
    per_year_lines,
    sep72
  )

  outputs_block <- if (length(copied_paths)) {
    c("Outputs saved:", paste0("  ", copied_paths))
  } else {
    "Outputs saved: none"
  }

  summary_block <- c(
    "",
    sep72,
    "Processing summary (stage_occurrences)",
    sprintf("%-20s %s", "Timestamp:", timestamp),
    sprintf("%-20s %s", "Alpha code:", alpha_code),
    sprintf("%-20s %d", "Total bins:", total_years),
    sprintf("%-20s %d", "Valid files:", valid_files),
    sprintf("%-20s %d", "Copied files:", copied_files),
    sprintf("%-20s %d", "Missing files:", missing_files),
    outputs_block,
    sprintf("%-20s %.2f sec", "Total elapsed:", elapsed),
    sprintf("%-20s %s", "Log file:", log_file)
  )

  ## ---- write log ------------------------------------------------------------
  cat(paste(c(yearly_block, summary_block, ""), collapse = "\n"),
      file = log_file, append = TRUE)

  ## ---- console echo ---------------------------------------------------------
  cat(paste(c(yearly_block, summary_block, ""), collapse = "\n"))

  invisible(list(
    alpha_code      = alpha_code,
    years           = years,
    copied_paths    = copied_paths,
    counts          = list(
      total_years   = total_years,
      valid_files   = valid_files,
      copied_files  = copied_files,
      missing_files = missing_files
    ),
    elapsed_seconds = elapsed,
    log_file        = log_file
  ))
}
