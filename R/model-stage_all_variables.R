#' Stage environmental variables into TimeSeries bins
#'
#' Copies raster variables from run-level storage into TimeSeries bins,
#' organizing predictors for downstream modeling while reporting progress
#' and recording a processing summary in the run log.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' Copies raster variables ('.tif', '.asc') from the run-level '_vars'
#' directory into the corresponding TimeSeries directories:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_vars/<year>/*.(tif|asc)}
#'
#' \strong{Outputs}:
#' Files are staged into:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{TimeSeries/<year>/vars/*}
#'
#' Existing files are overwritten.
#'
#' \strong{Methods}:
#' \itemize{
#'   \item Recursively scans source directories for '.tif' and '.asc' files
#'   \item Copies files into the target directory (flat structure)
#'   \item Overwrites existing files
#'   \item Warns on missing source directories
#'   \item Warns on duplicate basenames (last file wins)
#' }
#'
#' \strong{Logging}:
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' including:
#' \itemize{
#'   \item Timestamp
#'   \item Alpha code
#'   \item Number of files transferred
#'   \item File list
#'   \item Total elapsed time
#'   \item Output location
#' }
#'
#' @param alpha_code Character. Species alpha code (e.g., \code{"CASP"}).
#' @param year Integer. Vector of 5-year bins. Default:
#'   \code{seq(1980, 2020, by = 5)}.
#' @param project_dir Character, NULL. Optional project root directory.
#'   If NULL, resolved using \code{rENM_project_dir()}.
#'
#' @return Data frame (invisibly returned) with columns:
#' \itemize{
#'   \item \code{year}
#'   \item \code{src_dir}
#'   \item \code{dest_dir}
#'   \item \code{files_found}
#'   \item \code{files_copied}
#'   \item \code{duplicates}
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Copies raster files into TimeSeries directories
#'   \item Overwrites existing files in target directories
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @examples
#' \dontrun{
#'   stage_all_variables("CASP")
#'   stage_all_variables("CASP", year = 2005)
#'   stage_all_variables("CASP", year = c(1990, 2010, 2020))
#'   stage_all_variables("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @export
stage_all_variables <- function(alpha_code,
                                year = seq(1980, 2020, by = 5),
                                project_dir = NULL) {

  ## ---- validate inputs ------------------------------------------------------
  start_time <- Sys.time()

  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }

  if (!is.numeric(year) || any(!is.finite(year)) || any(year %% 1 != 0)) {
    stop("`year` must be an integer vector (e.g., seq(1980, 2020, by = 5)).", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))
  years <- sort(unique(as.integer(year)))

  ## ---- resolve project directory --------------------------------------------
  project_root <- rENM_project_dir(project_dir)

  runs_path <- function(...) file.path(project_root, "runs", alpha_code, ...)
  log_file_path <- runs_path("_log.txt")

  ## ---- helpers --------------------------------------------------------------
  fmt_elapsed <- function(dt) {
    secs <- as.numeric(dt, units = "secs")
    hrs  <- floor(secs / 3600); secs <- secs - 3600 * hrs
    mins <- floor(secs / 60);   secs <- round(secs - 60 * mins)
    sprintf("%d:%02d:%02d", hrs, mins, secs)
  }

  append_log_block <- function(logfile, alpha_code, outputs_saved,
                               output_file, elapsed_seconds, files_list) {

    sep_line <- paste(rep("-", 72), collapse = "")
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    elapsed   <- fmt_elapsed(as.difftime(elapsed_seconds, units = "secs"))

    files_header <- sprintf("%-34s %d", "Files transferred:", length(files_list))
    files_lines  <- if (length(files_list)) {
      paste0("\n", paste0("  - ", files_list, collapse = "\n"))
    } else ""

    block <- paste0(
      sep_line, "\n",
      "Processing summary (stage_all_variables)\n",
      sprintf("%-34s %s\n", "Timestamp:", timestamp),
      sprintf("%-34s %s\n", "Alpha code:", alpha_code),
      sprintf("%-34s %s\n", "Raster source:", "variable rasters (*.tif, *.asc)"),
      files_header, files_lines, "\n",
      sprintf("%-34s %s\n", "Outputs saved:", outputs_saved),
      sprintf("%-34s %s\n", "Total elapsed:", elapsed),
      sprintf("%-34s %s\n", "Output file:", ifelse(is.na(output_file) || !nzchar(output_file), "-", output_file))
    )

    con <- file(logfile, open = "a", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(block, con = con)
  }

  ## ---- processing -----------------------------------------------------------
  message(sprintf("-> Staging variables for '%s' into TimeSeries/<year>/vars ...", alpha_code))

  results <- data.frame(
    year = integer(0),
    src_dir = character(0),
    dest_dir = character(0),
    files_found = integer(0),
    files_copied = integer(0),
    duplicates = character(0),
    stringsAsFactors = FALSE
  )

  total_copied     <- 0L
  last_target      <- NA_character_
  transferred_list <- character(0)

  for (yy in years) {

    src  <- runs_path("_vars", yy)
    dest <- runs_path("TimeSeries", yy, "vars")

    message(sprintf("  * Year %d", yy))
    message(sprintf("    - Source: %s", src))
    message(sprintf("    - Target: %s", dest))

    if (!dir.exists(src)) {
      warning(sprintf("    - Skipped: source directory not found: %s", src))
      results <- rbind(results, data.frame(
        year = yy, src_dir = src, dest_dir = dest,
        files_found = 0L, files_copied = 0L, duplicates = "",
        stringsAsFactors = FALSE
      ))
      next
    }

    if (!dir.exists(dest)) {
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      message("    - Created target directory.")
    }

    files <- list.files(
      path = src,
      pattern = "\\.(tif|asc)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    n_found <- length(files)

    if (n_found == 0L) {
      message("    - No raster files found.")
      results <- rbind(results, data.frame(
        year = yy, src_dir = src, dest_dir = dest,
        files_found = 0L, files_copied = 0L, duplicates = "",
        stringsAsFactors = FALSE
      ))
      next
    }

    base_tbl <- sort(table(basename(files)), decreasing = TRUE)
    dups <- names(base_tbl[base_tbl > 1L])
    dup_msg <- if (length(dups)) paste(dups, collapse = ", ") else ""

    if (length(dups)) {
      warning(sprintf("    - Duplicate basenames detected (last one wins): %s", dup_msg))
    }

    copied_vec <- vapply(files, function(f) {
      to_path <- file.path(dest, basename(f))
      ok <- file.copy(f, to_path, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
      if (ok) transferred_list <<- c(transferred_list, to_path)
      ok
    }, logical(1))

    n_copied <- sum(copied_vec, na.rm = TRUE)
    total_copied <- total_copied + n_copied
    last_target  <- dest

    message(sprintf("    - Copied %d / %d files.", n_copied, n_found))

    results <- rbind(results, data.frame(
      year = yy, src_dir = src, dest_dir = dest,
      files_found = n_found, files_copied = n_copied, duplicates = dup_msg,
      stringsAsFactors = FALSE
    ))
  }

  ## ---- logging --------------------------------------------------------------
  elapsed <- difftime(Sys.time(), start_time, units = "secs")

  append_log_block(
    logfile         = log_file_path,
    alpha_code      = alpha_code,
    outputs_saved   = total_copied,
    output_file     = ifelse(is.na(last_target), "-", last_target),
    elapsed_seconds = as.numeric(elapsed, units = "secs"),
    files_list      = transferred_list
  )

  message(sprintf("Done. Total files copied: %d", total_copied))
  message(sprintf("Log updated: %s", log_file_path))

  invisible(results)
}
