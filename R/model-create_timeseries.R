#' Build an ensemble rENM time series for a species
#'
#' Schedules and executes create_ensemble_model() for each 5-year time
#' bin between 1980 and 2020 in parallel, producing a full rENM time
#' series for a given species.
#'
#' @details
#' This function is part of the rENM.model framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' A four-letter species code is required. The project directory is
#' resolved either from the provided argument or via
#' \code{rENM_project_dir()}.
#'
#' \strong{Processing steps}:
#' \itemize{
#'   \item Validates \code{alpha_code} and prepares the run directory
#'   \item Resolves the base project directory
#'   \item Spawns a PSOCK cluster with size:
#'         \code{min(number of years, detectCores())}
#'   \item Ensures each worker can call \code{create_ensemble_model()}
#'         within the rENM.model namespace
#'   \item Launches one job per 5-year bin using
#'         \code{parallel::parLapplyLB()} for load-balanced execution
#'   \item Collates per-year results including status, elapsed time,
#'         and output summaries
#'   \item Appends a formatted summary to the run log
#' }
#'
#' \strong{Worker function availability}:
#' Workers in a PSOCK cluster are fresh R sessions. The function
#' \code{create_ensemble_model()} is explicitly referenced from the
#' rENM.model namespace and made available to workers via closure.
#'
#' This supports both development (devtools::load_all()) and installed
#' package workflows.
#'
#' \strong{Fault tolerance}:
#' Each year is processed independently. Errors for individual years are
#' captured and logged without aborting the full time series. The cluster
#' is always stopped on exit.
#'
#' \strong{Outputs}:
#' A processing summary is appended to:
#'
#' \code{file.path(project_dir, "runs", <alpha_code>, "_log.txt")}
#'
#' @param alpha_code Character. Four-letter banding code, e.g., "CASP".
#' @param project_dir Character, NULL. Base directory for outputs and logs.
#'   If NULL, the directory is resolved using \code{rENM_project_dir()}.
#'
#' @return Invisibly returns a named List keyed by year. Each element is:
#' \itemize{
#'   \item ok: Logical indicating successful completion
#'   \item year: Integer representing the 5-year bin
#'   \item elapsed: Numeric elapsed seconds for the task
#'   \item result: Output from \code{create_ensemble_model()} or NULL
#'   \item notes: Character containing "ok" or an error message
#' }
#'
#' Attributes:
#' \itemize{
#'   \item total_elapsed_sec: Numeric total wall-clock time in seconds
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Executes parallel jobs using a PSOCK cluster
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom parallel detectCores makeCluster stopCluster clusterSetRNGStream parLapplyLB
#' @importFrom rENM.model create_ensemble_model
#'
#' @examples
#' \dontrun{
#'   create_timeseries("CASP")
#'   create_timeseries("CASP", project_dir = "/my/project")
#' }
#'
#' @seealso
#' [create_ensemble_model()], \pkg{parallel}
#'
#' @export
create_timeseries <- function(alpha_code, project_dir = NULL) {

  # -----------------------------
  # Project directory resolution
  # -----------------------------
  if (is.null(project_dir)) {
    project_dir <- rENM_project_dir()
  }
  if (!dir.exists(project_dir)) {
    stop("`project_dir` does not exist: ", project_dir)
  }

  # -----------------------------
  # small helpers
  # -----------------------------
  stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  say <- function(...) message(paste0("[", stamp(), "] ", paste(..., collapse = "")))

  kv_builder <- local({
    width <- NULL
    function(label, value = NULL) {
      if (is.null(value)) return(label)
      if (is.null(width)) {
        width <<- 2L + max(nchar(c(
          "Timestamp", "Alpha code", "Years processed",
          "Parallel workers", "Completed", "Failed",
          "Total elapsed"
        )))
      }
      sprintf("%-*s %s", width + 1L, paste0(label, ":"), value)
    }
  })

  sep_line <- strrep("-", 72L)

  # -----------------------------
  # validate inputs
  # -----------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L)
    stop("alpha_code must be a single character value", call. = FALSE)

  alpha_code <- toupper(trimws(alpha_code))

  if (!grepl("^[A-Z]{4}$", alpha_code))
    stop("alpha_code must be exactly four letters (A-Z)", call. = FALSE)

  # -----------------------------
  # timing
  # -----------------------------
  t_start <- Sys.time()

  # -----------------------------
  # paths & config
  # -----------------------------
  base_dir <- file.path(project_dir, "runs", alpha_code)
  log_fn   <- file.path(base_dir, "_log.txt")

  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }

  years <- seq(1980L, 2020L, by = 5L)

  # -----------------------------
  # parallel setup
  # -----------------------------
  say("Preparing parallel workers ...")

  n_cores <- max(1L, min(length(years), parallel::detectCores(logical = TRUE)))

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)

  parallel::clusterSetRNGStream(cl, as.integer(Sys.time()))

  say("Launched ", n_cores, " workers. Scheduling ", length(years), " runs ...")

  # -----------------------------
  # worker job
  # -----------------------------
  job_fun <- function(yr, alpha_code) {

    t0 <- proc.time()[["elapsed"]]

    out <- tryCatch({

      res <- create_ensemble_model(alpha_code, yr)

      list(
        ok      = TRUE,
        year    = yr,
        elapsed = proc.time()[["elapsed"]] - t0,
        result  = res,
        notes   = "ok"
      )

    }, error = function(e) {

      list(
        ok      = FALSE,
        year    = yr,
        elapsed = proc.time()[["elapsed"]] - t0,
        result  = NULL,
        notes   = conditionMessage(e)
      )
    })

    out
  }

  results <- parallel::parLapplyLB(cl, years, job_fun, alpha_code = alpha_code)

  names(results) <- as.character(years)

  say("All jobs completed. Collating results ...")

  # -----------------------------
  # logging
  # -----------------------------
  n_ok  <- sum(vapply(results, function(x) isTRUE(x$ok), logical(1)))
  n_err <- length(results) - n_ok

  fmt_hms <- function(secs) {
    if (is.na(secs) || !is.finite(secs)) return("-")
    s <- as.integer(round(secs))
    sprintf("%02d:%02d:%02d", s %/% 3600L, (s %% 3600L) %/% 60L, s %% 60L)
  }

  per_year_lines <- vapply(results, function(x) {

    status  <- if (x$ok) "completed" else "failed"
    elapsed <- fmt_hms(x$elapsed)
    outputs <- if (x$ok) "tif+asc+stats+pi" else "-"

    sprintf(
      "  %d  status=%-9s elapsed=%s   outputs=%s   notes=%s",
      as.integer(x$year), status, elapsed, outputs,
      if (x$ok) "-" else x$notes
    )

  }, character(1))

  total_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  header_lines <- c(
    sep_line,
    "Processing summary (create_timeseries)",
    kv_builder("Timestamp",        stamp()),
    kv_builder("Alpha code",       alpha_code),
    kv_builder("Bins processed",   length(years)),
    kv_builder("Parallel workers", n_cores),
    kv_builder("Completed",        n_ok),
    kv_builder("Failed",           n_err),
    kv_builder("Total elapsed",    sprintf("%s (%.2f s)", fmt_hms(total_elapsed), total_elapsed)),
    "Per-year status:",
    per_year_lines,
    ""
  )

  cat(paste0(paste(header_lines, collapse = "\n"), "\n"), file = log_fn, append = TRUE)

  say("Summary appended to ", log_fn)

  attr(results, "total_elapsed_sec") <- total_elapsed

  invisible(results)
}
