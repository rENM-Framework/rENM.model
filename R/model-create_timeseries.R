#' Build a full rENM time series for a species
#'
#' Schedules and executes \code{\link{create_ensemble_model}} for each 5-year
#' bin from 1980 to 2020 in parallel, producing a complete rENM time series
#' for a given species.
#'
#' @details
#' Workers are launched in a PSOCK cluster of size
#' \code{min(number_of_years, detectCores())}. Each year is processed
#' independently so errors for individual years are captured and logged without
#' aborting the full time series. The cluster is always stopped on exit.
#'
#' \code{create_ensemble_model()} is referenced via namespace closure and is
#' available to workers in both installed-package and \code{devtools::load_all()}
#' workflows.
#'
#' A processing summary is appended to
#' \code{<project_dir>/runs/<alpha_code>/_log.txt}.
#'
#' @param alpha_code Character. Four-letter banding code (e.g., \code{"CASP"}).
#' @param project_dir Character. Path to the rENM project root. If NULL,
#'   resolved via \code{\link{rENM_project_dir}}.
#'
#' @return Invisible named list keyed by year. Each element contains:
#'   \code{ok} (logical), \code{year} (integer), \code{elapsed} (seconds),
#'   \code{result} (output from \code{create_ensemble_model()} or NULL),
#'   \code{notes} (character: \code{"ok"} or an error message).
#'   The attribute \code{total_elapsed_sec} records total wall-clock time.
#'   Side effects are parallel job execution, file outputs per year, and a
#'   summary log entry.
#'
#' @seealso \code{\link{create_ensemble_model}}, \code{\link{rENM_project_dir}}
#' @importFrom parallel detectCores makeCluster stopCluster clusterSetRNGStream parLapplyLB
#'
#' @examples
#' \dontrun{
#'   create_timeseries("CASP")
#'   create_timeseries("CASP", project_dir = "/my/project")
#' }
#'
#' @export
create_timeseries <- function(alpha_code, project_dir = NULL) {

  if (is.null(project_dir)) {
    project_dir <- rENM_project_dir()
  }
  if (!dir.exists(project_dir)) {
    stop("`project_dir` does not exist: ", project_dir)
  }

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

  if (!is.character(alpha_code) || length(alpha_code) != 1L)
    stop("alpha_code must be a single character value", call. = FALSE)

  alpha_code <- toupper(trimws(alpha_code))

  if (!grepl("^[A-Z]{4}$", alpha_code))
    stop("alpha_code must be exactly four letters (A-Z)", call. = FALSE)

  t_start <- Sys.time()

  base_dir <- file.path(project_dir, "runs", alpha_code)
  log_fn   <- file.path(base_dir, "_log.txt")

  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }

  years <- seq(1980L, 2020L, by = 5L)

  say("Preparing parallel workers ...")

  n_cores <- max(1L, min(length(years), parallel::detectCores(logical = TRUE)))

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)

  parallel::clusterSetRNGStream(cl, as.integer(Sys.time()))

  say("Launched ", n_cores, " workers. Scheduling ", length(years), " runs ...")

  # Worker job; create_ensemble_model() is resolved via closure
  cem_fn <- create_ensemble_model

  job_fun <- function(yr, alpha_code) {

    t0 <- proc.time()[["elapsed"]]

    out <- tryCatch({

      res <- cem_fn(alpha_code, yr)

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
