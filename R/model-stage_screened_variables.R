#' Stage screened environmental variables into the time series
#'
#' Stage (copy) ranked raster predictor layers into per-year TimeSeries
#' directories based on a master variable-ranking table. Supports several
#' ranking metrics and an optional top-N filter.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Master ranking CSV, one per year:
#'     \code{<proj_root>/runs/<alpha_code>/TimeSeries/<year>/vars/}
#'     \code{_<alpha_code>-<year>-Variable-Ranking.csv}
#'   \item Source rasters:
#'     \code{<proj_root>/runs/<alpha_code>/_vars/<year>/}
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Ranked rasters copied to:
#'     \code{<proj_root>/runs/<alpha_code>/TimeSeries/<year>/vars/}
#'   \item Per-year status line and a 72-dash summary appended to:
#'     \code{<proj_root>/runs/<alpha_code>/_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Ranking columns:
#'     \code{rank_ranksum}, \code{rank_medianPI},
#'     \code{rank_ratio}, \code{rank_wz}
#'   \item Tie-breakers use the associated score column
#'     (higher or lower depending on metric).
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Years must be multiples of five in the range
#'     \code{\{1980, 1985, ..., 2020\}}.
#'   \item Raster filenames must match the \code{var} column
#'     in the ranking CSV.
#' }
#'
#' @param alpha_code Character. Species alpha code (e.g., "CASP"). Used to
#'   resolve paths under the project \code{runs} directory.
#' @param year Integer. Vector of years or single year. Defaults to
#'   \code{seq(1980, 2020, 5)}. Each value must be a 5-year step within
#'   1980-2020.
#' @param top_n Integer or NULL. Number of top-ranked variables to stage per
#'   year. Default = 10; use \code{NULL} to copy all variables.
#' @param rank_method Character. One of "ranksum","medianPI","ratio","wz".
#'   Selects the ranking column used to order variables. Default = "ratio".
#'
#' @return Invisible named list with:
#' \itemize{
#'   \item \code{alpha_code} Character
#'   \item \code{years} Integer vector
#'   \item \code{rank_method} Character
#'   \item \code{top_n} Integer or NULL
#'   \item \code{copied_paths} Character vector
#'   \item \code{counts} List: \code{total_years}, \code{attempted_files},
#'     \code{copied_files}, \code{missing_files},
#'     \code{years_missing_ranking}
#'   \item \code{elapsed_seconds} Numeric
#'   \item \code{log_file} Character
#' }
#'
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#' stage_screened_variables("CASP", top_n = 5, rank_method = "medianPI")
#' }
#'
#' @export
stage_screened_variables <- function(alpha_code,
                                     year = seq(1980, 2020, 5),
                                     top_n = 10,
                                     rank_method = "ratio") {
  start_time <- Sys.time()

  # ---- validate rank_method -> rank/score column mapping -------------------
  valid_methods <- c("ranksum", "medianPI", "ratio", "wz")
  if (!rank_method %in% valid_methods) {
    stop("Invalid rank_method: must be one of ", paste(valid_methods, collapse = ", "))
  }
  rank_col_map <- c(
    ranksum  = "rank_ranksum",
    medianPI = "rank_medianPI",
    ratio    = "rank_ratio",
    wz       = "rank_wz"
  )
  score_col_map <- c(
    ranksum  = "ranksum_score",  # NOTE: lower is better for tie-breaks
    medianPI = "median_PI",
    ratio    = "ratio_score",
    wz       = "wz_score"
  )
  chosen_rank_col  <- unname(rank_col_map[[rank_method]])
  chosen_score_col <- unname(score_col_map[[rank_method]])

  # ---- validate year argument ----------------------------------------------
  valid_steps <- seq(1980, 2020, 5)
  if (length(year) == 1 && year %% 5 == 0 && year >= 1980 && year <= 2020) {
    years <- as.integer(year)
  } else if (length(year) >= 1 && all(year %in% valid_steps)) {
    years <- as.integer(year)
  } else {
    stop("Invalid 'year'. Use a single 5-year step within 1980 - 2020 ",
         "or a vector subset of 1980, 1985, ..., 2020.")
  }

  # ---- validate top_n ------------------------------------------------------
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n <= 0) {
      stop("`top_n` must be a single positive integer or NULL.")
    }
    top_n <- as.integer(top_n)
  }

  # ---- paths ---------------------------------------------------------------
  proj_root    <- rENM_project_dir()
  base_dir     <- file.path(proj_root, "runs", alpha_code)
  src_vars_root <- file.path(base_dir, "_vars")
  log_file     <- file.path(base_dir, "_log.txt")

  total_years <- length(years)
  attempted_files <- copied_files <- missing_files <- years_missing_ranking <- 0L
  copied_paths <- character(0)
  per_year_lines <- character(0)
  sep72 <- paste(rep("-", 72), collapse = "")

  for (yr in years) {
    dest_dir <- file.path(base_dir, "TimeSeries", as.character(yr), "vars")
    rnk_file <- file.path(dest_dir, sprintf("_%s-%s-Variable-Ranking.csv", alpha_code, yr))
    src_year_dir <- file.path(src_vars_root, as.character(yr))

    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

    if (!file.exists(rnk_file)) {
      warning(sprintf("Missing ranking CSV for %s %s: %s", alpha_code, yr, rnk_file))
      years_missing_ranking <- years_missing_ranking + 1L
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6s %-10s %s\n", yr, "Missing", "-"))
      next
    }

    rnk_df <- tryCatch(read.csv(rnk_file, stringsAsFactors = FALSE),
                       error = function(e) NULL)
    if (is.null(rnk_df)) {
      warning(sprintf("Failed to read ranking CSV for %s %s: %s", alpha_code, yr, rnk_file))
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6s %-10s %s\n", yr, "Error", "-"))
      next
    }

    # Ensure required columns exist
    var_col <- if ("var" %in% names(rnk_df)) "var" else {
      idx <- match("var", tolower(names(rnk_df)))
      if (is.na(idx)) NA_character_ else names(rnk_df)[idx]
    }
    if (is.na(var_col)) {
      warning(sprintf("Ranking CSV lacks 'var' column for %s %s: %s", alpha_code, yr, rnk_file))
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6s %-10s %s\n", yr, "Error", dest_dir))
      next
    }
    if (!chosen_rank_col %in% names(rnk_df)) {
      warning(sprintf("Ranking CSV lacks '%s' column for %s %s: %s",
                      chosen_rank_col, alpha_code, yr, rnk_file))
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6s %-10s %s\n", yr, "Error", dest_dir))
      next
    }

    # ---- Fixed ordering: build full-length tie-breaker & correct direction ---
    score_vec <- if (chosen_score_col %in% names(rnk_df)) rnk_df[[chosen_score_col]] else rep(0, nrow(rnk_df))
    decreasing_score <- rank_method %in% c("medianPI", "wz", "ratio")
    tie_vec <- if (decreasing_score) -score_vec else score_vec
    ord <- order(rnk_df[[chosen_rank_col]], tie_vec, na.last = NA)
    ordered_df <- rnk_df[ord, , drop = FALSE]

    vars_list <- unique(na.omit(trimws(as.character(ordered_df[[var_col]]))))
    if (length(vars_list) == 0) {
      warning(sprintf("Ranking CSV has no variables after ordering for %s %s: %s", alpha_code, yr, rnk_file))
      per_year_lines <- c(per_year_lines,
                          sprintf("%-6s %-10s %s\n", yr, "Empty", dest_dir))
      next
    }

    # Restrict to top_n if provided
    if (!is.null(top_n) && length(vars_list) > top_n) {
      vars_list <- vars_list[seq_len(top_n)]
    }

    exts <- c(".tif", ".asc")
    attempted_this_year <- copied_this_year <- 0L

    for (v in vars_list) {
      for (ext in exts) {
        src_file  <- file.path(src_year_dir, paste0(v, ext))
        dest_file <- file.path(dest_dir, paste0(v, ext))
        attempted_this_year <- attempted_this_year + 1L
        attempted_files <- attempted_files + 1L

        if (file.exists(src_file)) {
          ok <- tryCatch(file.copy(src_file, dest_file, overwrite = TRUE),
                         error = function(e) FALSE)
          if (isTRUE(ok)) {
            copied_this_year <- copied_this_year + 1L
            copied_files <- copied_files + 1L
            copied_paths <- c(copied_paths, dest_file)
            message(paste("Copied:", src_file, "->", dest_file))
          } else {
            missing_files <- missing_files + 1L
            warning(sprintf("Copy failed: %s -> %s", src_file, dest_file))
          }
        } else {
          missing_files <- missing_files + 1L
          warning(sprintf("Missing source raster: %s", src_file))
        }
      }
    }

    status <- if (attempted_this_year == 0L) {
      "Empty"
    } else if (copied_this_year == 0L) {
      "Missing"
    } else if (copied_this_year < attempted_this_year) {
      sprintf("Partial %d/%d", copied_this_year, attempted_this_year)
    } else {
      sprintf("Copied %d/%d", copied_this_year, attempted_this_year)
    }

    per_year_lines <- c(per_year_lines,
                        sprintf("%-6s %-10s %s\n", yr, status, dest_dir))
  }

  # ---- timing & summary -----------------------------------------------------
  elapsed   <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  outputs_saved_block <- if (length(copied_paths) > 0) {
    paste0("\n", paste(copied_paths, collapse = "\n"))
  } else {
    "None"
  }

  yearly_text <- paste0(
    "\n", sep72, "\n",
    "Yearly staging details (stage_screened_variables)\n",
    sprintf("%-6s %-10s %s\n", "Year", "Status", "Destination"),
    paste(per_year_lines, collapse = ""),
    sep72, "\n"
  )

  summary_text <- paste0(
    "\n", sep72, "\n",
    "Processing summary (stage_screened_variables)\n",
    sprintf("%-20s %s\n", "Timestamp:", timestamp),
    sprintf("%-20s %s\n", "Alpha code:", alpha_code),
    sprintf("%-20s %s\n", "Rank method:", rank_method),
    sprintf("%-20s %s\n", "Rank column:", chosen_rank_col),
    sprintf("%-20s %s\n", "Top N:", if (is.null(top_n)) "all" else as.character(top_n)),
    sprintf("%-20s %s\n", "Year:", paste(years, collapse = ", ")),
    sprintf("%-20s %d\n", "Total years:", total_years),
    sprintf("%-20s %d\n", "Attempted files:", attempted_files),
    sprintf("%-20s %d\n", "Copied files:", copied_files),
    sprintf("%-20s %d\n", "Missing files:", missing_files),
    sprintf("%-20s %d\n", "Years missing ranking:", years_missing_ranking),
    sprintf("%-20s %s\n", "Outputs saved:", outputs_saved_block),
    sprintf("%-20s %s seconds\n", "Total elapsed:", elapsed)
  )

  # Append summaries to the run log
  # cat(yearly_text,  file = log_file, append = TRUE)
  cat(summary_text, file = log_file, append = TRUE)

  # Echo to console as well (optional)
  # cat(yearly_text)
  cat(summary_text)

  invisible(list(
    alpha_code = alpha_code,
    years = years,
    rank_method = rank_method,
    top_n = top_n,
    copied_paths = copied_paths,
    counts = list(
      total_years = total_years,
      attempted_files = attempted_files,
      copied_files = copied_files,
      missing_files = missing_files,
      years_missing_ranking = years_missing_ranking
    ),
    elapsed_seconds = elapsed,
    log_file = log_file
  ))
}
