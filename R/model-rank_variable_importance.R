#' Parse and rank relative variable importance from an SDM text report
#'
#' Parse a "Relative Variable Importance" style text report and return a
#' tidy data.frame. Supports multiple metric sections and optional
#' aggregation across metrics for variable ranking.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' The parser expects one or more metric sections of the form:
#' \itemize{
#'   \item Based on <Metric> metric:
#'   \item var1 (40\%)
#'   \item var2 (20\%)
#'   \item ====
#' }
#'
#' Input forms:
#' \itemize{
#'   \item If input is a path to a readable file, the file is read
#'         line-by-line.
#'   \item If input is a single character string containing the entire
#'         report, it is split on \\n.
#'   \item If input is a character vector, it is treated as pre-split
#'         lines.
#' }
#'
#' \strong{Variable line parsing}
#' Lines like varname (12.3\%). Variable names may include letters,
#' digits, underscores, hyphens, and dots. If a line does not match
#' the strict pattern, a more permissive fallback (.*)\\((number)\%) is
#' attempted. Percentages are parsed as numeric and are not constrained
#' to (0, 100).
#'
#' \strong{Metric section headers}
#' Lines of the form Based on <Metric> metric: start sections. A section
#' is terminated by the next header or a separator line of = characters.
#'
#' \strong{Combining metrics}
#' When combine_metrics = TRUE, values are averaged using mean(...,
#' na.rm = TRUE). Sorting is descending by percent, then variable.
#'
#' \strong{Logging behavior}
#' When log = TRUE, a summary block is appended to:
#' \code{<rENM_project_dir()>/runs/<ALPHA_CODE>/}
#' \code{_log.txt}
#'
#' @param input Character. File path or text input.
#' @param combine_metrics Logical. Average across metrics if TRUE.
#' @param log Logical. Write summary to log.
#' @param log_path Character. Optional log file path.
#' @param alpha_code Character. Optional 4-letter code.
#'
#' @return Data frame with ranked variable importance.
#'
#' @importFrom stats aggregate
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' df <- rank_variable_importance("file.txt")
#' }
#'
#' @export
rank_variable_importance <- function(
    input,
    combine_metrics = TRUE,
    log = FALSE,
    log_path = NULL,
    alpha_code = NULL
) {

  resolve_log_path <- function() {
    if (!is.null(log_path) && is.character(log_path) && nzchar(log_path)) {
      return(log_path)
    }
    if (!is.null(alpha_code) && is.character(alpha_code) &&
        nchar(alpha_code) == 4L) {
      return(file.path(
        rENM_project_dir(),
        "runs",
        toupper(alpha_code),
        "_log.txt"
      ))
    }
    NA_character_
  }

  safe_dir_create <- function(path) {
    dir <- dirname(path)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(TRUE)
  }

  # -------------------------------------------------------------------------
  # Internal logging utilities (rENM standard)
  # -------------------------------------------------------------------------

  #' @keywords internal
  .sep_line <- strrep("-", 72)

  #' @keywords internal
  .pad_kv <- function(key, value) {
    sprintf("%-22s : %s", key, value)
  }

  #' @keywords internal
  .stamp <- function() {
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  }

  ## ---- ingest --------------------------------------------------------------
  t0 <- Sys.time()

  if (length(input) == 1L && is.character(input) && file.exists(input)) {
    lines <- readLines(input, warn = FALSE)
    input_type <- sprintf(
      "file (%s)",
      normalizePath(input, winslash = "/", mustWork = FALSE)
    )
  } else if (is.character(input)) {
    if (length(input) == 1L) {
      lines <- unlist(strsplit(input, "\n", fixed = TRUE),
                      use.names = FALSE)
      input_type <- "character (single string)"
    } else {
      lines <- input
      input_type <- "character (vector of lines)"
    }
  } else {
    stop("`input` must be a file path or character vector.",
         call. = FALSE)
  }

  if (!length(lines)) {
    stop("No content to parse.", call. = FALSE)
  }

  ## ---- parse ---------------------------------------------------------------
  metric_start_idx <- grep(
    "^\\s*Based\\s+on\\s+.+\\s+metric:\\s*$",
    lines, perl = TRUE, ignore.case = TRUE
  )

  if (!length(metric_start_idx)) {
    stop("No metric sections found.", call. = FALSE)
  }

  is_separator <- function(x) grepl("^=+$", trimws(x))

  rx_strict <- "^\\s*([A-Za-z0-9_\\.-]+)\\s+.*\\(([-+]?[0-9]*\\.?[0-9]+)\\s*%\\)"
  rx_loose  <- "^\\s*(.*?)\\s*\\(([-+]?[0-9]*\\.?[0-9]+)\\s*%\\)"

  out <- list()
  rows_parsed <- 0L

  for (i in seq_along(metric_start_idx)) {
    start <- metric_start_idx[i]
    metric_line <- lines[start]

    metric <- sub(
      "^\\s*Based\\s+on\\s+(.+)\\s+metric:\\s*$",
      "\\1",
      metric_line,
      perl = TRUE,
      ignore.case = TRUE
    )

    next_metric <- if (i < length(metric_start_idx))
      metric_start_idx[i + 1L] else length(lines) + 1L

    sep_candidates <- which(is_separator(lines))
    sep_candidates <- sep_candidates[
      sep_candidates > start & sep_candidates < next_metric
    ]

    end <- if (length(sep_candidates))
      min(sep_candidates) - 1L else next_metric - 1L

    if (end <= start) next

    block <- lines[(start + 1L):end]

    hits_strict <- grepl(rx_strict, block, perl = TRUE)
    hits_loose  <- !hits_strict & grepl(rx_loose, block, perl = TRUE)

    any_hits <- hits_strict | hits_loose
    if (!any(any_hits)) next

    vars <- character(sum(any_hits))
    pcnt <- numeric(sum(any_hits))

    if (any(hits_strict)) {
      idx <- which(hits_strict)
      vars[seq_along(idx)] <- sub(rx_strict, "\\1", block[idx])
      pcnt[seq_along(idx)] <- as.numeric(sub(rx_strict, "\\2", block[idx]))
    }

    if (any(hits_loose)) {
      idx_l <- which(hits_loose)
      offset <- sum(hits_strict)
      vars[(offset + seq_along(idx_l))] <-
        trimws(sub(rx_loose, "\\1", block[idx_l]))
      pcnt[(offset + seq_along(idx_l))] <-
        as.numeric(sub(rx_loose, "\\2", block[idx_l]))
    }

    ok <- is.finite(pcnt) & nzchar(vars)
    vars <- vars[ok]
    pcnt <- pcnt[ok]
    if (!length(vars)) next

    df <- data.frame(
      variable = vars,
      percent  = pcnt,
      metric   = rep(metric, length(vars)),
      stringsAsFactors = FALSE
    )

    rows_parsed <- rows_parsed + nrow(df)
    out[[length(out) + 1L]] <- df
  }

  if (!length(out)) {
    stop("No variable lines parsed.", call. = FALSE)
  }

  all_df <- do.call(rbind, out)

  ## ---- combine -------------------------------------------------------------
  if (isTRUE(combine_metrics)) {
    agg <- stats::aggregate(
      percent ~ variable,
      data = all_df,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    agg <- agg[order(-agg$percent, agg$variable), ]
    rownames(agg) <- NULL
    result <- agg
  } else {
    all_df <- all_df[
      with(all_df, order(metric, -percent, variable)),
    ]
    rownames(all_df) <- NULL
    result <- all_df
  }

  ## ---- logging -------------------------------------------------------------
  if (isTRUE(log)) {
    lp <- resolve_log_path()
    if (!is.na(lp)) {
      safe_dir_create(lp)

      topN <- utils::head(result, n = min(5L, nrow(result)))

      top_lines <- if (nrow(topN)) {
        paste0(
          "  - ",
          sprintf("%-24s %6.2f%%", topN$variable, topN$percent)
        )
      } else {
        "  (none)"
      }

      block <- c(
        "",
        .sep_line,
        "Processing summary (rank_variable_importance)",
        .pad_kv("Timestamp", .stamp()),
        .pad_kv("Alpha code",
                if (!is.null(alpha_code)) toupper(alpha_code) else "-"),
        .pad_kv("Input type", input_type),
        .pad_kv("Rows parsed", rows_parsed),
        "Top variables:",
        top_lines,
        ""
      )

      con <- file(lp, open = "a", encoding = "UTF-8")
      on.exit(close(con), add = TRUE)
      writeLines(block, con, sep = "\n", useBytes = TRUE)
    }
  }

  result
}
