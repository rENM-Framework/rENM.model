#' Reduce variable covariance
#'
#' Detects multicollinearity in the stagged predictor rasters through an
#' adaptive Variable Inflation Factor (VIF) routine, records results, and
#' retains only non-correlated variables for further processing.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' \itemize{
#'   \item Input rasters: \code{<project_dir>/runs/<alpha_code>/TimeSeries/
#'         <year>/vars}
#'   \item Log file: \code{<project_dir>/runs/<alpha_code>/_log.txt}
#' }
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Predictor rasters matching \code{file_patterns}
#'   \item Optional exclusions via \code{exclude_vars}
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Results CSV: \code{_<alpha_code>-<year>-Correlation-Results.csv}
#'   \item Directory of moved rasters (only if variables are dropped):
#'         \code{_<alpha_code>-<year>-Correlated-Variables}
#'   \item Log entry appended to \code{_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Adaptive sampling over \code{sample_schedule} sizes until the kept
#'         set is stable for \code{stability_patience} iterations
#'   \item Multicollinearity screening with \code{usdm::vifstep()}
#'         at threshold \code{vif_threshold}
#'   \item Reporting of high-correlation pairs where
#'         \eqn{|r| \ge corr_threshold}
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Raster layers must share identical geometry
#'   \item Layer names must correspond to predictor identifiers
#' }
#'
#' @param alpha_code Character. Species alpha code (for example, \code{"CASP"}).
#' @param year Integer vector, or \code{NULL}. Years to process; if
#'   \code{NULL} (default) the sequence 1980, 1985, ..., 2020 is used.
#' @param vif_threshold Numeric. VIF threshold passed to
#'   \code{usdm::vifstep()} (default 10).
#' @param sample_schedule Integer vector. Adaptive sample sizes tried in
#'   order; the routine stops early once the kept set stabilizes. Default
#'   \code{c(10000, 20000, 50000)}.
#' @param stability_patience Integer. Consecutive identical kept sets needed
#'   for convergence (default 1).
#' @param corr_threshold Numeric. Absolute Pearson correlation threshold used
#'   when listing high-correlation pairs (default 0.95).
#' @param exclude_vars Character vector. Predictor names that must always be
#'   retained. Case-sensitive match against sampled column names (default
#'   \code{character(0)}).
#' @param file_patterns Character vector. Regular expressions selecting raster
#'   filenames to include. Matching stems trigger companion \code{.asc} moves.
#'   Default \code{c("\\\\.tif$", "\\\\.grd$", "\\\\.img$", "\\\\.bil$")}.
#' @param overwrite_results Logical. Overwrite an existing results CSV
#'   (default \code{TRUE}).
#'
#' @return Invisible list keyed by year containing:
#' \itemize{
#'   \item \code{kept} – character vector of retained predictors
#'   \item \code{dropped} – character vector of removed predictors
#'   \item \code{vif_table} – data frame of VIF results for all variables
#'   \item \code{cor_pairs} – data frame of high-correlation pairs
#'   \item \code{results_csv} – file path to the written CSV
#'   \item \code{moved_to} – directory path for relocated rasters or
#'         \code{NA} when none were moved
#' }
#'
#' @importFrom raster stack nlayers sampleRandom
#' @importFrom stats complete.cases na.omit cor
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.csv write.table
#' @importFrom usdm vifstep
#'
#' @examples
#' \dontrun{
#'   reduce_covariance("CASP")
#'   reduce_covariance("CASP", year = 2005,
#'                     exclude_vars = c("bio1", "bio12"))
#'   reduce_covariance("CASP", year = 2000,
#'                     sample_schedule = c(8000, 16000, 32000, 64000))
#' }
#'
#' @export
reduce_covariance <- function(alpha_code,
                              year = NULL,
                              vif_threshold = 10,
                              sample_schedule = c(10000, 20000, 50000),
                              stability_patience = 1,
                              corr_threshold = 0.95,
                              exclude_vars = character(0),
                              file_patterns = c("\\.tif$", "\\.grd$", "\\.img$", "\\.bil$"),
                              overwrite_results = TRUE) {

  # ============================== Setup & checks ==============================
  start_time_all <- Sys.time()
  stopifnot(is.character(alpha_code), length(alpha_code) == 1, nchar(alpha_code) > 0)
  stopifnot(is.numeric(vif_threshold),  length(vif_threshold)  == 1)
  stopifnot(is.numeric(corr_threshold), length(corr_threshold) == 1)
  stopifnot(is.numeric(sample_schedule), all(sample_schedule > 0))
  stopifnot(is.numeric(stability_patience), stability_patience >= 1)

  years <- if (is.null(year)) seq(1980, 2020, by = 5) else as.integer(year)
  stopifnot(all(is.numeric(years)))

  # Required packages
  pkgs <- c("raster", "usdm", "tools")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing))
    stop(sprintf("Missing required package(s): %s", paste(missing, collapse = ", ")))

  # ---------------------------------------------------------------------------
  # Use package helper to locate the project root
  # ---------------------------------------------------------------------------
  project_dir <- rENM_project_dir()

  # ------------------------------- Helpers ------------------------------------
  rule72 <- paste(rep("-", 72), collapse = "")

  append_log <- function(alpha_code, raster_source, n_total, n_valid,
                         outputs_saved, elapsed, output_file) {
    log_dir  <- file.path(project_dir, "runs", alpha_code)
    log_path <- file.path(log_dir, "_log.txt")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    block <- sprintf(
      paste0(
        "\n%s\n",
        "Processing summary (reduce_covariance)\n",
        "Timestamp          : %s\n",
        "Alpha code         : %s\n",
        "Raster source      : %s\n",
        "Total cells        : %d\n",
        "Valid cells        : %d\n",
        "Outputs saved      : %s\n",
        "Total elapsed      : %s\n",
        "Output file        : %s\n"
      ),
      rule72,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      alpha_code,
      raster_source,
      as.integer(n_total),
      as.integer(n_valid),
      outputs_saved,
      sprintf("%.2f sec", as.numeric(elapsed, units = "secs")),
      output_file
    )
    cat(block, file = log_path, append = TRUE)
    invisible(log_path)
  }

  list_var_files <- function(dir_path, patterns) {
    files <- unlist(lapply(patterns, function(pat)
      list.files(dir_path, pattern = pat, full.names = TRUE, ignore.case = TRUE)))
    unique(files)
  }

  high_corr_pairs <- function(cor_mat, thr = 0.7) {
    if (is.null(cor_mat))
      return(data.frame(var1 = character(), var2 = character(), r = numeric(),
                        stringsAsFactors = FALSE))
    cm <- cor_mat
    cm[upper.tri(cm, diag = TRUE)] <- NA
    idx <- which(abs(cm) >= thr, arr.ind = TRUE)
    if (!nrow(idx))
      return(data.frame(var1 = character(), var2 = character(), r = numeric(),
                        stringsAsFactors = FALSE))
    data.frame(var1 = rownames(cm)[idx[, 1]],
               var2 = colnames(cm)[idx[, 2]],
               r    = cm[idx], stringsAsFactors = FALSE)
  }

  find_companion_files <- function(dir_path, stems) {
    if (!length(stems)) return(character(0))
    norm <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
    stems_norm <- unique(norm(stems))
    cand <- list.files(dir_path, pattern = "\\.(tif|tiff|asc)$",
                       full.names = TRUE, ignore.case = TRUE)
    if (!length(cand)) return(character(0))
    file_stems <- tolower(gsub("[^a-z0-9]+", "",
                               tools::file_path_sans_ext(basename(cand))))
    cand[file_stems %in% stems_norm]
  }

  # ================================ Main loop =================================
  results <- list()

  for (yr in years) {
    year_start <- Sys.time()
    cat(sprintf("\n[reduce_covariance] %s - Processing year %d\n", alpha_code, yr))

    base_dir      <- file.path(project_dir, "runs", alpha_code,
                               "TimeSeries", yr, "vars")
    base_dir_exp  <- normalizePath(base_dir, mustWork = FALSE)

    if (!dir.exists(base_dir_exp)) {
      warning(sprintf("Directory does not exist: %s - skipping %d", base_dir, yr))
      next
    }

    var_files <- list_var_files(base_dir_exp, file_patterns)
    if (!length(var_files)) {
      warning(sprintf("No raster files found in %s - skipping %d", base_dir, yr))
      next
    }
    cat(sprintf("[reduce_covariance] Found %d file(s)\n", length(var_files)))

    # Load and prep stack
    rs <- raster::stack(var_files)
    if (any(duplicated(names(rs)))) names(rs) <- make.unique(names(rs))
    cat(sprintf("[reduce_covariance] Stack has %d layer(s)\n", raster::nlayers(rs)))

    # --------------------------- Adaptive sampling ----------------------------
    prev_kept         <- NULL
    stable_hits       <- 0
    kept_final        <- NULL
    vif_tbl_final     <- NULL
    all_vars          <- NULL
    cor_pairs_final   <- NULL
    last_sample_n     <- 0
    last_valid_n      <- 0

    for (n_request in sample_schedule) {
      n_request      <- min(n_request, raster::ncell(rs))
      last_sample_n  <- n_request
      cat(sprintf("[reduce_covariance] Adaptive sample: up to %d cells\n", n_request))

      samp      <- raster::sampleRandom(rs, size = n_request, na.rm = TRUE,
                                        as.data.frame = TRUE, sp = FALSE)
      valid_n   <- sum(stats::complete.cases(samp))
      if (valid_n == 0) {
        warning(sprintf("All samples NA for %d at n = %d - trying next size", yr, n_request))
        next
      }
      samp           <- stats::na.omit(samp)
      last_valid_n   <- valid_n
      all_vars       <- colnames(samp)

      # Correlation (reporting only)
      cor_mat        <- tryCatch(stats::cor(samp, use = "pairwise.complete.obs"),
                                 error = function(e) NULL)
      cor_pairs      <- high_corr_pairs(cor_mat, thr = corr_threshold)
      if (nrow(cor_pairs)) {
        cat(sprintf("[reduce_covariance] %d high-|r| pairs (>= %.2f)\n",
                    nrow(cor_pairs), corr_threshold))
      } else {
        cat(sprintf("[reduce_covariance] No high-|r| pairs at >= %.2f\n",
                    corr_threshold))
      }

      # VIF selection
      cat(sprintf("[reduce_covariance] usdm::vifstep(th = %.2f)\n", vif_threshold))
      vif_obj <- tryCatch(
        suppressMessages(usdm::vifstep(samp, th = vif_threshold, trace = FALSE)),
        error = function(e) suppressMessages(usdm::vifstep(samp, th = vif_threshold))
      )

      vif_tbl <- tryCatch(
        as.data.frame(vif_obj@results, stringsAsFactors = FALSE),
        error = function(e) {
          if (is.data.frame(vif_obj)) vif_obj
          else if (!is.null(vif_obj$results))
            as.data.frame(vif_obj$results, stringsAsFactors = FALSE)
          else stop("Cannot extract VIF results from usdm::vifstep() return value")
        })

      if (!all(c("Variables", "VIF") %in% names(vif_tbl))) {
        vcol <- names(vif_tbl)[grepl("Var", names(vif_tbl), ignore.case = TRUE)][1]
        ccol <- names(vif_tbl)[grepl("VIF", names(vif_tbl), ignore.case = TRUE)][1]
        if (is.na(vcol) || is.na(ccol))
          stop("Cannot locate Variables/VIF columns in VIF results table")
        names(vif_tbl)[names(vif_tbl) == vcol] <- "Variables"
        names(vif_tbl)[names(vif_tbl) == ccol] <- "VIF"
      }

      kept_now <- unique(as.character(vif_tbl$Variables))

      # Enforce exclusion list
      if (length(exclude_vars)) {
        present_exclusions <- intersect(exclude_vars, all_vars)
        missing_exclusions <- setdiff(exclude_vars, all_vars)
        if (length(missing_exclusions))
          warning(sprintf("[reduce_covariance] %d exclusion(s) not found for %d: %s",
                          length(missing_exclusions), yr,
                          paste(missing_exclusions, collapse = ", ")))
        kept_now <- union(kept_now, present_exclusions)
      }

      # Stability check
      if (!is.null(prev_kept) && setequal(prev_kept, kept_now)) {
        stable_hits <- stable_hits + 1
        cat(sprintf("[reduce_covariance] Kept set stable (%d/%d)\n",
                    stable_hits, stability_patience))
      } else {
        stable_hits <- 0
        cat("[reduce_covariance] Kept set changed; continuing sampling\n")
      }
      prev_kept        <- kept_now
      kept_final       <- kept_now
      vif_tbl_final    <- vif_tbl
      cor_pairs_final  <- cor_pairs

      if (stable_hits >= stability_patience) {
        cat("[reduce_covariance] Adaptive sampling converged; stopping early\n")
        break
      }
    } # adaptive loop

    if (is.null(kept_final)) {
      warning(sprintf("Adaptive sampling failed for %d - skipping", yr))
      next
    }

    dropped_final <- setdiff(all_vars, kept_final)
    cat(sprintf("[reduce_covariance] Kept %d; Dropped %d\n",
                length(kept_final), length(dropped_final)))

    # --------------------------- Write results CSV ---------------------------
    out_csv_name <- sprintf("_%s-%d-Correlation-Results.csv", alpha_code, yr)
    out_csv_path <- file.path(base_dir_exp, out_csv_name)

    if (!file.exists(out_csv_path) || overwrite_results) {
      tmpfile <- tempfile(fileext = ".csv")

      all_tbl <- data.frame(
        Variables = all_vars,
        VIF       = NA_real_,
        Kept      = all_vars %in% kept_final,
        stringsAsFactors = FALSE
      )
      if (nrow(vif_tbl_final)) {
        m   <- match(all_tbl$Variables, vif_tbl_final$Variables)
        has <- !is.na(m)
        all_tbl$VIF[has] <- vif_tbl_final$VIF[m[has]]
      }
      utils::write.csv(all_tbl, tmpfile, row.names = FALSE)

      # Second list
      cat("\nHighCorrelationPairs\n", file = tmpfile, append = TRUE)
      cat("var1,var2,r\n",           file = tmpfile, append = TRUE)
      if (nrow(cor_pairs_final)) {
        utils::write.table(cor_pairs_final[, c("var1", "var2", "r")],
                           file      = tmpfile, sep = ",",
                           row.names = FALSE, col.names = FALSE, append = TRUE)
      }

      cat(sprintf("\nSummary,All,%d\nSummary,Kept,%d\nSummary,Dropped,%d\n",
                  length(all_vars), length(kept_final), length(dropped_final)),
          file = tmpfile, append = TRUE)

      file.rename(tmpfile, out_csv_path)
      cat(sprintf("[reduce_covariance] Wrote CSV: %s\n", out_csv_path))
    } else {
      warning(sprintf("CSV exists and overwrite_results = FALSE: %s", out_csv_path))
    }

    # ------------------- Move correlated variables (if any) -------------------
    moved_files    <- character(0)
    moved_dir_path <- NA_character_

    if (length(dropped_final) == 0) {
      cat("[reduce_covariance] No correlated variables detected by VIF - nothing to move\n")
    } else {
      moved_dir_name <- sprintf("_%s-%d-Correlated-Variables", alpha_code, yr)
      moved_dir_path <- file.path(base_dir_exp, moved_dir_name)
      if (!dir.exists(moved_dir_path))
        dir.create(moved_dir_path, recursive = TRUE, showWarnings = FALSE)

      files_to_move <- find_companion_files(base_dir_exp, dropped_final)

      if (length(files_to_move)) {
        for (f in files_to_move) {
          dest <- file.path(moved_dir_path, basename(f))
          if (file.exists(dest))
            dest <- file.path(moved_dir_path,
                              sprintf("%s__%s", format(Sys.time(), "%Y%m%d%H%M%S"),
                                      basename(f)))
          ok <- file.rename(f, dest)
          if (!ok) {
            okc <- file.copy(f, dest, overwrite = FALSE,
                             copy.mode = TRUE, copy.date = TRUE)
            if (okc) unlink(f)
            ok <- okc
          }
          if (ok) moved_files <- c(moved_files, dest)
        }
        cat(sprintf("[reduce_covariance] Moved %d file(s) to %s\n",
                    length(moved_files), moved_dir_path))
      } else {
        cat("[reduce_covariance] Dropped variables identified, but no matching files found to move\n")
      }
    }

    # ------------------------------- Logging ----------------------------------
    elapsed <- Sys.time() - year_start
    outputs_saved <- paste(c(basename(out_csv_path),
                             paste0(length(moved_files), " moved file(s)")),
                           collapse = " | ")
    append_log(alpha_code    = alpha_code,
               raster_source = base_dir,
               n_total       = last_sample_n,
               n_valid       = last_valid_n,
               outputs_saved = outputs_saved,
               elapsed       = elapsed,
               output_file   = out_csv_path)

    # ------------------------------ Collect out -------------------------------
    results[[as.character(yr)]] <- list(
      kept        = kept_final,
      dropped     = dropped_final,
      vif_table   = vif_tbl_final,
      cor_pairs   = cor_pairs_final,
      results_csv = out_csv_path,
      moved_to    = moved_dir_path
    )

    cat(sprintf("[reduce_covariance] Year %d complete (%.2f sec)\n",
                yr, as.numeric(elapsed, units = "secs")))
  }

  # ============================== Wrap-up & return =============================
  total_elapsed <- Sys.time() - start_time_all
  cat(sprintf("\n[reduce_covariance] Finished all requested year(s) in %.2f sec\n",
              as.numeric(total_elapsed, units = "secs")))
  invisible(results)
}
