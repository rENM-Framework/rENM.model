## Global vars for dplyr / NSE columns in this file
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "var", "perm_imp", "n_pairs", "sd_PI", "median_PI", "IQR_PI",
    "rank_medianPI", "rank_stability", "z_PI", "z_IQR",
    "ranksum_score", "wz_score", "ratio_score",
    "mean_PI", "se_mean_PI", "se_median_PI",
    "rank_ranksum", "rank_wz", "rank_ratio"
  ))
}

#' Iterative bivariate ENM screening with adaptive, stability-aware convergence
#'
#' Repeatedly fits bivariate \code{maxnet} ecological niche models, aggregates
#' permutation importance (PI) scores, and stops when both the mean PI and
#' membership of the adaptive top-k variable set stabilize. Designed for rapid,
#' reproducible predictor down-selection across 5-year intervals without the
#' Java dependency required by \code{dismo::maxent()}.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' \code{rENM_project_dir()}.
#'
#' \strong{Pipeline context}
#' \itemize{
#'   \item Reads occurrences from \code{runs/<alpha_code>/_occs/of-<year>.csv}
#'   \item Reads raster predictors from \code{runs/<alpha_code>/_vars/<year>/}
#'   \item Writes outputs to
#'         \code{runs/<alpha_code>/TimeSeries/<year>/vars/}
#' }
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Occurrence table with \code{longitude}, \code{latitude}
#'   \item At least two GeoTIFF predictor layers
#' }
#'
#' \strong{Methods}
#' The function repeatedly samples bivariate predictor pairs and fits
#' \code{maxnet} models using occurrence and background environmental values
#' extracted from the selected raster layers. Background points are drawn
#' directly from non-\code{NA} raster cells using
#' \code{raster::sampleRandom(xy = TRUE, na.rm = TRUE)}, eliminating the need
#' for \code{dismo::randomPoints()}.
#'
#' For each fitted pair, permutation importance is computed manually in a
#' MaxEnt-like way: each predictor is permuted in turn on the model evaluation
#' table, predictions are recomputed without refitting, the drop in training AUC
#' is measured, and the resulting drops are normalized to percentages summing to
#' 100 across the two predictors in that bivariate model. This produces a
#' permutation-based importance metric analogous in spirit to Java MaxEnt's
#' permutation importance, while remaining fully native to R.
#'
#' Convergence requires both tests to pass:
#' \itemize{
#'   \item \emph{Mean stability} - relative change in top-k mean PI below an
#'         adaptive threshold (\code{10\%} if \eqn{p \le 20},
#'         \code{7.5\%} if \eqn{21 \le p \le 40}, else \code{5\%}).
#'   \item \emph{Membership stability} - Jaccard overlap of consecutive
#'         top-k sets above \code{set_threshold}.
#' }
#' Adaptive top-k: \eqn{k=\min(12,\max(5,\mathrm{round}(\sqrt{p})))}. Batch
#' sampling defaults to \code{max(1, floor(1.5 * sqrt(p)))} appearances per
#' variable and is halved in \code{fast} mode. Parallel execution uses
#' deterministic RNG streams.
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Variable ranking CSV
#'   \item Convergence trace CSV
#'   \item Log entry appended to \code{_log.txt}
#' }
#'
#' If \code{year = NULL}, all 5-year intervals 1980-2020 are processed.
#'
#' @param alpha_code Character. Species or run code (e.g., "CASP").
#' @param year Integer, Character, or NULL. Year to process; NULL runs all
#'   5-year intervals 1980-2020.
#' @param max_iter Integer. Maximum iterations safeguard; default 50.
#' @param batch_samples Integer or NULL. Target appearances per variable per
#'   iteration; NULL triggers an adaptive fallback.
#' @param w_stability Numeric. Instability weight in \code{wz_score};
#'   default 0.5.
#' @param background_n Integer. Background points per model; default 2500.
#' @param seed Integer or NULL. RNG seed; NULL draws a random seed that is
#'   returned.
#' @param ncores Integer. CPU cores for parallel fits; default
#'   \code{max(1, parallel::detectCores() - 1)}.
#' @param fast Logical. Enable fast triage mode; default FALSE.
#' @param set_threshold Numeric. Jaccard threshold for membership stability;
#'   default 0.90. Must lie in \code{(0, 1)}.
#' @param passes_required Integer or NULL. Consecutive passes needed for
#'   convergence; NULL defaults to 1 in fast mode, else 2.
#'
#' @return A list (invisible) with:
#' \describe{
#'   \item{\code{importance_summary}}{Data frame of variable-wise PI statistics
#'         and rank diagnostics.}
#'   \item{\code{convergence_trace}}{Data frame with per-iteration diagnostics.}
#'   \item{\code{output_files}}{Named vector with paths to the two CSV files.}
#'   \item{\code{seed_used}}{Integer seed actually used in this run.}
#' }
#' Side effects: two CSV files and a log block are written to disk inside the
#' project directory.
#'
#' @importFrom maxnet maxnet maxnet.formula
#' @importFrom raster stack subset extract nlayers extent sampleRandom
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterSetRNGStream
#' @importFrom dplyr %>% group_by summarise mutate arrange select if_else
#' @importFrom stats median IQR sd quantile predict
#' @importFrom utils read.csv write.csv flush.console
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' \dontrun{
#' ## Quick triage across all years
#' res <- screen_by_convergence2("CASP", fast = TRUE)
#'
#' ## Standard fidelity for a single year
#' res <- screen_by_convergence2("CASP", 1990, seed = 42)
#' }
#'
#' @seealso
#' \code{\link[maxnet]{maxnet}},
#' \code{\link[raster]{stack}}
#'
#' @export
screen_by_convergence2 <- function(alpha_code,
                                  year = NULL,
                                  max_iter = 50,
                                  batch_samples = NULL,
                                  w_stability = 0.5,
                                  background_n = 2500,
                                  seed = NULL,
                                  ncores = max(1, parallel::detectCores() - 1),
                                  fast = FALSE,
                                  set_threshold = 0.90,
                                  passes_required = NULL) {

  ## ---- multi-year wrapper (year = NULL) ------------------------------------
  if (is.null(year)) {
    years <- seq(1980, 2020, by = 5)
    res <- lapply(years, function(y) {
      message(sprintf("\n>> Processing %s (%d)%s",
                      alpha_code, y, if (fast) " [fast]" else ""))
      screen_by_convergence2(
        alpha_code, y, max_iter, batch_samples,
        w_stability, background_n, seed, ncores, fast, set_threshold, passes_required
      )
    })
    return(invisible(list(
      all_years = years,
      output_files = vapply(res, function(x) x$output_files$ranking, character(1L))
    )))
  }

  t_start <- Sys.time()

  ## --- seed handling --------------------------------------------------------
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    message(sprintf("[screen_by_convergence2] Using randomly chosen seed = %d", seed))
  } else {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be NULL or a finite numeric scalar.", call. = FALSE)
    }
    seed <- as.integer(seed)
    message(sprintf("[screen_by_convergence2] Using user-supplied seed = %d", seed))
  }

  ## --- require packages (no library calls) ---------------------------------
  req_pkgs <- c("maxnet", "raster", "doParallel", "foreach", "dplyr")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         ". Please install them before running.")
  }

  ## --- validate inputs ------------------------------------------------------
  if (missing(alpha_code) || !is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code))
    stop("`alpha_code` must be a non-empty character scalar.")
  if (!is.finite(max_iter) || max_iter < 1) stop("`max_iter` must be >= 1.")
  if (!is.null(batch_samples) && (!is.finite(batch_samples) || batch_samples < 1))
    stop("`batch_samples` must be NULL or a positive integer.")
  if (!is.finite(w_stability)) stop("`w_stability` must be finite.")
  if (!is.finite(background_n) || background_n < 1) stop("`background_n` must be a positive integer.")
  if (!is.numeric(set_threshold) || set_threshold < 0 || set_threshold > 1)
    stop("`set_threshold` must be in [0, 1].")
  if (!is.null(passes_required) && (!is.finite(passes_required) || passes_required < 1))
    stop("`passes_required` must be NULL or a positive integer.")

  ## --- paths ----------------------------------------------------------------
  project_dir <- rENM_project_dir()
  occ_path <- file.path(project_dir, "runs", alpha_code, "_occs", sprintf("of-%s.csv", year))
  vars_dir <- file.path(project_dir, "runs", alpha_code, "_vars", as.character(year))
  out_dir  <- file.path(project_dir, "runs", alpha_code, "TimeSeries", as.character(year), "vars")
  ranking_fn <- file.path(out_dir, paste0("_", alpha_code, "-", year, "-Variable-Ranking.csv"))
  conv_trace_fn <- file.path(out_dir, paste0("_", alpha_code, "-", year, "-Convergence-Trace.csv"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(occ_path)) stop("Occurrences file not found: ", occ_path)
  if (!dir.exists(vars_dir))  stop("Variables directory not found: ", vars_dir)

  ## --- data -----------------------------------------------------------------
  set.seed(seed)
  occ_df <- read.csv(occ_path)
  if (!all(c("longitude", "latitude") %in% names(occ_df)))
    stop("Occurrence file must contain columns: longitude, latitude.")
  occ <- occ_df[, c("longitude", "latitude")]

  env_files <- list.files(vars_dir, pattern = "\\.tif$", full.names = TRUE)
  if (length(env_files) < 2L) stop("Need at least 2 environmental layers (.tif) in: ", vars_dir)
  env_stack <- stack(env_files)
  names(env_stack) <- tools::file_path_sans_ext(basename(env_files))
  vars <- names(env_stack)
  p <- length(vars)

  ## --- adaptive top-k and threshold ----------------------------------------
  k_top <- min(12L, max(5L, as.integer(round(sqrt(p)))))
  conv_threshold <- if (p <= 20) 0.10 else if (p <= 40) 0.075 else 0.05

  ## --- FAST MODE tuning & passes_required (adaptive default) ---------------
  if (is.null(passes_required)) passes_required <- if (isTRUE(fast)) 1L else 2L
  if (isTRUE(fast)) {
    max_iter       <- min(max_iter, 10L)
    conv_threshold <- max(conv_threshold, 0.12)
    background_n   <- min(background_n, 1500L)
    if (missing(set_threshold) || isTRUE(all.equal(set_threshold, 0.90))) {
      set_threshold <- 0.80
    }
  }

  ## --- background (reused) --------------------------------------------------
  background_raw <- raster::sampleRandom(
    env_stack,
    size = background_n,
    xy = TRUE,
    na.rm = TRUE
  )

  if (is.null(background_raw) || nrow(background_raw) < 1L) {
    stop("Unable to sample background points from predictor stack.")
  }

  background <- background_raw[, c("x", "y"), drop = FALSE]
  colnames(background) <- c("x", "y")

  if (nrow(background) < background_n) {
    warning(sprintf(
      "Requested %d background points, but only %d non-NA cells were available.",
      background_n, nrow(background)
    ))
  }

  ## --- batch_samples fallback ----------------------------------------------
  compute_fallback <- function(p, fast) {
    base <- max(1L, as.integer(round(1.5 * sqrt(p))))
    if (isTRUE(fast)) base <- max(1L, base %/% 2)
    base
  }
  if (is.null(batch_samples)) {
    batch_samples <- compute_fallback(p, fast)
    message(sprintf("[screen_by_convergence2] using fallback batch_samples = %d (p = %d%s)",
                    batch_samples, p, if (fast) ", fast" else ""))
  } else {
    batch_samples <- as.integer(batch_samples)
    if (!is.finite(batch_samples) || batch_samples < 1L) {
      batch_samples <- compute_fallback(p, fast)
      message(sprintf("[screen_by_convergence2] provided batch_samples invalid; using fallback = %d", batch_samples))
    }
  }

  ## --- parallel (deterministic RNG for workers) ----------------------------
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, seed)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)

  ## --- helpers --------------------------------------------------------------
  # Fast AUC from ranks; y must be binary 0/1 and score numeric
  calc_auc <- function(y, score) {
    ok <- is.finite(y) & is.finite(score)
    y <- y[ok]
    score <- score[ok]

    n_pos <- sum(y == 1L)
    n_neg <- sum(y == 0L)

    if (n_pos < 1L || n_neg < 1L) return(NA_real_)

    r <- rank(score, ties.method = "average")
    sum_r_pos <- sum(r[y == 1L])
    (sum_r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  }

  # Compute MaxEnt-like permutation importance from maxnet using drop in AUC
  get_maxnet_perm_importance <- function(sub_env, occ, background, pair_vars, fast = FALSE) {
    occ_vals <- raster::extract(sub_env, occ)
    bg_vals  <- raster::extract(sub_env, background)

    occ_vals <- as.data.frame(occ_vals)
    bg_vals  <- as.data.frame(bg_vals)

    names(occ_vals) <- pair_vars
    names(bg_vals)  <- pair_vars

    occ_keep <- complete.cases(occ_vals)
    bg_keep  <- complete.cases(bg_vals)

    occ_vals <- occ_vals[occ_keep, , drop = FALSE]
    bg_vals  <- bg_vals[bg_keep, , drop = FALSE]

    if (nrow(occ_vals) < 5L || nrow(bg_vals) < 5L) {
      return(data.frame(var = pair_vars, perm_imp = NA_real_))
    }

    pres_bg <- rbind(occ_vals, bg_vals)
    pa <- c(rep(1L, nrow(occ_vals)), rep(0L, nrow(bg_vals)))

    classes <- if (isTRUE(fast)) "lq" else "lqph"
    regmult <- 1

    form <- maxnet::maxnet.formula(pa, pres_bg, classes = classes)

    model <- try(
      maxnet::maxnet(p = pa, data = pres_bg, f = form, regmult = regmult),
      silent = TRUE
    )
    if (inherits(model, "try-error")) {
      return(data.frame(var = pair_vars, perm_imp = NA_real_))
    }

    baseline_pred <- try(
      stats::predict(model, pres_bg, type = "cloglog"),
      silent = TRUE
    )
    if (inherits(baseline_pred, "try-error")) {
      return(data.frame(var = pair_vars, perm_imp = NA_real_))
    }

    baseline_auc <- calc_auc(pa, baseline_pred)
    if (!is.finite(baseline_auc)) {
      return(data.frame(var = pair_vars, perm_imp = NA_real_))
    }

    auc_drop <- vapply(pair_vars, function(v) {
      perm_data <- pres_bg
      perm_data[[v]] <- sample(perm_data[[v]], length(perm_data[[v]]), replace = FALSE)

      perm_pred <- try(
        stats::predict(model, perm_data, type = "cloglog"),
        silent = TRUE
      )
      if (inherits(perm_pred, "try-error")) return(NA_real_)

      perm_auc <- calc_auc(pa, perm_pred)
      if (!is.finite(perm_auc)) return(NA_real_)

      max(0, baseline_auc - perm_auc)
    }, numeric(1))

    if (all(!is.finite(auc_drop))) {
      perm_imp <- rep(NA_real_, length(pair_vars))
    } else {
      auc_drop[!is.finite(auc_drop)] <- 0
      s <- sum(auc_drop)
      if (s <= 0) {
        perm_imp <- rep(0, length(pair_vars))
      } else {
        perm_imp <- 100 * auc_drop / s
      }
    }

    data.frame(var = pair_vars, perm_imp = as.numeric(perm_imp))
  }

  # Summarize variable-wise PI statistics with rank diagnostics
  summarize_importance <- function(df, w_stability = 0.5) {
    if (is.null(df) || !nrow(df)) return(data.frame())
    k_median <- sqrt(pi / 2)

    summ <- df %>%
      group_by(var) %>%
      summarise(
        n_pairs   = sum(!is.na(perm_imp)),
        median_PI = median(perm_imp, na.rm = TRUE),
        IQR_PI    = IQR(perm_imp, na.rm = TRUE),
        mean_PI   = mean(perm_imp, na.rm = TRUE),
        sd_PI     = sd(perm_imp, na.rm = TRUE),
        .groups   = "drop"
      )

    m_mean <- mean(summ$median_PI, na.rm = TRUE)
    m_sd   <- sd(summ$median_PI,   na.rm = TRUE)
    i_mean <- mean(summ$IQR_PI,    na.rm = TRUE)
    i_sd   <- sd(summ$IQR_PI,      na.rm = TRUE)

    z_pi_vec <- if (is.finite(m_sd) && m_sd > 0) {
      (summ$median_PI - m_mean) / m_sd
    } else {
      rep(0, nrow(summ))
    }
    z_iqr_vec <- if (is.finite(i_sd) && i_sd > 0) {
      (summ$IQR_PI - i_mean) / i_sd
    } else {
      rep(0, nrow(summ))
    }

    summ %>%
      mutate(
        se_mean_PI   = if_else(n_pairs > 0, sd_PI / sqrt(pmax(1, n_pairs)), NA_real_),
        se_median_PI = if_else(n_pairs > 0, k_median * sd_PI / sqrt(pmax(1, n_pairs)), NA_real_),
        rank_medianPI = rank(-median_PI, ties.method = "first"),
        rank_stability = rank(IQR_PI, ties.method = "first"),
        ranksum_score  = rank_medianPI + rank_stability,
        z_PI  = z_pi_vec,
        z_IQR = z_iqr_vec,
        wz_score    = z_PI - w_stability * z_IQR,
        ratio_score = median_PI / (IQR_PI + 1e-6),
        rank_ranksum = rank(ranksum_score, ties.method = "first"),
        rank_wz      = rank(-wz_score,     ties.method = "first"),
        rank_ratio   = rank(-ratio_score,  ties.method = "first")
      ) %>%
      arrange(rank_medianPI) %>%
      select(
        var, n_pairs, median_PI, IQR_PI, mean_PI, se_mean_PI, se_median_PI,
        rank_medianPI, rank_ranksum, rank_wz, rank_ratio, z_PI, z_IQR, wz_score, ratio_score
      )
  }

  # Mean of top-k median_PI values in a summary table ordered by rank_medianPI
  topk_mean <- function(summary_df, k = 5L) {
    if (!nrow(summary_df)) return(NA_real_)
    kk <- min(k, nrow(summary_df))
    mean(summary_df$median_PI[seq_len(kk)], na.rm = TRUE)
  }

  pct_change <- function(prev, curr) {
    if (is.na(prev) || prev == 0) return(Inf)
    abs(curr - prev) / abs(prev)
  }

  # Balanced pairing sampler: each var appears ~appearances times with diverse
  # partners and self-pairs are not allowed
  balanced_pairs <- function(vars, appearances) {
    appearances <- max(1L, as.integer(appearances))
    vars <- as.character(vars)

    if (length(vars) < 2L) return(list())

    pool <- sample(rep(vars, each = appearances))
    pairs <- vector("list", 0L)

    while (length(pool) >= 2L) {
      v1 <- pool[1]
      pool <- pool[-1]

      idx <- which(pool != v1)

      if (!length(idx)) break

      j <- sample(idx, 1L)
      v2 <- pool[j]
      pool <- pool[-j]

      pairs[[length(pairs) + 1L]] <- c(v1, v2)
    }

    pairs
  }

  ## --- iterative convergence loop (mean + membership) -----------------------
  results_all <- NULL
  prev_topk_mean <- NA_real_
  prev_topk_set  <- character(0)
  passes_mean <- 0L
  passes_set  <- 0L
  conv_trace <- data.frame()
  cumulative_pairs <- 0L
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    appear_per_var <- max(1L, as.integer(ceiling(batch_samples)))
    if (isTRUE(fast)) appear_per_var <- max(1L, as.integer(ceiling(0.5 * appear_per_var)))

    set.seed(seed + iter)
    if (isTRUE(fast)) {
      pairings  <- balanced_pairs(vars, appearances = appear_per_var)
      replicates <- length(pairings)
      if (length(pairings) == 0L) {
        pairings <- list(sample(vars, 2))
        replicates <- 1L
      }
    } else {
      replicates <- max(1L, ceiling((appear_per_var * p) / 2))
      pairings <- replicate(replicates, sample(vars, 2), simplify = FALSE)
    }

    # --- fit & extract permutation importance -------------------------------
    batch_res <- foreach(
      p = pairings,
      .combine = rbind,
      .packages = c("maxnet", "raster", "stats")
    ) %dopar% {
      sub_env <- raster::subset(env_stack, p)
      get_maxnet_perm_importance(
        sub_env = sub_env,
        occ = occ,
        background = background,
        pair_vars = p,
        fast = fast
      )
    }

    # accumulate
    results_all <- if (is.null(results_all)) batch_res else rbind(results_all, batch_res)
    cumulative_pairs <- cumulative_pairs + replicates

    # summarize cumulatively and check convergence
    summary_now <- summarize_importance(results_all, w_stability = w_stability)

    ord <- order(summary_now$rank_medianPI)
    ordered_summary <- summary_now[ord, , drop = FALSE]
    curr_topk_mean <- topk_mean(ordered_summary, k = k_top)
    top_vars <- ordered_summary$var[seq_len(min(k_top, nrow(ordered_summary)))]

    # mean-stability pass
    change <- pct_change(prev_topk_mean, curr_topk_mean)
    if (iter >= 2 && is.finite(change) && change < conv_threshold) {
      passes_mean <- passes_mean + 1L
    } else if (iter >= 2) {
      passes_mean <- 0L
    }

    # membership-stability pass (Jaccard of consecutive top-k sets)
    jaccard <- if (length(prev_topk_set) && length(top_vars)) {
      inter <- length(intersect(top_vars, prev_topk_set))
      union <- length(union(top_vars, prev_topk_set))
      if (union > 0) inter / union else 0
    } else 0
    if (iter >= 2 && jaccard >= set_threshold) {
      passes_set <- passes_set + 1L
    } else if (iter >= 2) {
      passes_set <- 0L
    }

    passes <- min(passes_mean, passes_set)
    converged <- (iter >= 2) && (passes >= passes_required)

    conv_trace <- rbind(conv_trace, data.frame(
      iteration = iter,
      batch_pairs = replicates,
      cumulative_pairs = cumulative_pairs,
      topk = k_top,
      topk_mean_PI = curr_topk_mean,
      pct_change = if (is.finite(change)) change else NA_real_,
      jaccard = jaccard,
      passes_mean = passes_mean,
      passes_set = passes_set,
      passes_required = passes_required,
      threshold = conv_threshold,
      set_threshold = set_threshold,
      batch_samples = batch_samples,
      fast = fast,
      converged = converged
    ))

    if (converged) break
    prev_topk_mean <- curr_topk_mean
    prev_topk_set  <- top_vars
  }

  importance_summary <- summarize_importance(results_all, w_stability = w_stability)

  ## --- write ranking CSV & convergence trace -------------------------------
  message("[screen_by_convergence2] Writing updated files:")
  message("  Ranking CSV  -> ", normalizePath(ranking_fn, mustWork = FALSE))
  message("  Trace   CSV  -> ", normalizePath(conv_trace_fn, mustWork = FALSE))
  if (file.exists(ranking_fn)) file.remove(ranking_fn)
  if (file.exists(conv_trace_fn)) file.remove(conv_trace_fn)

  dir.create(dirname(ranking_fn), recursive = TRUE, showWarnings = FALSE)

  if (!nrow(conv_trace)) {
    conv_trace <- data.frame(
      iteration = integer(0), batch_pairs = integer(0), cumulative_pairs = integer(0),
      topk = integer(0), topk_mean_PI = numeric(0), pct_change = numeric(0),
      jaccard = numeric(0), passes_mean = integer(0), passes_set = integer(0),
      passes_required = integer(0), threshold = numeric(0), set_threshold = numeric(0),
      batch_samples = integer(0), fast = logical(0), converged = logical(0)
    )
  }

  write.csv(importance_summary, ranking_fn, row.names = FALSE)
  write.csv(conv_trace,      conv_trace_fn, row.names = FALSE)

  flush.console()
  Sys.sleep(0.2)

  ## --- eBird-style processing summary --------------------------------------
  t_end <- Sys.time()
  s <- as.numeric(difftime(t_end, t_start, units = "secs"))
  fmt_elapsed <- sprintf("%02d:%02d:%02d", s %/% 3600, (s %% 3600) %/% 60, round(s %% 60))
  log_file <- file.path(project_dir, "runs", alpha_code, "_log.txt")
  sep_line <- paste(rep("-", 72), collapse = "")
  ts_str   <- format(t_end, "%Y-%m-%d %H:%M:%S %Z")
  header   <- "Processing summary (screen_by_convergence2)"
  f <- function(key, val) sprintf("%-18s : %s", key, val)

  last <- if (nrow(conv_trace)) conv_trace[nrow(conv_trace), ] else NULL
  conv_line <- if (!is.null(last)) {
    sprintf("Iter %d | pairs(cum)=%d | topk=%d | topk_mean=%.4f | d%%=%.3f | J=%.2f | passes(mean/set)=%d/%d (need %d) | thr=%.3f | set_thr=%.2f | batch=%d | %s%s",
            last$iteration, last$cumulative_pairs, k_top, last$topk_mean_PI,
            ifelse(is.na(last$pct_change), NA_real_, last$pct_change),
            last$jaccard, last$passes_mean, last$passes_set, last$passes_required,
            conv_threshold, set_threshold, batch_samples,
            if (isTRUE(last$converged)) "CONVERGED" else "NOT CONVERGED",
            if (fast) " [fast]" else "")
  } else "No iterations recorded"

  if (nrow(conv_trace)) {
    ct_log <- conv_trace
    num_cols <- c("topk_mean_PI", "pct_change", "threshold", "jaccard", "set_threshold")
    ct_log[, num_cols] <- lapply(ct_log[, num_cols, drop = FALSE], function(x) round(x, 4))
    ct_log$fast <- ifelse(ct_log$fast, "Y", "N")
    fmt_row <- function(r) sprintf(
      "%3d | %5d | %6d | %3d | %7.4f | %7.4f | %4.2f | %1d/%1d | %1d | %7.4f | %4.2f | %5d | %1s | %s",
      r$iteration, r$batch_pairs, r$cumulative_pairs, r$topk,
      r$topk_mean_PI, r$pct_change, r$jaccard,
      r$passes_mean, r$passes_set, r$passes_required,
      r$threshold, r$set_threshold, r$batch_samples, r$fast, r$converged
    )
    ct_lines <- vapply(split(ct_log, seq_len(nrow(ct_log))), fmt_row, character(1))
  } else {
    ct_lines <- character(0)
  }

  outputs <- c(ranking_fn, conv_trace_fn)

  log_block <- c(
    "",
    sep_line,
    header,
    f("Timestamp",         ts_str),
    f("Alpha code",        alpha_code),
    f("Year",              year),
    f("Mode",              if (fast) "FAST" else "STANDARD"),
    f("Seed used",         seed),
    f("Top-k (adaptive)",  k_top),
    f("Threshold",         sprintf("%.3f", conv_threshold)),
    f("Set threshold",     sprintf("%.2f", set_threshold)),
    f("Passes required",   passes_required),
    f("Max iterations",    max_iter),
    f("Batch samples",     batch_samples),
    f("w_stability",       w_stability),
    f("Background_N",      background_n),
    f("Cores used",        ncores),
    f("Convergence",       conv_line),
    f("Outputs saved",     length(outputs)),
    vapply(outputs, function(x) paste0(" - ", x), character(1)),
    "Convergence trace (iter | pairs | cum_pairs | topk | topk_mean | d% | J | passes m/s | need | thr | set_thr | batch | fast | done):",
    if (length(ct_lines)) paste0("  ", ct_lines) else "  (no iterations recorded)",
    f("Total elapsed",     fmt_elapsed),
    ""
  )

  try({
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(paste0(log_block, collapse = "\n"), file = log_file, append = TRUE)
  }, silent = TRUE)

  invisible(list(
    importance_summary = importance_summary,
    convergence_trace  = conv_trace,
    output_files       = list(
      ranking = ranking_fn,
      convergence_trace = conv_trace_fn
    ),
    seed_used = seed
  ))
}
