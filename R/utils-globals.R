utils::globalVariables(c(
  ## common ggplot / terra column names
  "x", "y", "value",
  ## dplyr / tidy eval columns (screen_by_convergence1/2)
  "var", "perm_imp",
  "n_pairs", "sd_PI", "median_PI", "IQR_PI",
  "rank_medianPI", "rank_stability",
  "z_PI", "z_IQR",
  "ranksum_score", "wz_score", "ratio_score",
  "mean_PI", "se_mean_PI", "se_median_PI",
  "rank_ranksum", "rank_wz", "rank_ratio",
  ## misc
  "trend", ".data"
))

#' @importFrom rENM.core rENM_project_dir get_species_info show_species show_variables
NULL
