# rENM.model 0.1.0

* Initial release.
* Added `stage_occurrences()` to copy occurrence CSVs into TimeSeries bins.
* Added `stage_all_variables()` to copy predictor rasters into TimeSeries bins.
* Added `screen_by_convergence1()` for convergence-based variable screening via
  dismo MaxEnt (requires Java).
* Added `screen_by_convergence2()` for convergence-based variable screening via
  native R maxnet (no Java dependency).
* Added `reduce_covariance()` to remove collinear predictors via adaptive VIF
  screening.
* Added `stage_screened_variables()` to copy ranked predictors into TimeSeries
  bins.
* Added `create_ensemble_model()` to fit an ensemble ENM for a single species
  and time bin.
* Added `create_timeseries()` to run `create_ensemble_model()` across all time
  bins in parallel.
* Added `create_range_map()` to produce a binary presence-absence range map.
* Added `plot_suitability()` to plot a continuous climatic suitability raster.
* Added `save_suitability_plot()` to save a suitability ggplot to disk.
* Added `rank_variable_importance()` to parse and rank variable importance from
  an SDM report.
* Applied compatibility patch in `create_ensemble_model()` to fix an
  `"invalid 'scipen'"` error triggered by R 4.6.0's tightened option validation.
  The bug originates in `raster::writeValues()`, which saves `scipen` as a named
  list instead of a scalar; the patch replaces the broken `options("scipen")`
  save with `getOption("scipen")` so the restore call is valid under R 4.6.0.
