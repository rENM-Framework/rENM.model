# rENM.model 0.2.0.9000

* Fixed `screen_by_convergence2()` selecting different variables on repeated
  runs given the same seed. `clusterSetRNGStream()` makes each worker's RNG
  stream reproducible, but `%dopar%` does not guarantee which worker runs
  which task, so a given variable pairing could draw from a different stream
  between runs. Permutation importance is sampled inside the workers, so this
  changed the importance scores, the variable ranking, and ultimately which
  predictors entered the model — two seeded Pinyon Jay runs selected `bio8`
  versus `lai` for 1980, and differed in every other year as well. The loop
  now registers `doRNG`, which ties the stream to the iteration rather than
  the worker, making results independent of task scheduling and of `ncores`.
  Adds `doRNG` to Imports. `screen_by_convergence1()` shares the pattern but
  is not reached by the pipeline and is unchanged; it would also need its
  MaxEnt `randomseed=true` argument addressed, since that randomization lives
  in Java and is outside R's control.
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
