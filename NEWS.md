# rENM.model NEWS

## rENM.model 0.1.1.9000 (development)

### Bug fixes

- Fixed `"invalid 'scipen'"` error that caused all `create_timeseries()` workers to fail under R 4.6.0. R 4.6.0 tightened validation of the `scipen` option: `options(scipen = <list>)` now throws an error instead of silently coercing (see R 4.6.0 NEWS). The bug is in the `writeValues` S4 method for `signature("RasterLayer", "vector")` in the `raster` package: `opsci = options("scipen")` saves scipen as a named list (`list(scipen=0)`), not a scalar, so the subsequent `options(scipen = opsci)` restore call fails in R 4.6.0. This code path is hit for every raster write in `.asc` format, which `sdm::sdm()` triggers internally when writing model predictions during training. A compatibility patch is applied inline at the top of `create_ensemble_model()`, immediately after `raster` is loaded via `must_have_pkg()`, using `methods::setMethod()` to replace the broken method. Placing it there ensures it runs in both the main R process and the PSOCK workers spawned by `create_timeseries()` — workers load `raster` directly via `library()` and would not benefit from a package `.onLoad` hook. The fix replaces `opsci = options("scipen")` with `opsci = getOption("scipen")` so the saved value is a scalar and the restore call is valid. The patch is guarded and will silently do nothing once `raster` is fixed upstream. (#R46-scipen)

- `Hmisc` removed from `Imports`. It is an internal dependency of `sdm` and is never called directly by `rENM.model`; declaring it caused `"Namespace in Imports field not imported from: 'Hmisc'"` in `R CMD check`.

- `randomForest`, `gbm`, and `earth` moved from `Imports` to `Suggests`. These are `sdm` method-driver packages loaded at runtime via `must_have_pkg()`. Because they are loaded with `library()` rather than referenced via `::` or `@importFrom`, declaring them in `Imports` caused `"Namespaces in Imports field not imported from"` in `R CMD check`. `must_have_pkg()` already checks for their presence and stops with a clear message if any are missing.
