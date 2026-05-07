# rENM.model

![rENM](https://img.shields.io/badge/rENM-framework-blue)
![module](https://img.shields.io/badge/module-model-informational)

**Modeling and reconstruction engine for the rENM Framework**

## Overview

`rENM.model` implements the core ecological niche modeling and historical
reconstruction workflows within the rENM Framework. It transforms standardized
occurrence and environmental data into time-resolved estimates of climatic
suitability.

This package depends on `rENM.core` for project-directory resolution and
species metadata access. All functions accept an optional `project_dir`
argument; see `?rENM_project_dir` for configuration options.

## Key functions

| Function | Description |
|---|---|
| `stage_occurrences()` | Copy occurrence CSVs into TimeSeries bins |
| `stage_all_variables()` | Copy predictor rasters into TimeSeries bins |
| `screen_by_convergence1()` | Convergence-based variable screening via dismo MaxEnt (requires Java) |
| `screen_by_convergence2()` | Convergence-based variable screening via native R maxnet (no Java dependency) |
| `reduce_covariance()` | Remove collinear predictors via adaptive VIF screening |
| `stage_screened_variables()` | Copy ranked predictors into TimeSeries bins |
| `create_ensemble_model()` | Fit an ensemble ENM for a single species and time bin |
| `create_timeseries()` | Run `create_ensemble_model()` across all time bins in parallel |
| `create_range_map()` | Produce a binary presence-absence range map |
| `plot_suitability()` | Plot a continuous climatic suitability raster |
| `save_suitability_plot()` | Save a suitability ggplot to disk |
| `rank_variable_importance()` | Parse and rank variable importance from an SDM report |

## Installation

```r
# From GitHub
devtools::install_github("rENM-Framework/rENM.model")

# From a local source directory
devtools::install_local("rENM.model")
```

## Getting started

Set up a project directory and preprocess occurrence and predictor data first
(see `rENM.data`), then run the modeling pipeline in order:

```r
library(rENM.model)

proj <- "/path/to/your/rENM/project"

# 1. Stage occurrence and predictor data into TimeSeries bins
stage_occurrences("CASP", project_dir = proj)
stage_all_variables("CASP", project_dir = proj)

# 2. Screen variables (use convergence2 for a Java-free workflow)
screen_by_convergence2("CASP", project_dir = proj)

# 3. Remove collinear predictors and finalize variable sets
reduce_covariance("CASP", project_dir = proj)
stage_screened_variables("CASP", project_dir = proj)

# 4. Fit models and generate the full time series
create_ensemble_model("CASP", year = 2000, project_dir = proj)
create_timeseries("CASP", project_dir = proj)
```

For interactive work, configure the project directory once per session to
avoid passing it to every function:

```r
options(rENM.project_dir = "/path/to/your/rENM/project")

stage_occurrences("CASP")
stage_all_variables("CASP")
# ...
```

## Modeling pipeline

```
stage_occurrences()
stage_all_variables()
        ↓
screen_by_convergence1()  or  screen_by_convergence2()  (Java-free)
        ↓
reduce_covariance()
        ↓
stage_screened_variables()
        ↓
create_ensemble_model()   ← single time bin
        ↓
create_timeseries()       ← all bins in parallel
```

Variable screening reads from the run-level `_occs/` and `_vars/` directories.
Staging functions copy files into `TimeSeries/<year>/occs/` and
`TimeSeries/<year>/vars/`. Model outputs are written to
`TimeSeries/<year>/model/`.

## Role in the rENM framework

`rENM.model` is the third stage in the pipeline:

```
rENM.core → rENM.data → rENM.model → rENM.analysis → rENM.ai → rENM.reports
```

It consumes the run directory structure and preprocessed data produced by
`rENM.data` and generates the modeled suitability surfaces and range maps
consumed by `rENM.analysis`.

## License

See `LICENSE` for details.

---

**rENM Framework** — A modular system for reconstructing and analyzing
long-term ecological niche dynamics.
