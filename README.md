# rENM.model

![rENM](https://img.shields.io/badge/rENM-framework-blue) ![module](https://img.shields.io/badge/module-model-informational)[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20797840.svg)](https://doi.org/10.5281/zenodo.20797840)

**Modeling and reconstruction engine for the rENM Framework**

## Overview

`rENM.model` implements the core ecological niche modeling and historical reconstruction workflows within the rENM Framework. It transforms standardized occurrence and environmental data into time-resolved estimates of climatic suitability.

This package depends on `rENM.core` for project-directory resolution and species metadata access. Functions find the project directory through `rENM.core::rENM_project_dir()`; see `?rENM_project_dir` for configuration options.

## Key functions

| Function | Description |
|------------------------------------|------------------------------------|
| `stage_occurrences()` | Copy occurrence CSVs into TimeSeries bins |
| `stage_all_variables()` | Copy predictor rasters into TimeSeries bins |
| `screen_by_convergence1()` | Convergence-based variable screening via dismo MaxEnt (requires Java) |
| `screen_by_convergence2()` | Convergence-based variable screening via native R maxnet (no Java dependency) |
| `reduce_covariance()` | Remove collinear predictors via adaptive VIF screening |
| `stage_screened_variables()` | Copy ranked predictors into TimeSeries bins |
| `create_ensemble_model()` | Fit an ensemble ENM on land cells for a single species and time bin |
| `create_timeseries()` | Run `create_ensemble_model()` across all time bins in parallel |
| `create_range_map()` | Produce a binary presence-absence range map |
| `plot_suitability()` | Plot a continuous climatic suitability raster |
| `save_suitability_plot()` | Save a suitability ggplot to disk |
| `rank_variable_importance()` | Parse and rank variable importance from an SDM report |

## Installation

``` r
# From GitHub
remotes::install_github("rENM-Framework/rENM.model")

# From a local source directory
remotes::install_local("rENM.model")
```

## Getting started

Set up a project directory and preprocess occurrence and predictor data first (see `rENM.data`), then run the modeling pipeline in order:

``` r
library(rENM.model)

# set once per session, or set RENM_PROJECT_DIR in ~/.Renviron
options(rENM.project_dir = "/path/to/your/rENM/project")

# 1. Stage occurrence records into TimeSeries bins
stage_occurrences("CASP")

# 2. Screen variables and stage the selected set
#    (or stage_all_variables("CASP") to stage every variable unscreened)
screen_by_convergence2("CASP", seed = 42)
stage_screened_variables("CASP")

# 3. Optionally remove strongly collinear predictors
# reduce_covariance("CASP")

# 4. Fit the ensemble models for all nine bins
create_timeseries("CASP", seed = 42)
```

## Modeling pipeline

```         
stage_occurrences()
        ↓
screen_by_convergence2()  or  screen_by_convergence1()  (requires Java)
        ↓
stage_screened_variables()      or  stage_all_variables()  (no screening)
        ↓
reduce_covariance()             (optional)
        ↓
create_timeseries()             ← runs create_ensemble_model() for every bin
```

Variable screening reads from the run-level `_occs/` and `_vars/` directories. Staging functions copy files into `TimeSeries/<year>/occs/` and `TimeSeries/<year>/vars/`. Model outputs are written to `TimeSeries/<year>/model/`.

## Role in the rENM Framework

`rENM.model` is the third stage in the pipeline:

```         
rENM.core → rENM.data → rENM.model → rENM.analysis → rENM.ai → rENM.reports
```

It consumes the run directory structure and preprocessed data produced by `rENM.data` and generates the modeled suitability surfaces and range maps consumed by `rENM.analysis`.

## License

See `LICENSE` for details.

------------------------------------------------------------------------

**rENM Framework** — A modular system for reconstructing and analyzing long-term ecological niche dynamics.
