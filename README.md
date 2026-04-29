# rENM.model

**Modeling and reconstruction engine for the rENM Framework**

## Overview

`rENM.model` implements the core ecological niche modeling and historical reconstruction workflows within the rENM Framework. It transforms standardized occurrence and environmental data into time-resolved estimates of climatic suitability.

This package focuses on **model construction, variable screening, and temporal reconstruction**.

## Role in the rENM Framework

Within the modular rENM ecosystem, `rENM.model`: - Prepares staged inputs for modeling workflows - Performs **variable screening and selection** - Builds **ensemble ecological niche models** - Generates **time series of climatic suitability** - Produces spatial outputs (e.g., range maps, suitability surfaces)

It is the computational core that converts prepared data into modeled niche dynamics.

## Key Functions

-   `stage_occurrences()` — Prepare occurrence data for modeling
-   `stage_all_variables()` — Assemble full predictor sets
-   `screen_by_convergence1()` / `screen_by_convergence2()` — Monte Carlo-based variable screening
-   `reduce_covariance()` — Remove collinearity among predictors
-   `stage_screened_variables()` — Finalize selected variable sets
-   `create_ensemble_model()` — Build ensemble ENM models
-   `create_timeseries()` — Generate temporal suitability outputs
-   `create_range_map()` — Produce spatial range estimates
-   `plot_suitability()` / `save_suitability_plot()` — Visualization utilities
-   `rank_variable_importance()` — Evaluate predictor contributions

## Installation

``` r
devtools::install_local("rENM.model")
```

## Example

``` r
library(rENM.model)

# stage inputs
stage_occurrences("CASP")
stage_all_variables("CASP")

# screen variables
screen_by_convergence1("CASP")
screen_by_convergence2("CASP")

# build model
create_ensemble_model("CASP")

# generate time series
create_timeseries("CASP")
```

## Relationship to Other Packages

`rENM.model` provides the modeled outputs used by all downstream analysis and reporting layers.

## License

See `LICENSE` for details.

------------------------------------------------------------------------

**rENM Framework**\
A modular system for reconstructing and analyzing long-term ecological niche dynamics.
