# Load packages

## ---- general_packages_loaded ----
# Tidyverse packages, including those from dplyr and ggplot2
library(tidyverse)
library(tools)
library(parallel)
library(epiomics)
library(cowplot) # For combining figures
library(ComplexHeatmap) # Additional Plotting Functions
library(ggh4x) # For faceted plots
library(data.table) # For function setnames
library(colorspace)
library(tidytext)
library(table1)

# Not loaded
# library(fs)
# library(here)


## ---- hima_packages ----
# For running High Dimensional Mediation Analysis:
library(qvalue)
library(HIMA)
# For group lasso in high dimensional mediation/intermediate integration:
library(xtune)
# Mediated effect for intermediate integration: 
library(RMediation) # devtools::install_github("cran/RMediation")
# Bootstrap confidence intervals: 
library(boot) 
library(glmnet)


## ---- med_lf_packages ----
# Functions for JIVE 
library(r.jive)


## ---- lucid_packages ----
# For running LUCID analysis:
library(LUCIDus)
# Functions for LUCID in parallel:
library(mclust)
# For Plotting Sankey Diagrams:
library(networkD3)
library(plotly)
library(htmlwidgets)
library(glmnet)
library(glasso)
library(nnet)
library(boot)
library(progress)
library(jsonlite)

## ---- set_theme ----
ggplot2::theme_set(cowplot::theme_cowplot())