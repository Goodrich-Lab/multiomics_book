# Load packages

## ---- packages_loaded ----
# General Packages:
library(tidyverse)
library(tools)
library(parallel)
library(boot) 
library(table1)
# Packages for Plotting:
library(ggplot2)
library(cowplot) 
library(ComplexHeatmap) 
library(ggh4x)
# Packages for High Dimensional Mediation:
library(HIMA)
library(xtune)
library(RMediation)
library(glmnet)
# Packages for Mediation with Latent Factors:
library(r.jive)
# Packages for Quasi-mediation:
library(LUCIDus)
library(mclust)
library(networkD3)
library(plotly)
library(htmlwidgets)
library(glasso)
library(nnet)
library(progress)
library(jsonlite)

## ---- set_theme ----
ggplot2::theme_set(cowplot::theme_cowplot())