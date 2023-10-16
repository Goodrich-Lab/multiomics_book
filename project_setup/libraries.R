# Load packages

## ---- packages_loaded ----
# General Packages:
library(tidyverse)
library(tools)
library(parallel)
library(boot) 
library(table1)

# Plotting Packages:
library(ggplot2)
library(cowplot) 
library(ComplexHeatmap) 

# High Dimensional Mediation Packages:
library(HIMA)
library(xtune)
library(RMediation)
library(glmnet)

# Mediation with Latent Factors Packages:
library(r.jive)

# Quasi-mediation Packages:
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