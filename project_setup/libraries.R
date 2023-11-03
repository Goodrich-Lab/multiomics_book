## ---- install_packages ----
# Download EnvirOmix package if not already installed:
if(!requireNamespace("EnvirOmix", quietly = TRUE)){
  if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("Goodrich-Lab/EnvirOmix")
}

# List all CRAN packages used in this book:
cran_packages <- c("tidyverse",
                   "tools",
                   "parallel",
                   "boot",
                   "table1",
                   "BiocManager",
                   "ggplot2",
                   "cowplot",
                   "ggh4x",
                   "LUCIDus",
                   "networkD3",
                   "plotly",
                   "htmlwidgets",
                   "jsonlite")

# Install cran packages if not already installed:
for(package_name in cran_packages){
  if(!requireNamespace(package_name, quietly = TRUE)){
    install.packages(package_name, repos = "http://cran.r-project.org/")
  }
}

# List all Bioconductor packages used in this book:
bioconductor_packages <- c("BiocGenerics", 
                           "ComplexHeatmap", 
                           "IRanges",
                           "S4Vectors")

# Install bioconductor packages if not already installed:
for(package_name in bioconductor_packages){
  if(!requireNamespace(package_name, quietly = TRUE)){
    BiocManager::install(package_name)
  }
}

## ---- load_packages ----
# General Packages:
library(EnvirOmix)
library(tidyverse)
library(table1)
# Load plotting packages:
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(ggh4x)
# Packages for Quasi-mediation:
#LUCIDus Version 3.0.1 is required for this analysis:
library(LUCIDus)
library(networkD3)
library(plotly)
library(htmlwidgets)
library(jsonlite)

## ---- set_theme ----
ggplot2::theme_set(cowplot::theme_cowplot())
