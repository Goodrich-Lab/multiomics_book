## ---- install_packages ----
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
                   "HIMA",
                   "xtune",
                   "RMediation",
                   "glmnet",
                   "r.jive",
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
                           "qvalue", 
                           "S4Vectors")

# Install bioconductor packages if not already installed:
for(package_name in bioconductor_packages){
  if(!requireNamespace(package_name, quietly = TRUE)){
    BiocManager::install(package_name)
  }
}

## ---- load_packages ----

# General Packages:
library(tidyverse) 
# library(tools) 
library(parallel) 
library(boot)  
library(table1) 
# Load plotting packages:
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
library(networkD3)
library(plotly)
library(htmlwidgets)
library(jsonlite)

## ---- set_theme ----
ggplot2::theme_set(cowplot::theme_cowplot())
