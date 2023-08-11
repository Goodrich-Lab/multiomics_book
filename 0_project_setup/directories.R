# directories

# home directory for project
dir_home <- here::here() 

# Original Data directory
dir_data_og <- dir_home %>%
  dirname() %>%
  fs::path("multiomics_methods_in_HELIX", "0_data")

# Mercury Data directory
dir_data_hg <- fs::path(dir_home, "0_data")

# results folder
dir_results <- fs::path(dir_home, "2_results")

# figures 
dir_figs <- fs::path(dir_home, "3_figures")

# functions
dir_fxn <- fs::path(here::here(), 
                    "0_project_setup",
                    "functions")
rm(dir_home)
