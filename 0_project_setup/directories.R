# directories

# home directory for project
dir_home <- here::here() 

# Mercury Data directory
dir_data_hg <- fs::path(dir_home, "0_data")

# results folder
dir_results <- fs::path(dir_home, "2_results")

# figures 
dir_figs <- fs::path(dir_home, "figures")

# functions
dir_fxn <- fs::path(here::here(), 
                    "0_project_setup",
                    "functions")
rm(dir_home)
