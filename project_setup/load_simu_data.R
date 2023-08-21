## ---- load_data ----

# Load simulated data
simulated_data <- read_rds(fs::path(dir_data_hg, "simulated_HELIX_data_2.RDS")) 

# Define exposure and outcome name
covars <- c("e3_sex_None", "hs_child_age_yrs_None")

# Extract exposure and outcome data
# outcomes <- simulated_data[["phenotype"]]
exposure <- simulated_data[["phenotype"]]$hs_hg_m_scaled
outcome  <- simulated_data[["phenotype"]]$ck18_scaled

# Get numeric matrix of covariates 
covs <- simulated_data[["phenotype"]][covars] 

# create list of omics data 
omics_lst <- simulated_data[-which(names(simulated_data) == "phenotype")]

# Create data frame of omics data
omics_df <- omics_lst %>% 
  purrr::map(~as_tibble(.x, rownames = "name")) %>%
  purrr::reduce(left_join, by = "name") %>%
  column_to_rownames("name")


## ---- set_color_pal ----

# Set Color Palettes 
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

# Color pallet for sankey diagrams
sankey_colors <- matrix(c("exposure", col_pal[6],
                          "lc1",      col_pal[1],
                          "lc2",      col_pal[2],
                          "lc3",      col_pal[3],
                          "lc4",      col_pal[4],
                          "layer1",   col_pal[1],
                          "layer2",   col_pal[2],
                          "layer3",   col_pal[3],
                          "layer4",   col_pal[4],
                          "Outcome",  col_pal[8],
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red", 
                          "neg_clus_to_out", "#e4e5f2"), 
                        byrow = TRUE, nrow = 14)

# Change to dataframe
colnames(sankey_colors) <- c("domain", "range")
sankey_colors <- as_tibble(sankey_colors)

# Assign colors to omics layers
annotation_colors <- list(Type = c(Methylation = sankey_colors$range[2], 
                                   Transcriptome = sankey_colors$range[5], 
                                   miRNA = sankey_colors$range[4], 
                                   Proteins = sankey_colors$range[3],
                                   Metabolome = "blue"))