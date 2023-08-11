# Load Data

# Read in reduced full scaled data
helix_dat_reduced <- read_rds(fs::path(dir_data_hg,
                                       "simulate_Hg_ck18_screened_scaled_omics.RDS")) 

# Get exposure, outcome and covariate data ----

## Define exposure and outcome name ----
exposure_name <- "hs_hg_m_scaled"
outcome_name  <- "ck18_scaled"
covars <- c("e3_sex_None", 
            "hs_child_age_yrs_None", 
            "h_fish_preg_Ter")

## Extract exposure and outcome data----
outcomes <- helix_dat_reduced[["phenotype"]]
exposure <- as.numeric(helix_dat_reduced[["phenotype"]][[exposure_name]]) 
outcome <- as.numeric(helix_dat_reduced[["phenotype"]][[outcome_name]]) 

## Get matrix of covariates ----
covs <- helix_dat_reduced[["phenotype"]][covars] %>% 
  mutate(h_fish_preg_Ter = as.numeric(h_fish_preg_Ter)) %>%
  droplevels() %>% 
  fastDummies::dummy_cols(remove_selected_columns = TRUE,
                          remove_first_dummy = TRUE) 

# Omics data-----
## create list of omics data------
omics_name <- c("methylome", "transcriptome","proteome", "miRNA","metabolome")
omics_lst <- helix_dat_reduced[which(names(helix_dat_reduced) %in% omics_name)]

# Change omics list elements to dataframes
# Originally "omics"
omics_lst_df <- purrr::map(omics_lst, ~as_tibble(.x, rownames = "name"))

# Create data frame of omics data
# Originally "omics_comb"
omics_df <- omics_lst_df  %>%
  purrr::reduce(left_join, by = "name") %>%
  column_to_rownames("name")

## Omics annotations -------------
omics_names <- readRDS(fs::path(dir_data_hg, "feature_metadata_v2.RDS"))

# Set Color Palettes ----
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

# Color pallet for 
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

