## ---- load_data ----
# Load R package
library(EnvirOmix)

# Load simulated data
data(simulated_data)

# Define exposure and outcome name
exposure <- simulated_data[["phenotype"]]$hs_hg_m_scaled
outcome  <- simulated_data[["phenotype"]]$ck18_scaled

# Get numeric matrix of covariates 
covs <- simulated_data[["phenotype"]][c("e3_sex_None", "hs_child_age_yrs_None")] 
covs$e3_sex_None <- ifelse(covs$e3_sex_None == "male", 1, 0)

# create list of omics data 
omics_lst <- simulated_data[-which(names(simulated_data) == "phenotype")]

# Create data frame of omics data
omics_df <- omics_lst |> 
  purrr::map(~as_tibble(.x, rownames = "name")) |>
  purrr::reduce(left_join, by = "name") |>
  column_to_rownames("name")