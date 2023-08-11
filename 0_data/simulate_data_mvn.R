set.seed(28894)
# Load Data
# Read in reduced scaled omics data
dir_home_1 <- here::here() 
helix_dat_reduced <- read_rds(fs::path(dir_data_hg,
                                       "Hg_ck18_screened_scaled_omics.RDS"))
#extract the omics
methylome_d <- data.frame(helix_dat_reduced$methylome)
rownames(methylome_d) <- c(1:420)
sum(sapply(methylome_d,is.factor))



transcriptome_d <- data.frame(helix_dat_reduced$transcriptome) 
rownames(transcriptome_d) <- c(1:420)
sum(sapply(transcriptome_d,is.factor))

miRNA_d <- data.frame(helix_dat_reduced$miRNA) 
rownames(miRNA_d) <- c(1:420)
sum(sapply(miRNA_d,is.factor))

proteome_d <- data.frame(helix_dat_reduced$proteome)
rownames(proteome_d) <- c(1:420)
sum(sapply(proteome_d,is.factor))

metabolome_d <- data.frame(helix_dat_reduced$metabolome)
rownames(metabolome_d) <- c(1:420)
sum(sapply(metabolome_d,is.factor))

pheno_d <- data.frame(helix_dat_reduced$phenotype) 
rownames(pheno_d) <- c(1:420)
pheno_d$helixid <- c(1:420)
sum(sapply(pheno_d,is.factor))
#extract the covariates and outcomes only
index = sapply(pheno_d,is.factor)
pheno_d_covariate = subset(pheno_d,select = c(#"h_cohort",
  "e3_sex_None", 
  "h_fish_preg_Ter"))

#randomize it 
pheno_d_covariate <- pheno_d_covariate[sample(1:nrow(pheno_d_covariate)), ] 


omics_d = cbind(methylome_d,transcriptome_d,miRNA_d,proteome_d,metabolome_d,
                pheno_d$hs_child_age_yrs_None,
                pheno_d$hs_hg_m_scaled, pheno_d$ck18_scaled)

sum(sapply(omics_d,is.factor))
#use MVN to simulate data
mean_omics_d = colMeans(omics_d)
cov_omics_d = cov(omics_d)
omics_d_sim = mvrnorm(n=420, mu = mean_omics_d, Sigma = cov_omics_d)
#visualize the correlation pattern
library(corrplot)
corr = cor(omics_d_sim)
corrplot(corr)
corr_r = cor(omics_d)
corrplot(corr_r)

#extract simulated omics layers
s_methylome_d <- omics_d_sim[,1:28]
s_transcriptome_d <- omics_d_sim[,29:56]
s_miRNA_d <- omics_d_sim[,57:84]
s_proteome_d <- omics_d_sim[,85:112]
s_metabolome_d <- omics_d_sim[,113:140]
s_pheno_d <- cbind(omics_d_sim[,141:143],pheno_d_covariate)
colnames(s_pheno_d)[1] ="hs_child_age_yrs_None"
colnames(s_pheno_d)[2] ="hs_hg_m_scaled"
colnames(s_pheno_d)[3] ="ck18_scaled"
rownames(s_pheno_d) <- c(1:420)


#make the data-analysis ready simulated reduced helix dataset 
s_helix_dat_reduced <- helix_dat_reduced
s_helix_dat_reduced$methylome <- s_methylome_d
s_helix_dat_reduced$transcriptome <- s_transcriptome_d
s_helix_dat_reduced$miRNA <- s_miRNA_d
s_helix_dat_reduced$proteome <- s_proteome_d
s_helix_dat_reduced$metabolome <- s_metabolome_d
s_helix_dat_reduced$phenotype <- s_pheno_d

#save the simulated data
write_rds(s_helix_dat_reduced, file = fs::path(dir_data_hg,
                                               "simulate_Hg_ck18_screened_scaled_omics.RDS"))


# Reduce annotation dataset to only include simulated data
s_helix_dat_reduced <- readRDS(file = fs::path(dir_data_hg,
                                               "simulate_Hg_ck18_screened_scaled_omics.RDS"))
# Get all 
# ftrnms <- lapply(s_helix_dat_reduced[-6], colnames) %>% unlist() %>% as.character()
# anodat <- readRDS(file = fs::path(dir_data_hg, "all_omics_annotation_v2.RDS"))
# # Get name to match the column names 
# anodat$ftr_name_match = str_replace_all(anodat$ftr_name, "-", ".")
# # Reduce omic metadata to only features in the data
# anodat_reduced = anodat %>%
#   dplyr::filter(ftr_name_match %in% ftrnms)
# anodat_reduced <- anodat_reduced %>%
#   mutate(ftr_name = ftr_name_match)  %>%
#   select(-ftr_name_match)
# # Save file
# saveRDS(anodat_reduced, file = fs::path(dir_data_hg, "feature_metadata.RDS"))

## Read back in
anodat_reduced <- readRDS(file = fs::path(dir_data_hg, "feature_metadata.RDS"))

# Create omic names for plotting
anodat_reduced <- anodat_reduced %>%
  mutate(ftr_name_for_plots = case_when(
    omic_layer == "miRNA" ~ str_remove(ftr_name, "hsa-"),
    omic_layer == "metabolome" ~ common_name,
    omic_layer == "methylome" ~ if_else(is.na(gene_names),
                                        ftr_name,
                                        gene_names),
    omic_layer == "proteome" ~ str_remove(ftr_name, "pro_"),
    omic_layer == "transcriptome" ~ gene_assignment)) %>%
  mutate(ftr_name_for_plots = case_when(
    ftr_name_for_plots == "---" ~ paste0(ftr_name, " (NONCODE)"),
    is.na(ftr_name_for_plots)  ~ ftr_name,
    str_detect(ftr_name_for_plots, ";") ~ str_split_fixed(ftr_name_for_plots, ";", 2)[,1],
    str_detect(ftr_name_for_plots, "//") ~ str_split_fixed(ftr_name_for_plots, "//", 3)[,2],
    TRUE ~ ftr_name_for_plots)) %>%
  group_by(ftr_name_for_plots) %>% 
  mutate(ftr_name_for_plots = str_c(ftr_name_for_plots, row_number(), sep = "___")) %>%
  ungroup() %>%
  mutate(ftr_name_for_plots = str_replace(ftr_name_for_plots, " ___2", ".2") %>% 
           str_replace("___2", ".2") %>%
           str_remove("___1"))


anodat_reduced$ftr_name_for_plots

anodat_reduced2 <- anodat_reduced %>%
  dplyr::select(omic_layer, ftr_name, ftr_name_for_plots)

saveRDS(anodat_reduced2, file = fs::path(dir_data_hg, "feature_metadata_v2.RDS"))

