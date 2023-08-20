# Load functions

# get name of omics layers based on feature names
get_omic_layer <- function(x){ 
  omic = case_when(str_detect(x, "pro_") ~ "Proteins", 
                   str_detect(x, "miR") ~ "miRNA", 
                   str_detect(x, "cg") ~ "Methylation", 
                   str_detect(x, "TC") ~ "Transcriptome", 
                   str_detect(x, "met") ~ "Metabolome")
  omic = if_else(is.na(omic), "Metabolome", omic)
  return(omic)
}

# get name of omics layers based on feature names
get_omic_layer_lowercase <- function(x){ 
  omic = case_when(str_detect(x, "pro_") ~ "proteins", 
                   str_detect(x, "miR") ~ "mirna", 
                   str_detect(x, "cg") ~ "methylome", 
                   str_detect(x, "TC") ~ "transcriptome", 
                   str_detect(x, "met") ~ "metabolome")
  omic = if_else(is.na(omic), "Metabolome", omic)
  return(omic)
}


# get order of omics layers based on feature names
get_omic_layer_numeric <- function(x){ 
  x <- tolower(x)
  if(str_detect(paste(x, collapse = " "), "cg")){
  omic_num = case_when(str_detect(x, "cg") ~ 1, 
                       str_detect(x, "tc") ~ 2, 
                       str_detect(x, "miR") ~ 3,
                       str_detect(x, "pro_") ~ 4, 
                       str_detect(x, "met") ~ 5)
  } else{
    omic_num = case_when(str_detect(x, "meth") ~ 1, 
                         str_detect(x, "transc") ~ 2, 
                         str_detect(x, "miR") ~ 3,
                         str_detect(x, "pro") ~ 4, 
                         str_detect(x, "met") ~ 5)
  }
}

# rename features
rename_ftrs <- function(x){ 
  x %>% 
    str_remove("met_") %>%
    str_remove("hsa-") %>%
    str_remove("pro_")
}
