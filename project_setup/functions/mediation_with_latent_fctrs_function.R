## ---- Mediation with latent fctrs early ----
#' Conduct principal component analysis (PCA) as a dimensional reduction step 
#' and selected the top i principal components which explained >80% of the variance. 
#' Following the joint dimensional reduction step, conduct mediation analysis.

#' Given exposure, outcome and multiple omics data,
#' function runs HIMA mediation analysis with latent factors in early integration
#' and returns a list including two tidy dataframes for the results
#'
#' @param exposure A numeric vector for the exposure variable
#' @param outcome A numeric vector for the outcome variable
#' @param omics_lst A list of numeric matrices representing omics data
#' @param covs A numeric matrix representing the covariates
#'
#' @return A list including two tidy dataframes summarizing the results of HIMA analysis
#' one dataframe including the significant result of HIMA analysis with PCs as mediators.
#' another dataframe including the result of feature correlation with significant PCs.
#' 
#' @import dplyr
#' @importFrom tidyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom purrr map
#' @importFrom stats prcomp cor
#' @importFrom HIMA hima
#' @importFrom base cumsum min sum scale apply
#' @importFrom stringr str_replace
#'
med_with_latent_fctrs_early_integration <- function(exposure, 
                                   outcome, 
                                   omics_lst, 
                                   covs = NULL, 
                                   Y.family = "gaussian",
                                   M.family = "gaussian",
                                   fdr.level = 0.05) {
  
  # Combines omics data into one dataframe
  omics_lst_df <- purrr::map(omics_lst, ~as_tibble(.x, rownames = "name"))
  
  # Meta data
  meta_df <- imap_dfr(omics_lst_df, ~tibble(omic_layer = .y, ftr_name = names(.x)))%>%
    filter(ftr_name != "name") %>%
    mutate(omic_num = case_when(str_detect(omic_layer, "meth") ~ 1, 
                                str_detect(omic_layer, "transc") ~ 2, 
                                str_detect(omic_layer, "miR") ~ 3,
                                str_detect(omic_layer, "pro") ~ 4, 
                                str_detect(omic_layer, "met") ~ 5))
  
  # i. Obtain PCs----
  omics_df_pca <- prcomp(omics_df, center = TRUE, scale. = TRUE)
  
  # Calculate variance and proportion of variance explained by each PC
  vars <- apply(omics_df_pca$x, 2, var)
  props <- vars / sum(vars)
  cum_props <- cumsum(props)
  
  # Determine the number of PCs needed to explain > 80% of the total variance
  n_80_pct <- min((1:length(cum_props))[cum_props > 0.8])
  
  # The first PCs which explain >80% of the variance are used as latent mediators
  PCs <- omics_df_pca$x[, 1:n_80_pct] %>% scale()
  
  # ii. Perform HIMA with PCs as mediator----
  result_hima_comb_pc <- hima(X = exposure,
                              Y = outcome,
                              M = PCs,
                              COV.MY = covs,
                              COV.XM = covs,
                              Y.family =  Y.family,
                              M.family = M.family, 
                              scale = FALSE)
  
  # Change to tibble and select significant PCs
  result_hima_pca_early <- as_tibble(result_hima_comb_pc, rownames = "pc_num") %>% 
    filter(BH.FDR < fdr.level) 
  
  # Filter Significant PCs, Create scaled %TE variable
  result_hima_pca_early_sig <- result_hima_pca_early %>% 
    mutate(
      te_direction = if_else(beta<0, -1*`% total effect`, `% total effect`), 
      `% Total Effect scaled` = 100*`% total effect`/sum(`% total effect`) %>%
        round(1), 
      multiomic_mthd = "Early") %>%
    mutate(pc_named = str_replace(pc_num, "PC", "Joint Comp. "),
           pcs_ordered = forcats::fct_reorder(pc_named, `te_direction`)) %>%
    rename(Alpha = alpha, 
           Beta = beta,
           `TE (%)` = `% Total Effect scaled`) 
  
  # ii correlation of features vs PC's --------------
  # Extract variable correlation with principal components
  var.cor <- cor(omics_df, omics_df_pca$x)
  
  # Select only significant PCs
  ftr_cor_sig_pcs <- var.cor[,(colnames(var.cor) %in% 
                                 result_hima_pca_early_sig$pc_num)]
  
  ftr_cor_sig_pcs_df <- ftr_cor_sig_pcs %>%
    as_tibble(rownames = "feature") %>%
    left_join(meta_df, by = c("feature" = "ftr_name"))
    
  
  res = list(result_hima_pca_early_sig = result_hima_pca_early_sig, 
             result_ftr_cor_sig_pcs_early = ftr_cor_sig_pcs_df)
  return(res)
}

## ---- Mediation with latent fctrs intermediate ----
#' Perform a joint dimensionality reduction step using Joint and Individual Variance Explained (JIVE). 
#' Following the joint dimensionality reduction step, conduct mediation analysis.


#' Given exposure, outcome and multiple omics data,
#' function runs HIMA mediation analysis with latent factors in intermediate integration
#' and returns a list including two tidy dataframes for the results
#'
#' @param exposure A numeric vector for the exposure variable
#' @param outcome A numeric vector for the outcome variable
#' @param omics_lst A list of numeric matrices representing omics data
#' @param covs A numeric matrix representing the covariates
#'
#' @return A list including two tidy dataframes summarizing the results of HIMA analysis
#' one dataframe including the significant result of HIMA analysis with latent factors as mediators.
#' another dataframe including the result of feature correlation with significant PCs.
#' 
#' @import dplyr
#' @importFrom tidyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom r.jive jive
#' @importFrom stats prcomp cor
#' @importFrom HIMA hima
#' @importFrom base min sum scale apply svd diag
#' @importFrom stringr str_replace
#'
med_with_latent_fctrs_intermediate_integration <- function(exposure, 
                                                    outcome, 
                                                    omics_lst, 
                                                    covs = NULL, 
                                                    Y.family = "gaussian",
                                                    M.family = "gaussian",
                                                    fdr.level = 0.05) {
  # Combines omics data into one dataframe
  omics_lst_df <- purrr::map(omics_lst, ~as_tibble(.x, rownames = "name"))
  
  # Meta data
  meta_df <- imap_dfr(omics_lst_df, ~tibble(omic_layer = .y, ftr_name = names(.x)))%>%
    filter(ftr_name != "name") %>%
    mutate(omic_num = case_when(str_detect(omic_layer, "meth") ~ 1, 
                                str_detect(omic_layer, "transc") ~ 2, 
                                str_detect(omic_layer, "miR") ~ 3,
                                str_detect(omic_layer, "pro") ~ 4, 
                                str_detect(omic_layer, "met") ~ 5))
  # Transpose omics matrices
  omics_t <- lapply(omics_lst, t)
  
  # Rename omics datasets
  names(omics_t) = names(omics_lst)
  
  # Run JIVE with the optimal ranks provided:----- 
  result_jive2 <- jive(data = omics_t,
                       rankJ = 22,
                       rankA = c(6, 9, 5, 5, 8),
                       method = "given",
                       conv = 1e-04,
                       maxiter = 100,
                       showProgress = FALSE)
  
  # Get the components from JIVE
  # Function to extract joint and individual PCs ----
  get_PCs_jive <- function (result){
    # Number of joint factors 
    n_joint = result$rankJ
    # Number of individual factors
    n_indiv = result$rankA
    # Get number of data matrices in result
    l <- length(result$data)    
    # Calculate total number of PCs to compute
    nPCs = n_joint + sum(n_indiv)    
    # Initialize matrix to hold PC scores
    PCs = matrix(nrow = nPCs, ncol = dim(result$data[[1]])[2])    
    # Initialize vector to hold PC names
    PC_names = rep("", nPCs)    
    # If joint structure is present
    if (n_joint > 0) {   
      # Compute SVD on joint structure
      SVD = svd(do.call(rbind, result$joint), nu = n_joint, nv = n_joint)
      # Compute PC scores for joint structure
      PCs[1:n_joint,] = diag(SVD$d)[1:n_joint, 1:n_joint] %*% t(SVD$v[, 1:n_joint])    
      # Assign names to joint PCs
      PC_names[1:n_joint] = paste("Joint ", 1:n_joint)    
    }
    # Loop over data matrices
    for (i in 1:l) {    
      # If individual structure is present for this matrix
      if (n_indiv[i] > 0) {    
        # Compute SVD on individual structure
        SVD = svd(result$individual[[i]], nu = n_indiv[i], 
                  nv = n_indiv[i])    
        # Get indices for PCs corresponding to this data matrix
        indices = (n_joint + sum(n_indiv[0:(i - 1)]) + 1):(n_joint + 
                                                             sum(n_indiv[0:i]))    
        # Compute PC scores for individual structure
        PCs[indices, ] = diag(SVD$d)[1:n_indiv[i], 1:n_indiv[i]] %*% 
          t(SVD$v[, 1:n_indiv[i]])   
        # Assign names to individual PCs
        PC_names[indices] = paste0(names(result$data)[i], 
                                   "_", 1:n_indiv[i]) 
      }
    }
    # Rename PCs
    rownames(PCs) <- PC_names %>%
      str_replace("  ", "_")
    # Transpose and change to data.frame
    out <- as.data.frame(t(PCs))
    # Return output
    return(out)
  }
  
  factors_jive <- get_PCs_jive(result_jive2) %>% 
    dplyr::mutate(across(everything(), ~as.vector(scale(.))))
  
  # run mediation analysis------
  result_hima_jive <- hima(X = exposure,
                           Y = outcome,
                           M = factors_jive,
                           COV.MY = covs,
                           COV.XM = covs,
                           Y.family = c("gaussian"),
                           M.family = c("gaussian"), 
                           verbose = TRUE,
                           scale = FALSE)
  
  # Modify and filter significant                
  result_hima_jive_1 <- result_hima_jive %>% 
    rownames_to_column("component") %>% 
    mutate(multiomic_mthd = "Intermediate",
           ind_joint = str_split_fixed(component, fixed("_"), 2)[,1], 
           ind_joint_num = case_when(ind_joint == "Joint" ~ 1, 
                                     ind_joint == "methylome" ~ 2,
                                     ind_joint == "transcriptome" ~ 3)) %>% 
    dplyr::select(multiomic_mthd, everything())
  
  # Filter significant components and create scaled %TE variable
  result_hima_jive_sig <-  result_hima_jive_1 %>% 
    filter(BH.FDR<fdr.level) %>%   
    mutate(`% Total Effect scaled` =
             round(100*`% total effect`/sum(`% total effect`),1),
           te_direction = if_else(beta < 0, 
                                  -1 * `% total effect`, 
                                  `% total effect`)) %>%
    mutate(components_ord = forcats::fct_reorder(component, `te_direction`)) %>%
    rename(Alpha = alpha, 
           Beta = beta,
           `TE (%)` = `% Total Effect scaled`)
  
  # correlation of features vs JIVE factors --------------
  # Extract variable correlation with JIVE Factors
  var.cor <- cor(omics_df, factors_jive) 
  
  # Select only significant PCs
  ftr_cor_sig_pcs_jive <- var.cor[,(colnames(var.cor) %in% 
                                          result_hima_jive_sig$component)] 
  
  ftr_cor_sig_pcs_jive_df <- ftr_cor_sig_pcs_jive %>%
    as_tibble(rownames = "feature") %>%
    left_join(meta_df, by = c("feature" = "ftr_name"))
  
  res = list(result_hima_jive_sig = result_hima_jive_sig, 
             result_ftr_cor_sig_pcs_jive = ftr_cor_sig_pcs_jive_df)
  return(res)
}


## ---- Mediation with latent fctrs late ----
#' Conduct principal component analysis (PCA) as a dimensionality reduction step  
#' on each omics layer separately 
#' and selected the top i principal components which explained >80% of the variance. 
#' Following the joint dimensionality reduction step, conduct mediation analysis.

#' Given exposure, outcome and multiple omics data,
#' function runs HIMA mediation analysis with latent factors in late integration
#' and returns a list including two tidy dataframes for the results
#'
#' @param exposure A numeric vector for the exposure variable
#' @param outcome A numeric vector for the outcome variable
#' @param omics_lst A list of numeric matrices representing omics data
#' @param covs A numeric matrix representing the covariates
#'
#' @return A list including two tidy dataframes summarizing the results of HIMA analysis
#' one dataframe including the significant result of HIMA analysis with PCs as mediators for each omics layer.
#' another dataframe including the result of feature correlation with significant PCs.
#' 
#' @import dplyr
#' @importFrom tidyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom purrr map2 mapreduce
#' @importFrom stats prcomp cor
#' @importFrom HIMA hima
#' @importFrom base cumsum min sum scale apply
#' @importFrom stringr str_replace
#'
med_with_latent_fctrs_late_integration <- function(exposure, 
                                                    outcome, 
                                                    omics_lst, 
                                                    covs = NULL, 
                                                    Y.family = "gaussian",
                                                    M.family = "gaussian",
                                                    fdr.level = 0.05) {
  
  # Combines omics data into one dataframe
  omics_lst_df <- purrr::map(omics_lst, ~as_tibble(.x, rownames = "name"))
  
  # Meta data
  meta_df <- imap_dfr(omics_lst_df, ~tibble(omic_layer = .y, ftr_name = names(.x)))%>%
    filter(ftr_name != "name") %>%
    mutate(omic_num = case_when(str_detect(omic_layer, "meth") ~ 1, 
                                str_detect(omic_layer, "transc") ~ 2, 
                                str_detect(omic_layer, "miR") ~ 3,
                                str_detect(omic_layer, "pro") ~ 4, 
                                str_detect(omic_layer, "met") ~ 5))
  
  # Define function to run PCA on a matrix and return loadings and scores
  run_pca <- function(mat, omic_name) {
    # perform PCA on input matrix with scaling
    pca <- prcomp(mat, scale. = TRUE)
    # calculate variance explained by each principal component
    vars <- apply(pca$x, 2, var)
    props <- vars / sum(vars)
    # calculate cumulative proportion of variance explained
    cum_props <- cumsum(props)
    
    # determine the number of principal components needed to explain 80% of the variance
    n_80_pct <- min((1:length(cum_props))[cum_props > 0.8])
    # create a dataframe of scores for the principal components and scale them
    PCs <- pca$x[, 1:n_80_pct] %>% scale() %>% as.data.frame()
    colnames(PCs) <- paste0(omic_name, "_", colnames(PCs))
    # create a dataframe of loadings for the principal components and scale them
    loadings_df <- pca$rotation[, 1:n_80_pct] %>% scale() %>% as.data.frame()
    colnames(loadings_df) <- paste0(omic_name, "_", colnames(loadings_df))
    loadings_df <- rownames_to_column(loadings_df, "feature")
    
    # Including all PCs 
    # create a dataframe of scores for the principal components and scale them
    PCs_full <- pca$x %>% scale() %>% as.data.frame()
    colnames(PCs_full) <- paste0(omic_name, "_", colnames(PCs_full))
    
    # create a dataframe of proportion of variance explained by each principal component
    props_df <- data.frame(pc_num = paste0(omic_name, "_",
                                           names(props)),
                           pc_var_explained = props)
    rownames(props_df) <- NULL
    # return a list of results
    return(list(loadings = loadings_df, scores = PCs, scores_full = PCs_full,
                n_pcs_80_pct = n_80_pct, pc_var_explained = props_df))
  }
  # i. get PCs and HIMA -----
  # PCA scores: Apply function to each matrix in the list and collect 
  scores_list_late_int <-  map2(omics_lst, 
                                names(omics_lst), 
                                ~run_pca(.x, .y)$scores)
  scores_df <- purrr::reduce(scores_list_late_int, cbind) %>% as.data.frame()
  
  # Loadings: Apply function to each matrix in the list and collect PC 
  loadings_list_late_int <- map2(omics_lst, 
                                 names(omics_lst),
                                 ~run_pca(.x, .y)$loadings)
  loadings_df <- purrr::reduce(loadings_list_late_int, full_join)
  
  # Number of PCs explain >80%: Get number of PCs for each omic
  late_int_pcs_80 <- purrr::map(omics_lst, ~run_pca(.x, NULL)$n_pcs_80_pct) %>%
    bind_rows()
  
  # ii. run HIMA with the current dataset and principal components ----
  result_hima_late_integration <- hima(X = exposure,
                                       Y = outcome,
                                       M = scores_df, 
                                       COV.MY = covs,
                                       COV.XM = covs,
                                       Y.family = c("gaussian"),
                                       M.family = c("gaussian"), 
                                       scale = FALSE) 
  
  # Filter significant pcs, create scaled %TE variable
  result_hima_late_sig <- result_hima_late_integration %>%
    as_tibble(rownames = "pc_num") %>%
    filter(BH.FDR < 0.05) %>%
    mutate(
      te_direction = if_else(beta<0, -1*`% total effect`, `% total effect`), 
      `% Total Effect scaled` = 100*`% total effect`/sum(`% total effect`)) 
  
  result_hima_late_sig <- result_hima_late_sig  %>% 
    mutate(pc_numeric = str_split_fixed(pc_num, fixed("_"), 2)[,2],
           omic_layer = str_split_fixed(pc_num, fixed("_"), 2)[,1] %>%
             str_to_sentence(),
           multiomic_mthd = "Late")%>%
    mutate(omic_layer = str_replace(omic_layer, "Mirna", "miRNA" ),
           omic_pc = str_c(omic_layer, " ", pc_numeric) %>%
             str_replace("PC", "Comp. ")) %>%
    mutate(pcs_ordered = forcats::fct_reorder(omic_pc, `te_direction`)) %>%
    dplyr::select(multiomic_mthd, omic_pc, omic_layer, pc_num, pcs_ordered,
                  alpha, beta, `% Total Effect scaled`) %>%
    rename(Alpha = alpha, 
           Beta = beta,
           `TE (%)` = `% Total Effect scaled`)%>%
    mutate(omic_num = case_when(str_detect(pc_num, "meth") ~ 1, 
                                str_detect(pc_num, "transc") ~ 2, 
                                str_detect(pc_num, "miR") ~ 3,
                                str_detect(pc_num, "pro") ~ 4, 
                                str_detect(pc_num, "met") ~ 5))
  
  
  # Extract variable correlation with principal components
  scores_full_list_late_int <- map2(omics_lst, 
                                    names(omics_lst),
                                    ~run_pca(.x, .y)$scores_full)
  
  scores_df_full <- purrr::reduce(scores_full_list_late_int, cbind) %>% as.data.frame()
  
  # Get correlation of omics and PCs by omic layer
  var.cor <- map_df(names(omics_lst), function(type) {
    cor(omics_lst[[type]], scores_full_list_late_int[[type]]) %>%
      as.data.frame()
  })
  
  # Select only significant PCs
  ftr_cor_sig_pcs_late <- var.cor %>%
    dplyr::select(result_hima_late_sig$pc_num) %>%
    as_tibble(rownames = "feature") %>%
    left_join(meta_df, by = c("feature" = "ftr_name"))
  
  res = list(result_hima_late_sig = result_hima_late_sig, 
             result_ftr_cor_sig_pcs_late = ftr_cor_sig_pcs_late)
  return(res)
  
}
