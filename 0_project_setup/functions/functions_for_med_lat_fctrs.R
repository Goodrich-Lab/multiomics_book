
# create function to reorder data
cluster_reorder <- function(df){
  df <-  column_to_rownames(df, "feature")
  # perform hierarchical clustering on the rows of the dataframe
  dist_mat <- dist(df)
  hc <- hclust(dist_mat)
  # get the numeric order from the hierarchical clustering results
  order <- hc$order
  # order the rows of the dataframe by the clustering results
  df_ordered <- df[order, ]
  
  df_ordered <-  rownames_to_column(df_ordered, "feature")
}

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
  PCs <- pca$x[, 1:n_80_pct] |> scale() |> as.data.frame()
  colnames(PCs) <- paste0(omic_name, "_", colnames(PCs))
  # create a dataframe of loadings for the principal components and scale them
  loadings_df <- pca$rotation[, 1:n_80_pct] |> scale() |> as.data.frame()
  colnames(loadings_df) <- paste0(omic_name, "_", colnames(loadings_df))
  loadings_df <- rownames_to_column(loadings_df, "feature")
  
  # Including all PCs 
  # create a dataframe of scores for the principal components and scale them
  PCs_full <- pca$x |> scale() |> as.data.frame()
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
  rownames(PCs) <- PC_names |>
    str_replace("  ", "_")
  # Transpose and change to data.frame
  out <- as.data.frame(t(PCs))
  # Return output
  return(out)
}
