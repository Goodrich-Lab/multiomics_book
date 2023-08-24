
# rearrange cluster order
#
# for continuous outcome - use the cluster combination corresponding to smallest
# mean as the reference cluster
get_ref_cluster <- function(Delta) {
  K <- Delta$K
  mu <- Delta$mu
  mu_matrix <- vec_to_array(K = K, mu = mu)
  ref_index <- which(mu_matrix == min(mu_matrix))
  ref <- arrayInd(ref_index, .dim = K)
  return(ref)
}


# re-arrange parameters for Delta
reorder_Delta <- function(ref, Delta) {
  K <- Delta$K
  mu_matrix <- vec_to_array(K = K, mu = Delta$mu)
  mu <- mu_matrix[ref]
  
  # if 2 omics layers
  if(length(K) == 2) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])
    
    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i] - mu_matrix[1, ref[2]]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])
    
    K_order <- list(K1 = k1_order,
                    K2 = k2_order)
  }
  
  
  # if 3 omics layers
  if(length(K) == 3) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1] - mu_matrix[ref[1], 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])
    
    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1] - mu_matrix[1, ref[2], 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])
    
    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i] - mu_matrix[1, 1, ref[3]]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])
    
    
    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order)
  }
  
  
  # if 4 omics layers
  
  # if 5 omics layers
  if(length(K) == 5) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])
    
    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])
    
    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])
    
    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1] - mu_matrix[1, 1, 1, ref[4], 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])
    
    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, ref[5]]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])
    
    
    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order)
  }
  
  
  Delta$mu <- mu
  return(list(Delta = Delta,
              K_order = K_order))
}



reorder_Mu_Sigma <- function(Mu_Sigma, K_order) {
  for(i in 1:length(K_order)) {
    temp_Mu <- Mu_Sigma$Mu[[i]]
    temp_Sigma <- Mu_Sigma$Sigma[[i]]
    # reorder Mu
    Mu_Sigma$Mu[[i]] <- temp_Mu[, K_order[[i]]]
    Mu_Sigma$Sigma[[i]] <- temp_Sigma[, , K_order[[i]]]
  }
  return(Mu_Sigma)
}



reorder_Beta <- function(Beta, K_order) {
  for(i in 1:length(K_order)) {
    temp_Beta <- Beta[[i]]
    temp_Beta <- rbind(rep(0, ncol(temp_Beta)),
                       temp_Beta)
    temp_Beta_reorder <- temp_Beta[K_order[[i]], ]
    ref <- temp_Beta_reorder[1, ]
    for(j in 1:nrow(temp_Beta_reorder)) {
      temp_Beta_reorder[j, ] <- temp_Beta_reorder[j, ] - ref
    }
    Beta[[i]] <- temp_Beta_reorder[-1, ]
  }
  
  return(Beta)
}



reorder_z <- function(z, K_order) {
  if(length(K_order) == 2) {
    z <- z[K_order[[1]], K_order[[2]], ]
  }
  return(z)
}



# Reorder lucid M 

# ----------- reorder LUCIDusM model by specifying reference cluster ------------
# note: only works for K = 2 in each omic layer
# reference = c(1,1,2)
# lucidus_fit <- fit_reordered
reorder_lucidM <- function(lucidus_fit,
                           reference = NULL) {
  if(is.null(reference)) {
    warning("no reference specified, return the original model")
    return(lucidus_fit)
  }
  
  n_omic <- length(reference)
  
  # reorder beta
  GtoX <- lucidus_fit$res_Beta$Beta
  lucidus_fit$res_Beta$Beta <- lapply(1:n_omic, function(i) {
    (-1)^(reference[i] - 1) * GtoX[[i]] # if reference = 1, no changes; 
    # if reference = 2, flip the reference and negate the estimates
  })
  # reorder mu
  XtoZ <- lucidus_fit$res_Mu_Sigma$Mu
  lucidus_fit$res_Mu_Sigma$Mu <- lapply(1:n_omic, function(i) {
    x <- c(1, 2) # order of clusters
    if(reference[i] == 2) {
      x <- c(2, 1)
      XtoZ[[i]][, x]
    } else{
      XtoZ[[i]][, x]
    }
  }) 
  # reorder gamma
  XtoY <- lucidus_fit$res_Delta$Gamma$mu
  XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # reference level using the new reference
  XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if reference = 2, flip the estimates
  lucidus_fit$res_Delta$Gamma$mu <- XtoY
  lucidus_fit$res_Delta$fit$coefficients <- XtoY
  
  # return the object using the new reference
  return(lucidus_fit)
}




