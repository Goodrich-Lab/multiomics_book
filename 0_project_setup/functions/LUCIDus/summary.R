#' @title Summarize results of LUCID model
#'
#' @param object A LUCID model fitted by \code{\link{estimate_lucid}}
#' @param boot.se An object returned by \code{\link{boot_lucid}},
#' which contains the bootstrap confidence intervals
#'
#' @export
#' @examples
#' \dontrun{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#'
#' # fit lucid model
#' fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2,
#' seed = 1008)
#'
#' # conduct bootstrap resampling
#' boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100)
#'
#' # summarize lucid model
#' summary_lucid(fit1)
#'
#' # summarize lucid model with bootstrap CIs
#' summary_lucid(fit1, boot.se = boot1)
#' }

summary_lucid <- function(object, boot.se = NULL){
  if (class(object) == "early_lucid"){
    s1 <- object$select$selectG
    s2 <- object$select$selectZ
    nG <- sum(s1)
    nZ <- sum(s2)
    K <- object$K
    gamma <- object$res_Gamma$beta
    #obtain number of parameters
    if(object$family == "normal"){
      nY <- length(object$res_Gamma$beta) + length(object$res_Gamma$sigma)
    }
    if(object$family == "binary"){
      nY <- length(object$res_Gamma$beta)
    }
    npars <- (nG + 1) * (K - 1) + (nZ * K + nZ^2 * K) + nY
    BIC <- -2 * object$likelihood + npars * log(nrow(object$post.p))
    results <- list(beta = object$res_Beta[, c(TRUE, s1)],
                    mu = object$res_Mu[, s2],
                    gamma = object$res_Gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                    loglik = object$likelihood,
                    boot.se = boot.se)
    class(results) <- "sumlucid"
    return(results)
  }
  if (class(object) == "lucid_parallel"){
    ##not having regularity yet, to be added
    nG <- ncol(object$res_Beta) -1
    K <- object$K
    Z <- object$Z
    #obtain number of parameters
    if(object$family == "gaussian"){
      nY <- length(object$res_Gamma$Gamma$mu) + length(object$res_Gamma$Gamma$sd)
    }
    if(object$family == "binomial"){
      #binary summary res_Gamma$Gamma$mu is unclear, use object$res_Gamma$fit$coefficients instead
      nY <- length(object$res_Gamma$fit$coefficients)
    }
    npars = 0
    for (i in 1:length(K)){
      nZ = ncol(Z[[i]])
      npars_new = (nG + 1) * (K[i] - 1) + (nZ * K[i] + nZ * nZ * K[i])
      npars = npars + npars_new
    }
    npars = npars + nY
    BIC <- -2 * object$likelihood + npars * log(nrow(object$post.p[[1]]))

    results <- list(beta = object$res_Beta,
                    mu = object$res_Mu,
                    Gamma = object$res_Gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                    loglik = object$likelihood
                    #boot.se = boot.se
    )
    class(results) <- "sumlucid"
    return(results)

  }
  if (class(object) == "lucid_serial"){

  }
}


#' Print the output of LUCID in a nicer table
#'
#' @param x An object returned by \code{summary_lucid}
#' @param ... Other parameters to be passed to \code{print}
#' @export
#'
print.sumlucid <- function(x, ...){
  K <- x$K
  beta <- as.data.frame(x$beta)
  dim1 <- ncol(beta) - 1
  z.mean <- as.data.frame(t(x$mu))
  cat("----------Summary of the LUCID model---------- \n \n")
  cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")
  y <- switch(x$family, normal = f.normal,
              binary = f.binary)
  y(x$gamma, K, se = x$boot.se$gamma)
  cat("\n")
  cat("(2) Z: mean of omics data for each latent cluster \n")
  if(is.null(x$boot.se)){
    colnames(z.mean) <- paste0("mu_cluster", 1:K)
    print(z.mean)
  } else{
    print(x$boot.se$mu)
  }
  cat("\n")
  cat("(3) E: odds ratio of being assigned to each latent cluster for each exposure \n")
  if(is.null(ncol(beta))) {
    cat("no exposure is selected given current penalty Rho_G, please use a smaller penalty")
  } else {
    dd <- as.matrix(as.data.frame(beta)[2:K, 2:ncol(beta)])
    g.or <- data.frame(beta = unlist(split(dd, row(dd))))
    rownames(g.or) <- paste0(colnames(beta)[-1], ".cluster", sapply(2:K, function(x) return(rep(x, dim1))))
    if(is.null(x$boot.se)){
      g.or$OR <- exp(g.or$beta)
      print(g.or)
    } else{
      print(x$boot.se$beta)
    }
  }
}


# summarize output of normal outcome
f.normal <- function(x, K, se){

  cat("(1) Y (continuous outcome): mean of Y for each latent cluster (and effect of covariates if included) \n")

  if(!is.null(se)){
    y <- se
  } else {
    beta <- x$beta
    y <- as.data.frame(beta)
    row.names(y)[1:K] <- paste0("cluster", 1:K)
    colnames(y) <- "Gamma"
  }
  print(y)
}


# summarize output of binary outcome
f.binary <- function(x, K, se){
  cat("(1) Y (binary outcome): log odds of Y for cluster 1 (reference) and log OR for rest cluster (and log OR of covariate if included)\n")
  gamma <- as.data.frame(x$beta)
  colnames(gamma) <- "gamma"
  if(is.null(se)){
    gamma$`exp(gamma)` <- exp(gamma$gamma)
  } else{
    gamma <- cbind(gamma, se[, -1])
  }
  print(gamma)
}



##########functions for LUCID in parallel##########




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

  # if 1 omics layers
  if(length(K) == 1) {
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


    K_order <- list(K1 = k1_order)
  }

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
  if(length(K) == 4) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1,1] - mu_matrix[ref[1], 1, 1,1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1,1] - mu_matrix[1, ref[2], 1,1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i,1] - mu_matrix[1, 1, ref[3],1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i] - mu_matrix[1, 1, 1, ref[4]]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order)
  }

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


  # if 6 omics layers
  if(length(K) == 6) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1, 1] - mu_matrix[1, 1, 1, ref[4], 1, 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i, 1] - mu_matrix[1, 1, 1, 1, ref[5], 1]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])

    # reorder K6
    mu_k6 <- rep(0, K[6])
    for(i in 1:K[6]) {
      mu_k6[i] <- mu_matrix[1, 1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, 1, ref[6]]
    }
    mu_k6_sort <- sort(mu_k6)
    mu <- c(mu, mu_k6_sort[mu_k6_sort != 0])
    # order of re-arranged cluster for omics 6
    k6 <- order(mu_k6)
    k6_order <- c(ref[6], k6[k6 != ref[6]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order,
                    K6 = k6_order)
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


#' function to reorder all model parameters
#'
#' @param model A model returned by EM_lucid
#'
#' @return A LUCID model reordered by effect size of outcome
#' @export
#'
reorder_lucid <- function(model) {
  ref <- get_ref_cluster(Delta = model$res_Delta$Delta)
  r_Delta <- reorder_Delta(ref = ref,
                           Delta = model$res_Delta$Delta)
  r_Mu_Sigma <- reorder_Mu_Sigma(model$res_Mu_Sigma,
                                 K_order = r_Delta$K_order)
  r_Beta <- reorder_Beta(Beta = model$res_Beta$Beta,
                         K_order = r_Delta$K_order)
  model$res_Delta$Delta <- r_Delta
  model$res_Mu_Sigma$Mu <- r_Mu_Sigma$Mu
  model$res_Mu_Sigma$Sigma <- r_Mu_Sigma$Sigma
  model$res_Beta$Beta <- r_Beta
  model$z <- reorder_z(model$z, K_order = r_Delta$K_order)
  return(model)
}


# function to calculate BIC
cal_bic <- function(model) {
  nOmics <- length(model$K)
  Beta <- model$res_Beta$Beta
  Mu <- model$res_Mu_Sigma$Mu
  Sigma <- model$res_Mu_Sigma$Sigma
  Delta <- model$res_Delta$Delta

  # calculate number of parameters
  n_Beta <- sum(sapply(1:nOmics, function(i) {
    length(Beta[[i]])
  }))
  n_Mu <- sum(sapply(1:nOmics, function(i) {
    length(Mu[[i]])
  }))
  n_Sigma <- sum(sapply(1:nOmics, function(i) {
    length(Sigma[[i]])
  }))
  n_Delta <- length(Delta$mu) + 1
  n_par <- n_Beta + n_Mu + n_Sigma + n_Delta

  # calculate BIC
  bic <- -2 * model$loglik + n_par * log(model$N)
  return(bic)
}

