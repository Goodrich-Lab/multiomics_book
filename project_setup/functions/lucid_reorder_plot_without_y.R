

#----- reorder clusters from LUCID ------------------
#' reorder cluster estimated from LUCID
#'
#' @param model a model returned the function lucid
#' @param order the desired order of the cluster label. For example, if 2 clusters
#' are estiamted and you want to flip the cluster label, you should input c(2, 1).
#' The first element will be used as the reference cluster
#'
reorder_lucid <- function(model,
                          order) {
  # record parameters 
  pars <- model$pars
  beta <- pars$beta
  mu <- pars$mu
  sigma <- pars$sigma
  gamma <- pars$gamma
  K <- model$K
  
  # 1 - reorder exposure effect
  # use the estimate from the 
  ref_cluster <- order[1]
  beta <- t(t(beta) - beta[ref_cluster, ])[order, ]
  
  # 2 - reorder omic effect
  mu <- mu[order, ]
  var <- vector(mode = "list", length = K)
  for (i in 1:K) {
    var[[i]] <- sigma[[order[i]]]
  }
  
  # 3 - reorder outcome effect
  gamma$beta[1:K] <- gamma$beta[order]
  gamma$sigma <- gamma$sigma[order]
  
  model$pars <- list(beta = beta,
                     mu = mu,
                     sigma = sigma,
                     gamma = gamma)
  
  return(model)
}

# example ----
# G <- sim_data$G
# Z <- sim_data$Z
# Y_normal <- sim_data$Y_normal
# cov <- sim_data$Covariate
# 
# # In this model, the exposure effect should be around 2
# fit <- lucid(G = G, Z = Z, Y = Y_normal, CoY = cov)
# summary_lucid(fit)
# flip the label of the 2 clusters, 
# which is equivalent to use the second cluster as the reference cluster
# fit_reorder <- reorder_lucid(model = fit,
#                              order = c(2, 1))
# the exposure effect should be around 1/2


#' 2. change variable name in the sankey diagram ------------------
#' @title Visualize LUCID model through a Sankey diagram
#' @description In the Sankey diagram, each node either represents a variable (exposure,
#' omics or outcome) or a latent cluster. Each line represents an association. The
#' color of the node represents variable type, either exposure, omics or outcome.
#' The width of the line represents the effect size of a certain association; the
#' color of the line represents the direction of a certain association. 
#' 
#' @param x A LUCID model fitted by \code{\link{est_lucid}}
#' @param G_color Color of node for exposure
#' @param X_color Color of node for latent cluster
#' @param Z_color Color of node for omics data
#' @param G_name variable name for latent exposure variables
#' @param Z_name variable name for omic features
#' @param pos_link_color Color of link corresponds to positive association
#' @param neg_link_color Color of link corresponds to negative association
#' @param fontsize Font size for annotation
#' 
#' @return A DAG graph created by \code{\link{sankeyNetwork}}
#' 

plot_lucid_without_outcome <- function(x,
                                       G_color = "dimgray",
                                       X_color = "#eb8c30",
                                       Z_color = "#2fa4da",
                                       G_name = NULL,
                                       Z_name = NULL,
                                       pos_link_color = "#67928b",
                                       neg_link_color = "#d1e5eb",
                                       fontsize = 7
) {
  K <- x$K
  var.names <- x$var.names
  pars <- x$pars
  dimG <- length(var.names$Gnames)
  dimZ <- length(var.names$Znames)
  valueGtoX <- as.vector(t(x$pars$beta[, -1]))
  valueXtoZ <- as.vector(t(x$pars$mu))
  # valueXtoY <- as.vector(x$pars$gamma$beta)[1:K]
  if(is.null(G_name)) {
    G_name = x$var.names$Gnames
  }
  GtoX <- data.frame(
    source = rep(G_name, K),
    target = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, function(x) rep(x, dimG)))),
    value = abs(valueGtoX),
    group = as.factor(valueGtoX > 0))
  if(is.null(Z_name)) {
    Z_name = var.names$Znames
  }
  XtoZ <- data.frame(
    source = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, function(x) rep(x, dimZ)))),
    target = rep(Z_name, K),
    value = abs(valueXtoZ),
    group = as.factor(valueXtoZ > 0))
  # if(is.null(Y_name)) {
  #   Y_name = var.names$Ynames
  # }
  # XtoY <- data.frame(source = paste0("Latent Cluster", 1:K),
  #                    target = rep(Y_name, K),
  #                    value = abs(valueXtoY),
  #                    group = as.factor(valueXtoY > 0))
  
  links <- rbind(GtoX, XtoZ) #, XtoY
  nodes <- data.frame(name = unique(c(as.character(links$source), 
                                      as.character(links$target))),
                      group = as.factor(c(rep("exposure", dimG), 
                                          rep("lc", K), 
                                          rep("biomarker", dimZ))))
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1 
  color_scale <- data.frame(
    domain = c("exposure", "lc", "biomarker", "TRUE", "FALSE"),
    range = c(G_color, X_color, Z_color, pos_link_color, neg_link_color))
  
  p <- sankeyNetwork(Links = links, 
                     Nodes = nodes,
                     Source = "IDsource", 
                     Target = "IDtarget",
                     Value = "value", 
                     NodeID = "name",
                     colourScale = JS(
                       sprintf(
                         'd3.scaleOrdinal()
                        .domain(%s)
                        .range(%s)
                       ',
                       jsonlite::toJSON(color_scale$domain),
                       jsonlite::toJSON(color_scale$range)
                       )), 
                     LinkGroup ="group", 
                     NodeGroup ="group",
                     sinksRight = FALSE, 
                     fontSize = fontsize)
  p
}

# example
# plot_lucid(fit)
# plot_lucid_without_outcome(fit)
# # you can also re-name exposure names and the omic feature name
# plot_lucid_without_outcome(fit,
#            G_name = paste0("PFAS_", 1:10))


#' reorder LUCID in Parallel model by specifying reference cluster ------------
# note: only works for K = 2 in each omic layer
# reference = c(1,1,2)
# lucidus_fit <- fit_reordered
reorder_lucid_parallel <- function(lucidus_fit,
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
  XtoZ <- lucidus_fit$res_Mu
  lucidus_fit$res_Mu <- lapply(1:n_omic, function(i) {
    x <- c(1, 2) # order of clusters
    if(reference[i] == 2) {
      x <- c(2, 1)
      XtoZ[[i]][, x]
    } else{
      XtoZ[[i]][, x]
    }
  }) 
  # reorder gamma
  XtoY <- lucidus_fit$res_Gamma$Gamma$mu
  XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # level w/ new reference
  XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if ref = 2, flip estimates
  lucidus_fit$res_Gamma$Gamma$mu <- XtoY
  lucidus_fit$res_Gamma$fit$coefficients <- XtoY
  
  # return the object using the new reference
  return(lucidus_fit)
}


