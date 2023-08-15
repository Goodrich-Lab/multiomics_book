#' @title A wrapper function to perform model selection for LUCID
#'
#' @description Given a grid of K and L1 penalties (incluing Rho_G, Rho_Z_mu and
#' Rho_Z_Cov), fit LUCID model over all combinations of K and L1 penalties to
#' determine the optimal penalty.
#'
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable
#' should be transformed into dummy variables. If a matrix or data frame, rows
#' represent observations and columns correspond to variables.
#' @param Z Omics data, a numeric matrix or data frame. Rows correspond to observations
#' and columns correspond to variables.
#' @param Y Outcome, a numeric vector. Categorical variable is not allowed. Binary
#' outcome should be coded as 0 and 1.
#' @param CoG Optional, covariates to be adjusted for estimating the latent cluster.
#' A numeric vector, matrix or data frame. Categorical variable should be transformed
#' into dummy variables.
#' @param CoY Optional, covariates to be adjusted for estimating the association
#' between latent cluster and the outcome. A numeric vector, matrix or data frame.
#' Categorical variable should be transformed into dummy variables.
#' @param K Number of latent clusters. An integer greater or equal to 2. If K is
#' a vector, model selection on K is performed
#' @param family Distribution of outcome. For continuous outcome, use "normal";
#' for binary outcome, use "binary". Default is "normal".
#' @param Rho_G A scalar or a vector. This parameter is the LASSO penalty to regularize
#' exposures. If it is a vector, \code{tune_lucid} will conduct
#' model selection and variable selection. User can try penalties from 0 to 1.
#' @param Rho_Z_Mu A scalar or a vector. This parameter is the LASSO penalty to
#' regularize cluster-specific means for omics data (Z). If it is a vector,
#' \code{tune_lucid} will conduct model selection and
#' variable selection. User can try penalties from 1 to 100.
#' @param Rho_Z_Cov A scalar or a vector. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics
#' data (Z). If it is a vector, \code{tune_lucid} will conduct
#' model selection and variable selection. User can try penalties from 0 to 1.
#' @param ... Other parameters passed to \code{est_lucid}
#'
#' @export
#'
#' @return A list:
#' \item{best_model}{the best model over different combination of tuning parameters}
#' \item{tune_list}{a data frame contains combination of tuning parameters and c
#' orresponding BIC}
#' \item{res_model}{a list of LUCID models corresponding to each combination of
#' tuning parameters}
#'
#' @examples
#' \dontrun{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#'
#' # find the optimal model over the grid of K
#' tune_K <- tune_lucid(G = G, Z = Z, Y = Y_normal, useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2:5)
#'
#' # tune penalties
#' tune_Rho_G <- tune_lucid(G = G, Z = Z, Y = Y_normal, useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2, Rho_G = c(0.1, 0.2, 0.3, 0.4))
#' tune_Rho_Z_Mu <- tune_lucid(G = G, Z = Z, Y = Y_normal, useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2, Rho_Z_Mu = c(10, 20, 30, 40))
#' tune_Rho_Z_Cov <- tune_lucid(G = G, Z = Z, Y = Y_normal, useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2, Rho_Z_Cov = c(0.1, 0.2, 0.3))
#'
#' }
tune_lucid <- function(G,
                       Z,
                       Y,
                       CoG = NULL,
                       CoY = NULL,
                       family = "normal",
                       K = 2:5,
                       lucid_model = c("early", "parallel"),
                       Rho_G = 0,
                       Rho_Z_Mu = 0,
                       Rho_Z_Cov = 0,
                       ...){

  if (match.arg(lucid_model) == "early"){
  # combinations of tuning parameters
  tune_list <- expand.grid(K, Rho_G, Rho_Z_Mu, Rho_Z_Cov)
  colnames(tune_list) <- c("K", "Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")
  m <- nrow(tune_list)
  tune_list$BIC <- rep(0, m)
  if(m > 1) {
    cat("Tuning LUCID model \n \n")
  }

  # fit models for each combination
  res_model <- vector(mode = "list",
                      length = m)
  for(i in 1:m) {
    fit <- try(est_lucid(G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG,
                         CoY = CoY,
                         family = family,
                         lucid_model = "early",
                         K = tune_list[i, 1],
                         Rho_G = tune_list[i, 2],
                         Rho_Z_Mu = tune_list[i, 3],
                         Rho_Z_Cov = tune_list[i, 4],
                         ...))
    if("try-error" %in% class(fit)) {
      tune_list[i, 5] <- NA
    } else {
      tune_list[i, 5] <- summary_lucid(fit)$BIC
    }
    res_model[[i]] <- fit
  }
  x <- min(tune_list[, 5], na.rm = TRUE)
  if(is.na(x)) {
    stop("LUCID model fails to converge given current tuning parameters")
  }
  best_model <- res_model[[which(tune_list[, 5]== x)]]
  return(list(best_model = best_model,
              tune_list = tune_list,
              res_model = res_model))
  }

######LUCID in parallel, needs work (from model_select)....######
  else{
    # use bic to conduct model selection
    tune_lucid <- function(K_list, G, Z, Y, ...) {
      # change K_list into a matrix
      K_matrix <- as.matrix(expand.grid(K_list))
      if(min(K_matrix) < 2) {
        stop("minimum K should be 2")
      }
      bic <- rep(0, nrow(K_matrix))
      model_list <- vector(mode = "list", length = nrow(K_matrix))
      for(i in 1:nrow(K_matrix)) {
        model_list[[i]] <- EM_lucid(G = G, Z = Z, Y = Y,
                                    K = K_matrix[i, ],
                                    ...)
        bic[i] <- cal_bic(model_list[[i]])
      }
      model_opt_index <- which(bic == min(bic))
      K_matrix <- cbind(K_matrix, bic)
      return(list(tune_K = K_matrix,
                  model_list = model_list,
                  model_opt = model_list[[model_opt_index]]))
    }
  }
}
