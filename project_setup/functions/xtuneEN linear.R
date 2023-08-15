library(glmnet)
library(crayon)
library(selectiveInference)

# xtune function ------

#' Estimate Elastic-net Regression
#'
#' The function \code{xtuneEN()} estimates the elastic-net regression with the given parameters.
#'
#' @param X A \code{matrix} of the independent variables.
#' @param Y A \code{vector} of the dependent variable.
#' @param U A \code{matrix} of the additional independent variables.
#' @param Z A \code{matrix} of the covariates for the penalty vector.
#' @param c A \code{numeric} of the elastic-net mixing parameter.
#' @param sigma.square A \code{numeric} of the variance of the errors.
#' @param alpha.est.init A \code{vector} of the initial value of the penalty parameter.
#' @param maxstep A \code{numeric} of the maximum number of iterations for the outer loop.
#' @param tolerance A \code{numeric} of the tolerance for the outer loop.
#' @param maxstep_inner A \code{numeric} of the maximum number of iterations for the inner loop.
#' @param tolerance_inner A \code{numeric} of the tolerance for the inner loop.
#' @param compute.likelihood A \code{logical} of whether to compute the likelihood score.
#' @param verbosity A \code{logical} of whether to print the iteration progress.
#' @param standardize A \code{logical} of whether to standardize the input data.
#' @param intercept A \code{logical} of whether to include an intercept term in the model.
#' @param epsilon A \code{numeric} of the parameter for the inner loop.
#'
#' @return A \code{list} containing the estimated coefficients, penalty vector, lambda, alpha, number of iterations, variance and likelihood score.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100), ncol=5)
#' Y <- rnorm(100)
#' xtuneEN(X, Y)
#'
#' @importFrom stats glm
#' 
xtuneEN <- function(X, 
                    Y, 
                    U = NULL,
                    Z = as.matrix(rep(1, ncol(X))),
                    c = 0.5,
                    sigma.square = estimateVariance(X,Y), 
                    alpha.est.init = rep(0,ncol(Z)+1), 
                    maxstep = 100, 
                    tolerance = 0.001,
                    maxstep_inner = 50, 
                    tolerance_inner = 0.1, 
                    compute.likelihood = FALSE, 
                    verbosity = FALSE,
                    standardize = TRUE, 
                    intercept = TRUE,
                    epsilon = 3){ 
  
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)-1
  
  ##------------ Elastic-net regression
  ## Initialize 
  alpha.old = alpha.est.init 
  likelihood.score = c() 
  k = 1 
  
  ## calculate alpha.max
  alpha.max = max(abs(colSums(X*Y)))/ (c * n)
  
  ## reparameterize Z
  Z.original = Z
  Z = sweep(Z,2,colMeans(Z))
  Z = cbind(rep(1,p), Z)
  
  while(k < maxstep){
    # Given alpha, update theta 
    lambda = exp(Z%*%alpha.old)
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
    Sigma_y = sigma.square * diag(n) + (t(t(X) * c(gamma))) %*% t(X)
    theta = colSums(X * solve(Sigma_y, X))
    
    # Compute likelihood
    if (compute.likelihood == TRUE) {
      likelihood.score = c(likelihood.score, 
                           approx_likelihood.EN(alpha.old,c,
                                                X, Y, Z, sigma.square))
    }
    
    # Given theta, update alpha
    cat(yellow$italic$bold("Start estimating alpha:\n"))
    update.result <- update_alpha.EN(X, Y, Z,c=c, 
                                     alpha.old = alpha.old, 
                                     alpha.max = alpha.max, 
                                     epsilon = epsilon,
                                     sigma.square = sigma.square, 
                                     theta = theta, 
                                     maxstep_inner = maxstep_inner,
                                     tolerance_inner = tolerance_inner)
    
    alpha.new <- update.result$alpha.est
    
    # Check convergence
    if (sum(abs(alpha.new - alpha.old)) < tolerance) {
      cat(red$bold("Done!\n"))
      break
    }
    cat("Difference between alpha_old and alpha_new:",sum(abs(alpha.new - alpha.old)),"\n")
    alpha.old <- alpha.new
    
    # Track iteration progress
    if (verbosity == TRUE) {
      #cat("#-----------------Iteration ", k, " Done -----------------#\n",
      #    sep = "")
      cat(green$italic("#---"),green$bold("Outer loop Iteration",k,"Done"),green$italic("---#\n"),sep = "")
    }
    k <- k + 1
  }
  
  tauEst = exp(Z%*%alpha.old)
  pen_vec = tauEst * sigma.square/n
  
  pen_vec[pen_vec>1e6] <- 1e6
  if(is.null(U)) {
    pen_vec_cov = pen_vec
  } else {
      pen_vec_cov = c(pen_vec, rep(0,ncol(U)))
  }
  C = sum(pen_vec_cov)/p
  
  obj <- glmnet(cbind(X, U), 
                Y, 
                alpha = c,
                family="gaussian",
                lambda = C, 
                penalty.factor = pen_vec_cov,
                standardize = standardize,
                intercept = intercept)
  
  cus.coef <- tryCatch(
    coef(obj,x=cbind(X, U),y=Y,alpha = c, exact=TRUE,s=C, 
         penalty.factor = pen_vec_cov,
         standardize=standardize,intercept = intercept), 
    error = function(e) e
  )
  
  if(!inherits(cus.coef, "error")) {
    return(list(beta.est = cus.coef, 
                penalty.vector = pen_vec_cov, 
                lambda = C, 
                alpha.est = alpha.old, 
                n_iter = k - 1, 
                sigma.square = sigma.square, 
                likelihood.score = likelihood.score))
  } else {
    cus.coef <- data.frame(s1 = rep(NA, length(colnames(X))))
    rownames(cus.coef) <- colnames(X)
    
    return(list(beta.est = cus.coef, 
                penalty.vector = pen_vec_cov, 
                lambda = C, 
                alpha.est = alpha.old, 
                n_iter = k - 1, 
                sigma.square = sigma.square, 
                likelihood.score = likelihood.score))
  }
  
  # return(list(beta.est = cus.coef, penalty.vector = pen_vec_cov, lambda = C, alpha.est = alpha.old, 
  #             n_iter = k - 1, sigma.square = sigma.square, likelihood.score = likelihood.score))
  
}


# Estimate Variance Function ------
estimateVariance <- function(X, Y, n_rep = 5) {
  Y <- as.double(drop(Y))
  dimY = dim(Y)
  nrowY = ifelse(is.null(dimY), length(Y), dimY[1])
  if (nrowY < 10) {
    stop("Need at least 10 observations to estimate variance")
  }
  
  temp = array(NA, n_rep)
  for (i in 1:n_rep) {
    c = suppressWarnings(estimateSigma(X, Y)$sigmahat^2)
    temp[i] = ifelse(is.infinite(c), NA, c)
  }
  return(mean(temp, na.rm = T))
}


# approx_likelihood function ----
approx_likelihood.EN <- function(to_estimate,c, X, Y, Z, sigma.square.est) {
  n = nrow(X)
  lambda = exp(Z %*% to_estimate) 
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  K = sigma.square.est * diag(n) + X %*% diag(c(gamma)) %*% t(X)
  logdetK = determinant(K)$modulus[1]
  part1 = t(Y) %*% solve(K, Y)
  normapprox = 1/2 * (part1 + logdetK)
  return(-as.numeric(normapprox))
}

# update_alpha.EN function ----
update_alpha.EN <- function(X, Y, Z,c,alpha.old, alpha.max, 
                            epsilon, sigma.square, theta, maxstep_inner,
                            tolerance_inner) {
  ## initial
  alpha.inner.old = alpha.old
  k_inner = 1
  n = nrow(X)
  p = ncol(X)
  while (k_inner < maxstep_inner) {
    # given alpha update delta
    lambda = exp(Z %*% alpha.inner.old) 
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
    
    sd_y <- sqrt(var(Y) * (n - 1)/n)
    C = sum(1/gamma)/p * sd_y * sigma.square/n
    delta.est = coef(glmnet(X, Y, alpha = 0, penalty.factor = 1/gamma, lambda = C,
                            standardize = F, intercept = FALSE))[-1]
    
    ## given delta update alpha
    alpha.inner.new <- optim(alpha.old, likelihood.alpha.theta.EN,likelihood.alpha.theta.gradient.EN,
                             c =c,Z = Z, theta = theta, delta = delta.est,method = "L-BFGS-B", upper = c(alpha.max*epsilon, rep(Inf, length(alpha.old)-1)))$par
    if (sum(abs(alpha.inner.new - alpha.inner.old)) < tolerance_inner) {
      break
    }
    cat(blue$italic("#-----------------"),magenta("Inner loop Iteration",k_inner,"Done"),blue$italic("-----------------#\n"),sep = "")
    k_inner = k_inner + 1
    alpha.inner.old <- alpha.inner.new
  }
  return(list(alpha.est = alpha.inner.old, inner_iter = k_inner))
}

# likelihood.alpha.theta.EN function ----
likelihood.alpha.theta.EN <- function(Z,c, alpha, theta, delta) {
  lambda = exp(Z %*% alpha) 
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  
  return(as.numeric(t(theta) %*% gamma + delta^2 %*% (1/gamma)))
}

# likelihood.alpha.theta.gradient.EN function ----
likelihood.alpha.theta.gradient.EN <- function(Z,c, alpha, theta, delta) {
  lambda = exp(Z %*% alpha) 
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  
  dev_gamma = theta - delta^2/(gamma^2)
  dev_gamma_alpha = as.vector((-2*(2*(1-c) + 2*c^2*lambda))/(2*lambda*(1-c) + (c*lambda)^2)^2 *lambda) *Z
  return(crossprod(dev_gamma, dev_gamma_alpha))
}

relative_change<-function(a,b){
  return((a-b)/b)
}