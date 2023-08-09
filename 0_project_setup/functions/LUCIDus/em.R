#' @title Fit LUCID models with one or multiple omics layers
#' @description EM algorithm to estimate LUCID with one or multiple omics layers
#' @param lucid_model Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel
#' @param G an N by P matrix representing exposures
#' @param Z If "early", an N by M matrix; If "parallel", a list, each element is a matrix with N rows
#' @param Y a length N vector
#' @param CoG an N by V matrix representing covariates to be adjusted for G -> X
#' @param CoY an N by K matrix representing covariates to be adjusted for X -> Y
#' @param K If "early", an integer;If "parallel",an integer vector, same length as Z
#' @param modelNames a vector of strings specifies the geometric model of omics
#' data. If NULL, See more in ?mclust::mclustModelNames
#' @param useY logical, if TRUE, EM algorithm fits a supervised LUCID; otherwise
#' unsupervised LUCID.
#' @param max_itr Maximum iterations of the EM algorithm. If the EM algorithm iterates
#' more than max_itr without converging, the EM algorithm is forced to stop.
#' @param max_tot.itr Max number of total iterations for \code{est_lucid} function.
#' \code{est_lucid} may conduct EM algorithm for multiple times if the algorithm
#' fails to converge.(to be done)
#' @param Rho_G A scalar. This parameter is the LASSO penalty to regularize
#' exposures. If user wants to tune the penalty, use the wrapper
#' function \code{lucid}
#' @param Rho_Z_Mu A scalar. This parameter is the LASSO penalty to
#' regularize cluster-specific means for omics data (Z). If user wants to tune the
#' penalty, use the wrapper function \code{lucid}
#' @param Rho_Z_Cov A scalar. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics
#' data (Z). If user wants to tune the penalty, use the wrapper function \code{lucid}
#' @param family The distribution of the outcome
#' @param seed Random seed to initialize the EM algorithm
#' @param init_impute Method to initialize the imputation of missing values in
#' LUCID. "mclust" will use \code{mclust:imputeData} to implement EM Algorithm
#' for Unrestricted General Location Model to impute the missing values in omics
#' data; \code{lod} will initialize the imputation via relacing missing values by
#' LOD / sqrt(2). LOD is determined by the minimum of each variable in omics data.
#' @param init_par For "early", an interface to initialize EM algorithm (to be done);
#' For "parallel", mclust is the default for quick convergence
#' @param verbose A flag indicates whether detailed information for each iteration
#' of EM algorithm is printed in console. Default is FALSE.(to be done)
#'
#' @import mclust
#'
#' @return A list contains the object below:
#' 1. res_Beta: estimation for G->X associations
#' 2. res_Mu_Sigma: estimation for X->Z associations
#' 3. res_Gamma: estimation for X->Y associations
#' 4. loglik: log likelihood of LUCID model
#' 5. z: posterior cluster assignment for each observation
#' 6. K: number of latent clusters
#' 7. N: number of observations
#'
#' @examples
#' i <- 1008
#' set.seed(i)
#' G <- matrix(rnorm(500), nrow = 100)
#' Z1 <- matrix(rnorm(1000),nrow = 100)
#' Z2 <- matrix(rnorm(1000), nrow = 100)
#' Z3 <- matrix(rnorm(1000), nrow = 100)
#' Z <- list(Z1 = Z1, Z2 = Z2, Z2 = Z3)
#' CoY <- matrix(rnorm(200), nrow = 100)
#' CoG <- matrix(rnorm(200), nrow = 100)
#' Y <- rnorm(100)
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2), CoG = CoG, CoY = CoY,
#' lucid_model = "parallel",
#' family = "normal",
#' modelName = "VVV",
#' seed = i,
#' init_impute = "mclust",
#' init_par = "mclust",
#' useY = TRUE)

#'
#' @export
est_lucid <- function(lucid_model = c("early", "parallel"),
                      G, Z, Y, CoG = NULL, CoY = NULL, K,
                      modelName = "VVV",
                      useY = TRUE,
                      tol = 1e-3,
                      max_itr = 1e3,
                      max_tot.itr = 1e4,
                      Rho_G = 0,
                      Rho_Z_Mu = 0,
                      Rho_Z_Cov = 0,
                      family = c("normal", "binary"),
                      seed = 123,
                      init_impute = c("mclust", "lod"),
                      init_par = c("mclust", "random"),
                      verbose = FALSE) {
  ## 1.0 basic setup, apply for both lucid models ====
  init_impute <- match.arg(init_impute)
  init_par <- match.arg(init_par)
  Select_G <- FALSE
  Select_Z <- FALSE
  if(Rho_G != 0) {
    Select_G <- TRUE
  }
  if(Rho_Z_Mu != 0 | Rho_Z_Cov != 0) {
    Select_Z <- TRUE
  }

  ## 1.1 check data format, apply for both lucid models ====

  #check for G and get Gnames
  if(is.null(G)) {
    stop("Input data 'G' is missing")
  } else {
    if(!is.matrix(G)) {
      G <- as.matrix(G)
      if(!is.numeric(G)) {
        stop("Input data 'G' should be numeric; categorical variables should be transformed into dummies")
      }
    }
  }
  if(is.null(colnames(G))){
    Gnames <- paste0("G", 1:ncol(G))
  } else {
    Gnames <- colnames(G)
  }
  colnames(G) <- Gnames

  #check for Y and get Ynames
  if(is.null(Y)) {
    stop("Input data 'Y' is missing")
  } else {
    if(!is.matrix(Y)) {
      Y <- as.matrix(Y)
      if(!is.numeric(Y)) {
        stop("Input data 'Y' should be numeric; binary outcome should be transformed them into dummies")
      }
      if(ncol(Y) > 1) {
        stop("Only continuous 'Y' or binary 'Y' is accepted")
      }
    }
  }
  if(is.null(colnames(Y))) {
    Ynames <- "outcome"
  } else {
    Ynames <- colnames(Y)
  }
  colnames(Y) <- Ynames
  if(family == "binary") {
    if(!(all(Y %in% c(0, 1)))) {
      stop("Binary outcome should be coded as 0 and 1")
    }
  }

  #check for CoG and CoY and get their names
  CoGnames <- NULL
  if(!is.null(CoG)) {
    if(!is.matrix(CoG)) {
      CoG <- as.matrix(CoG)
      if(!is.numeric(CoG)) {
        stop("Input data 'CoG' should be numeric; categroical variables should be transformed into dummies")
      }
    }
    if(is.null(colnames(CoG))) {
      CoGnames <- paste0("CoG", 1:ncol(CoG))
    } else {
      CoGnames <- colnames(CoG)
    }
    colnames(CoG) <- CoGnames
  }

  CoYnames <- NULL
  if(!is.null(CoY)) {
    if(!is.matrix(CoY)) {
      CoY <- as.matrix(CoY)
      if(!is.numeric(CoY)) {
        stop("Input data 'CoY' should be numeric; categorical variables should be transformed into dummies")
      }
    }
    if(is.null(colnames(CoY))) {
      CoYnames <- paste0("CoY", 1:ncol(CoY))
    } else {
      CoYnames <- colnames(CoY)
    }
    colnames(CoY) <- CoYnames
  }

  if (match.arg(lucid_model) == "early"){

    # ========================== Early Integration ==========================

    # 1. basic setup for estimation function under early =============
    family <- match.arg(family)
    modelName = match.arg(modelName)

    ## 1.1 check data format ====  special for Z under early
    if(is.null(Z)) {
      stop("Input data 'Z' is missing")
    }else {
      if(is.list(Z)) {
        stop("Input data 'Z' should not be a list for early integration LUCID but a N by M matrix!")
      }
      if(!is.matrix(Z)) {
        Z <- as.matrix(Z)
        if(!is.numeric(Z)) {
          stop("Input data 'Z' should be numeric")
        }
      }
    }

    if(is.null(colnames(Z))){
      Znames <- paste0("Z", 1:ncol(Z))
    } else {
      Znames <- colnames(Z)
    }


    ## 1.2 record input dimensions, family function ====
    N <- nrow(Y)
    dimG <- ncol(G)
    dimZ <- ncol(Z);
    dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
    dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
    G <- cbind(G, CoG)
    Gnames <- c(Gnames, CoGnames)
    family.list <- switch(family, normal = normal(K = K, dimCoY),
                          binary = binary(K = K, dimCoY))
    Mstep_Y <- family.list$f.maxY
    switch_Y <- family.list$f.switch


    ## 1.3. check missing pattern ====
    na_pattern <- check_na(Z)
    if(na_pattern$impute_flag) {
      # initialize imputation
      if(init_impute == "mclust") {
        cat("Intializing imputation of missing values in 'Z' via the mix package \n\n")
        invisible(capture.output(Z <- mclust::imputeData(Z, seed = seed)))
        Z[na_pattern$indicator_na == 3, ] <- NA
      }
      if(init_impute == "lod") {
        cat("Intializing imputation of missing values in 'Z' via LOD / sqrt(2) \n\n")
        Z <- apply(Z, 2, fill_data_lod)
        colnames(Z) <- Znames
      }

    }


    # 2. EM algorithm for LUCID ================
    tot.itr <- 0
    convergence <- FALSE
    while(!convergence && tot.itr <= max_tot.itr) {
      if(tot.itr > 0) {
        seed <- seed + 10
      }
      set.seed(seed)

      ## 2.1 initialize model parameters ====

      # initialize beta
      res.beta <- matrix(data = runif(K * (dimG + dimCoG + 1)), nrow = K)
      res.beta[1, ] <- 0

      # initialize mu and sigma
      # initialize by mclust
      if(init_par == "mclust") {
        cat("Initialize LUCID with mclust \n")
        invisible(capture.output(mclust.fit <- Mclust(Z[na_pattern$indicator_na != 3, ],
                                                      G = K,
                                                      modelNames = modelName)))
        if(is.null(mclust.fit)) {
          stop("mclust failed for specified model - please set modelName to `NULL` to conduct automatic model selection ")
        }
        if(is.null(modelName)){
          model.best <- mclust.fit$modelName
        } else{
          model.best <- modelName
        }
        res.mu <- t(mclust.fit$parameters$mean)
        res.sigma <- mclust.fit$parameters$variance$sigma
        # browser()

      } else { # initialize by random guess
        cat("Initialize LUCID with random values from uniform distribution \n")
        if(is.null(modelName)){
          model.best <- "VVV"
          cat("GMM model for LUCID is not specified, 'VVV' model is used by default \n")
        } else{
          model.best <- modelName
        }
        res.mu <- matrix(runif(dimZ * K, min = -0.5, max = 0.5),
                         nrow = K)

        res.sigma <- gen_cov_matrices(dimZ = dimZ, K = K)
      }

      # initialize family specific parameters gamma
      res.gamma <- family.list$initial.gamma(K, dimCoY)


      # start EM algorithm
      cat("Fitting Early Integration LUCID model",
          paste0("(", "K = ", K, ", Rho_G = ", Rho_G, ", Rho_Z_Mu = ", Rho_Z_Mu, ", Rho_Z_Cov = ", Rho_Z_Cov, ")"),
          "\n")
      res.loglik <- -Inf
      itr <- 0
      while(!convergence && itr <= max_itr){
        itr <- itr + 1
        tot.itr <- tot.itr + 1
        check.gamma <-  TRUE

        # 2.2 E-step ====
        # calculate log-likelihood for observation i being assigned to cluster j
        new.likelihood <- Estep_early(beta = res.beta,
                                mu = res.mu,
                                sigma = res.sigma,
                                gamma = res.gamma,
                                G = G,
                                Z = Z,
                                Y = Y,
                                CoY = CoY,
                                N = N,
                                K = K,
                                family.list = family.list,
                                itr = itr,
                                useY = useY,
                                dimCoY = dimCoY,
                                ind.na = na_pattern$indicator_na)
        # normalize the log-likelihood to probability
        res.r <- t(apply(new.likelihood, 1, lse_vec))

        if(!all(is.finite(res.r))){
          cat("iteration", itr,": EM algorithm collapsed: invalid estiamtes due to over/underflow, try LUCID with another seed \n")
          break
        } else{
          if(isTRUE(verbose)) {
            cat("iteration", itr,": E-step finished.\n")
          }
        }


        # 2.3 M-step - parameters ====
        # update model parameters to maximize the expected likelihood
        invisible(capture.output(new.beta <- Mstep_G(G = G,
                                                     r = res.r,
                                                     selectG = Select_G,
                                                     penalty = Rho_G,
                                                     dimG = dimG,
                                                     dimCoG = dimCoG,
                                                     K = K)))
        new.mu.sigma <- Mstep_Z(Z = Z,
                                r = res.r,
                                selectZ = Select_Z,
                                penalty.mu = Rho_Z_Mu,
                                penalty.cov = Rho_Z_Cov,
                                model.name = model.best,
                                K = K,
                                ind.na = na_pattern$indicator_na,
                                mu = res.mu)
        if(is.null(new.mu.sigma$mu)){
          cat("variable selection failed, try LUCID with another seed \n")
          break
        }
        if(useY){
          new.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames)
          check.gamma <- is.finite(unlist(new.gamma))
        }


        # 2.4 M step - impute missing values ====
        if(na_pattern$impute_flag){
          Z <- Istep_Z(Z = Z,
                       p = res.r,
                       mu = res.mu,
                       sigma = res.sigma,
                       index = na_pattern$index, lucid_model = "early")
        }


        # 2.5 control step ====
        check.value <- all(is.finite(new.beta),
                           is.finite(unlist(new.mu.sigma)),
                           check.gamma)

        if(!check.value){
          cat("iteration", itr,": Invalid estimates, try LUCID with another seed \n")
          break
        } else{
          res.beta <- new.beta
          res.mu <- new.mu.sigma$mu
          res.sigma <- new.mu.sigma$sigma
          if(useY){
            res.gamma <- new.gamma
          }

          new.loglik <- sum(rowSums(res.r * new.likelihood))

          if(Select_G) {
            new.loglik <- new.loglik - Rho_G * sum(abs(res.beta))
          }
          if(Select_Z) {
            new.loglik <- new.loglik - Rho_Z_Mu * sum(abs(res.mu)) - Rho_Z_Cov * sum(abs(res.sigma))
          }
          if(isTRUE(verbose)) {
            if(Select_G | Select_Z) {
              cat("iteration", itr,": M-step finished, ", "penalized loglike = ", sprintf("%.3f", new.loglik), "\n")
            } else{
              cat("iteration", itr,": M-step finished, ", "loglike = ", sprintf("%.3f", new.loglik), "\n")
            }
          } else {
            cat(".")
          }

          if(abs(res.loglik - new.loglik) < tol){
            convergence <- TRUE
            cat("Success: Early Integration LUCID converges!", "\n\n")
          }
          res.loglik <- new.loglik
        }
      }
    }

    # 3. summarize results ===============
    if(!useY){
      res.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, K = K, CoYnames = CoYnames)
    }

    res.likelihood <- Estep_early(beta = res.beta,
                            mu = res.mu,
                            sigma = res.sigma,
                            gamma = res.gamma,
                            G = G,
                            Z = Z,
                            Y = Y,
                            family.list = family.list,
                            itr = itr,
                            CoY = CoY,
                            N = N,
                            K = K,
                            dimCoY = dimCoY,
                            useY = useY,
                            ind.na = na_pattern$indicator_na)
    res.r <- t(apply(res.likelihood, 1, lse_vec))


    res.loglik <- sum(rowSums(res.r * res.likelihood))

    if(Select_G) {
      res.loglik <- res.loglik - Rho_G * sum(abs(res.beta))
    }
    if(Select_Z) {
      res.loglik <- res.loglik - Rho_Z_Mu * sum(abs(res.mu)) - Rho_Z_Cov * sum(abs(res.sigma))
    }
    # browser()
    pars <- switch_Y(beta = res.beta, mu = res.mu, sigma = res.sigma, gamma = res.gamma, K = K)
    res.r <- res.r[, pars$index]
    colnames(pars$beta) <- c("intercept", Gnames)
    colnames(pars$mu) <- Znames
    if(Select_G){
      tt1 <- apply(pars$beta[, -1], 2, range)
      selectG <- abs(tt1[2, ] - tt1[1, ]) > 0.001
    } else{
      selectG <- rep(TRUE, dimG)
    }
    if(Select_Z){
      tt2 <- apply(pars$mu, 2, range)
      selectZ <- abs(tt2[2, ] - tt2[1, ]) > 0.001
    } else{
      selectZ <- rep(TRUE, dimZ)
    }
    results <- list(pars = list(beta = pars$beta, 
                                mu = pars$mu, 
                                sigma = pars$sigma, 
                                gamma = pars$gamma),
                    K = K, 
                    var.names =list(Gnames = Gnames, 
                                    Znames = Znames, 
                                    Ynames = Ynames), 
                    modelName = model.best, 
                    likelihood = res.loglik,
                    post.p = res.r, 
                    family = family,
                    select = list(selectG = selectG, selectZ = selectZ), 
                    useY = useY,
                    Z = Z,
                    init_impute = init_impute,
                    init_par = init_par,
                    Rho = list(Rho_G = Rho_G,
                               Rho_Z_Mu = Rho_Z_Mu,
                               Rho_Z_Cov = Rho_Z_Cov)
    )
    class(results) <- c("early_lucid")
    return(results)
  }
  else{
    # ========================== LUCID IN PARALLEL ==========================


    ## 1.1 check data format ==== special for Z  under parallel
    if(is.null(Z)) {
      stop("Input data 'Z' is missing")
    }
    if(!is.list(Z)) {
      stop("Input data 'Z' should be a list for LUCID in Parallel!")
    } else {
    for(i in 1:length(Z)) {
      if(!is.matrix(Z[[i]])) {
        Z[[i]] <- as.matrix(Z[[i]])
        if(!is.numeric(Z[[i]])) {
          stop("Input data 'Z' should be numeric")
        }
      }
    }}

    Znames <- vector("list", length(Z))
    for(i in 1:length(Z)) {
      if(is.null(colnames(Z[[i]]))){
        Znames[[i]] <- paste0("Z_", i, "_", 1:ncol(Z[[i]]))
      } else {
        Znames[[i]] <- colnames(Z[[i]])
      }}

    ## 1.2 basic setup ====
    N <- nrow(G)
    nOmics <- length(Z)
    nG <- ncol(G)
    nZ <- as.integer(sapply(Z, ncol))


    if (match.arg(family) == "normal"){
      family = "gaussian"
    }else{
      family = "binomial"
    }
    modelNames =  rep(modelName, length(K))

    dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
    dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
    # combine G and CoG to adjust for CoG
    if(!is.null(CoG)) {
      G <- cbind(G, CoG)
      Gnames <- c(Gnames, CoGnames)
    }



    ## 1.3. check missing pattern & initial imputation ====
    na_pattern <- vector("list", nOmics)
    for(i in 1:nOmics) {
      na_pattern[[i]] <- check_na(Z[[i]])
      if(na_pattern[[i]]$impute_flag) {
        # initialize imputation
        if(init_impute == "mclust") {
          cat("Intializing imputation of missing values in 'Z' via the mix package \n\n")
          invisible(capture.output(Z[[i]] <- mclust::imputeData(Z[[i]], seed = seed)))
          Z[[i]][na_pattern[[i]]$indicator_na == 3, ] <- NA
        }
        if(init_impute == "lod") {
          cat("Intializing imputation of missing values in 'Z' via LOD / sqrt(2) \n\n")
          Z[[i]] <- apply(Z[[i]], 2, fill_data_lod)
          colnames(Z[[i]]) <- Znames[[i]]
        }

      }}

    tot.itr <- 0
    flag_converge <- FALSE
    while(!flag_converge && tot.itr <= max_tot.itr) {
      if(tot.itr > 0) {
        seed <- seed + 10
      }
      set.seed(seed)

    ## 1.4 initialize model parameters ====
    Mu_Sigma <- initialize_Mu_Sigma(K = K, Z = Z, modelNames = modelNames, na_pattern = na_pattern)
    Mu <- Mu_Sigma$Mu
    Sigma <- Mu_Sigma$Sigma
    Beta <- vector(mode = "list", length = nOmics)
    for(i in 1:nOmics) {
      invisible(capture.output(temp_fit <- nnet::multinom(Mu_Sigma$z[[i]] ~ G)))
      Beta[[i]] <- coef(temp_fit)
    }
    # Beta <- initialize_Beta(K = K, nG = nG)
    Gamma <- initialize_Delta(K = K, CoY = CoY, family = family,
                              z = Mu_Sigma$z, Y = Y)
    loglik <- -Inf



    ## 2.1 start EM algorithm ====
    cat("Fitting LUCID in Parallel model",
        paste0("(", "K = ", K, ", Rho_G = ", Rho_G, ", Rho_Z_Mu = ", Rho_Z_Mu, ", Rho_Z_Cov = ", Rho_Z_Cov, ")"),
        "\n")
    itr <- 0
    while(!flag_converge & itr < max_itr) {
      itr <- itr + 1
      tot.itr <- tot.itr + 1

      # E-step
      Estep_array <- Estep(G = G, Z = Z, Y = Y,
                          Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Gamma,
                          family = family, useY = useY, na_pattern = na_pattern)
      Estep_r <- Estep_to_r(Estep_array = Estep_array,
                            K = K,
                            N = N)

      if(!all(is.finite(Estep_r))){
        cat("iteration", itr,": EM algorithm collapsed: invalid estiamtes due to over/underflow, try LUCID with another seed \n")
        break
      } else{
        if(isTRUE(verbose)) {
          cat("iteration", itr,": E-step finished.\n")
        }
      }

      # M-step 1
      #Mstep_GtoX added penalty, but how is G exluded?
      res_Beta <- Mstep_GtoX(G = G, r = Estep_r, selectG = Select_G, penalty = Rho_G, K = K, N = N)
      res_Mu_Sigma <- Mstep_XtoZ(Z = Z, r = Estep_r, K = K,
                                modelNames = modelNames, N = N, na_pattern = na_pattern)
      if(useY) {
        res_Gamma <- Mstep_XtoY(Y = Y, CoY = CoY,r = Estep_r, K = K, N = N,
                              family = family)
      }

      if(is.null(res_Mu_Sigma$Mu)){
        cat("variable selection failed, try LUCID with another seed \n")
        break
      }

      #control step ====
      check.value <- all(is.finite(unlist(res_Beta$Beta)),
                        is.finite(unlist(res_Mu_Sigma$Mu)),
                        if(useY){is.finite(unlist(res_Gamma$Gamma))})
      if(!check.value){
        cat("iteration", itr,": Invalid estimates, try LUCID with another seed \n")
        break
      }else{
        # update parameters
        Beta <- res_Beta$Beta
        Mu <- res_Mu_Sigma$Mu
        Sigma <- res_Mu_Sigma$Sigma
        if(useY) {
          Gamma <- res_Gamma$Gamma
        }
      }

      post.p <- vector(mode = "list", length = nOmics)
      for(i in 1:nOmics) {
        post.p[[i]] = compute_res_r(r = Estep_r, N = N,layer = i)
      }

      # M step 2 - impute missing values ====
      for(i in 1:nOmics) {
        if(na_pattern[[i]]$impute_flag){
          Z[[i]] <- Istep_Z(Z = Z[[i]],
                            p = post.p[[i]],
                            mu = Mu[[i]],
                            sigma = Sigma[[i]],
                            index = na_pattern[[i]]$index, lucid_model = "parallel")
        }}


      # check convergence
      loglik_update <- cal_loglik(Estep_array = Estep_array,
                                  Estep_r = Estep_r)
      ### regularity to be added
      ###  if(Select_G) {
      ### new.loglik <- new.loglik - Rho_G * sum(abs(res.beta))
      ### }
      ### if(Select_Z) {
      ### new.loglik <- new.loglik - Rho_Z_Mu * sum(abs(res.mu)) - Rho_Z_Cov * sum(abs(res.sigma))
      ### }

      if(abs(loglik - loglik_update) < tol) {
        flag_converge <- TRUE
        cat("Success: LUCID in parallel converges!", "\n\n")
      } else {
        loglik <- loglik_update
        cat(paste0("iteration ", itr, ": log-likelihood = ", loglik_update, "\n"))
      }
    }
    }



    # 3. summarize results ===============
    if(!useY) {
      res_Gamma <- Mstep_XtoY(Y = Y, CoY = CoY, r = Estep_r, K = K, N = N,
                              family = family)
      Gamma <- res_Gamma$Gamma
    }

    results <- list(res_Beta = res_Beta,
                     res_Mu_Sigma = res_Mu_Sigma,
                     res_Delta = res_Gamma,
                     loglik = loglik_update,
                     z = Estep_r,
                     K = K,
                     N = N)
    
    class(results) <- c("lucid_parallel")
    return(results)
  }}
