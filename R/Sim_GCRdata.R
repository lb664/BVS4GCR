# 135 characters #####################################################################################################################
#' @title Sim_GCRdata
#'
#' @description Simulate data from the Gaussian Copula Regression model
#'
#' @param m Integer, number of responses
#' @param n Integer, number of samples
#' @param p Integer, number of predictors
#' @param pFixed Integer, number of covariates (including the intercept). If \code{pFixed = 1} only the intercept will be simulated
#' @param responseType \code{m}-dimensional character vector specifying the response's type. Response's types currently supported are 
#' \code{c("Gaussian", "binary", "ordinal", "binomial", "negative binomial")}
#' @param extras List specifying response-specific arguments. These are the one-lag correlation \code{rhoY} in the autoregressive 
#' model used to simulate the responses, the correlation \code{rhoX} between the predictors (for both, see details in 
#' \insertCite{Rothman2010;textual}{BVS4GCR}), \code{fixedEffect} the vector of the regression coefficients for the fixed 
#' predictors (if simulated) where the first element is assumed to be the intercept, the mean \code{muEffect} and variance 
#' \code{varEffect} of the Gaussian distribution used to simulate the regression coefficients, the probabilities \code{s1Lev} and 
#' \code{s2Lev} utilized to induce sparsity on the regression coefficients, the standard deviation \code{sdResponse} of the response 
#' variables, the number of categories \code{nCat} of each ordinal categorical response (if simulated), the cut-off points 
#' \code{cutPoint} used to simulate the ordinal variables (if simulated), \code{nTrial} and \code{negBinomPar} the numer of trials in 
#' the binomial experiment and the over-dispersion parameter of the negative binomial distribution. See also the examples below and 
#' details in \insertCite{Alexopoulos2021;textual}{BVS4GCR}
#' @param std Logical parameter to standardise the predictors before the analysis. Default value is set at \code{TRUE} 
#' @param seed Seed used to initialise the simulation
#'
#' @details
#' The parameters \code{s1Lev} and \code{s2Lev} in the \code{extras} list can be specified following the fact that 
#' ((1 - \code{s2Lev}) * \code{p}) predictors are irrelevant for all \code{m} responses and that each relevant predictor is common 
#' across (\code{s1Lev} * \code{m}) predictors \insertCite{Rothman2010}{BVS4GCR}. Type of responses currently supported: (continuous) 
#' Gaussian, (discrete) binary, ordinal categorical and count (binomial and negative binomial distributions)
#'
#' @export
#'
#' @return The value returned is a list object \code{list(Y, X, B, R, par)}
#' \itemize{
#'  \item{\code{Y}}{ \code{n} times \code{m} matrix of responses}
#'  \item{\code{X}}{ \code{n} times (\code{pFixed} + \code{p}) matrix of covariates (including the intercept if simulated) and 
#'  predictors on which Bayesian Variable Selection is performed}
#'  \item{\code{B}}{ \code{m} times (\code{pFixed} + \code{p}) matrix of regression coefficients}
#'  \item{\code{R}}{ \code{m} times \code{m} correlation matrix used to simulate the responses}
#'  \item{\code{extra}}{ List of all parameters used in the simulation \code{list(rhoY, rhoX, fixedEffect, muEffect, varEffect, 
#'        s1Lev, s2Lev, sdResponse, nCat, cutPoint, nTrial, negBinomPar, std, seed)} } }
#'
#' @references
#' \insertAllCited{}
#'
### 100 characters ################################################################################
#' @examples
#' # Example 1: Simulate a combination of Gaussian, binary and ordinal responses
#'
#' d <- Sim_GCRdata(m = 4, n = 500, p = 30, pFixed = 1, responseType = c("Gaussian", "binary", 
#'                  "ordinal", "ordinal"), 
#'                  extras = list(rhoY = 0.8, rhoX = 0.7, fixedEffect = 1, muEffect = 
#'                                c(1, rep(0.5, 3)), varEffect = c(1, rep(0.2, 3)), s1Lev = 0.15, 
#'                                s2Lev = 0.95, sdResponse = rep(1, 4), nCat = c(3, 4), cutPoint = 
#'                                cbind(c(0.5, NA, NA), c(0.5, 1.2, 2))), 
#'                  seed = 28061971)
#'
#' par(mfrow = c(2, 2))
#' hist(d$Y[,1], main = "Gaussian", xlab = expression(Y[1]), prob = TRUE, ylab = "")
#' barplot(prop.table(table(d$Y[, 2])), main = "binary", xlab = expression(Y[2]))
#' barplot(prop.table(table(d$Y[, 3])), main = "ordinal (3 categories)", xlab = expression(Y[3]))
#' barplot(prop.table(table(d$Y[, 4])), main = "ordinal (4 categories)", xlab = expression(Y[4]))
#' 
#' par(mfrow = c(2, 2))
#' plot(d$B[, 1], ylab=expression(beta[1]), main=expression(Y[1]), xlab="Pred.", pch=16, cex=.7)
#' plot(d$B[, 2], ylab=expression(beta[2]), main=expression(Y[2]), xlab="Pred.", pch=16, cex=.7)
#' plot(d$B[, 3], ylab=expression(beta[3]), main=expression(Y[3]), xlab="Pred.", pch=16, cex=.7)
#' plot(d$B[, 4], ylab=expression(beta[4]), main=expression(Y[4]), xlab="Pred.", pch=16, cex=.7)
#'
#'
#' # Example 2: As in Example 1 with more categories for the forth response 
#'
#' d <- Sim_GCRdata(m = 4, n = 500, p = 30, pFixed = 1, responseType = c("Gaussian", "binary", 
#'                  "ordinal", "ordinal"), 
#'                  extras = list(rhoY = 0.8, rhoX = 0.7, fixedEffect = 1, muEffect = 
#'                                c(1, rep(0.5, 3)), varEffect = c(1, rep(0.2, 3)), s1Lev = 0.15, 
#'                                s2Lev = 0.95, sdResponse = rep(1, 4), nCat = c(3, 5), cutPoint = 
#'                                cbind(c(0.5, 1, NA, NA), c(0.5, 0.7, 1.5, 2))), 
#'                  seed = 28061971)
#'
#' par(mfrow = c(2, 2))
#' hist(d$Y[, 1], main = "Gaussian", xlab = expression(Y[1]), prob = TRUE, ylab = "")
#' barplot(prop.table(table(d$Y[, 2])), main = "binary", xlab = expression(Y[2]))
#' barplot(prop.table(table(d$Y[, 3])), main = "ordinal (3 categories)", xlab = expression(Y[3]))
#' barplot(prop.table(table(d$Y[, 4])), main = "ordinal (5 categories)", xlab = expression(Y[4]))
#'
#'
#' # Example 3: Simulate a combination of Gaussian and count responses
#'
#' d <- Sim_GCRdata(m = 4, n = 1000, p = 100, pFixed = 2, responseType = c("Gaussian", "binomial", 
#'                  "negative binomial", "negative binomial"), 
#'                  extras = list(rhoY = 0.8, rhoX = 0.7, fixedEffect = c(-0.5, 0.5), muEffect = 
#'                  c(1, rep(0.5, 3)), varEffect = c(1, rep(0.2, 3)), s1Lev = 0.05, s2Lev = 0.95, 
#'                  sdResponse = rep(1, 4), nTrial = 10, negBinomPar = c(0.5, 0.75)), 
#'                  seed = 28061971)
#'
#' par(mfrow = c(2, 2))
#' hist(d$Y[,1], main = "Gaussian", xlab = expression(Y[1]), prob = TRUE, ylab = "")
#' hist(table(d$Y[, 2]), main = "binomial", xlab = expression(Y[2]))
#' hist(table(d$Y[, 3]), main = "negative binomial", xlab = expression(Y[3]))
#' hist(table(d$Y[, 4]), main = "negative binomial", xlab = expression(Y[4]))
#'
#' @importFrom MASS mvrnorm


Sim_GCRdata <- function(m, n, p, pFixed, responseType, 
                        extras = list(rhoY = 0.8, rhoX = 0.7, 
                                      fixedEffect = rep(1, pFixed), muEffect = rep(0, m), varEffect = rep(1, m), 
                                      s1Lev = 0.10, s2Lev = 1, sdResponse = rep(1, m), 
                                      nCat = NULL, cutPoint = NULL, nTrial = NULL, negBinomPar = NULL), 
                        std = TRUE, seed = 31122021)
{
  
  if (std == FALSE)
  {
    sdXY <- c(FALSE, FALSE, FALSE, FALSE)
  } else {
    sdXY <- c(FALSE, TRUE, FALSE, FALSE)
  }
  
  p <- p + pFixed
  mC <- mO <- mOrdinal <- mBinary <- mBinom <- mNegBinom <- 0
  
  for (k in 1 : m)
  {
    if (responseType[k] == "Gaussian")
    {
      mC <- mC + 1
    }
    if (responseType[k] == "binary")
    {
      mBinary <- mBinary + 1
      mO <- mO + 1
    }
    if (responseType[k] == "ordinal")
    {
      mOrdinal <- mOrdinal + 1
      mO <- mO + 1
    }
    if (responseType[k] == "binomial")
    {
      mBinom <- mBinom + 1
      mO <- mO + 1
    }
    if (responseType[k] == "negative binomial")
    {
      mNegBinom <- mNegBinom + 1
      mO <- mO + 1
    }
  }
  
  YSim <- Z <- matrix(NA, n, m)
  
  ########## Creating covariance matrix for responses ##########
  
  Sigma <- matrix(NA, m, m)
  
  for (k in 1 : m)
  {
    for (kk in 1 : m)
    {
      Sigma[k, kk] <- extras$rhoY ^(abs(k - kk))
    }
  }
  
  ########## Creating band covariance for predictors ##########
  
  corX <- matrix(0, p - pFixed, p - pFixed)
  
  for (j in 1 : (p - pFixed))
  {
    for (jj in 1 : (p - pFixed))
    {
      corX[j, jj] <- extras$rhoX ^(abs(j - jj))
    }
  }
  
  ########## Simulating predictors ##########
  
  XSim <- matrix(NA, n, p)
  
  if (pFixed >= 1)
  {
    XSim[, 1] <- rep(1, n)
    if (pFixed >= 2)
    {
      XSim[, 2 : pFixed] <- rnorm(n * length(2 : pFixed))
    }
  }
  
  for (i in 1 : n)
  {
    XSim[i, (pFixed + 1) : p] <- rmvn(1, rep(0, p - pFixed), corX)
  }
  
  ########## Simulating effects ##########
  
  pStar <- p - pFixed
  K <- matrix(0, pStar, m)
  W <- matrix(0, pStar, m)
  
  for (k in 1 : m)
  {
    K[, k] <- rbinom(pStar, 1, extras$s1Lev)
    W[, k] <- rnorm(pStar, extras$muEffect[k], sqrt(extras$varEffect[k]))
  }
  
  Q <- matrix(0, pStar, m)
  rowInds <- rbinom(pStar, 1, extras$s2Lev)
  Q[which(rowInds != 0), ] <- 1
  BInds <- t(W * K * Q)
  BSim <- matrix(0, m, p)
  
  if (pFixed > 0)
  {
    fixedEffect <- extras$fixedEffect
    for (k in 1 : m)
    {
      BSim[k, c(1 : pFixed, which(BInds[k, ] != 0) + pFixed)] <- c(fixedEffect, BInds[k, c(which(BInds[k, ] != 0))])
    }
  } else {
    for (k in 1 : m)
    {
      BSim[k, c(which(BInds[k, ] != 0))] <- c(BInds[k, c(which(BInds[k, ] != 0))])
    }
  }
  
  BSimVec <- as.vector(t(BSim))
  R <- cov2cor(Sigma)
  D <- diag(extras$sdResponse)
  R <- D %*% R %*% D
  XB <- matrix(NA, n, m)
  
  for (i in 1 : n)
  {
    Z[i, ] <- mvrnorm(1, rep(0, m), R)
    for (k in 1 : m)
    {
      XB[i, k] <- sum(XSim[i, ] * BSim[k, ])
    }
  }
  
  ########## Simulating Gaussian responses ##########
  
  if (mC > 0)
  {
    YSim[, which(responseType == "Gaussian")] <- Z[, which(responseType == "Gaussian")] + XB[, which(responseType == "Gaussian")]
  }
  
  ########## Simulating binary, categorical, count responses by specifying suitable cdfs ##########
    
  if (mO > 0)
  {
    if (mBinary > 0)
    {
      for (i in 1 : n)
      {
        for (j in which(responseType == "binary"))
        {
          if ((Z[i, j] > -Inf) & (Z[i, j] <= -XB[i, j]))
          {
            YSim[i, j] <- 0
          }
          if ((Z[i, j] > -XB[i, j]) & (Z[i, j] <= Inf))
          {
            YSim[i, j] <- 1
          }
        }
      }
    }
    
    ########## Ordinal categorical responses ##########
    
    if (mOrdinal > 0)
    {
      countCat <- 0
      for (k in which(responseType == "ordinal"))
      {
        countCat <- countCat + 1
        nCat <- extras$nCat[countCat]
        for (i in 1 : n)
        {
          if ((Z[i, k] <= -XB[i, k]) & (Z[i, k] >- Inf))
          {
            YSim[i, k] <- 0
          }
          if ((Z[i, k] > -XB[i, k]) & (Z[i, k] <= extras$cutPoint[1, countCat] - XB[i, k]))
          {
            YSim[i, k] <- 1
          }
          if (nCat == 3)
          {
            if ((Z[i, k] > extras$cutPoint[1, countCat] - XB[i, k]))
            {
              YSim[i, k] <- 2
            }
          } else {
            if ((Z[i, k] > extras$cutPoint[1, countCat] - XB[i, k]) & (Z[i, k] <= extras$cutPoint[2, countCat] - XB[i, k]))
            {
              YSim[i, k] <- 2
            }
            if (nCat == 4)
            {
              if ((Z[i, k] > extras$cutPoint[2, countCat] - XB[i, k]))
              {
                YSim[i, k] <- 3
              }
            } else {
              if ((Z[i, k] > extras$cutPoint[2, countCat] - XB[i, k] ) & (Z[i, k] <= extras$cutPoint[3, countCat] - XB[i,k]))
              {
                YSim[i, k] <- 3
              }
              if (nCat == 5)
              {
                if ((Z[i, k] > extras$cutPoint[3, countCat] - XB[i, k]))
                {
                  YSim[i, k] <- 4
                }
              }
            }
          }
        }
      }
    }
    
    ########## Binomial responses ##########
    
    if (mBinom > 0)
    {
      nTrial <- extras$nTrial
      countBinom <- 0
      for (k in which(responseType == "binomial"))
      {
        countBinom <- countBinom + 1
        for (i in 1 : n)
        {
          c <- 0
          while((Z[i, k] > qnorm(pbinom(c, nTrial[countBinom], exp(XB[i, k]) / (1 + exp(XB[i,k]))))) | 
                (Z[i, k] < qnorm(pbinom(c - 1, nTrial[countBinom], exp(XB[i, k]) / (1 + exp(XB[i, k]))))))
          {
            c <- c + 1
            # cat(c, "\n")
          }
          YSim[i, k] <- c
        }
      }
    }
    
    ########## Negative binomial responses ##########
    
    if (mNegBinom > 0)
    {
      negBinomPar <- extras$negBinomPar
      countNegBinom <- 0
      for (k in which(responseType == "negative binomial"))
      {
        countNegBinom <- countNegBinom + 1
        for (i in 1 : n)
        {
          c <- 0
          while ((Z[i, k] > qnorm(pnbinom(c, negBinomPar[countNegBinom], 1 - (exp(XB[i, k]) / (1 + exp(XB[i, k])))))) | 
                 (Z[i, k] < qnorm(pnbinom(c - 1, negBinomPar[countNegBinom], 1 - (exp(XB[i, k]) / (1 + exp(XB[i, k])))))))
          {
            c <- c + 1
            # cat(c, "\n")
          }
          YSim[i, k] <- c
        }
      }
    }
  }
  
  colnames_Y <- responseType
  rownames_Y <- paste0("i", seq(1 : n))
  colnames_X <- paste0("P", seq(1 : p))
  rownames_X <- paste0("i", seq(1 : n))
  
  rownames(YSim) <- rownames_Y
  colnames(YSim) <- colnames_Y
  rownames(XSim) <- rownames_X
  colnames(XSim) <- colnames_X
  rownames(BSim) <- colnames_Y
  colnames(BSim) <- colnames_X
  rownames(Sigma) <- colnames_Y
  colnames(Sigma) <- colnames_Y
  rownames(R) <- colnames_Y
  colnames(R) <- colnames_Y
  
  ########## Data transformation ##########
  
  XSim <- scale(XSim, center = sdXY[1], scale = sdXY[2])
  YSim <- scale(YSim, center = sdXY[3], scale = sdXY[4])
  
  extra <- c(extras, std = std, seed = seed)
  
  output <- list(Y = YSim, X = XSim, B = t(BSim), Sigma = Sigma, R = R, extra = extra)
  
  return(output)
}
