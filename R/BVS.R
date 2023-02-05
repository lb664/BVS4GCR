# 135 characters #####################################################################################################################
#' @title Bayesian Variable Selection
#' @description Internal function that performs one Markov chain Monte Carlo iteration to sample from the posterior distribution of 
#' all unknowns involved in the Bayesian Variable Selection step. For details, see \insertCite{Alexopoulos2021;textual}{BVS4GCR}
#'
#' @references
#' \insertAllCited{}


########## BVS ##########

BVS <- function(Y, NA_Y, X, XB, GammaCurr, BCurr, ZCurr, RCurr, invRCurr, d, 
                  responseType, pFixed, 
                  nCat, cutPoint, negBinomPar, 
                  hyperPar, samplerPar)
{
  mlogpval <- samplerPar$mlogpval
  
  aPi <- hyperPar$aPi
  bPi <- hyperPar$aPi
  tau <- hyperPar$tau
  nTrial <- hyperPar$nTrial
  
  GammaProp <- rep(0, length(GammaCurr))
  GammaPropList <- list()
  countBinom <- 0
  countNegBinom <- 0
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  pStar <- p - pFixed
  ind <- 1 : m
  
  mu <- Sigma <- TU <- TL <- matrix(NA, n, m)
  TUProp <- TLProp <- rep(NA, n)
  ZTilde <- ZCurr + XB
  
  for (k in 1 : m)
  {
    NAInd <- NA_Y[[k]] 
    start <- (k - 1) * p + 1
    end <- k * p
    GammaInd <- which(GammaCurr[start : end] > 0)
    
    ########## Proposing Gamma ##########
    
    GammaPropList <- Sample_gamma(GammaCurr[(start + pFixed) : end], pStar, mlogpval[, k], samplerPar)
    logGammaPropListRatio <- GammaPropList[[2]]
     
    if (pFixed > 0)
    {
      GammaPropList[[1]] <- c(1 : pFixed, GammaPropList[[1]] + pFixed)
    }
     
    if (length(GammaPropList[[1]] > 0))
    {
      setOne <- GammaPropList[[1]] + (k - 1) * p
      setZero <- ((start : end)[-GammaPropList[[1]]])
    } else {
      setOne <- NULL
      setZero <- (start : end)
    }
    
    if (pFixed > 0)
    {
      GammaStarInit <- GammaPropList[[1]][-c(1 : pFixed)]
    } else {
      GammaStarInit <- GammaPropList[[1]]
    }
    
    logpriorGammaRatio <- lbeta(length(GammaStarInit) + aPi, pStar - length(GammaStarInit) + bPi) - 
                          lbeta(sum(GammaCurr[start : end]) - pFixed + aPi, pStar - sum(GammaCurr[start : end]) + bPi)
    
    if (length(GammaInd) > 0)
    {
      XGamma <- as.matrix(X[, GammaInd])
      B_k <- BCurr[start : end]
      B_kGamma <- B_k[GammaInd]
      XBGamma <- XGamma %*% B_kGamma
      tXGamma <- t(XGamma)
    } else {
      XGamma <- matrix(0, n, 1)
      B_k <- rep(0, pStar)
      B_kGamma <- 0
      XBGamma <- rep(0, n)
      tXGamma <- t(XGamma)
    }
    
    if (length(GammaPropList[[1]]) > 0)
    {
      if (pFixed > 0)
      {
        GammaStar <- c(1 : pFixed, GammaStarInit)
      } else {
        GammaStar <- c(GammaStarInit)
      }
      XGammaStar <- as.matrix(X[, GammaStar])
      tXGammaStar <- t(XGammaStar)
    } else {
      if (pFixed > 0)
      {
        GammaStar <- c(1 : pFixed)
        XGammaStar <- as.matrix(X[, GammaStar])
        tXGammaStar <- t(XGammaStar)
      } else {
        GammaStar <- NULL
        XGammaStar <- matrix(0, n, 1)
        tXGammaStar <- t(XGammaStar)
      }
    }
    
    ########## Proposing Beta ##########
    
    RR <- RCurr[k, -k] %*% solve(RCurr[-k, -k])
    mu[, k] <- RR %*% t(ZCurr[, -k])
    Sigma[, k] <- RCurr[k, k] - RR %*% RCurr[-k, k]
    if (responseType[k] == "Gaussian")
    {
      ZTilde[, k] <- d[k] * ZCurr[, k] + XB[, k]
      const <- invRCurr[k, k] / (d[k] * d[k])
      T <- rep(0, n)
      for (kk in ind[-k])
      {
        T <- T + ZCurr[, kk] * invRCurr[k, kk] / d[k]
      }
      if (length(GammaStar) > 1) 
      {
        Q <- diag(rep(1 /(tau), length(GammaStar))) + const * tXGammaStar %*% XGammaStar
        invQ <- solve(Q)
      } else {
        Q <- as.matrix(1 /(tau) ) + const * tXGammaStar %*% XGammaStar
        invQ <- solve(Q)
      }
      muProp <- invQ %*% tXGammaStar %*% (const * ZTilde[, k] + T)
      if (length(GammaInd) > 1)
      {
        Q_OLD <- diag(rep(1 /(tau), length(GammaInd))) + const * tXGamma %*% XGamma
        invQ_OLD <- solve(Q_OLD)
      } else {
        Q_OLD <- as.matrix(1 /(tau)) + const * tXGamma %*% XGamma
        invQ_OLD <- solve(as.matrix(1 /(tau)) + const * tXGamma %*% XGamma)
      }
      muProp_OLD <- invQ_OLD %*% tXGamma %*% (const * ZTilde[, k] + T)
      BProp <- mvrnorm(1, muProp, invQ)
    } else {
      ZTilde[, k] <- ZCurr[, k] + XB[, k]
      const <- invRCurr[k, k]
      T <- rep(0, n)
      for (kk in ind[-k])
      {
        T <- T + ZCurr[, kk] * invRCurr[k, kk]
      }
      if (length(GammaStar) > 1)
      {
        Q <- diag(rep(1 /(tau), length(GammaStar))) + const * tXGammaStar %*% XGammaStar
        invQ <- solve(Q)
      } else {
        Q <- as.matrix(1 /(tau)) + const*tXGammaStar%*%XGammaStar
        invQ <- solve(Q)
      }
      muProp <- invQ %*% tXGammaStar %*% (const * ZTilde[, k] + T)
      if (length(GammaInd) > 1)
      {
        Q_OLD <- diag(rep(1 /(tau), length(GammaInd))) + const * tXGamma %*% XGamma
        invQ_OLD <- solve(Q_OLD)
      } else {
        Q_OLD <- as.matrix(1 /(tau)) + const * tXGamma %*% XGamma
        invQ_OLD <- solve(as.matrix(1 /(tau))  + const * tXGamma %*% XGamma)
      }
      muProp_OLD <- invQ_OLD %*% tXGamma %*% (const * ZTilde[, k] + T)
      BProp <- mvrnorm(1, muProp, invQ)
    }
    if (length(GammaStar) > 0)
    {
      XBProp<-(as.matrix(X[, GammaStar]) %*% BProp)[, 1]
    } else {
      XBProp <- rep(0,n)
    }
    
    ########## Updating Beta and Gamma ##########
    
    GammaProp[setOne] <- 1
    GammaProp[setZero] <- 0
    
    if (responseType[k] == "Gaussian")
    {
      likeHolmes <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus + 0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) - 0.5 * length(GammaInd) * log(tau)
      likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muProp) %*% Q %*% as.matrix(muProp) - 0.5 * length(GammaStar) * log(tau)
      logAccRatio <- likeHolmesProp - likeHolmes + logGammaPropListRatio + logpriorGammaRatio
      if (log(runif(1)) < logAccRatio)
      {
        if (length(GammaStar) > 0)
        {
          BCurr[setOne] <- BProp
          GammaCurr[setOne] <- 1
        }
        GammaCurr[setZero] <- 0
        BCurr[setZero] <- 0
        XBGamma <- XBProp
        XB[, k] <- XBProp
      }
      for (i in 1 : n)
      {
        if (is.na(Y[i, k]))
        {
          ZCurr[i, k] <- rnorm(1, RR %*% (ZCurr[i, -k]), sqrt(RCurr[k, k] - RR %*% RCurr[-k, k]))
        } else {
          ZCurr[i, k] <- (Y[i, k] - XBGamma[i]) / d[k]
        }
      }
    } else {
      if (responseType[k] == "binary")
      {
        likeHolmes <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus + 0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) - 0.5*length(GammaInd) * log(tau)
        likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muProp) %*% Q %*% as.matrix(muProp) - 0.5 * length(GammaStar) * log(tau)
        logAccRatio <- likeHolmesProp - likeHolmes + logGammaPropListRatio + logpriorGammaRatio
        if (log(runif(1)) < logAccRatio)
        {
          if (length(GammaStar) > 0)
          {
            BCurr[setOne] <- BProp
            GammaCurr[setOne] <- 1
          }
          GammaCurr[setZero] <- 0
          BCurr[setZero] <- 0
          XBGamma <- XBProp
          XB[, k] <- XBProp
        }
        ZCurr[, k] <- ZTilde[, k] - XB[, k]
        TU[, k] <- cutPoint[[k]][Y[, k] + 2] - XBGamma
        TL[, k] <- cutPoint[[k]][Y[, k] + 1] - XBGamma
      }
      
      if (responseType[k] == "ordinal") 
      {
        likeHolmes <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus + 0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) - 0.5 * length(GammaInd) * log(tau)
        likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muProp) %*% Q %*% as.matrix(muProp) - 0.5 * length(GammaStar) * log(tau)
        logAccRatio <- likeHolmesProp - likeHolmes + logGammaPropListRatio + logpriorGammaRatio
        if (log(runif(1)) < logAccRatio)
        {
          if (length(GammaStar) > 0)
          {
            BCurr[setOne] <- BProp
            GammaCurr[setOne] <- 1
          }
          GammaCurr[setZero] <- 0
          BCurr[setZero] <- 0
          XBGamma <- XBProp
          XB[, k] <- XBProp
        }
        ZCurr[, k] <- ZTilde[, k] - XB[, k]
        TU[, k] <- cutPoint[[k]][Y[, k] + 2] - XBGamma
        TL[, k] <- cutPoint[[k]][Y[, k] + 1] - XBGamma
      }
      
      if (responseType[k] == "binomial")
      {
        countBinom <- countBinom + 1
        TUProp <- qnorm(pbinom(Y[, k], nTrial[countBinom], 1 / (1 + exp(-XBProp))))
        TU[, k] <- qnorm(pbinom(Y[, k], nTrial[countBinom], 1 / (1 + exp(-XBGamma))))
        TLProp <- qnorm(pbinom(Y[, k] - 1, nTrial[countBinom], 1 / (1 + exp(-XBProp))))
        TL[, k] <- qnorm(pbinom(Y[, k] - 1, nTrial[countBinom], 1 / (1 + exp(-XBGamma))))
      }
      
      if (responseType[k] == "negative binomial")
      {
        countNegBinom <- countNegBinom + 1
        TUProp <- qnorm(pnbinom(Y[, k], negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBProp)))))
        TU[, k] <- qnorm(pnbinom(Y[, k], negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBGamma)))))
        TLProp <- qnorm(pnbinom(Y[, k] - 1, negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBProp)))))
        TL[, k]<- qnorm( pnbinom(Y[, k] - 1, negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBGamma)))))
        approxInd <- union(which(TU[, k] == TL[, k]), which(TL[, k] == Inf))
        if ((length(approxInd) > 0))
        {
          TL[, k][approxInd] <- qnorm(0.9999999999999)
          TU[, k][approxInd] <- qnorm(0.99999999999999)          
        }
        approxInd <- union(which(TUProp == TLProp), which(TLProp == Inf))
        if ((length(approxInd) > 0))
        {
          TLProp[approxInd] <- qnorm(0.9999999999999)
          TUProp[approxInd] <- qnorm(0.99999999999999)          
        }
      }
      
      if (responseType[k] == "binomial" || responseType[k] == "negative binomial")
      {
        diffLikProp <- pnorm((TUProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TLProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))
        diffLikProp[which(diffLikProp <= 0)] <- pnorm(Inf) - pnorm(8.2)
        logLikProp <- sum(log(diffLikProp))
        diffLik <- ((pnorm((TU[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TL[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))))
        diffLik[which(diffLik <= 0)] <- pnorm(Inf) - pnorm(8.2)
        logLik_OLD <- sum(log(diffLik))
        logPrior <- sum(dnorm(BProp, 0, sqrt(tau), log = TRUE)) - sum(dnorm(B_kGamma, 0, sqrt(tau), log = TRUE))
        logAccRatio <- logPrior + dmvnorm(B_kGamma, muProp_OLD, invQ_OLD, log = TRUE) - dmvnorm(BProp, muProp, invQ, log = TRUE) + logpriorGammaRatio
        logAccRatio <- logAccRatio + logLikProp - logLik_OLD + logGammaPropListRatio
        if (log(runif(1)) < logAccRatio)
        {
          if (length(GammaStar) > 0)
          {
            BCurr[setOne] <- BProp
            GammaCurr[setOne] <- 1
            GammaCurr[setZero] <- 0
          }
          BCurr[setZero] <- 0
          XBGamma <- XBProp
          XB[, k] <- XBProp
          TU[, k] <- TUProp
          TL[, k] <- TLProp
        }
      }
      
      ########## Updating cut-points ##########
      
      if (responseType[k] == "ordinal")
      {
        GammaInd <- which(GammaCurr[start : end] > 0)
        if (length(GammaInd) > 0)
        {
          XGamma <- as.matrix(X[, GammaInd])
          B_k <- BCurr[start : end]
          B_kGamma <- B_k[GammaInd]
          XBGamma <- XGamma %*% B_kGamma
          tXGamma <- t(XGamma)
        } else {
          XGamma <- matrix(0, n, 1)
          B_k <- rep(0, pStar)
          B_kGamma <- 0
          XBGamma <- rep(0, n)
          tXGamma <- t(XGamma)
        }
        cutPoint_TMP <- rep(NA, length(cutPoint[[k]]) - 2)
        cutPoint_TMP[1] <- 0
        cutPoint_TMP[2 : length(cutPoint_TMP)] <- log(cutPoint[[k]][3 : (length(cutPoint[[k]]) - 1)] - cutPoint[[k]][2 : (length(cutPoint[[k]]) - 2)])
        cutPointProp_TMP <- cutPoint_TMP[2 : length(cutPoint_TMP)] + rnorm(length(2 : length(cutPoint_TMP)), 0, sqrt(0.1))
        cutPointProp <- cutPoint[[k]]
        cutPointProp[3] <- exp(cutPointProp_TMP)[1]
        if (nCat[k] >= 4)
        {
          for (jOrd in 4 : nCat[k])
          {
            cutPointProp[jOrd] <- sum(exp(cutPointProp_TMP)[1 : (jOrd - 2)])
          }
        }
        TUProp <- cutPointProp[Y[, k] + 2] - XBGamma
        TU[, k] <- cutPoint[[k]][Y[, k] + 2] - XBGamma
        TLProp <- cutPointProp[Y[, k] + 1] - XBGamma
        TL[, k] <- cutPoint[[k]][Y[, k] + 1] - XBGamma
        priorRatio <- sum(dnorm(cutPointProp_TMP, 0, sqrt(10), log = TRUE)) - sum(dnorm(cutPoint_TMP[2 : length(cutPoint_TMP)], 0, sqrt(10), log = TRUE))
        diffLikProp <- pnorm((TUProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TLProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))
        diffLikProp[which(diffLikProp == 0)] < pnorm(Inf) - pnorm(8.2)
        logLikProp <- sum(log(diffLikProp))
        diffLik <- ((pnorm((TU[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TL[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))))
        diffLik[which(diffLik == 0)] <- pnorm(Inf) - pnorm(8.2)
        logLik_OLD <- sum(log(diffLik))
        likRatio <- logLikProp - logLik_OLD
        
        if (log(runif(1)) < (likRatio + priorRatio))
        {
          cutPoint[[k]][3] <- cutPointProp[3]
          if (nCat[k] >= 4)
          {
            for (jOrd in 4 : nCat[k])
            {
              cutPoint[[k]][jOrd] <- cutPointProp[jOrd]
            }
          }
          TU[, k] <- TUProp
          TL[, k] <- TLProp
        }
      }
      
      ########## Updating negative binomial overdispersion parameter ##########
      
      if (responseType[k] == "negative binomial")
      {
        GammaInd <- which(GammaCurr[start : end] > 0)
        if (length(GammaInd) > 0)
        {
          XGamma <- as.matrix(X[, GammaInd])
          B_k <-BCurr[start : end]
          B_kGamma <- B_k[GammaInd]
          XBGamma <- XGamma %*% B_kGamma
          tXGamma <-t(XGamma)
        } else {
          XGamma <- matrix(0, n, 1)
          B_k <- rep(0, pStar)
          B_kGamma <- 0
          XBGamma <- rep(0, n)
          tXGamma <- t(XGamma)
        }
        RR <- RCurr[k, -k] %*% solve(RCurr[-k, -k])
        mu[, k] <- RR %*% t(ZCurr[, -k])
        Sigma[, k] <- RCurr[k, k] - RR %*% RCurr[-k, k]
        rhoProp <- log(negBinomPar[countNegBinom]) + rnorm(1, 0, sqrt(0.05))
        TUProp <- qnorm(pnbinom(Y[, k], exp(rhoProp), 1 - (1 / (1 + exp(-XBGamma)))))
        TLProp <- qnorm(pnbinom(Y[, k] - 1, exp(rhoProp), 1 - (1 / (1 + exp(-XBGamma)))))
        approxInd <- union(which(TUProp == TLProp), which(TLProp == Inf))
        if ((length(approxInd) > 0))
        {
          TLProp[approxInd] <- qnorm(0.9999999999999)
          TUProp[approxInd] <- qnorm(0.99999999999999)
        }
        priorRatio <- dgamma(exp(rhoProp), 2, rate = 1, log = TRUE) - dgamma(negBinomPar[countNegBinom], 2, rate = 1, log = TRUE)
        diffLikProp <- pnorm((TUProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TLProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))
        diffLikProp[which(diffLikProp == 0)] <- pnorm(Inf) - pnorm(8.2)
        logLikProp <- sum(log(diffLikProp))
        diffLik <- ((pnorm((TU[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TL[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))))
        diffLik[which(diffLik == 0)] <- pnorm(Inf) - pnorm(8.2)
        logLik_OLD <- sum(log(diffLik))
        likRatio <- logLikProp - logLik_OLD
        if (log(runif(1)) < (likRatio + priorRatio + rhoProp - log(negBinomPar[countNegBinom])))
        {
          negBinomPar[countNegBinom] <- exp(rhoProp)
          TU[, k] <- TUProp
          TL[, k] <- TLProp
        }
      }
      
      ########## Updating Z latent variables ##########
      
      for (i in 1 : n)
      {
        if (!is.na(Y[i, k]))
        {
          ZCurr[i, k] <- rtruncnorm(1, a = TL[i, k], b = TU[i, k], RR %*% (ZCurr[i, -k]), sqrt(RCurr[k, k] - RR %*% RCurr[-k, k]))
          if (is.na(ZCurr[i, k]))
          {
            cat(c(TL[i, k], TU[i, k]), "\n")
          }
        } else {
          ZCurr[i, k] <- rnorm(1, RR %*% (ZCurr[i, -k]), sqrt(RCurr[k, k] - RR %*% RCurr[-k, k]))
        }
      }
    }
  }
  
  output <- list(GammaProp, GammaCurr, BCurr, ZCurr, XB, cutPoint, negBinomPar)
  
  return(output)
}
