# 135 characters #####################################################################################################################
#' @title Sampler of the binary latent variable
#' @description Internal function for the proposal distribution of the binary latent vector for each response based on 
#' \insertCite{Guan2011;textual}{BVS4GCR}. For details, see \insertCite{Alexopoulos2021;textual}{BVS4GCR}
#'
#' @param gammaCurr Current vector of the predictor inclusion indicators of one response
#' @param p Number of predictors subject to variable selection
#' @param mlogpval Vector of minus log10 p-values of the one-at-a-time regressions, used to build the proposal distribution
#' @param samplerPar List of the parameters of the Markov chain Monte Carlo algorithm
#' @keywords internal
#'
#' @references
#' \insertAllCited{}


# CORRECTED (2026-11-06): the original proposal had a 0-based draw used as a 1-based index (no-op proposals), a wrong reverse position for deletions, missing factors for the
# forced moves and an ill-defined swap, so the reported proposal ratio was not the true ratio. See Departures_from_BVS4GCR (BayesCopReg, Article_2/Documents).

########## .Sample_gamma ##########

# Proposal for the inclusion vector of one response (add / delete / swap moves; the position of the predictor to add is drawn from a mixture of a uniform and a
# truncated geometric distribution on the ranking of the excluded predictors). Returns the proposed set, the log of q(gamma | gamma*) / q(gamma* | gamma)
# and the type of move. The proposal probabilities are computed exactly, including the forced moves at the boundaries (empty and full model) and the position of
# the reverse move in the ranking of the NEW model, so that the Metropolis-Hastings ratio is exact.

.Sample_gamma <- function(gammaCurr, p, mlogpval, samplerPar)
{
  probAddDel <- samplerPar$probAddDel
  probAdd <- samplerPar$probAdd
  probMix <- samplerPar$probMix
  geomMean <- samplerPar$geomMean

  ord <- sort(mlogpval, decreasing = TRUE, index.return = TRUE)$ix      # predictors from the strongest to the weakest univariable evidence
  setOne <- which(gammaCurr != 0)
  pgam <- length(setOne)

  # probability of each type of move at a model with q included predictors
  pAdd <- function(q) if (q >= p) 0 else if (q == 0) 1 else probAddDel * probAdd
  pDel <- function(q) if (q == 0) 0 else if (q >= p) 1 else probAddDel * (1 - probAdd)
  pSwap <- function(q) if (q == 0 || q >= p) 0 else 1 - probAddDel
  # position (1-based) of predictor j among the excluded predictors of the model 'set', in the order 'ord'
  posExcl <- function(j, set) which(setdiff(ord, set) == j)

  u <- runif(1)

  if (pgam == 0)
  {
    moveType <- "add"
  } else if (pgam == p) {
    moveType <- "delete"
  } else if (u < probAddDel) {
    moveType <- if (runif(1) < probAdd) "add" else "delete"
  } else {
    moveType <- "swap"
  }

  if (moveType == "add")
  {
    r <- Qrnd(p, pgam, probMix, geomMean)
    addJ <- setdiff(ord, setOne)[r]
    setOneProp <- c(setOne, addJ)
    qFwd <- log(pAdd(pgam)) + log(Qpdf(r, p, pgam, probMix, geomMean))
    qRev <- log(pDel(pgam + 1)) - log(pgam + 1)
  }

  if (moveType == "delete")
  {
    delJ <- setOne[sample.int(pgam, 1)]
    setOneProp <- setdiff(setOne, delJ)
    qFwd <- log(pDel(pgam)) - log(pgam)
    qRev <- log(pAdd(pgam - 1)) + log(Qpdf(posExcl(delJ, setOneProp), p, pgam - 1, probMix, geomMean))
  }

  if (moveType == "swap")
  {
    r <- Qrnd(p, pgam, probMix, geomMean)
    addJ <- setdiff(ord, setOne)[r]
    delJ <- setOne[sample.int(pgam, 1)]
    setOneProp <- c(setdiff(setOne, delJ), addJ)
    qFwd <- log(pSwap(pgam)) + log(Qpdf(r, p, pgam, probMix, geomMean)) - log(pgam)
    qRev <- log(pSwap(pgam)) + log(Qpdf(posExcl(delJ, setOneProp), p, pgam, probMix, geomMean)) - log(pgam)
  }

  setOneProp <- sort(setOneProp)
  output <- list(setOneProp, qRev - qFwd, moveType)

  return(output)
}

########## Qrnd ##########

# Position (1, ..., p - pgam) of the excluded predictor to add: with probability u uniform, otherwise a geometric distribution (mean geomMean) truncated to 1, ..., p - pgam
Qrnd <- function(p, pgam, u, geomMean)
{
  m <- p - pgam

  if (runif(1) < u)
  {
    r <- sample.int(m, 1)
  } else {
    r <- 1 + qgeom(runif(1) * pgeom(m - 1, 1 / geomMean), 1 / geomMean)
  }

  return(r)
}

########## Qpdf ##########

Qpdf <- function(r, p, pgam, u, geomMean)
{
  m <- p - pgam
  ff <- u * (1 / m) + (1 - u) * (dgeom(r - 1, 1 / geomMean) / pgeom(m - 1, 1 / geomMean))

  return(ff)
}
