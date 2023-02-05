# 135 characters #####################################################################################################################
#' @title Sampler of Hyper-Inverse Wishart
#' @description Internal functions to simulate from the Hyper Inverse-Wishart distribution given a decomposable graph. The R code of 
#' this function is a modified version of the Matlab code of \insertCite{Carvalho2007;textual}{BVS4GCR} and 
#' \insertCite{Talhouk2012;textual}{BVS4GCR}. For details, see details, see \insertCite{Alexopoulos2021;textual}{BVS4GCR}
#'
#' @references
#' \insertAllCited{}


########## Sim_HIW ##########

Sim_HIW <- function(G, S, n)
{
  TT <- max(dim(G)[1], dim(G)[2])
  cc <- check_chordal(G)
  check <- cc[[1]]
  order <- cc[[2]]
  if (check == 0)
  {
    stop('graph not decomp')
  }
  cliques <- chordal_to_ripcliques_zo(G, order)
  ss <- seps_resids_hists_zo(cliques)
  Tost_Thi_hat <- g_constrain_zo(S, order, cliques, ss)[[1]]
  sigma_identity <- generate_HIW_g_delta_identity_zo(G, cliques, n, ss)[[1]]
  Sigma <- transform_g_conditional_HIW_no_mcs_zo(sigma_identity, G, cliques, n, diag(rep(1, TT)), Tost_Thi_hat, ss)
  output <- list(Sigma[[1]], Sigma[[2]])
  
  return(output)
}


########## check_chordal ##########

check_chordal <- function(g)
{
  diag(g) <- 1
  p <- dim(g)[1]
  order <- rep(0, p)
  chordal <- 1
  order[1] <- 1
  capU <- 2 : p
  output <- cpp_check_chordal(g, capU, p, chordal, order)   # cpp
  return(output)
}

########## chordal_to_ripcliques_zo ##########

chordal_to_ripcliques_zo <- function(G, order)
{
  p <- dim(G)[1]
  pa <- matrix(0, p, p)
  num_pa <- rep(0, p)
  ladder <- rep(0, p)
  cliques <- matrix(0, p, p)
  index_non_zero_ladder <- 0
  pre_v <- rep(0, p)
  for (i in 2 : p)
  {
    v <- order[i]
    ns <- G[, v]
    sum <- 0
    pa_i <- rep(0, p)
    for (k in 1 : (i - 1))
    {
      aux <- order[k]
      if (ns[aux] != 0)
      {
        sum <- sum + 1
        pa_i[aux] <- 1
      }
    }
    num_pa[i] <- sum
    if (i == 2)
    {
      if ( num_pa[i - 1] >= num_pa[i])
      {
        ladder[i - 1] <- order[i - 1]
        if (ladder[i - 1] != 0)
        {
          index_non_zero_ladder <- index_non_zero_ladder + 1
          reserve <- rep(0, p)
          reserve[ladder[i - 1]] <- 1
          cliques[, index_non_zero_ladder] <- reserve
        }
      }
    }
    if (i != 2)
    {
      if ( num_pa[i - 1] >= num_pa[i])
      {
        ladder[i - 1] <- order[i - 1]
        if (ladder[i - 1] != 0)
        {
          index_non_zero_ladder <- index_non_zero_ladder + 1
          reserve[ladder[i - 1]] <- 1
          cliques[, index_non_zero_ladder] <- reserve
        }
      }
    }
    if (i == p)
    {
      ladder[p] <- order[p]
      if (ladder[p] != 0)
        {
          index_non_zero_ladder <- index_non_zero_ladder + 1
          pa_i[ladder[p]] <- 1
          cliques[, index_non_zero_ladder] <- pa_i
        }
    }
    reserve <- pa_i
  }
  output <- cliques
  
  return(output)
}

########## seps_resids_hists_zo ##########

seps_resids_hists_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  num_cliques <- max(which(colSums(cliques) != 0))
  seps <- matrix(0, p, p)
  resids <- hists <- matrix(0, p, p)
  hists[, 1] <- cliques[, 1]
  if (num_cliques >= 2)
  {
    for (index in 2 : num_cliques)
    {
      hists[, index] <- cpp_union_zo((cliques[, index]) + (hists[, index - 1]))   # cpp
      seps[, index] <- (cliques[, index]) * (hists[, index - 1])
      resids[, index] <- cpp_setdiff_zo((cliques[, index]), (cliques[, index]) - (hists[, index - 1]))   # cpp
    }
  }
  output <- list(seps, resids, hists)
  
  return(output)
}

########## generate_HIW_g_delta_identity_zo ##########

generate_HIW_g_delta_identity_zo <- function(g, cliques, delta,ss)
{
  index_finish <- 0
  index_start <- 0
  num_cliques <- max(which(colSums(cliques) != 0))
  seps <- ss[[1]]
  residuals <- ss[[2]]
  histories <- ss[[3]]
  num_Rj <- rep(0, num_cliques)
  p <- max(dim(g)[1], dim(g)[2])
  perfect_order <- rep(0, p)
  c_1 <- which(cliques[, 1] != 0)
  perfect_order[1 : length(c_1)] <- c_1
  index_finish <- length(c_1)
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      Rj <- which(residuals[, j] != 0)
      num_Rj[j] <- length(Rj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start : index_finish] <- Rj
    }
  }
  rev_perf <- rep(0, p)
  for (i in 1 : p)
  {
    rev_perf[i] <- perfect_order[p - i + 1]
  }
  g_rev_perf <- matrix(0, p, p)
  g_rev_perf <- g[rev_perf, rev_perf]
  Psi <- matrix(0, p, p)
  for (i in 1 : p)
  {
    if (i == p)
    {
      Psi[p, p] <- rchisq(1, df = delta) ^0.5
    } else {
      nu_i <- sum(g_rev_perf[i, (i + 1) : p])
      Psi[i, i] <- rchisq(1, df = delta + nu_i) ^0.5
    }
  }
  for (i in 1 : (p - 1))
  {
    for (j in (i + 1) : p)
    {
      if (g_rev_perf[i, j] == 0)
      {
        Psi[i, j] <- 0
      } else {
        Psi[i, j] <- rnorm(1)
      }
    }
  }
  K_rev_perf <- sigma_id_rev_perf <- matrix(0, p, p)
  K_rev_perf <- t(Psi) %*% Psi
  sigma_id_rev_perf <- solve(K_rev_perf)
  inverse_permute <- rep(0, p)
  for (j in 1 : p)
  {
    inverse_permute[j] <- which(rev_perf == j)
  }
  K_id <- sigma_id <- matrix(0, p, p)
  K_id <- K_rev_perf[inverse_permute, inverse_permute]
  sigma_id <- sigma_id_rev_perf[inverse_permute, inverse_permute]
  output <- list(sigma_id, K_id, sigma_id_rev_perf, K_rev_perf)
  
  return(output)
}

########## g_constrain_zo ##########

g_constrain_zo <- function(B, order, cliques, ss)
{
  p <- dim(B)[1]
  seps <- ss[[1]]
  num_cliques <- max(which(colSums(cliques) != 0))
  sum_big_K_cliques <- sum_big_K_seps <- matrix(0, p, p)
  if (num_cliques >= 1)
  {
    for (j in 1 : num_cliques)
    {
      cj <- which(cliques[, j] != 0)
      B_cj <- B[cj, cj]
      K_cj <- solve(B_cj)
      big_Kj <- matrix(0, p, p)
      big_Kj[cj, cj] <- K_cj
      sum_big_K_cliques <- sum_big_K_cliques + big_Kj
    }
  }
  if (num_cliques >= 1)
  {
    for (j in 1 : num_cliques)
    {
      sj <- which(seps[, j] != 0)
      B_sj <- B[sj, sj]
      if (length(sj) > 0)
      {
        K_sj <- solve(B_sj)
      }
      big_Kj <- matrix(0, p, p)
      if (length(sj) > 0)
      {
        big_Kj[sj, sj] <- K_sj
      }
      sum_big_K_seps <- sum_big_K_seps + big_Kj
    }
  }
  K_hat <- sum_big_K_cliques - sum_big_K_seps
  B_hat <- solve(K_hat + diag(rep(0.00000001, dim(K_hat)[2])))
  output <- list(B_hat, K_hat)
  
  return(output)
}

########## transform_g_conditional_HIW_no_mcs_zo ##########

transform_g_conditional_HIW_no_mcs_zo <- function(sigma_B, g, cliques, delta, B, D, ss)
{
  p <- dim(g)[1]
  num_cliques <- max(which(colSums(cliques) != 0))
  seps <- ss[[1]]
  residuals <- ss[[2]]
  histories <- ss[[3]]
  num_Cj <- rep(0, num_cliques)
  num_Sj <- rep(0, num_cliques)
  num_Rj <- rep(0, num_cliques)
  num_Cj[1] <- sum(cliques[, 1])
  perfect_order <- rep(0, p)
  perfect_order[1 : num_Cj[1]] <- which(cliques[, 1] !=0 )
  index_finish <- num_Cj[1]
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      Cj <- which(cliques[, j] != 0)
      Rj <- which(residuals[, j] != 0)
      Sj <- which(seps[, j] != 0)
      num_Cj[j] <- length(Cj)
      num_Rj[j] <- length(Rj)
      num_Sj[j] <- length(Sj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start : index_finish] <- Rj
    }
  }
  rev_perf <- rep(0, p)
  index <- 0
  for (i in 1 : p)
  {
    rev_perf[i] <- perfect_order[p - index]
    index <- index + 1
  }
  g_rev_perf <- g[rev_perf, rev_perf]
  B_rev_perf <- B[rev_perf, rev_perf]
  D_rev_perf <- D[rev_perf, rev_perf]
  sigma_B_rev_perf <- sigma_B[rev_perf, rev_perf]
  K_B <- solve(sigma_B_rev_perf)
  choleskyK_B <- chol(K_B)
  c1 <- num_Cj[1]
  index_start <- p - c1 + 1
  index_finish <- p
  indexC1_in_Upsilon_D <- rep(0, c1)
  indexC1_in_Upsilon_D <- index_start : index_finish
  Upsilon_D <- matrix(0, p, p)
  B_1 <- D_1 <- Q_1 <- P_1 <- O_1 <- matrix(0, c1, c1)
  B_1 <- B_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  Q_1 <- chol(solve(B_1))
  D_1 <- D_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  P_1 <- chol(solve(D_1))
  O_1 <- solve(Q_1) %*% P_1
  Upsilon_D[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] <- choleskyK_B[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] %*% O_1
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      cj <- num_Cj[j]
      indexCj_in_Upsilon_D <- rep(0, cj)
      Rj <- which(residuals[, j] != 0)
      rj <- num_Rj[j]
      indexRj_in_Upsilon_D <- rep(0, rj)
      indexRj_inOj <- rep(0, rj)
      unsort_indexRj_in_Upsilon_D <- rep(0, rj)
      Sj <- which(seps[, j] != 0)
      sj <- num_Sj[j]
      unsort_indexSj_in_Upsilon_D <- rep(0, sj)
      indexSj_in_Upsilon_D <- rep(0, sj)
      B_j <- D_j <- Q_j <- P_j <- O_j <- matrix(0, cj, cj)
      if (rj > 0)
      {
        for (k in 1 : rj)
        {
          unsort_indexRj_in_Upsilon_D[k] <- which(rev_perf == Rj[k])
        }
        indexRj_in_Upsilon_D <- sort(unsort_indexRj_in_Upsilon_D)
      }
      if (sj > 0)
      {
        for (k in 1 : sj)
        {
          unsort_indexSj_in_Upsilon_D[k] <- which(rev_perf == Sj[k])
        }
      indexSj_in_Upsilon_D <- sort(unsort_indexSj_in_Upsilon_D)
      }
    indexCj_in_Upsilon_D <- c(indexRj_in_Upsilon_D, indexSj_in_Upsilon_D)
    B_j <- (B_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D])
    Q_j <- chol(solve(B_j))
    D_j <- D_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D]
    P_j <- chol(solve(D_j))
    O_j <- solve(Q_j) %*% P_j
    if (rj > 0)
    {
      indexRj_inOj <- 1 : rj
    } else {
      indexRj_inOj <- NULL
    }
    if ((rj + 1) <= cj)
    {
      indexSj_inOj <- (rj + 1) : cj
    } else {
      indexSj_inOj <- NULL
    }
    Upsilon_D[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D] <- as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% as.matrix(O_j[indexRj_inOj, indexRj_inOj])
    Upsilon_D[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] <- as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% O_j[indexRj_inOj,indexSj_inOj] + choleskyK_B[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] %*% as.matrix(O_j[indexSj_inOj, indexSj_inOj])
    }
  }
  K_D_rev_perf <- sigma_D_rev_perf <- matrix(0, p, p)
  K_D_rev_perf <- t(Upsilon_D) %*% Upsilon_D
  sigma_D_rev_perf <- solve(K_D_rev_perf)
  inverse_permute <- rep(0, p)
  for (j in 1 : p)
  {
    inverse_permute[j] <- which(rev_perf == j)
  }
  K_D <- sigma_D <- matrix(0, p, p)
  K_D <- K_D_rev_perf[inverse_permute, inverse_permute]
  sigma_D <- sigma_D_rev_perf[inverse_permute, inverse_permute]
  output <- list(sigma_D, K_D)
  
  return(output)
}
