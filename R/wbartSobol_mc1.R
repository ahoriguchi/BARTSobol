GetVHR <- function(tree, cutpt_list, p) {
  #' Get nodes represented as valued hyperrectangles (VHR).
  #'
  #' Helper function used in \code{\link{GetVE}}.
  #'
  #' @param tree nodes of one tree represented as a vector of strings.
  #' @param cutpt_list list of cutpoint vectors.
  #' @param p input space dimension.
  #'
  #' @return \code{list(mu, am, bm)}
  #' \item{mu}{A vector of terminal node parameter values.}
  #' \item{am}{
  #' Left endpoints of the VHR. A matrix with n_b rows and p columns.
  #' Each row corresponds to a terminal node of tree.
  #' Each columns corresponds to a dimension.
  #' The (i, j) value is the left endpoint of dimension j of the hyperrectangle corresponding to the ith terminal node.
  #' }
  #' \item{bm}{Right endpoints of the VHR. Otherwise, same as am.}

  nl <- lapply(tree, function(s) unlist(strsplit(s, split = ' ')))  # node list
  nid_v <- as.integer(sapply(nl, `[[`, 1))
  b_inds <- !((2*nid_v) %in% nid_v)
  nid_v_b <- nid_v[b_inds]
  n_b <- length(nid_v_b)

  mu <- as.numeric(sapply(nl, `[[`, 4))[b_inds]
  am <- matrix(0, nrow = n_b, ncol = p)
  bm <- matrix(1, nrow = n_b, ncol = p)

  for (i in 1:n_b) {
    nid <- nid_v_b[i]
    tmp <- nid
    while (tmp > 1) {
      is_right <- as.logical(tmp %% 2)
      tmp <- floor(tmp/2)
      ind <- which(tmp == nid_v)
      v_b <- as.integer(nl[[ind]][2]) + 1
      c_b_ind <- as.integer(nl[[ind]][3]) + 1
      c_b <- cutpt_list[[v_b]][c_b_ind]
      if (!is_right) {  # look at Eq 5 of Matt's 2016 paper.
        bm[i, v_b] <- min(bm[i, v_b], c_b)
      } else {
        am[i, v_b] <- max(am[i, v_b], c_b)
      }
    }
  }

  return(list(mu, am, bm))
}


GetVE <- function(P, mu, am, bm, to_ign = TRUE) {
  #' Get variance expression of ensemble.
  #'
  #' Ensemble is expressed as VHRs.
  #' Each row of vector mu and matrices am and bm represents a terminal node.
  #' Each terminal node has a mean parameter (\code{mu[i]})
  #' and a hyperrectangle (\code{am[i, ]}, \code{bm[i, ]}).
  #' Helper function used in \code{\link[BARTSobol]{GetVA1E}}.
  #'
  #' @param P variable index set, represented as a vector.
  #' @param mu a vector of mean parameters of terminal nodes.
  #' @param am a length(mu) x p matrix of 'left end points' of hyperrectangle.
  #' @param bm a length(mu) x p matrix of 'right end points' of hyperrectangle.
  #' @param to_ign logical to ignore nodes that is proven to contribute
  #' nothing to variance expression calculation.
  #'
  #' @return Variance of conditional expectation of ensemble.
  #'
  #' @seealso \code{\link{GetVHR}}

  # We can ignore nodes that don't split on any variable in P.
  # This will save us some computation.
  if (to_ign) {
    left_same <- (am[, P, drop = FALSE] == 0)
    right_same <- (bm[, P, drop = FALSE] == 1)
    nodes_ignore <- apply(left_same & right_same, 1, function(r) Reduce(`&`, r))

    mu <- mu[!nodes_ignore]
    am <- am[!nodes_ignore, , drop = FALSE]
    bm <- bm[!nodes_ignore, , drop = FALSE]

    if (length(mu) == 0) {
      return(0)
    }
  }

  # GetVol <- function(hr) apply(as.matrix(hr), 1, prod)
  GetVol <- function(hr) matrixStats::rowProds(as.matrix(hr))
  vhr_vol_nP <- mu * GetVol((bm-am)[, -P])
  term1 <- outer(vhr_vol_nP, vhr_vol_nP)  # Coefficients d_k^p

  hr_vol_P <- GetVol((bm-am)[, P])
  term3 <- outer(hr_vol_P, hr_vol_P)  # Product terms P_i(I_k^i)P_i(I_l^i)

  term2 <- lapply(P, function(i) outer(bm[, i], bm[, i], pmin) - outer(am[, i], am[, i], pmax))
  term2 <- Reduce("*", term2)
  term2 <- term2 * (term2 > 0)  # Intersection terms P_i(I_k^i \cap I_l^i)

  return(sum(term1 * (term2 - term3)))  # a single number

}



GetVA1E <- function(tree_list_vhr, p) {
  #' Compute Sobol' sensitivity indices for one ensemble (i.e. one posterior sample).
  #'
  #' Helper function used in \code{\link{wbartSobol}}.
  #'
  #' @param tree_list_vhr ensemble as expressed as list of VHRs.
  #' @param p input dimension
  #'
  #' @return \code{list(Si, Sij, TSi)}.
  #' \item{Si}{
  #' Main-effects Sobol' indices (normalized).
  #' A matrix with p rows and 2 columns.
  #' Each row corresponds to one of the p$variables.
  #' The first column contains the variable indices.
  #' The second column corresponds to the p main-effects Sobol' indices of the posterior draw of the BART ensemble.
  #' }
  #' \item{Sij}{
  #' Two-way Sobol' indices (normalized).
  #' A data frame with p rows and 3 columns.
  #' Each row corresponds to a pair of the p$variables.
  #' The first two columns contain the variable indices.
  #' The third column corresponds to the two-way Sobol' indices of the posterior draw of the BART ensemble.
  #' }
  #' \item{TSi}{Total-effects Sobol' indices (normalized). Otherwise, same as Si.}
  #'

  # Combine everything
  mu <- Reduce(c, lapply(tree_list_vhr, `[[`, 1))
  am <- Reduce(rbind, lapply(tree_list_vhr, `[[`, 2))
  bm <- Reduce(rbind, lapply(tree_list_vhr, `[[`, 3))

  # Compute Sobol' indices.
  v_1w <- sapply(1:p, function(i) GetVE(i, mu, am, bm))
  # v_2w <- apply(expand.grid(1:p, 1:p), 1, function(P) GetVE(P, mu, am, bm)) - apply(expand.grid(v_1w, v_1w), 1, sum)
  v_2w <- apply(data.matrix(expand.grid(1:p, 1:p)), 1, function(P) GetVE(P, mu, am, bm)) -
    matrixStats::rowSums2(data.matrix(expand.grid(v_1w, v_1w)))
  v_pw <- GetVE(1:p, mu, am, bm)  # variance of ensemble function
  tsi <- 1 - sapply(1:p, function(i) GetVE((1:p)[-i], mu, am, bm))/v_pw
  rm(mu, am, bm)

  # Cosmetics
  v_1w <- cbind(Var1 = 1:p, v_1w)
  v_2w <- cbind(expand.grid(1:p, 1:p), v_2w)
  v_2w <- v_2w[v_2w[, 1] < v_2w[, 2], ]
  tsi <- cbind(Var1 = 1:p, tsi)

  # Normalize
  v_1w[, 2] <- v_1w[, 2]/v_pw
  v_2w[, 3] <- v_2w[, 3]/v_pw

  return(list(Si = v_1w, Sij = v_2w, TSi = tsi))
}



wbartSobol <- function(wb_obj, p, n_tree, nd_post) {
  #' Get Sobol' sensitivity indices from BART.
  #'
  #' Compute main-effects, two-way, and total-effects Sobol' indices for BART::wbart() object.
  #' For each VA measure, capture all across ndpost ensembles.
  #' Prior to calling wbartSobol(), x_train from wb_obj should have been scaled to unit hypercube.
  #'
  #' @param wb_obj BART::wbart() object.
  #' @param p input dimension.
  #' @param n_tree number of trees per ensemble.
  #' @param nd_post number of posterior samples (i.e. number of ensembles).
  #'
  #' @return \code{list(Si, Sij, TSi)}.
  #' \item{Si}{
  #' Main-effects Sobol' indices (normalized).
  #' A matrix with p rows and (ndpost + 1) columns.
  #' Each row corresponds to one of the p$variables.
  #' The first column contains the variable indices.
  #' Each remaining column corresponds to the p main-effects Sobol' indices of the posterior draw of the BART ensemble.
  #' The (i, j+1) value is S_i for the jth kept draw of the ensemble.
  #' }
  #' \item{Sij}{Two-way Sobol' indices (normalized). A data frame with p rows and (ndpost + 2) columns. Each row corresponds to a pair of the p$variables. The first two columns contain the variable indices. Each remaining column corresponds to the two-way Sobol' indices of the posterior draw of the BART ensemble.}
  #' \item{TSi}{Total-effects Sobol' indices (normalized). Otherwise, same as Si.}
  #'
  #' @examples
  #' p <- 5L
  #' n_obs <- 100L
  #' n_tree <- 20L
  #' nd_post <- 100L
  #'
  #' # Create training data
  #' x_train <- sapply(1:p, function(i) rnorm(n_obs))
  #' ScaleToUnitHypercube <- function(x) (x - min(x)) / (max(x) - min(x))
  #' x_train <- apply(x_train, 2, ScaleToUnitHypercube)
  #' DataGenFun <- function(x) (x[1] - 0.5) * (x[2] - 0.5) + (x[3] - 0.5)^2
  #' y_train <- apply(x_train, 1, DataGenFun) + rnorm(n_obs)
  #'
  #' # Train BART model; call wbartSobol()
  #' wb_obj <- BART::wbart(x_train, y_train, ntree = n_tree, ndpost = nd_post)
  #' wbartSobol(wb_obj, p, n_tree, nd_post)

  to_debug <- FALSE

  wb_info <- unlist(strsplit(wb_obj$treedraws$trees, split = '\n'))
  tree_marks <- which(sapply(wb_info[-1], function(s) nchar(s)) < 7) + 1  # Hacky way of finding tree marks
  tree_list <- mapply(function(i, j) i:(i+j-1), tree_marks+1, as.integer(wb_info[tree_marks]))  # Express tree list as node indices.
  tree_list <- lapply(tree_list, function(tree) wb_info[tree])  # Express tree list as wbart() output
  tree_list <- lapply(tree_list, function(tree) GetVHR(tree, wb_obj$treedraws$cutpoints, p))  # Express tree list as VHRs
  rm(wb_info, tree_marks)
  gc()
  ens_list_vhr <- lapply(1:nd_post, function(i) tree_list[1:n_tree + (i-1)*n_tree]) # List of ndpost ensembles whose trees are expressed as VHRs
  rm(tree_list)
  gc()

  # Compute means over posterior samples
  va_list <- lapply(ens_list_vhr, GetVA1E, p = p)
  rm(ens_list_vhr)
  gc()

  Si_list <- lapply(va_list, `[[`, 1)
  Sij_list <- lapply(va_list, `[[`, 2)
  TSi_list <- lapply(va_list, `[[`, 3)
  rm(va_list)
  gc()

  # If I just want the keep all the posterior samples
  tmp <- Reduce(cbind, TSi_list)
  rm(TSi_list)
  gc()
  TSi_allpost <- tmp[, c(1, 2*(1:nd_post))]
  tmp <- Reduce(cbind, Si_list)
  rm(Si_list)
  gc()
  Si_allpost <- tmp[, c(1, 2*(1:nd_post))]
  tmp <- Reduce(cbind, Sij_list)
  rm(Sij_list)
  gc()
  Sij_allpost <- tmp[, c(1, 2, 3*(1:nd_post))]

  my_vals_l <- list(Si = Si_allpost, Sij = Sij_allpost, TSi = TSi_allpost)

  return(my_vals_l)

}



#' wbartSobolText <- function(wb_obj_text, wb_obj_cutpoints, p, n_tree, nd_post) {
#'   #' Get variable activity measures (Si, Sij, cnt, TSi) from BART.
#'   #' For each VA measure, take the sample mean across ndpost ensembles.
#'   #' Return the four sample means.
#'   #'
#'   #' @param wb_obj_text wb_obj$treedraws$trees
#'   #' @param wb_obj_cutpoints wb_obj$treedraws$cutpoints
#'   #' @param p input dimension.
#'   #' @param n_tree number of trees per ensemble.
#'   #' @param nd_post number of posterior samples (i.e. number of ensembles).
#'   #'
#'   #' @return List of Si, Sij, TSi.
#'
#'   to_debug <- FALSE
#'
#'   wb_info <- unlist(strsplit(wb_obj_text, split = '\n'))
#'   tree_marks <- which(sapply(wb_info[-1], function(s) nchar(s)) < 7) + 1  # Hacky way of finding tree marks
#'   tree_list <- mapply(function(i, j) i:(i+j-1), tree_marks+1, as.integer(wb_info[tree_marks]))  # Express tree list as node indices.
#'   tree_list <- lapply(tree_list, function(tree) wb_info[tree])  # Express tree list as wbart() output
#'   tree_list <- lapply(tree_list, function(tree) GetVHR(tree, wb_obj_cutpoints, p))  # Express tree list as VHRs
#'   rm(wb_info, tree_marks)
#'   gc()
#'   ens_list_vhr <- lapply(1:nd_post, function(i) tree_list[1:n_tree + (i-1)*n_tree]) # List of ndpost ensembles whose trees are expressed as VHRs
#'   rm(tree_list)
#'   gc()
#'   if (to_debug) {
#'     save(ens_list_vhr, file = "ens_list_vhr.RData")
#'   }
#'
#'
#'
#'   # Compute means over posterior samples
#'   # va_list <- lapply(tree_list_vhr, GetVA1E, p = p)
#'   # va_list <- parallel::mclapply(tree_list_vhr, GetVA1E, p = p, mc.cores = parallel::detectCores())
#'   # bs <- floor(nd_post/5000)
#'   # bs_remain <- nd_post %% (5000 * bs)
#'   bs <- ceiling(nd_post/5000)
#'   va_list_list <- lapply(1:bs, function(i) {
#'     a <- min((i-1)*5000 + 1, nd_post)
#'     b <- min((i-1)*5000 + 5000, nd_post)
#'     parallel::mclapply(ens_list_vhr[a:b], GetVA1E, p = p, mc.cores = 32L)
#'   })
#'   va_list <- Reduce(c, va_list_list)
#'   # va_list <- parallel::mclapply(ens_list_vhr, GetVA1E, p = p, mc.cores = 32L)  # List of ndpost lists of VA measures.
#'   rm(ens_list_vhr)
#'   gc()
#'
#'   if (to_debug) {
#'     save(va_list, file = "va_list.RData")
#'   }
#'
#'   Si_list <- lapply(va_list, `[[`, 1)
#'   Sij_list <- lapply(va_list, `[[`, 2)
#'   TSi_list <- lapply(va_list, `[[`, 3)
#'   rm(va_list)
#'   gc()
#'
#'   # # If I just want the point estimates over the posterior samples
#'   # Si_shell <- Si_list[[1]]
#'   # tmp <- Reduce(cbind, Si_list)
#'   # Si_shell[, 2] <- matrixStats::rowMeans2(tmp[, 2*(1:nd_post)])
#'   # Sij_list <- lapply(va_list, `[[`, 2)
#'   # Sij_shell <- Sij_list[[1]]
#'   # Sij_shell[, 3] <- apply(sapply(Sij_list, function(x) x[, 3]), 1, mean)
#'   # TSi_list <- lapply(va_list, `[[`, 3)
#'   # TSi_shell <- TSi_list[[1]]
#'   # TSi_shell[, 2] <- apply(sapply(TSi_list, function(x) x[, 2]), 1, mean)
#'
#'   # If I just want the keep all the posterior samples
#'   tmp <- Reduce(cbind, TSi_list)
#'   print(paste("dim(Reduce(cbind, TSi_list)):", dim(tmp)))
#'   TSi_allpost <- tmp[, c(1, 2*(1:nd_post))]
#'   tmp <- Reduce(cbind, Si_list)
#'   print(paste("dim(Reduce(cbind, Si_list)):", dim(tmp)))
#'   Si_allpost <- tmp[, c(1, 2*(1:nd_post))]
#'   tmp <- Reduce(cbind, Sij_list)
#'   print(paste("dim(Reduce(cbind, Sij_list)):", dim(tmp)))
#'   Sij_allpost <- tmp[, c(1, 2, 3*(1:nd_post))]
#'   rm(Si_list, Sij_list, TSi_list)
#'   gc()
#'
#'   my_vals_l <- list(Si = Si_allpost, Sij = Sij_allpost, TSi = TSi_allpost)
#'
#'   return(my_vals_l)
#'
#' }








