competing.fun <- function(W) {

  # res.BMDD <- bmdd(W, type = 'count', trace = TRUE)
  # beta <- res.BMDD$beta
  # Xhat.BMDD <- t(t(beta) / colSums(beta))
  # 
  # if(!is.null(Z) & is.null(U)) {
  #   res.BMDD <- bmdd(W, type = 'count', Z = Z, alp.eta = TRUE, pi.xi = TRUE, trace = TRUE)
  # } else if(is.null(Z) & !is.null(U)) {
  #   res.BMDD <- bmdd(W, type = 'count', U = U, alp.kap = TRUE, pi.zeta = TRUE, trace = TRUE)
  # }
  # if(!is.null(Z) | !is.null(U)) {
  #   beta <- res.BMDD$beta
  #   Xhat.BMDD.2 <- t(t(beta) / colSums(beta))
  # }
  # 
  # res.mbDenoise <- ZIPPCApn(X = t(W), V = NULL, family = "negative.binomial",
  #                           n.factors = 2, rank = FALSE,
  #                           trace = FALSE, maxit = 100, parallel = FALSE)
  # Xhat.mbDenoise <- t(res.mbDenoise$Q)
  # 
  # if(!is.null(Z)) {
  #   res.mbDenoise <- ZIPPCApn(X = t(W), V = Z[, 1], family = "negative.binomial",
  #                             n.factors = 2, rank = FALSE,
  #                             trace = FALSE, maxit = 100, parallel = FALSE)
  #   Xhat.mbDenoise.2 <- t(res.mbDenoise$Q)
  # }
  # 
  # res.mbImpute <- mbImpute(condition = NULL, otu_tab = t(W), metadata = NULL, 
  #                          D = NULL, k = 5, parallel = FALSE, ncores = 1, 
  #                          unnormalized = TRUE)
  # What.mbImpute <- res.mbImpute$imp_count_mat_norm
  # Xhat.mbImpute <- t(What.mbImpute / rowSums(What.mbImpute))
  # 
  # if(!is.null(Z)) {
  #   res.mbImpute <- mbImpute(condition = NULL, otu_tab = t(W), metadata = Z, 
  #                            D = NULL, k = 5, parallel = FALSE, ncores = 1, 
  #                            unnormalized = TRUE)
  #   What.mbImpute <- res.mbImpute$imp_count_mat_norm
  #   Xhat.mbImpute.2 <- t(What.mbImpute / rowSums(What.mbImpute))
  # }
  # 
  # res.SAVER <- saver(as.data.frame(W), estimates.only = TRUE)
  # Xhat.SAVER <- t(t(res.SAVER) / colSums(res.SAVER))
  # 
  # out_dir <- paste0(scImpute.folder, ii, '/')
  # dir.create(out_dir)
  # count_path <- paste0(out_dir, 'W.csv')
  # write.csv(W, count_path)
  # res.scImpute <- scimpute(count_path = count_path, out_dir = out_dir, 
  #                          Kcluster = 1, ncores = 1)
  # What.scImpute <- read.csv(paste0(out_dir, 'scimpute_count.csv'))[, -1]
  # Xhat.scImpute <- t(t(What.scImpute) / colSums(What.scImpute))
  # 
  # m <- nrow(W); n <- ncol(W)
  # W_norm = normalize_data(t(W))
  # k_chosen = choose_k(W_norm, K = min(m, n), noise_start = min(m, n) - 20)$k 
  # k <- max(2, k_chosen)
  # res.ALRA <- alra(W_norm, k)
  # What.ALRA <- res.ALRA$A_norm_rank_k_cor_sc
  # Xhat.ALRA <- t(What.ALRA / rowSums(What.ALRA))
  
  m <- nrow(W); n <- ncol(W)
  rownames(W) <- paste0('taxon', 1 : m)
  Xhat.DMM <- array(NaN, c(m, n))
  while(is.nan(Xhat.DMM[1, 1])) {
    res.DMM <- dmn(t(W), 1) 
    theta <- fitted(res.DMM)
    pi <- mixture(res.DMM)
    Xhat.DMM <- t(t(theta) / colSums(theta)) %*% t(pi)
  }
  
  # Xhat.naive.1 <- t(t(W) / colSums(W))
  # What.naive <- W + 1
  # Xhat.naive.2 <- t(t(What.naive) / colSums(What.naive))
  # What.naive <- W * (W > 0) + 0.5 * (W == 0)
  # Xhat.naive.3 <- t(t(What.naive) / colSums(What.naive))
  # 
  # N <- colSums(W)
  # N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
  # N.mat[W > 0] <- 0
  # tmp <- N[max.col(N.mat)]
  # What.naive <- W + N.mat / tmp
  # Xhat.naive.4 <- t(t(What.naive) / colSums(What.naive))
  # 
  # if(is.null(Z) & is.null(U)) {
  #   Xhat.list <- list(Xhat.BMDD = Xhat.BMDD, Xhat.mbDenoise = Xhat.mbDenoise, 
  #                     Xhat.mbImpute = Xhat.mbImpute, Xhat.SAVER = Xhat.SAVER, 
  #                     Xhat.scImpute = Xhat.scImpute, Xhat.ALRA = Xhat.ALRA, Xhat.DMM = Xhat.DMM,
  #                     Xhat.naive.1 = Xhat.naive.1, Xhat.naive.2 = Xhat.naive.2, 
  #                     Xhat.naive.3 = Xhat.naive.3, Xhat.naive.4 = Xhat.naive.4)
  # } else if(!is.null(Z) & is.null(U)) {
  #   Xhat.list <- list(Xhat.BMDD = Xhat.BMDD, Xhat.mbDenoise = Xhat.mbDenoise, 
  #                     Xhat.mbImpute = Xhat.mbImpute, Xhat.SAVER = Xhat.SAVER, 
  #                     Xhat.scImpute = Xhat.scImpute, Xhat.ALRA = Xhat.ALRA, Xhat.DMM = Xhat.DMM,
  #                     Xhat.naive.1 = Xhat.naive.1, Xhat.naive.2 = Xhat.naive.2, 
  #                     Xhat.naive.3 = Xhat.naive.3, Xhat.naive.4 = Xhat.naive.4, 
  #                     Xhat.BMDD.2 = Xhat.BMDD.2,
  #                     Xhat.mbDenoise.2 = Xhat.mbDenoise.2, Xhat.mbImpute.2 = Xhat.mbImpute.2)
  # } else if(is.null(Z) & !is.null(U)) {
  #   Xhat.list <- list(Xhat.BMDD = Xhat.BMDD, Xhat.mbDenoise = Xhat.mbDenoise, 
  #                     Xhat.mbImpute = Xhat.mbImpute, Xhat.SAVER = Xhat.SAVER, 
  #                     Xhat.scImpute = Xhat.scImpute, Xhat.ALRA = Xhat.ALRA, Xhat.DMM = Xhat.DMM,
  #                     Xhat.naive.1 = Xhat.naive.1, Xhat.naive.2 = Xhat.naive.2, 
  #                     Xhat.naive.3 = Xhat.naive.3, Xhat.naive.4 = Xhat.naive.4, 
  #                     Xhat.BMDD.2 = Xhat.BMDD.2)
  # }
  
  Xhat.list <- list(Xhat.DMM = Xhat.DMM)
  return(Xhat.list)
}



