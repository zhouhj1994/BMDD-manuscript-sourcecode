run.fun <- function(m, n, q, q1, p, p1, r, model, scImpute.folder) {
  data <- parametric.simu(m, n, q, q1, p, p1, r, model, covariate = NULL)
  X <- data$X
  W <- data$W
 
  res.bmdd <- bmdd(W, type = 'count')
  beta <- res.bmdd$beta
  Xhat.bmdd <- t(t(beta) / colSums(beta))
  Xhat.bmdd.2 <- array(0, c(m, n))
  for(j in 1 : 20) {
    Xhat.j <- matrix(rgamma(m * n, beta, 1), m)
    Xhat.j <- t(t(Xhat.j) / colSums(Xhat.j))
    Xhat.bmdd.2 <- Xhat.bmdd.2 + Xhat.j / 20
  }
  
  family <- 'negative.binomial'
  if(model == 'Poisson') family <- 'poisson'
  res.mbDenoise <- ZIPPCApn(X = t(W), V = NULL, family = family,
                            n.factors = 2, rank = FALSE,
                            trace = FALSE, maxit = 100, parallel = FALSE)
  Xhat.mbDenoise <- t(res.mbDenoise$Q)

  res.mbImpute <- mbImpute(condition = NULL, otu_tab = t(W), metadata = NULL, 
                           D = NULL, k = 5, parallel = FALSE, ncores = 1, 
                           unnormalized = TRUE)
  What.mbImpute <- res.mbImpute$imp_count_mat_norm
  Xhat.mbImpute <- t(What.mbImpute / rowSums(What.mbImpute))
  
  err <- try(res.pmr <- autoTuneProxGradient(t(W), n_grid = 5), silent = TRUE)
  if(class(err) == 'try-error') {
    Xhat.pmr <- NA
  } else {
    Xhat.pmr <- t(res.pmr$X_hat)
  }
  
  res.SAVER <- saver(as.data.frame(W), estimates.only = TRUE)
  Xhat.SAVER <- t(t(res.SAVER) / colSums(res.SAVER))
  
  out_dir <- paste0(scImpute.folder, ii, '_', jj, '/')
  dir.create(out_dir)
  count_path <- paste0(out_dir, 'W.csv')
  write.csv(W, count_path)
  res.scImpute <- scimpute(count_path = count_path, out_dir = out_dir, 
                           Kcluster = 1, ncores = 1)
  What.scImpute <- read.csv(paste0(out_dir, 'scimpute_count.csv'))[, -1]
  Xhat.scImpute <- t(t(What.scImpute) / colSums(What.scImpute))
  
  W_norm = normalize_data(t(W))
  k_chosen = choose_k(W_norm, K = min(m, n), noise_start = min(m, n) - 20)$k 
  k <- max(2, k_chosen)
  res.ALRA <- alra(W_norm, k)
  What.ALRA <- res.ALRA$A_norm_rank_k_cor_sc
  Xhat.ALRA <- t(What.ALRA / rowSums(What.ALRA))
  
  rownames(W) <- paste0('taxon', 1 : m)
  fit <- list()
  ks <- 1 : 7
  for(k in ks){
    fit[[k]] <- dmn(t(W), k) 
  }
  lplc <- sapply(fit, laplace)
  k <- ks[which.min(lplc)]
  res.DMM <- fit[[k]]
  theta <- fitted(res.DMM)
  pi <- mixture(res.DMM)
  Xhat.DMM <- t(t(theta) / colSums(theta)) %*% t(pi)
  Xhat.DMM.2 <- array(0, c(m, n))
  for(j in 1 : 20) {
    beta <- array(NA, c(m, n))
    for(i in 1 : n) beta[, i] <- theta[, sample(1 : k, 1, prob = pi[i, ])]
    Xhat.j <- matrix(rgamma(m * n, beta, 1), m)
    Xhat.j <- t(t(Xhat.j) / colSums(Xhat.j))
    Xhat.DMM.2 <- Xhat.DMM.2 + Xhat.j / 20
  }
  
  Xhat.naive <- t(t(W) / colSums(W))
  What.naive <- W + 1
  Xhat.naive.2 <- t(t(What.naive) / colSums(What.naive))
  What.naive <- W * (W > 0) + 0.5 * (W == 0)
  Xhat.naive.3 <- t(t(What.naive) / colSums(What.naive))
  
  Xhat.list <- list(Xhat.bmdd, Xhat.bmdd.2, Xhat.mbDenoise, Xhat.mbImpute,
                    Xhat.pmr, Xhat.SAVER, Xhat.scImpute, Xhat.ALRA, 
                    Xhat.DMM, Xhat.DMM.2, Xhat.naive, Xhat.naive.2, Xhat.naive.3)
  
  res <- array(NA, c(9, 13))
  for(i in 1 : 13) {
    if(!all(is.na(Xhat.list[[i]]))){
      res[, i] <- evaluate.fun(X, Xhat.list[[i]])
    }
  }
  return(as.vector(res))
}

###############################################################
## Simulation runs
###############################################################
library(doSNOW)
cl <- makeCluster(20, type = "SOCK") 
registerDoSNOW(cl)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}

scImpute.folder <- paste0('scImpute/m_', m, '_n_', n, '_q_', q, 
                          '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, '/')
output.folder <- paste0('output/m_', m, '_n_', n, '_q_', q, 
                        '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, '/')
if(!dir.exists(scImpute.folder)) dir.create(scImpute.folder)
if(!dir.exists(output.folder)) dir.create(output.folder)

model.all <- c('gamma', 'neg-binomial', 'Poisson', 'log-normal', 'log-multinormal')
s <- 5

output.raw <- foreach (ii = 1 : nsim, .combine = 'rbind') %dopar% {
  library(BMDD)
  library(mbDenoise)
  library(mbImpute)
  source('pmr_tune_proximal_gradient.R')
  library(SAVER)
  library(scImpute)
  source("alra.R")
  library(DirichletMultinomial)
  library(foreach)
  library(Matrix)
  source('parametric_simu.R')
  source('evaluate_fun.R')
  
  result <- foreach(jj = 1 : s, .combine = 'rbind') %do% {
    unregister_dopar()
    model <- model.all[jj]
    set.seed(ii)
    run.fun(m, n, q, q1, p, p1, r, model, scImpute.folder)
  }
  write.table(result, paste0(output.folder, ii, ".txt"))
  result
}

write.table(output.raw, paste0(output.folder, 'output_raw_m_', m, '_n_', n, '_q_', q, 
                        '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt"))

output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
  res <- output.raw[seq(i, nsim * s, s), ]
  est <- colMeans(res, na.rm = TRUE)
  num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
  est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
  return(c(est, est.sd))
}

write.table(output, paste0(output.folder, 'output_m_', m, '_n_', n, '_q_', q, 
                           '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt"))

stopCluster(cl)










