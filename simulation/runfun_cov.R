run.fun <- function(m, n, q, q1, p, p1, r, model, covariate) {
  data <- parametric.simu(m, n, q, q1, p, p1, r, model, covariate)
  X <- data$X
  W <- data$W
  Z <- data$Z
  Z <- as.data.frame(Z)
 
  res.bmdd <- bmdd(W, type = 'count')
  beta <- res.bmdd$beta
  Xhat.bmdd <- t(t(beta) / colSums(beta))
  Xhat.bmdd.2 <- array(0, c(m, n))
  for(j in 1 : 20) {
    Xhat.j <- matrix(rgamma(m * n, beta, 1), m)
    Xhat.j <- t(t(Xhat.j) / colSums(Xhat.j))
    Xhat.bmdd.2 <- Xhat.bmdd.2 + Xhat.j / 20
  }
  
  res.bmdd <- bmdd(W, type = 'count', Z, formula.Z = '~Z', alp.eta = TRUE, pi.xi = TRUE)
  beta <- res.bmdd$beta
  Xhat.bmdd.3 <- t(t(beta) / colSums(beta))
  Xhat.bmdd.4 <- array(0, c(m, n))
  for(j in 1 : 20) {
    Xhat.j <- matrix(rgamma(m * n, beta, 1), m)
    Xhat.j <- t(t(Xhat.j) / colSums(Xhat.j))
    Xhat.bmdd.4 <- Xhat.bmdd.4 + Xhat.j / 20
  }
  
  family <- 'negative.binomial'
  if(model == 'Poisson') family <- 'poisson'
  res.mbDenoise <- ZIPPCApn(X = t(W), V = NULL, family = family,
                            n.factors = 2, rank = FALSE,
                            trace = FALSE, maxit = 100, parallel = FALSE)
  Xhat.mbDenoise <- t(res.mbDenoise$Q)
  
  res.mbDenoise <- ZIPPCApn(X = t(W), V = Z$Z, family = family,
                            n.factors = 2, rank = FALSE,
                            trace = FALSE, maxit = 100, parallel = FALSE)
  Xhat.mbDenoise.2 <- t(res.mbDenoise$Q)

  res.mbImpute <- mbImpute(condition = NULL, otu_tab = t(W), metadata = NULL, 
                           D = NULL, k = 5, parallel = FALSE, ncores = 1, 
                           unnormalized = TRUE)
  What.mbImpute <- res.mbImpute$imp_count_mat_norm
  Xhat.mbImpute <- t(What.mbImpute / rowSums(What.mbImpute))
  
  res.mbImpute <- mbImpute(condition = NULL, otu_tab = t(W), metadata = Z, 
                           D = NULL, k = 5, parallel = FALSE, ncores = 1, 
                           unnormalized = TRUE)
  What.mbImpute <- res.mbImpute$imp_count_mat_norm
  Xhat.mbImpute.2 <- t(What.mbImpute / rowSums(What.mbImpute))
  
  Xhat.list <- list(Xhat.bmdd, Xhat.bmdd.2, Xhat.bmdd.3, Xhat.bmdd.4,
                    Xhat.mbDenoise, Xhat.mbDenoise.2, Xhat.mbImpute, Xhat.mbImpute.2)
  
  res <- array(NA, c(9, 8))
  for(i in 1 : 8) {
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

output.folder <- paste0('output/cov_m_', m, '_n_', n, '_q_', q, 
                        '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, '/')
if(!dir.exists(output.folder)) dir.create(output.folder)

model.all <- c('gamma', 'neg-binomial', 'Poisson', 'log-normal', 'log-multinormal')
covariate.all <- c('Bernoulli', 'normal')
setting <- cbind(rep(model.all, each = 2), rep(covariate.all, 5))
s <- 10

output.raw <- foreach (ii = 1 : nsim, .combine = 'rbind') %dopar% {
  library(BMDD)
  library(mbDenoise)
  library(mbImpute)
  library(foreach)
  library(Matrix)
  source('parametric_simu.R')
  source('evaluate_fun.R')
  
  result <- foreach(jj = 1 : s, .combine = 'rbind') %do% {
    unregister_dopar()
    model <- setting[jj, 1]
    covariate <- setting[jj, 2]
    set.seed(ii)
    run.fun(m, n, q, q1, p, p1, r, model, covariate)
  }
  write.table(result, paste0(output.folder, ii, ".txt"))
  result
}

write.table(output.raw, paste0(output.folder, 'output_raw_cov_m_', m, '_n_', n, '_q_', q, 
                        '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt"))

output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
  res <- output.raw[seq(i, nsim * s, s), ]
  est <- colMeans(res, na.rm = TRUE)
  num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
  est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
  return(c(est, est.sd))
}

write.table(output, paste0(output.folder, 'output_cov_m_', m, '_n_', n, '_q_', q, 
                           '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt"))

stopCluster(cl)










