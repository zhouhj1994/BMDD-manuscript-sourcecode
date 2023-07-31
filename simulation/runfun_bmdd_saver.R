run.fun <- function(m, n, q, q1, p, p1, r, model) {
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
  
  res.SAVER <- saver(as.data.frame(W), estimates.only = TRUE)
  Xhat.SAVER <- t(t(res.SAVER) / colSums(res.SAVER))
  
  Xhat.naive <- t(t(W) / colSums(W))
  What.naive <- W + 1
  Xhat.naive.2 <- t(t(What.naive) / colSums(What.naive))
  What.naive <- W * (W > 0) + 0.5 * (W == 0)
  Xhat.naive.3 <- t(t(What.naive) / colSums(What.naive))
  
  Xhat.list <- list(Xhat.bmdd, Xhat.bmdd.2, Xhat.SAVER, Xhat.naive, Xhat.naive.2, Xhat.naive.3)
  
  res <- array(NA, c(9, 6))
  for(i in 1 : 6) {
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

output.folder <- paste0('output/bmdd_saver_m_', m, '_n_', n, '_q_', q, 
                        '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, '/')
if(!dir.exists(output.folder)) dir.create(output.folder)

model.all <- c('gamma', 'neg-binomial', 'Poisson', 'log-normal', 'log-multinormal')
s <- 5

output.raw <- foreach (ii = 1 : nsim, .combine = 'rbind') %dopar% {
  library(BMDD)
  library(SAVER)
  library(foreach)
  library(Matrix)
  source('parametric_simu.R')
  source('evaluate_fun.R')
  
  result <- foreach(jj = 1 : s, .combine = 'rbind') %do% {
    unregister_dopar()
    model <- model.all[jj]
    set.seed(ii)
    run.fun(m, n, q, q1, p, p1, r, model)
  }
  write.table(result, paste0(output.folder, ii, ".txt"))
  result
}

write.table(output.raw, paste0(output.folder, 'output_raw_bmdd_saver_m_', m, '_n_', n, '_q_', q, 
                        '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt"))

output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
  res <- output.raw[seq(i, nsim * s, s), ]
  est <- colMeans(res, na.rm = TRUE)
  num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
  est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
  return(c(est, est.sd))
}

write.table(output, paste0(output.folder, 'output_bmdd_saver_m_', m, '_n_', n, '_q_', q, 
                           '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt"))

stopCluster(cl)










