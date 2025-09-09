bmdd.run.time <- function(alp0, alp1, pi, m, n, a, b, q, r, Z, U, iterlim) {
  
  data <- data.generation(alp0, alp1, pi, m, n, a, b, q, r, Z, U)
  W <- data$W
  X <- data$X
  Z <- data$Z
  U <- data$U
  if(!is.null(Z)) Z <- data.frame(Z)
  if(!is.null(U)) U <- data.frame(U)
  
  if(is.null(Z) & is.null(U)) {
    t <- system.time(res.BMDD <- bmdd2(W, type = 'count', iterlim = iterlim))
  } else if(!is.null(Z) & is.null(U)) {
    t <- system.time(res.BMDD <- bmdd2(W, type = 'count', Z = Z, alp.eta = TRUE, pi.xi = TRUE, iterlim = iterlim))
  } else if(is.null(Z) & !is.null(U)) {
    t <- system.time(res.BMDD <- bmdd2(W, type = 'count', U = U, alp.kap = TRUE, pi.zeta = TRUE, iterlim = iterlim))
  }
  
  iter <- res.BMDD$iter
  if(iter == iterlim + 1) iter <- iterlim
  t <- as.vector(t[1] / iter)
  return(t)
}

bmdd_fit_on_combo <- readRDS("bmdd_fit_on_combo.rds")
fit <- bmdd_fit_on_combo$res[[s]]
pi <- fit$pi
alp0 <- fit$alpha$alp0
alp1 <- fit$alpha$alp1
Z <- NULL
U <- NULL
if(s == 2) {
  Z <- bmdd_fit_on_combo$cova$Z1
} else if(s == 3) {
  Z <- bmdd_fit_on_combo$cova$Z2
} else if(s == 4) {
  U <- bmdd_fit_on_combo$cova$U1
} else if(s == 5) {
  U <- bmdd_fit_on_combo$cova$U2
}

lm <- length(m.vec)
ln <- length(n.vec)
setting <- cbind(rep(m.vec, each = ln), rep(n.vec, lm))
l <- nrow(setting)
nsim <- 10

library(doSNOW)
cl <- makeCluster(10, type = "SOCK") 
registerDoSNOW(cl)

output <- foreach(ii = 1 : nsim, .combine = 'rbind') %dopar% {
  
  source('data_generation.R')
  source('bmdd2.R')
  library(foreach)
  sink(paste0('output/bmdd_run_time/time', s, '_', ii, '.txt'))
  t <- foreach(jj = 1 : l) %do% {
    para <- setting[jj, ]
    m <- para[1]
    n <- para[2]
    a <- 1000 * (m / 100) ^ (1 / 1.5)
    set.seed(ii)
    res <- bmdd.run.time(alp0, alp1, pi, m, n, a, b, q, r, Z, U, iterlim)
    print(jj)
    res
  }
  sink()
  unlist(t)
}

write.table(output, paste0('output/bmdd_run_time/time', s, '.txt'))


