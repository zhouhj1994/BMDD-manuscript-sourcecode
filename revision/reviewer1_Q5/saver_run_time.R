saver.run.time <- function(alp0, alp1, pi, m, n, a, b, q, r, Z, U, iterlim) {
  
  data <- data.generation(alp0, alp1, pi, m, n, a, b, q, r, Z, U)
  W <- data$W
  
  t <- system.time(res.saver <- saver(as.data.frame(W)))
  
  return(as.vector(t[1]))
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

output <- foreach(ii = 21 : (20 + nsim), .combine = 'rbind') %dopar% {
  
  library(SAVER)
  source('data_generation.R')
  #source('bmdd2.R')
  library(foreach)
  sink(paste0('output/saver_run_time_3/time', s, '_', ii, '.txt'))
  t <- foreach(jj = 1 : l) %do% {
    para <- setting[jj, ]
    m <- para[1]
    n <- para[2]
    a <- 1000 * (m / 100) ^ (1 / 1.5)
    set.seed(ii)
    res <- saver.run.time(alp0, alp1, pi, m, n, a, b, q, r, Z, U, iterlim)
    print(jj)
    res
  }
  sink()
  unlist(t)
}

write.table(output, paste0('output/saver_run_time_3/time', s, '.txt'))


