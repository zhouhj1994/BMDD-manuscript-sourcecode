zero.fun <- function(X) {
  if(any(X == 0)) {
    X <- t(apply(X, 1, function (x) {
      if(all(x == 0)) {
        x[x == 0] <- min(X[X != 0])
      } else {
        x[x == 0] <- min(x[x != 0]) 
      }
      return(x)
    }))
  }
  return(X)
}

run.fun <- function(m, n, gamma, beta0, sigma2, eta0, theta0, kappa0, lambda0, model, 
                    formula, alpha, K) {
  
  data <- simulate.data(m, n, gamma, beta0, sigma2, eta0, theta0, kappa0, lambda0, model)
  Y <- data$Y
  Z <- data$Z
  H <- data$H
  
  ## LinDA
  res <- linda(Y, Z, formula)
  pval <- res$output[[1]]$pvalue
  qval <- p.adjust(pval, method = 'BH')
  rej.linda <- which(qval <= alpha)

  ## BMDD
  res.bmdd <- bmdd(Y, type = 'count', trace = TRUE)
  beta <- res.bmdd$beta
  beta <- beta[, rep(1 : n, K)]
  X <- matrix(rgamma(m * n * K, beta, 1), m)
  X <- t(t(X) / colSums(X))
  X <- zero.fun(X)
  X <- t(t(X) / colSums(X))
  
  dens <- lgamma(colSums(beta)) - colSums(lgamma(beta)) + colSums((beta - 1) * log(X))
  dens <- matrix(dens, nrow = n)
  dens.ord <- t(apply(dens, 1, order, decreasing = TRUE))
  ind <- (dens.ord - 1) * n + 1 : n
  X.bmdd <- X[, ind]
  
  ## SAVER
  res.saver <- saver(as.data.frame(Y))
  X <- do.call(cbind, sample.saver(res.saver, rep = K))
  X[which(is.nan(X))] <- 0
  X <- t(t(X) / colSums(X))
  X <- zero.fun(X)
  X <- t(t(X) / colSums(X))
  
  ave <- res.saver$estimate
  se <- res.saver$se
  a <- ave ^ 2 / se ^ 2
  b <- ave / se ^ 2
  a <- a[, rep(1 : n, K)]
  b <- b[, rep(1 : n, K)]
  dens <- colSums(a * log(b) - lgamma(a) + (a - 1) * log(X) - b * X)
  dens <- matrix(dens, nrow = n)
  dens.ord <- t(apply(dens, 1, order, decreasing = TRUE))
  ind <- (dens.ord - 1) * n + 1 : n
  X.saver <- X[, ind]
  
  ## LMM
  lmm.fun <- function(X) {
    l <- ncol(X) / n
    id2 <- factor(rep(1 : n, l))
    Z1 <- as.data.frame(Z[rep(1 : n, l), ])
    colnames(Z1) <- colnames(Z)
    Z2 <- cbind(Z1, id2)
    res <- linda(X, Z2, formula = paste0(formula, "+(1|id2)"), type = "proportion")
    pval <- res$output[[1]]$pvalue
    qval <- p.adjust(pval, method = 'BH')
    rej <- which(qval <= alpha)
    return(rej)
  }
  
  ##
  l <- floor(K * c(0.2, 0.5, 1))
  rej.lmm.bmdd <- list()
  rej.lmm.saver <- list()
  
  for(i in 1 : length(l)) {
    rej.lmm.bmdd[[i]] <- lmm.fun(X.bmdd[, 1 : (n * l[i])])
    rej.lmm.saver[[i]] <- lmm.fun(X.saver[, 1 : (n * l[i])])
  }
  
  rej.list <- c(list(rej.linda), rej.lmm.bmdd, rej.lmm.saver)
  
  n.met <- length(rej.list)
  fdp <- rep(NA, n.met)
  power <- rep(NA, n.met)
  
  if (sum(H) == 0) {
    for (i in 1 : n.met) {
      rej <- rej.list[[i]]
      if(length(rej) == 0) {
        fdp[i] <- 0
      } else {
        fdp[i] <- 1
      }
    }
  } else {
    for (i in 1 : n.met) {
      rej <- rej.list[[i]]
      if(length(rej) == 0) {
        fdp[i] <- 0
        power[i] <- 0
      } else {
        fdp[i] <- sum(H[rej] == 0) / length(rej)
        power[i] <- sum(H[rej] == 1) / sum(H)
      }
    }
  }
  return(c(fdp, power))
}

###############################################################
## Simulation setups
###############################################################
source("simulate_data.R")
para0 <- readRDS("log.normal.para.rds")
para1 <- readRDS("Gamma.para.rds")
para2 <- readRDS("negBinom.para.rds")
para3 <- readRDS("Poisson.para.rds")

beta0 <- para0$beta0
sigma2 <- para0$sigma2 
eta0 <- para1$eta0
theta0 <- para2$theta0
kappa0 <- para2$kappa0
lambda0 <- para3$lambda0

gamma <- 0.1
alpha <- 0.05
K <- 100
formula <- '~u'
nsim <- 500

feature.size.vec <- c(50, 200, 500)
sample.size.vec <- c(50, 200)
s1 <- length(feature.size.vec); s2 <- length(sample.size.vec)

feature.size <- rep(feature.size.vec, each = s2)
sample.size <- rep(sample.size.vec, s1)
setting <- cbind(feature.size, sample.size)

s <- s1 * s2

###############################################################
## Simulation runs
###############################################################
library(doSNOW)
cl <- makeCluster(20, type = "SOCK") 
registerDoSNOW(cl)

output.raw <- foreach (i = 1 : nsim,.combine = 'rbind') %dopar% {
  library(LinDA)
  library(SAVER)
  library(BMDD)
  library(foreach)
  
  result <- foreach(j = 1 : s, .combine = 'rbind') %do% {
    para <- setting[j, ]
    m <- as.numeric(para[1])
    n <- as.numeric(para[2])
    set.seed((i - 1) * s + j)
    run.fun(m, n, gamma, beta0 = beta0, sigma2 = sigma2, eta0 = eta0, 
            theta0 = theta0, kappa0 = kappa0, lambda0 = lambda0,
            model = model, formula = formula, alpha = alpha, K = K) 
  }
  write.table(result, paste0("output/output_raw_", model, "_", i, ".txt"))
  result
}

write.table(output.raw, paste0("output/output_raw_", model, ".txt"))

n.met <- ncol(output.raw) / 2
output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
  res <- output.raw[seq(i, nsim * s, s), ]
  est <- colMeans(res, na.rm = TRUE)
  num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
  est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
  
  fdr.est <- est[1 : n.met]
  fdr.est.sd <- est.sd[1 : n.met]
  power.est <- est[(n.met + 1) : (2 * n.met)]
  power.est.sd <- est.sd[(n.met + 1) : (2 * n.met)]
  
  res <- c(fdr.est, power.est, fdr.est.sd, power.est.sd)
  res[which(is.nan(res))] <- NA
  res
}

write.table(output, paste0("output/output_", model, ".txt"))

stopCluster(cl)



