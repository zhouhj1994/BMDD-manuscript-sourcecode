library(doSNOW)
cl <- makeCluster(10, type = "SOCK") 
registerDoSNOW(cl)

if(type %in% c('Theoretical', 'Correlation', 'Theoretical2')) {
  # scImpute.folder <- paste0('scImpute/', setup, '/m_', m, '_n_', n, '_a_', a, 
  #                           '_b_', b, '_q_', q, '_r_', r, '/')
  output.folder <- paste0('output/', setup, '/m_', m, '_n_', n, '_a_', a, 
                          '_b_', b, '_q_', q, '_r_', r, '/')
} else if(type == 'Non-parametric') {
  # scImpute.folder <- paste0('scImpute/', setup, '/m_', m, '_n_', n, '/')
  output.folder <- paste0('output/', setup, '/m_', m, '_n_', n, '/')
}

# if(!dir.exists(scImpute.folder)) dir.create(scImpute.folder)
if(!dir.exists(output.folder)) dir.create(output.folder)

if(type %in% c('Theoretical', 'Theoretical2')) {
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
  if(type == 'Theoretical2') {
    fun <- function(s, alp) {
      if(s == 2) {
        Z1 <- as.integer(as.factor(Z)) - 1
        Y1 <- t(log(alp))
        X1 <- rbind(rep(1, n), Z1)
      } else if(s == 3) {
        Z1 <- Z
        Y1 <- t(log(alp))
        X1 <- rbind(rep(1, n), Z1)
      } else if(s == 4) {
        Z1 <- as.integer(as.factor(U)) - 1
        Y1 <- log(alp)
        X1 <- rbind(rep(1, m), Z1)
      } else if(s == 5) {
        Z1 <- U
        Y1 <- log(alp)
        X1 <- rbind(rep(1, m), Z1)
      }
      regr <- coef(lm(Y1 ~ Z1))
      regr[1, ] <- regr[1, ] / 2
      regr[2, ] <- regr[2, ] + abs(regr[1, ]) * (regr[2, ] > 0) -
        abs(regr[1, ]) * (regr[2, ] < 0)
      alp <- exp(t(regr) %*% X1)
      if(s %in% c(4, 5)) {
        alp <- t(alp)
      }
      return(alp)
    }
    alp0 <- fun(s, alp0)
    alp1 <- fun(s, alp1)
  }
} else if(type == 'Non-parametric') {
  rds <- readRDS('CDI_IBD_RA_SMOKE.rds')
  if (ref == 'CDI') {
    otu.tab <- rds$CDI.otu
  } else if(ref == 'IBD') {
    otu.tab <- rds$IBD.otu
  } else if(ref == 'RA') {
    otu.tab <- rds$RA.otu
  } else if(ref == 'SMOKE') {
    otu.tab <- rds$SMOKE.otu
  } 
}

output <- foreach(ii = 1 : nsim, .combine = 'rbind') %dopar% {
  library(BMDD)
  # library(mbDenoise)
  # library(mbImpute)
  # library(SAVER)
  # library(scImpute)
  # library(ALRA)
  # library(DirichletMultinomial)
  
  library(DescTools) 
  library(gdata) 
  
  source('competing_fun.R')
  source('evaluate_fun.R')
  
  set.seed(ii)
  if(type == 'Theoretical') {
    source('data_generation.R')
    data1 <- data.generation(alp0, alp1, pi, m, n, a, b, q, r, Z, U, rho = 10)
    data2 <- data.generation(alp0, alp1, pi, m, n, a, b, q, r, Z, U, rho = 100)
  } else if(type == 'Correlation') {
    source('data_generation2.R')
    data <- data.generation2(m, n, a, b, q, r, model)
  } else if(type == 'Non-parametric') {
    source('data_generation3.R')
    data <- data.generation3(otu.tab, m, n, a, b, q, r)
  } else if(type == 'Theoretical2') {
    source('data_generation.R')
    data <- data.generation(alp0, alp1, pi, m, n, a, b, q, r, Z, U)
  }
  
  W1 <- data1$W
  X1 <- data1$X
  # Z <- data$Z
  # U <- data$U
  # if(!is.null(Z)) Z <- data.frame(Z)
  # if(!is.null(U)) U <- data.frame(U)
  W2 <- data2$W
  X2 <- data2$X
  
  Xhat.list1 <- competing.fun(W1)
  res1 <- evaluate.fun(Xhat.list1, X1)
  
  Xhat.list2 <- competing.fun(W2)
  res2 <- evaluate.fun(Xhat.list2, X2)
  
  res <- cbind(res1, res2)
  
  saveRDS(list(res = res), 
          paste0(output.folder, ii, ".rds"))
}










