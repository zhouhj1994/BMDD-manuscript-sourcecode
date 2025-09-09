winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

preprocess.fun <- function(otu.tab, meta, prev.cut = 0.2, lib.cut = 1000, 
                           winsor.quan = 0.97) {
  keep.sam <- which(colSums(otu.tab) >= lib.cut)
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(meta[keep.sam, ])
  names(Z) <- names(meta)
  rownames(Z) <- rownames(meta)[keep.sam]
  keep.tax <- which(rowSums(Y > 0) / ncol(Y) >= prev.cut)
  Y <- Y[keep.tax, ]
  
  N <- colSums(Y)
  keep.sam1 <- which(colSums(Y) >= 1)
  Y1 <- Y[, keep.sam1]
  Z1 <- as.data.frame(Z[keep.sam1, ])
  names(Z1) <- names(Z)
  rownames(Z1) <- rownames(Z)[keep.sam1]
  
  Y1 <- winsor.fun(Y1, winsor.quan) 
  
  return(list(Y = Y1, Z = Z1, keep.sam = keep.sam[keep.sam1], keep.tax = keep.tax))
}

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

run.fun <- function(Y, Z, formula, K) {
  
  m <- nrow(Y)
  n <- ncol(Y)
  
  ## LinDA
  res <- linda(Y, Z, paste0('~', formula))
  pval.linda <- res$output[[1]]$pvalue
  print('-----LinDA-----')
  
  ## ANCOM-BC
  OTU = otu_table(Y, taxa_are_rows = TRUE)
  META = sample_data(Z)
  PHYSEQ = phyloseq(OTU, META)
  res <- try(ancombc2(data = PHYSEQ, fix_formula = formula))
  if(inherits(res, "try-error")) {
    pval.ancombc <- rep(NA, m)
  } else {
    pval.ancombc <- res$res[, 9]
  }
  print('-----ANCOMBC-----')
  
  ## BMDD
  res.bmdd <- bmdd(Y, type = 'count', trace = TRUE)
  beta <- res.bmdd$beta
  X <- t(t(beta) / colSums(beta))
  print('-----BMDD-----')
  
  ## LinDA-BMDD-Mean
  res <- linda(X, Z, paste0('~', formula), type = "proportion")
  pval.linda.bmdd.mean <- res$output[[1]]$pvalue
  print('-----LinDA-BMDD-Mean-----')
  
  ## ANCOMBC-BMDD-Mean
  OTU = otu_table(X, taxa_are_rows = TRUE)
  PHYSEQ = phyloseq(OTU, META)
  res <- try(ancombc2(data = PHYSEQ, fix_formula = formula))
  if(inherits(res, "try-error")) {
    pval.ancombc.bmdd.mean <- rep(NA, m)
  } else {
    pval.ancombc.bmdd.mean <- res$res[, 9]
  }
  print('-----ANCOMBC-BMDD-Mean-----')
  
  ## BMDD posterior samples
  beta <- beta[, rep(1 : n, K)]
  X <- matrix(rgamma(m * n * K, beta, 1), m)
  X <- t(t(X) / colSums(X))
  X <- zero.fun(X)
  X.bmdd <- t(t(X) / colSums(X))
  colnames(X.bmdd) <- paste0('sample', 1 : (n * K))
  rownames(X.bmdd) <- rownames(beta)
  
  ## SAVER
  res.saver <- saver(as.data.frame(Y))
  X <- t(t(res.saver$estimate) / colSums(res.saver$estimate))
  print('-----SAVER-----')
  
  ## LinDA-SAVER-Mean
  res <- linda(X, Z, paste0('~', formula), type = "proportion")
  pval.linda.saver.mean <- res$output[[1]]$pvalue
  print('-----LinDA-SAVER-Mean-----')
  
  ## ANCOMBC-SAVER-Mean
  OTU = otu_table(X, taxa_are_rows = TRUE)
  PHYSEQ = phyloseq(OTU, META)
  res <- try(ancombc2(data = PHYSEQ, fix_formula = formula))
  if(inherits(res, "try-error")) {
    pval.ancombc.saver.mean <- rep(NA, m)
  } else {
    pval.ancombc.saver.mean <- res$res[, 9]
  }
  print('-----ANCOMBC-SAVER-Mean-----')

  ## SAVER posterior samples
  X <- do.call(cbind, sample.saver(res.saver, rep = K))
  X[which(is.nan(X))] <- 0
  X <- t(t(X) / colSums(X))
  X <- zero.fun(X)
  X.saver <- t(t(X) / colSums(X))
  colnames(X.saver) <- paste0('sample', 1 : (n * K))
  
  ## LMM
  lmm.fun <- function(X) {
    l <- ncol(X) / n
    id2 <- factor(rep(1 : n, l))
    Z1 <- as.data.frame(Z[rep(1 : n, l), ])
    colnames(Z1) <- colnames(Z)
    Z2 <- cbind(Z1, id2)
    rownames(Z2) <- colnames(X)
    
    res <- linda(X, Z2, paste0('~', formula, "+(1|id2)"), type = "proportion")
    pval.linda <- res$output[[1]]$pvalue
    print('-----LinDA-LMM-----')
    
    OTU = otu_table(X[, 1 : (n * 20)], taxa_are_rows = TRUE)
    META = sample_data(Z2[1 : (n * 20), ])
    PHYSEQ = phyloseq(OTU, META)
    res <- try(ancombc2(data = PHYSEQ, fix_formula = formula, rand_formula = "(1|id2)"))
    if(inherits(res, "try-error")) {
      pval.ancombc <- rep(NA, m)
    } else {
      pval.ancombc <- res$res[, 9]
    }
    print('-----ANCOMBC-LMM-----')
    
    return(list(pval.linda, pval.ancombc))
  }
  
  ####
  res1 <- lmm.fun(X.bmdd)
  res2 <- lmm.fun(X.saver)
  
  pval.linda.bmdd <- res1[[1]]
  pval.ancombc.bmdd <- res1[[2]]
  pval.linda.saver <- res2[[1]]
  pval.ancombc.saver <- res2[[2]]
  
  pval.mat <- cbind(pval.linda, pval.ancombc, 
                    pval.linda.bmdd, pval.ancombc.bmdd,
                    pval.linda.saver, pval.ancombc.saver,
                    pval.linda.bmdd.mean, pval.ancombc.bmdd.mean, 
                    pval.linda.saver.mean, pval.ancombc.saver.mean)
  return(pval.mat)
}

###############################################################
## Run
###############################################################
filenames <- list.files("IBD")
K <- 100

data.list <- list()
for(i in 1 : 4) {
  load(paste0("IBD/", filenames[i]))
  res <- preprocess.fun(otu_table, metadata)
  data.list[[2 * (i - 1) + 1]] <- res$Y
  data.list[[2 * i]] <- res$Z$grp
}

library(doSNOW)
cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)

foreach(j = nshuf1 : nshuf2) %dopar% {
  library(LinDA)
  library(BMDD)
  library(SAVER)
  library(ANCOMBC)
  library(phyloseq)
  
  pval.mat.list <- list()
  sink(paste0("log/logfile", j, ".txt"))
  for(i in 1 : 4) {
    Y <- data.list[[2 * (i - 1) + 1]]
    Z <- data.list[[2 * i]]
    n <- length(Z)
    
    set.seed(j)
    Z <- Z[sample(1 : n, n)]
    Z <- as.data.frame(Z)
    names(Z) <- 'grp'
    rownames(Z) <- colnames(Y)
    pmat <- run.fun(Y, Z, 'grp', K)
    
    colnames(pmat) <- c('LinDA', 'ABCOMBC', 
                        'LinDA-BMDD', 'ANCOMBC-BMDD', 'LinDA-SAVER', 'ANCOMBC-SAVER',
                        'LinDA-BMDD-Mean', 'ANCOMBC-BMDD-Mean', 
                        'LinDA-SAVER-Mean', 'ANCOMBC-SAVER-Mean')
    rownames(pmat) <- rownames(Y)
    pval.mat.list[[i]] <- pmat
  }
  sink()
  saveRDS(pval.mat.list, paste0('output/pval.mat.list.', j, '.rds'))
}
stopCluster(cl)








