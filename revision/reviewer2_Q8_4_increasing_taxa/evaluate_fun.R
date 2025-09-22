evaluate.fun <- function(Xhat.list, X) {
  
  fun <- function(Xhat, X) {
    tmp <- (colSums(X ^ 2) - colSums(Xhat ^ 2)) ^ 2
    Msp <- mean(tmp)
    Msp.r <- median(tmp)
    
    X1 <- pmin(X, Xhat)
    tmp <- 1 - colSums(X1)
    Mbc <- mean(tmp)
    Mbc.r <- median(tmp)
    
    tmp <- sqrt(colSums((sqrt(X) - sqrt(Xhat)) ^ 2)) / sqrt(2)
    Mhd <- mean(tmp)
    Mhd.r <- median(tmp)
    
    tmp <- (X - Xhat) ^ 2
    Mse <- mean(tmp)
    Mse.r <- median(tmp)
    
    X1 <- t(t(X) - colMeans(X))
    X2 <- t(t(Xhat) - colMeans(Xhat))
    tmp <- colSums(X1 * X2) / sqrt(colSums(X1 ^ 2) * colSums(X2 ^ 2))
    Msc <- mean(tmp)
    Msc.r <- median(tmp)
    
    X1 <- X - rowMeans(X)
    X2 <- Xhat - rowMeans(Xhat)
    tmp <- rowSums(X1 * X2) / sqrt(rowSums(X1 ^ 2) * rowSums(X2 ^ 2))
    Mtc <- mean(tmp)
    Mtc.r <- median(tmp)
    
    tmp <- sqrt(colSums((X - Xhat) ^ 2))
    Msd <- mean(tmp)
    Msd.r <- median(tmp)
    
    tmp <- sqrt(rowSums((X - Xhat) ^ 2))
    Mtd <- mean(tmp)
    Mtd.r <- median(tmp)
    
    tmp <- (apply(Xhat, 1, Gini) - apply(X, 1, Gini)) ^ 2
    Mgi <- mean(tmp)
    Mgi.r <- median(tmp)
    
    ave1 <- rowMeans(Xhat); ave0 <- rowMeans(X)
    std1 <- apply(Xhat, 1, sd); std0 <- apply(X, 1, sd)
    tmp <- (ave1 - ave0) ^ 2 + (std1 - std0) ^ 2
    Mms <- mean(tmp)
    Mms.r <- median(tmp)
    
    tmp <- (std1 / ave1 - std0 / ave0) ^ 2
    Mcv <- mean(tmp)
    Mcv.r <- median(tmp)
    
    tmp <- (lowerTriangle(cor(t(Xhat))) - lowerTriangle(cor(t(X)))) ^ 2
    Mtt <- mean(tmp)
    Mtt.r <- median(tmp)
    
    m <- nrow(X)
    tmp <- sapply(1 : m, function(j)as.vector(ks.test(X[j, ], Xhat[j, ])$statistic))
    Mks <- mean(tmp)
    Mks.r <- median(tmp)
    
    tmp <- rowSums(abs(t(apply(Xhat, 1, sort)) - t(apply(X, 1, sort))))
    Mws <- mean(tmp)
    Mws.r <- median(tmp)
    
    if(any(X == 0)) X[X == 0] <- 1e-9
    if(any(Xhat == 0)) Xhat[Xhat == 0] <- 1e-9
    
    tmp <- colSums(X * log(X / Xhat))
    Mkl <- mean(tmp)
    Mkl.r <- median(tmp)
    
    tmp <- (-colSums(X * log(X)) + colSums(Xhat * log(Xhat))) ^ 2
    Msh <- mean(tmp)
    Msh.r <- median(tmp)
    
    X1 <- (X + Xhat) / 2
    tmp <- (colSums(X * log(X / X1)) + colSums(Xhat * log(Xhat / X1))) / 2
    Mjs <- mean(tmp)
    Mjs.r <- median(tmp)
    
    return(c(Mse, Msc, Mtc, Msd, Mtd, Msh, Msp, Mbc, Mkl, Mjs, Mhd, Mgi, Mms, Mcv, Mks, Mws, Mtt,
             Mse.r, Msc.r, Mtc.r, Msd.r, Mtd.r, Msh.r, Msp.r, Mbc.r, Mkl.r, Mjs.r, Mhd.r, 
             Mgi.r, Mms.r, Mcv.r, Mks.r, Mws.r, Mtt.r))
  }
  
  l <- length(Xhat.list)
  eval <- sapply(1 : l, function(i)fun(Xhat.list[[i]], X))
  return(eval)
}





