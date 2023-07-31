evaluate.fun <- function(X, Xhat) {
  X.cp <- X
  Xhat.cp <- Xhat
  if(any(X == 0)) {
    X <- t(apply(X, 1, function(x){
      if(any(x == 0)){
        if(all(x == 0)) {
          a <- 0.5 * min(X[X != 0])
          x[x == 0] <- ifelse(a > 0, a, min(X[X != 0]))
        } else {
          a <- 0.5 * min(x[x != 0])
          x[x == 0] <- ifelse(a > 0, a, min(x[x != 0]))
        }
      } 
      return(x)
    }))
    X <- t(t(X) / colSums(X))
  }
  
  if(any(Xhat == 0)) {
    Xhat <- t(apply(Xhat, 1, function(x){
      if(any(x == 0)){
        if(all(x == 0)) {
          a <- 0.5 * min(Xhat[Xhat != 0])
          x[x == 0] <- ifelse(a > 0, a, min(Xhat[Xhat != 0]))
        } else {
          a <- 0.5 * min(x[x != 0])
          x[x == 0] <- ifelse(a > 0, a, min(x[x != 0]))
        }
      } 
      return(x)
    }))
    Xhat <- t(t(Xhat) / colSums(Xhat))
  }
  
  Msh <- mean((-colSums(X * log(X)) + colSums(Xhat * log(Xhat))) ^ 2)

  Msp <- mean((colSums(X.cp ^ 2) - colSums(Xhat.cp ^ 2)) ^ 2)
  
  Mkl <- mean(colSums(X * log(X / Xhat)))
  
  Mhd <- mean(sqrt(colSums((sqrt(X.cp) - sqrt(Xhat.cp)) ^ 2)))
  
  Mfr <- sum((X.cp - Xhat.cp) ^ 2) / ncol(X)
  
  tmp <- (X + Xhat) / 2
  Mjs <- mean((colSums(X * log(X / tmp)) + colSums(Xhat * log(Xhat / tmp))) / 2)
  
  X1 <- as.matrix(dist(t(X.cp), diag = TRUE, upper = TRUE))
  Xhat1 <- as.matrix(dist(t(Xhat.cp), diag = TRUE, upper = TRUE))
  Mdd <- mean(sqrt(colSums((X1 - Xhat1) ^ 2)))
  
  X1 <- pmin(X.cp, Xhat.cp)
  Mbc <- mean(1 - 2 * colSums(X1) / colSums(X.cp + Xhat.cp))
  
  rank.fun <- function(x, y) {
    x1 <- rank(x)
    y1 <- rank(y)
    l <- length(x)
    k <- 0
    for(i in 1 : (l - 1)) {
      for(j in (i + 1) : l) {
        a1 <- x1[i] <= x1[j]
        a2 <- y1[i] <= y1[j]
        a3 <- x1[i] >= x1[j]
        a4 <- y1[i] >= y1[j]
        a <- (a1 & a2) | (a3 & a4)
        if(!a) {k <- k + 1}
      }
    }
    return(k*2/(l*(l-1)))
  }
  Mkt <- mean(sapply(1 : ncol(X.cp), function(i)rank.fun(X.cp[, i], Xhat.cp[, i])))
  
  return(c(Msh, Msp, Mkl, Mhd, Mfr, Mjs, Mdd, Mbc, Mkt))
}


