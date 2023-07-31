nonparametric.simu <- function(otu.tab, m, n, q = 0, q1 = 0, r = 0) {
  otu <- array(0, c(m, n))
  while(any(colSums(otu) == 0) | any(rowSums(otu) == 0)) {
    keep.samp <- sample(1 : ncol(otu.tab), n)
    keep.taxa <- sample(1 : nrow(otu.tab), m)
    otu <- otu.tab[keep.taxa, keep.samp]
  }
  
  res <- dmn(t(otu), 1) 
  gamma <- as.vector(fitted(res))
  
  W <- otu + gamma
  X <- t(t(W) / colSums(W))
  X.ori <- X
  
  ind.zero <- NULL
  if(q > 0) {
    k <- ceiling(m * q)
    ind.zero <- sample(m, k)
    bina.zero <- matrix(rbinom(k * n, 1, q1), k)
    X[ind.zero, ] <- X[ind.zero, ] * (1 - bina.zero)
  }
  
  if(r > 0) {
    bina.out <- matrix(rbinom(m * n, 1, r), m)
    X <- X * (1 - bina.out) + matrix(runif(m * n, 1.5, 2.5), m) * bina.out
  }
  
  if(q > 0 | r > 0) X <- t(t(X) / colSums(X))
  N <- rnbinom(n, size = 2, mu = 50 * log(m))
  while(any(N == 0)) {
    N <- rnbinom(n, size = 2, mu = 50 * log(m))
  }
  W <- sapply(1 : n, function(i)rmultinom(1, N[i], X[, i]))
  tmp <- rowSums(W)
  if(any(tmp == 0)){
    k <- sum(tmp == 0)
    W[tmp == 0, ] <- matrix(rep(round(N / max(N)), k), k, byrow = TRUE)
  }
  
  return(list(W = W, X = X, X.ori = X.ori, ind.zero = ind.zero))
}




