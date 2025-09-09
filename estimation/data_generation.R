data.generation <- function(alp0, alp1, pi, m, n, a, b, q, r, Z = NULL, U = NULL) {
  ind.taxa <- sample(1 : 100, m, replace = TRUE)
  ind.samp <- sample(1 : 80, n, replace = TRUE)
  if(is.null(Z) & is.null(U)) {
    pi <- pi[ind.taxa]
    alp0 <- alp0[ind.taxa]
    alp1 <- alp1[ind.taxa]
  } else {
    alp0 <- alp0[ind.taxa, ind.samp]
    alp1 <- alp1[ind.taxa, ind.samp]
    pi <- pi[ind.taxa, ind.samp]
    if(!is.null(Z) & is.null(U)) {
      Z <- Z[ind.samp]
    } else if(is.null(Z) & !is.null(U)) {
      U <- U[ind.taxa]
    }
  }
  
  del <- matrix(rbinom(m * n, 1, pi), m)
  beta <- alp0 * (1 - del) + alp1 * del
  
  X <- matrix(rgamma(m * n, beta, 1), m)
  X <- t(t(X) / colSums(X))
  
  if(q > 0) {
    X[X < 10 ^ (-q)] <- 0
    X <- t(t(X) / colSums(X))
  }
  
  N <- rnbinom(n, size = b, mu = a)
  while(any(N == 0)) {
    N <- rnbinom(n, size = b, mu = a)
  }
  W <- sapply(1 : n, function(i)rmultinom(1, N[i], X[, i]))
  tmp <- rowSums(W)
  if(any(tmp == 0)){
    k <- sum(tmp == 0)
    W[tmp == 0, ] <- matrix(rep(round(N / max(N)), k), k, byrow = TRUE)
  }
  
  if(r > 0) {
    ind <- sample(1 : n, r)
    for(i in ind) {
      W[, i] <- sample(W[, i], m)
    }
  }
  
  return(list(W = W, X = X, Z = Z, U = U, alp0 = alp0, alp1 = alp1, pi = pi))
}

