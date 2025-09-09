data.generation2 <- function(m = 100, n, a, b, q, r, model) {

  if(model == 'gamma') {
    alp0 <- runif(30, 0.1, 0.3); alp1 <- runif(30, 1.5, 2)
    alp2 <- c(runif(20, 0.1, 0.3), runif(50, 0.1, 2))
    alp2 <- sample(alp2, 70)
  } else if(model == 'log-normal') {
    alp0 <- runif(30, 0.1, 0.3); alp1 <- runif(30, 1.5, 2)
    alp2 <- c(runif(20, 0.1, 0.3), runif(50, 0.1, 2))
    alp2 <- sample(alp2, 70)
  } else if(model == 'Poisson') {
    alp0 <- runif(30, 1, 3); alp1 <- runif(30, 15, 2000)
    alp2 <- c(runif(20, 1, 3), runif(50, 1, 2000))
    alp2 <- sample(alp2, 70)
  } else if(model == 'neg-binomial') {
    alp0 <- runif(30, 10, 30); alp1 <- runif(30, 150, 2000)
    alp2 <- c(runif(20, 10, 30), runif(50, 10, 2000))
    alp2 <- sample(alp2, 70)
  }
  
  pi <- rep(0.5, 30)
  
  del <- matrix(rbinom(30 * n, 1, pi), 30)
  del2 <- t(del)
  del2[, 2 : 5] <- del2[, 1]
  del2[, 6 : 10] <- 1 - del2[, 1]
  del <- t(del2)
  
  ind <- colSums(del[c(11, 12), ]) == 2
  del[13 : 16, ind] <- 1
  del[17 : 20, ind] <- 0
  
  beta <- array(NA, c(m, n))
  beta[1 : 30, ] <- alp0 * (1 - del) + alp1 * del
  beta[31 : 100, ] <- alp2
  
  if(model == 'gamma') {
    X <- rgamma(m * n, beta, 1)
  } else if(model == 'log-normal') {
    X <- exp(rnorm(m * n, log(beta), 2))
  } else if(model == 'Poisson') {
    X <- rpois(m * n, beta)
  } else if(model == 'neg-binomial') {
    k <- 0.5
    X <- rnbinom(m * n, size = k, prob = k / (k + beta))
  }
  X <- matrix(X, m)
  
  u <- sum(ind)
  X[31 : 35, ind] <- t(sapply(31 : 35, function(i)runif(u, quantile(X[i, ], 0.75), max(X[i, ]))))
  X[36 : 40, ind] <- t(sapply(36 : 40, function(i)runif(u, min(X[i, ]), quantile(X[i, ], 0.25))))
  
  ind <- X[41, ] > quantile(X[41, ], 0.75)
  u <- sum(ind)
  X[42 : 45, ind] <- t(sapply(42 : 45, function(i)runif(u, quantile(X[i, ], 0.75), max(X[i, ]))))
  X[46 : 50, ind] <- t(sapply(46 : 50, function(i)runif(u, min(X[i, ]), quantile(X[i, ], 0.25))))
  
  X[51 : 60, ] <- t(apply(X[51 : 60, ], 1, sort))
  X[61 : 70, ] <- t(apply(X[61 : 70, ], 1, sort, decreasing = TRUE))
  
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
  
  return(list(W = W, X = X, Z = NULL, U = NULL))
}

