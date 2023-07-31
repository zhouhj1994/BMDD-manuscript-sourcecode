parametric.simu <- function(m, n, q = 0, q1 = 0, p = 0, p1 = 0, r = 0,
                            model = c('gamma', 'log-normal', 'Poisson',
                                      'neg-binomial', 'log-multinormal'),
                            covariate = c(NULL, 'Bernoulli', 'normal')) {
  if(model == 'gamma' | model == 'neg-binomial') {
    a0 <- 0.1; a1 <- 1.5
  } else if(model == 'log-normal' | model == 'log-multinormal') {
    a0 <- 0.1; a1 <- 1.5
  } else if(model == 'Poisson') {
    a0 <- 0.3; a1 <- 10 
  }
  alp0 <- rep(a0, m); alp1 <- rep(a1, m); Z <- NULL
  
  if(!is.null(covariate)) {
    if(covariate == 'Bernoulli') {
      Z <- 2 * rbinom(n, 1, 0.5) - 1
    } else if(covariate == 'normal') {
      Z <- rnorm(n, 0, 1)
    }
    alp0 <- alp0 * exp(rep(0.5, m) %*% t(Z))
    alp1 <- alp1 * exp(rep(1, m) %*% t(Z))
  }
  
  ind.zero <- NULL
  if(q > 0) {
    k <- ceiling(m * q)
    ind.zero <- sample(m, k)
    bina.zero <- matrix(rbinom(k * n, 1, q1), k)
  }
  ind.bimodal <- NULL
  if(p > 0) {
    k <- ceiling(m * (1 - q) * p)
    if(q == 0) {
      ind.bimodal <- sample(m, k)
    } else {
      ind.bimodal <- sample((1 : m)[-ind.zero], k)
    }
    if(!is.null(covariate)) {
      pi <- 1 / (1 + exp(-(log(p1 / (1 - p1)) + rep(1, m) %*% t(Z))))
      pi[-ind.bimodal, ] <- 0
    } else {
      pi <- rep(p1, m)
      pi[-ind.bimodal] <- 0
    }
  } else {
    if(!is.null(covariate)) {
      pi <- matrix(0, m, m)
    } else {
      pi <- rep(0, m)
    }
  }
  del <- matrix(rbinom(m * n, 1, pi), m)
  beta <- alp0 * (1 - del) + alp1 * del
  
  if(model == 'gamma') {
    X <- rgamma(m * n, beta, 1)
  } else if(model == 'log-normal') {
    X <- exp(rnorm(m * n, log(beta), 3))
  } else if(model == 'Poisson') {
    X <- rpois(m * n, beta)
  } else if(model == 'neg-binomial') {
    X <- rnbinom(m * n, size = beta, prob = 0.0001)
  } else if(model == 'log-multinormal') {
    tmp <- array(0.5, c(10, 10))
    rho.sub <- rbind(cbind(tmp, -tmp), cbind(-tmp, tmp))
    diag(rho.sub) <- 1
    rho.list <- list()
    for(i in 1 : (m / 20)) rho.list[[i]] <- rho.sub
    rho <- bdiag(rho.list)
    Sig <- 3 ^ 2 * rho
    eg <- eigen(Sig)
    Sig.sqrt <- eg$vectors %*% diag(sqrt(eg$values)) %*% t(eg$vectors)
    X <- exp(Sig.sqrt %*% matrix(rnorm(m * n), m) + log(beta))
  }
  X <- matrix(X, m)
  tmp <- colSums(X)
  if(any(tmp == 0)) {
    ind <- which(tmp > 0)
    X[, -ind] <- X[, ind[which.min(tmp[ind])]] / 2
  }
  X <- t(t(X) / colSums(X))
  X.ori <- X
  
  if(q > 0) {
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
  
  res <- list(W = W, X = X, X.ori = X.ori, Z = Z, ind.zero = ind.zero, ind.bimodal = ind.bimodal,
              delta = del, alpha = list(alp0 = alp0, alp1 = alp1), pi = pi)
  return(res)
}

