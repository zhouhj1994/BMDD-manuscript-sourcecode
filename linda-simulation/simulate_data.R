simulate.data <- function(m, n, gamma, beta0, sigma2, eta0, theta0, kappa0, lambda0, model) {
  if(m < 500) {
    ind <- sample(1 : 500, m)
    beta0 <- beta0[ind]
    sigma2 <- sigma2[ind]
    eta0 <- eta0[ind]
    theta0 <- theta0[ind]
    kappa0 <- kappa0[ind]
    lambda0 <- lambda0[ind]
  }
  
  if(model == 'S1') {
    if(n == 50) {
      mu <- 2
    } else if(n == 200) {
      mu <- 1.5
    }
  } else if(model == 'S2') {
    if(n == 50) {
      mu <- 3
    } else if(n == 200) {
      mu <- 1.5
    }
  } else if(model == 'S3') {
    if(n == 50) {
      mu <- 6
    } else if(n == 200) {
      mu <- 3
    }
  } else if(model == 'S4') {
    if(n == 50) {
      mu <- 10
    } else if(n == 200) {
      mu <- 5
    }
  }
    
  if(m == 50) {
    lib <- 1000
  } else if(m == 200) {
    lib <- 3000
  } else if(m == 500) {
    lib <- 7645
  }
  
  if(model == 'S1') {
    X0 <- matrix(rgamma(m * n, shape = eta0, rate = 1), nrow = m)
  } else if(model == 'S2') {
    X0 <- matrix(exp(rnorm(m * n, beta0, sqrt(sigma2))), nrow = m)
  } else if(model == 'S3') {
    X0 <- matrix(rpois(m * n, lambda = exp(lambda0 * lib)), nrow = m)
  } else if(model == 'S4') {
    X0 <- matrix(rnbinom(m * n, size = theta0, mu = exp(kappa0 * lib)), nrow = m)
  }
  
  pi0 <- t(t(X0) / colSums(X0))
  pi0.ave <- rowMeans(pi0)
  if (any(pi0.ave == 0)) {
    ind <- which(pi0.ave == 0)
    pi0.ave[ind] <- min(pi0.ave[-ind]) / 10
    pi0.ave <- pi0.ave / sum(pi0.ave)
  }
  tmp <- (pi0.ave > 0.005)
  mu.1 <- log(mu * tmp + mu * (0.005 / pi0.ave) ^ (1 / 3) * (1 - tmp))
  H <- rbinom(m, 1, gamma)
  alpha <- mu.1 * H
  
  u <- rbinom(n, 1, 0.5)
  while(length(unique(u)) == 1) {
    u <- rbinom(n, 1, 0.5)
  }
  Z <- cbind(u)
  beta <- alpha
  
  N <- rnbinom(n, size = 5.3, mu = lib)
  if(model == 'S1') {
    tmp <- exp(beta %*% t(Z)) * eta0
    X <- matrix(rgamma(m * n, shape = tmp, rate = 1), nrow = m)
    pi <- t(t(X) / colSums(X))
    Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pi[, s]))
  } else if(model == 'S2') {
    tmp <- beta %*% t(Z) + beta0
    X <- matrix(exp(rnorm(m * n, tmp, rep(sqrt(sigma2), n))), nrow = m)
    pi <- t(t(X) / colSums(X))
    Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pi[, s]))
  } else if(model == 'S3') {
    X <- exp(lambda0 %*% t(N) + beta %*% t(Z))
    Y <- matrix(rpois(m * n, lambda = X), nrow = m)
  } else if(model == 'S4') {
    X <- exp(kappa0 %*% t(N) + beta %*% t(Z))
    Y <- matrix(rnbinom(m * n, size = theta0, mu = X), nrow = m)
  }
  
  Z <- as.data.frame(Z)
  Y <- as.data.frame(Y)
  colnames(Y) <- rownames(Z) <- paste0('sample', 1 : n)
  rownames(Y) <- paste0('taxon', 1 : m)
  return(list(Y = Y, Z = Z, H = H, X = X, beta = beta))
}


