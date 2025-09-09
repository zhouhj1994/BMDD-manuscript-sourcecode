data.generation3 <- function(otu.tab, m, n) {
  otu <- array(0, c(m, n))
  if(n > ncol(otu.tab)) n <- ncol(otu.tab)
  
  if(sum(otu.tab == 0) / (dim(otu.tab)[1] * dim(otu.tab)[2]) > 0.8) {
    while(any(colSums(otu) == 0) | any(rowSums(otu) == 0)) { 
      keep.samp <- order(colSums(otu.tab > 0), decreasing = TRUE)[1 : n]
      otu <- otu.tab[, keep.samp]
      keep.taxa <- sample(order(rowSums(otu > 0), decreasing = TRUE)[1 : 400], m)
      otu <- otu[keep.taxa, ]
    }
  } else {
    while(any(colSums(otu) == 0) | any(rowSums(otu) == 0)) {
      keep.samp <- sample(1 : ncol(otu.tab), n)
      keep.taxa <- sample(1 : nrow(otu.tab), m)
      otu <- otu.tab[keep.taxa, keep.samp]
    }
  }

  W <- as.matrix(otu)
  X <- t(t(W) / colSums(W))
  
  ave <- sapply(1 : m, function(i) mean(log(W[i, W[i, ] > 0])))
  y <- cbind(rowSums(W == 0), rep(n, m))
  fit <- glm(y ~ ave, family = binomial)
  b <- coef(fit)
  v <- rowSums(W > 0)
  drop <- sapply(1 : m, function(i) {
    if(v[i] == 1) {
      return(NULL)
    } else {
      k <- ceiling(v[i] / 5)
      ind <- which(W[i, ] > 0)
      sample(ind, k, prob = 1 / (1 + exp(-b[1] - b[2] * log(W[i, ind]))))
    }
  })
  for(i in 1 : m) {
    W[i, drop[[i]]] <- 0
  }
  
  return(list(W = W, X = X, Z = NULL, U = NULL))
} 