mbDenoise.simu.fun <- function(k){
  n.n = 100
  n.w = 50
  n.factors = 20
  iter = 100
  
  {
    zero.prob <- rep(0,iter)
    zeroinfl.prob <- rep(0,iter)
    depth <- rep(0,iter)
  }
  
  #print(k)
  set.seed(k)
  
  Qn <- matrix(-1,n.n,n.w)
  Q.iter <-1
  while(sum(Qn<0)!=0 ){
    U <- matrix(0,nrow = n.n, ncol = n.factors,byrow = TRUE)
    for(j in 1:n.n){
      U[j,] <- abs(rnorm(n.factors, mean =0, sd = 1))
    }
    
    pr <- rep(0.3,n.factors)
    V1 <- matrix(0,n.w,n.factors,byrow = TRUE)
    for(j in 1:n.w){
      V1[j,] <- rbinom(n.factors, size=1, prob=pr)
    }
    diag(V1) <- 1
    V2 <- matrix(0,nrow = n.w, ncol = n.factors,byrow = TRUE)
    for(j in 1:n.w){
      V2[j,] <- rnorm(n.factors, mean =0, sd = 1e-3)
    }
    V <- V1+V2
    ZV <- U %*% t(V)
    Qn <- ZV/rowSums(ZV)
    X_ori <- ZV
    Q.iter <- Q.iter+1
  }
  
  Pi <- runif(n.n,1,10)
  Ri <- Pi/sum(Pi)
  Ni <- 3*n.n*n.w*Ri
  X <- matrix(0,n.n,n.w,byrow = TRUE)
  for(i in 1:n.n){
    X[i,] <- rmultinom(1, size = Ni[i], prob = Qn[i,])
  }
  colnames(X) <- c(1:ncol(X))
  
  X_cov<- c(rep(1,n.n/2),rep(0,n.n/2))
  
  Y <- log2(1+X)
  zerorow <- which(rowSums(X)==0)
  if(length(zerorow) >0 ){
    X <- X[-zerorow,];X_ori <- X_ori[-zerorow,];Y <- Y[-zerorow,];X_cov<-X_cov[-zerorow];
    Qn <- Qn[-zerorow,];
  }
  zerocol <- which(colSums(X)==0)
  if(length(zerocol) >0 ){
    X <- X[,-zerocol];X_ori <- X_ori[,-zerocol];Y <- Y[,-zerocol];
    Qn <- Qn[,-zerocol]
  }
  
  n.nn <- nrow(X)
  n.mm <- ncol(X)
  
  depth[k] <- summary(rowSums(X))[6]/summary(rowSums(X))[1]
  zero.prob[k]  <- sum(X==0)/(dim(X)[1]*dim(X)[2])
  return(list(W=X, X=Qn, Z=X_cov))
}


