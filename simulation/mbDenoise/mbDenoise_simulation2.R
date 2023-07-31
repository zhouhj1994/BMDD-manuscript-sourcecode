mbDenoise.simu.fun <- function(k){
  n.n = 60
  n.w = 100
  n.factors = 2
  iter = 100
  
  {
    zero.prob <- rep(0,iter)
    zeroinfl.prob <- rep(0,iter)
    depth <- rep(0,iter)
  }
  
  Q.iter <-1
  n.mm <- 0
  X <- matrix(0,n.n,n.w,byrow = TRUE)
  while(n.mm!=n.w | sum(is.na(X)) >0){
    if(Q.iter>1){
      s <- sample(100:2000,1)
      set.seed(k+s)
      print(k+s)
    }else{
      set.seed(k)
      #print(k)
    }
    f <- matrix(0,nrow = n.n, ncol = n.factors)
    for(i in 1:n.n){
      f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
    }
    betaj <- matrix(0,nrow = n.w, ncol = n.factors)
    for(j in 1:n.w){
      betaj[j,] <- runif(n.factors,-3,3)
    }
    alpha <- runif(n.n,-5,5)
    beta0 <- rep(0,n.w)
    g <- rep(0,n.w*0.5*0.5)
    ## in DA test, 0= 1,2,3,4,5, such as g <- rep(1,n.w*0.5*0.5)
    gamma <- c(g,-g,rep(0,(n.w-n.w*0.5)))
    X_cov<- c(rep(1,n.n/2),rep(0,n.n/2))
    
    ll <- f %*% t(betaj) + matrix(gamma,n.n,n.w,byrow=TRUE)*X_cov +matrix(alpha,n.n,n.w)+matrix(beta0,n.n,n.w,byrow=TRUE)
    exp_mat <- exp(ll)
    eta_mat <- matrix(0.25,n.n,n.w,byrow=TRUE)
    
    z <- matrix(0,n.n,n.w,byrow = TRUE)
    for(i in 1:n.n){
      z[i,] <- rbinom(n.w, size=1, prob=eta_mat[i,])
    }
    sum <- rowSums((1-z)*exp_mat)
    Qn_z <- (1-z)*exp_mat/sum
    
    sum <- rowSums(exp_mat)
    Qn <- exp_mat/sum
    
    X <- matrix(0,n.n,n.w,byrow = TRUE)
    for(i in 1:n.n){
      for(j in 1:n.w){
        X[i,j] <- rpois(n=1,lambda = exp_mat[i,j])
      }
    }
    X[z==1]=0
    colnames(X) <- c(1:ncol(X))
    
    Y <- log2(1+X)
    X_ori <- (1-z)*exp(f %*% t(betaj)+ matrix(gamma,n.n,n.w,byrow=TRUE)*X_cov +matrix(beta0,n.n,n.w,byrow=TRUE))
    
    zerorow <- which(rowSums(X)==0)
    if(length(zerorow) >0 ){
      X <- X[-zerorow,];X_ori <- X_ori[-zerorow,];Y <- Y[-zerorow,];X_cov<-X_cov[-zerorow];
      f <- f[-zerorow,];Qn <- Qn[-zerorow,];z <- z[-zerorow,]
    }
    zerocol <- which(colSums(X)==0)
    if(length(zerocol) >0 ){
      X <- X[,-zerocol];X_ori <- X_ori[,-zerocol];Y <- Y[,-zerocol];
      betaj <- t(t(betaj)[,-zerocol]);Qn <- Qn[,-zerocol];z <- z[,-zerocol]
    }
    
    n.nn <- nrow(X)
    n.mm <- ncol(X)
    
    depth[k] <- summary(rowSums(X))[6]/summary(rowSums(X))[1]
    zero.prob[k]  <- sum(X==0)/(dim(X)[1]*dim(X)[2])
    zeroinfl.prob[k] <- sum(z==1)/(sum(X==0))
    
    Q.iter <- Q.iter+1
  }
  return(list(W=X, X=Qn, Z=X_cov))
}


