mbImpute.simu.fun <- function() {
  folder <- 'mbImpute/Simulation1_Figure_1ABCD/'
  
  control_coeff <- readRDS(paste0(folder,"control_dat1_sim_add_filter_coef.rds"))
  T2D_coeff <- readRDS(paste0(folder,"T2D_dat2_sim_add_filter_coef.rds"))
  
  otu_real_data_T2D <- read.csv(paste0(folder,"otu_real_data_T2D.csv"), row.names = "X")
  otu_real_data_control <- read.csv(paste0(folder,"otu_real_data_control.csv"), row.names = "X")
  
  D <- read.csv(paste0(folder,"D.csv"), row.names = "X")
  
  meta_data_T2D <- read.csv(paste0(folder,"meta_data_T2D.csv"), row.names = "X")
  meta_data_control <- read.csv(paste0(folder,"meta_data_control.csv"), row.names = "X")
  
  #set.seed(seed)
  
  y_sim = otu_real_data_T2D
  x = meta_data_T2D
  k = 5
  c1 = T2D_coeff
  c1 = c1*5
  c1[1] = T2D_coeff[1]
  c1[length(c1)-6] = 0.2
  c1[length(c1)-11] = -0.3
  
  #loading used functions
  design_mat_row_gen2 <- function(count_mat, covariate_mat, row_index, col_index, close_taxa){
    n = dim(count_mat)[1]
    m = dim(count_mat)[2]
    k = length(close_taxa[[1]])
    if(is.vector(covariate_mat)){
      p = 1
      i = row_index
      j = col_index
      
      close_taxa_set <- close_taxa[[j]]
      #generate a row including response and a row of design matrix.
      row_gen <- rep(0, m*k + (n-1) * n + n*p)
      
      row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i])
    }else{
      p = dim(covariate_mat)[2]
      i = row_index
      j = col_index
      
      close_taxa_set <- close_taxa[[j]]
      #generate a row including response and a row of design matrix.
      row_gen <- rep(0, m*k + (n-1) * n + n*p)
      
      row_gen[((j-1)*k+1):(j*k)] = as.numeric(count_mat[i,close_taxa[[j]]])
      row_gen[(m*k + (i-1)*(n-1)+1):(k*m + i*(n-1))] = as.numeric(count_mat[-i,j])
      row_gen[(m*k + n*(n-1)+ p*(i-1) + 1):(m*k + n*(n-1) + i*p)] = as.numeric(covariate_mat[i,])
    }
    return(row_gen)
  }
  
  #identify the low abundance taxa by binomial test
  filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    n = length(y)
    nz <- sum(y <= (log10(1.01) + 1e-6))
    pz = 1 - nz/n
    test = pz - 1.96 * sqrt(pz * (1-pz)/n)
    if(nz == n || test <= 0){
      return(0)
    }else{
      return(1)
    }
  })
  
  #keep the original matrix
  #y_imp <- y_sim
  # print(dim(y_imp))
  # #keep the vectors that we neither impute nor borrow information
  # remain_vec <- which(unlist(filter) == 0)
  # y_rem <- y_sim[, remain_vec]
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]
  D = D[filter_vec,filter_vec]
  y_preserve <- y_sim
  #apply the imputation method on the simulated matrix
  m = dim(y_sim)[2]
  n = dim(y_sim)[1]
  # 53 * 193
  
  #x[1:50,]
  #X <- cbind(rep(1,n), x[1:n,])
  X <- as.matrix(x)
  
  #identifying the group needs imputation and group doesn't need imputation
  
  #find k closest taxa
  # k = 5
  #generate close set
  close_taxa <- list()
  for(j in 1:m){
    #close_taxa[[j-1]] = which(D[,j] %in% sort(unique(D[,j]))[2:(k+1)])
    close_dist <- D[D[,j] %in% sort(D[,j])[2:(k+1)],j]
    close_taxa_vec = which(D[,j] %in% close_dist[close_dist != max(close_dist)])
    if(length(close_taxa_vec) < k){
      close_taxa_vec <- c(close_taxa_vec, which(D[,j] == max(close_dist))[1:(k-length(close_taxa_vec))])
    }
    close_taxa[[j]] = close_taxa_vec
  }
  #generate design matrix
  p = dim(x)[2]
  if(is.null(p)){
    p = 1
  }
  #generate a row including response and a row of design matrix.
  row_length <- m * k + (n-1) * n + n*p
  idx_set <- matrix(NA, m*n, 2)
  for(i in 1:n){
    for(j in 1:m){
      idx_set[(i-1)*m+j, ] <- c(i, j)
    }
  }
  
  
  # for(itr in 1:5){
  #   print(itr)
  #perform imputation on the rest
  
  design_mat_gen <- matrix(0, nrow = dim(idx_set)[1], ncol = row_length)
  for(i in 1:dim(idx_set)[1]){
    design_mat_gen[i,] <- design_mat_row_gen2(y_sim, x[1:n,], idx_set[i,1], idx_set[i,2], close_taxa)
  }
  design_mat_gen <- cbind(rep(1, dim(design_mat_gen)[1]), design_mat_gen)
  imputed_value <- design_mat_gen %*% c1
  impute_mat <- y_sim
  for(i in 1:dim(idx_set)[1]){
    impute_mat[idx_set[i,1], idx_set[i,2]] = max(imputed_value[i], log10(1.01))
  }
  #print(sqrt(sum((impute_mat - y_sim)^2)))
  y_sim <- impute_mat
  # }
  
  y_sim <- y_sim - 1.5
  #write.csv(y_sim, "Karlsson_simulated.csv")
  
  gamma_norm_mix <- function(y, X){
    loglik <- function(p, alpha, beta, cov_par, var1, X, y){
      n = length(y)
      lkval <- 0
      fgam <- dgamma(y, shape = alpha, rate = beta)
      for(i in 1:n){
        if(!is.vector(X)){
          lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i,] %*% cov_par, sd = sqrt(var1)))
        }else{
          lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i] * cov_par, sd = sqrt(var1)))
        }
        
      }
      return(lkval)
    }
    n = length(y)
    alpha_init <- 1
    beta_init <- 10
    p_init <- 0.5
    cov_par_init <- solve(t(X) %*% X) %*% t(X) %*% y
    var_init <- t(y - X %*% cov_par_init) %*% (y - X %*% cov_par_init) / n
    
    alpha_t <- alpha_init
    beta_t <- beta_init
    cov_par_t <- cov_par_init
    var_t <- var_init
    p_t <- p_init
    
    #update gamam param
    #Wei's Method
    ### root-finding equation
    fn = function(alpha, target){
      log(alpha) - digamma(alpha) - target
    }
    update_gmm_pars = function(x, wt){
      if(max(wt) > 0.00001){
        tp_s = sum(wt)
        tp_t = sum(wt * x)
        tp_u = sum(wt * log(x))
        tp_v = -tp_u / tp_s - log(tp_s / tp_t)
        if (tp_v <= 0){
          alpha = 20
        }else{
          alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
          if (alpha0 >= 20){alpha = 20
          }else{
            alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v,
                            extendInt = "yes")$root
          }
        }
        ## need to solve log(x) - digamma(x) = tp_v
        ## We use this approximation to compute the initial value
        beta = tp_s / tp_t * alpha
      }else{
        alpha = 0.001
        beta = 1000
      }
      
      return(c(alpha, beta))
    }
    #convergence criteria
    flag = TRUE
    maxitr = 300
    itr = 0
    while(flag){
      #E_step
      mean_t <- X %*% cov_par_t
      n = length(y)
      a_hat_t <- rep(0, n)
      dg_t <- dgamma(y, shape = alpha_t, rate = beta_t)
      for(i in 1:n){
        if(dg_t[i] == 0){
          a_hat_t[i] = 0
        }else{
          a_hat_t[i] <- p_t * dg_t[i]/(p_t * dg_t[i]+ (1-p_t)*dnorm(y[i], mean = mean_t[i], sd = sqrt(var_t)))
        }
      }
      #maximization
      #fit p
      p_t1 <- sum(a_hat_t)/n
      X_tilta <- sqrt(1-a_hat_t) * X
      y_tilta <- sqrt(1-a_hat_t) * y
      #fit normal
      out <- tryCatch(
        {
          # Just to highlight: if you want to use more than one
          # R expression in the "try" part then you'll have to
          # use curly brackets.
          # 'tryCatch()' will return the last evaluated expression
          # in case the "try" part was completed successfully
          cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
          TRUE
        },
        error=function(cond) {
          FALSE
        }
      )
      if(!out){
        return(list("d" = y < log10(1.01) + 10^(-3) ))
      }
      cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
      var_t1 <- sum((1 - a_hat_t) * (y - X %*% cov_par_t)^2) / sum(1-a_hat_t)
      #fit gamma
      par_gamm <- update_gmm_pars(x = y, wt = a_hat_t)
      alpha_t1 = par_gamm[1]
      beta_t1 <- par_gamm[2]
      loglik1 <- loglik(p = p_t, alpha = alpha_t, beta = beta_t, cov_par = cov_par_t, var1 = var_t, X = X, y = y)
      loglik2 <- loglik(p = p_t1, alpha = alpha_t1, beta = beta_t1, cov_par = cov_par_t1, var1 = var_t1, X = X, y = y)
      if((abs(loglik1 - loglik2)) < 0.05 || itr > maxitr){
        flag = FALSE
      }else{
        alpha_t <- alpha_t1
        beta_t <- beta_t1
        cov_par_t <- cov_par_t1
        var_t <- var_t1
        p_t <- p_t1
        itr = itr + 1
      }
    }
    #fit a normal curve
    eta_hat <- cov_par_init
    omega_hat <- var_init
    #calculate new likelihood
    norm_log_lik = 0
    for(i in 1:n){
      if(!is.vector(X)){
        norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i,] %*% eta_hat, sd = sqrt(omega_hat)))
      }else{
        norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i] * eta_hat, sd = sqrt(omega_hat)))
      }
    }
    Dev = -2 * norm_log_lik - (-2 * loglik2)
    judgement = pchisq(Dev, df = 3, lower.tail = FALSE) < 0.05
    
    if(!judgement || (alpha_t/beta_t > 1)){
      p_t = 0
      a_hat_t <- rep(0,n)
    }
    return(list("p" = p_t, "alpha" = alpha_t, "beta" = beta_t, "cov_par" = cov_par_t, "var" = var_t, "d" = a_hat_t, "eta" = eta_hat, "omega" = omega_hat, "Deviance" = Dev))
  }
  
  ####### learn from learning missing rate empirically ############
  
  otable <- y_preserve
  meta_tab <- meta_data_T2D
  #D
  
  #dim(otable)
  # read in the gamma_norm_mix function
  mean_record <- c()
  percentage_record <- c()
  D_vale_record <- c()
  beta_record <- list()
  for(j in 1:dim(otable)[2]){
    result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab))
    beta_record[[j]] = result$cov_par
    mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
    percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
    D_vale_record <- c(D_vale_record, result$d)
  }
  #filter <- which(is.nan(mean_record) | mean_record < log10(1.01) + 0.001)
  
  # filter out the nan values.
  #plot(mean_record, percentage_record)
  
  # build a map between the mean and percentage of missing
  missing_rate <- function(mean, emp_mean, emp_miss){
    win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/3
    mean_up <- mean + win_len
    mean_lo <- mean - win_len
    sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
  }
  # missing_rate(1, mean_record, percentage_record)
  
  col_mean <- colMeans(y_sim)
  zero_rate <- unlist(lapply(col_mean, FUN = function(x){
    #print(x)
    return(missing_rate(x, mean_record, percentage_record))
  }))
  zero_mat <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:m){
    zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
  }
  sim_tab_zi = y_sim * zero_mat
  sim_tab_zi[sim_tab_zi < log10(1.01)+1e-6] = log10(1.01)
  
  #sqrt(sum((sim_tab_zi - y_preserve)^2))
  sim_tab_zi_mbImpute <- sim_tab_zi
  
  #write.csv(sim_tab_zi, "simulated_zi_matrix_mbImpute.csv")
  
  sim_tab_zi <- 10^(sim_tab_zi) - 1.01
  sim_tab_zi_trans <- t(sim_tab_zi)
  #write.csv(sim_tab_zi_trans, "simulated_zi_matrix.csv")
  
  y_sim_rec <- y_sim
  sim_tab_zi_rec <- sim_tab_zi
  
  #write.csv(y_sim_rec, "truth.csv")
  
  W <- sim_tab_zi_trans
  Y <- 10^y_sim_rec-1.01
  X <- t(Y / rowSums(Y))
  Z <- meta_data_T2D
  D <- D
  return(list(W=W, X=X, Z=Z, D=D))
}























