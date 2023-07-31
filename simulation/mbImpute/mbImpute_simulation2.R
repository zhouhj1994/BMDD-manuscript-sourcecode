mbImpute.simu.fun <- function() {
  folder <- 'mbImpute/Simulation2_Figure_1AB/'
  
  T2D_real <- read.csv(paste0(folder,"otu_real_data_T2D.csv"), row.names = "X")
  control_real <- read.csv(paste0(folder,"otu_real_data_control.csv"), row.names = "X")
  
  D <- read.csv(paste0(folder,"D.csv"), row.names = "X")
  
  meta_data_T2D <- read.csv(paste0(folder,"meta_data_T2D.csv"), row.names = "X")
  meta_data_control <- read.csv(paste0(folder,"meta_data_control.csv"), row.names = "X")
  
  real_data <- rbind(T2D_real, control_real)
  real_data_zi_rate <- apply(real_data, 2, FUN = function(x){
    sum(x < (log10(1.01) + 0.01))/length(x)
  })
  chosen_taxa <- which(real_data_zi_rate < 0.85)
  
  meta_real_data <- rbind(meta_data_T2D, meta_data_control)
  
  #set.seed(seed)
  
  meta_sim <- matrix(0, nrow = 50, ncol = 3)
  simulated2 <- matrix(NA, nrow = 50, ncol = length(chosen_taxa))
  # use majority vote to get the meta_data
  for(i in 1:length(chosen_taxa)){
    j = chosen_taxa[i]
    counts <- real_data[,j]
    sample_idx <- sample(which(counts > log10(1.01) + 1e-6), 50)
    #print(length(which(counts > log10(1.01) + 1e-6)))
    #print(length(sample_idx))
    meta_sim <- meta_sim + meta_real_data[sample_idx,]
    simulated2[,i] = real_data[sample_idx,j]
  }
  #dim(simulated2)
  # continuous, take mean
  # categorical, take majority vote
  meta_simulated2 <- meta_sim/length(chosen_taxa)
  
  D_sim <- D[chosen_taxa, chosen_taxa]
  #head(D_sim)
  
  ####### Simulation #######
  #setwd("~/Dropbox/mbimpute/code/Simulation2")
  
  # introduce zeros based on real data.
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
  
  otable <- real_data
  meta_tab <- meta_real_data
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
    win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/5
    mean_up <- mean + win_len 
    mean_lo <- mean - win_len
    sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
  }
  # missing_rate(1, mean_record, percentage_record)
  
  y_sim <- simulated2
  y_preserve <- y_sim
  col_mean <- colMeans(y_sim)
  zero_rate <- unlist(lapply(col_mean, FUN = function(x){
    #print(x)
    return(missing_rate(x, mean_record, percentage_record))
  }))
  #zero_rate[zero_rate > 0.9] = 0.9
  zero_rate[zero_rate > 0.8] = 0.8
  n = dim(y_sim)[1]
  m = dim(y_sim)[2]
  zero_mat <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:m){
    zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
  }
  sim_tab_zi = y_sim * zero_mat
  sim_tab_zi[sim_tab_zi < log10(1.01)+1e-6] = log10(1.01)
  
  #sqrt(sum((sim_tab_zi - y_sim)^2))
  sim_tab_zi_mbImpute <- sim_tab_zi
  
  #write.csv(sim_tab_zi, "simulated_zi_matrix_mbImpute.csv")
  #write.csv(meta_simulated2, "simulated_meta_data_mbImpute.csv")
  #write.csv(D_sim, "D_sim.csv")
  
  #dim(read.csv("simulated_zi_matrix_mbImpute.csv", row.names = "X"))
  
  sim_tab_zi <- 10^(sim_tab_zi) - 1.01
  sim_tab_zi_trans <- t(sim_tab_zi)
  #write.csv(sim_tab_zi_trans, "simulated_zi_matrix.csv")
  
  #write.csv(t(sim_tab_zi_trans), "simulated_zi_matrix_Magic.csv")
  
  y_sim_rec <- y_sim
  sim_tab_zi_rec <- sim_tab_zi
  # 50 * 145
  
  #write.csv(y_sim_rec, "truth.csv")
  
  W <- sim_tab_zi_trans
  Y <- 10^y_sim_rec-1.01
  X <- t(Y / rowSums(Y))
  Z <- meta_simulated2
  D <- D_sim
  return(list(W=W, X=X, Z=Z, D=D))
}



