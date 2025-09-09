library(foreach)
library(xtable)

parentheses <- function(x){
  sprintf("(%s)", x)
}

tab.fun <- function(setup, m = NULL, n = NULL, a = NULL, b = NULL, q = NULL, r = NULL, 
                    nsim, method, metric, type, robust = FALSE) {
  if(type == 'Non-parametric') {
    output.folder <- paste0('output/', setup, '/m_', m, '_n_', n, '/')
  } else {
    output.folder <- paste0('output/', setup, '/m_', m, '_n_', n, '_a_', a, 
                            '_b_', b, '_q_', q, '_r_', r, '/')
  }
  
  output.raw <- foreach(i = 1 : nsim, .combine = 'rbind') %do% {
    res <- readRDS(paste0(output.folder, i, ".rds"))
    res$res
  }
  
  s <- nrow(output.raw) / nsim
  k <- s / 2
  l <- ncol(output.raw)
  # 
  # if(l == 11) {
  #   ind <- c(1, 4, 8, 2, 3, 5 : 7, 9 : 11)
  # } else if(l == 14) {
  #   ind <- c(1, 12, 4, 8, 2, 13, 3, 14, 5 : 7, 9 : 11)
  # } else if(l == 12) {
  #   ind <- c(1, 12, 4, 8, 2, 3, 5 : 7, 9 : 11)
  # }
  # method <- method[ind]
  # output.raw <- output.raw[, ind]

  output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
    res <- output.raw[seq(i, nsim * s, s), ]
    est <- colMeans(res, na.rm = TRUE)
    num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
    est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
    return(c(est, est.sd))
  }
  
  M <- ifelse(output > .001, sprintf("%.4f", output), sprintf("%.2e", output))
  M[, (l + 1) : (2 * l)] <- parentheses(M[, (l + 1) : (2 * l)])
  
  v <- rep(0, s); v[c(2, 3, 2 + k, 3 + k)] <- 1
  ord <- t(sapply(1 : s, function(i) {
    order(output[i, 1 : l], decreasing = v[i])
  }))
  ind1 <- ord[, 1]
  ind2 <- ord[, 2]
  
  # for(i in 1 : s) {
  #   #M[i, ind1[i]] <- paste0(M[i, ind1[i]], '*')
  #   #M[i, ind2[i]] <- paste0(M[i, ind2[i]], '**')
  #   M[i, ind1[i]] <- paste0('{\\color{red}\\textbf{', M[i, ind1[i]], '}}')
  #   M[i, ind2[i]] <- paste0('{\\color{red}', M[i, ind2[i]], '}')
  # }

  M <- foreach(i = 1 : l, .combine = 'cbind') %do% {
    paste(M[, i], M[, i + l], sep = '') 
  }
  colnames(M) <- method
  rownames(M) <- rep(metric, 2)
  
  if(!robust) {
    table_latex <- xtable(M[1 : k, ])
  } else {
    table_latex <- xtable(M[(k + 1) : (2 * k), ])
  }
  print(table_latex, include.rownames = TRUE, comment = FALSE, 
        sanitize.text.function = function(x) x)
}

method <- c('BMDD1', 'BMDD2')

metric <- c('Mean squared error', 'Sample-wise correlation', 'Taxon-wise correlation',
            'Sample-wise distance', 'Taxon-wise distance', 
            'Shannon\'s index', 'Simpson\'s index', 'Bray-Curtis dissimilarity',
            'Kullback–Leibler divergence', 'Jensen–Shannon divergence', 'Hellinger distance',
            'Gini coefficient', '(Mean, Standard deviation)', 'Coefficient of variation',
            'Kolmogorov–Smirnov distance', 'Wasserstein distance', 'Pairwise taxon-to-taxon correlation')

nsim <- 10
m <- 100
n <- 80
a <- 1000
b <- 0.5
q <- 6
r <- 0

tab.fun('S1', m, n, a, b, q, r, nsim, method, metric, type = 'Theoretical')

