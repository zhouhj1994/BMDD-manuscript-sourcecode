library(foreach)
library(ggplot2)

fig.fun <- function(setup, m = NULL, n = NULL, a = NULL, b = NULL, q = NULL, r = NULL, 
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
  
  if(l == 11) {
    ind <- c(1, 4, 8, 2, 3, 5 : 7, 9 : 11)
  } else if(l == 14) {
    ind <- c(1, 12, 4, 8, 2, 13, 3, 14, 5 : 7, 9 : 11)
  } else if(l == 12) {
    ind <- c(1, 12, 4, 8, 2, 3, 5 : 7, 9 : 11)
  }
  method <- method[ind]
  output.raw <- output.raw[, ind]
  
  output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
    res <- output.raw[seq(i, nsim * s, s), ]
    est <- colMeans(res, na.rm = TRUE)
    num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
    est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
    return(c(est, est.sd))
  }
  
  
  M <- output[1 : k, 1 : l][-c(2, 3), ]
  tmp1 <- (M[, -1] - M[, 1]) / M[, -1]
  tmp2 <- (M[, 1] - M[, -1]) / M[, 1]
  
  A <- tmp1 * (tmp1 > 0) - tmp2 * (tmp2 > 0)
  colnames(A) <- method[-1]

  A[A < 0] <- A[A < 0] / 10
  
  method1 <- factor(1 : (l - 1), labels = method[-1])
  df <- data.frame(method = rep(method1, each = k - 2),
                   value = as.vector(A))
  
  custom_breaks <- c(1, 0.5, 0, -0.05, -0.1)
  custom_labels <- c("100%", "50%", "0%", "50%", "100%")
  
  pdf(paste0(setup,'_boxplot.pdf'), width = 12, height = 8)
  gp <- ggplot(df, aes(x = method, y = value))  +
    geom_boxplot(fill = "steelblue", color = "black", na.rm = TRUE) +
    theme_minimal(base_size = 25) +
    scale_y_continuous(
      breaks = custom_breaks,
      labels = custom_labels
    ) +
    xlab("") +
    ylab("BMDD Improvement Percentage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = unit(c(0.5,0.5,-0.5,0.5), "cm")) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1)
  print(gp)
  dev.off()
}


method <- c('BMDD', 'mbDenoise', 'mbImpute', 'SAVER', 'scImpute', 'ALRA', 'DMM',
            'naive-1', 'naive-2', 'naive-3', 'naive-4')

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

# setup = 'S1'
# robust = FALSE
# type = 'Theoretical'

fig.fun('S1', m, n, a, b, q, r, nsim, method, metric, type = 'Theoretical')

fig.fun('S6', m, n, a, b, q, r, nsim, method, metric, type = 'Correlation')
fig.fun('S7', m, n, a, b, q, r, nsim, method, metric, type = 'Correlation')
fig.fun('S8', m, n, a, b, q, r, nsim, method, metric, type = 'Correlation')
fig.fun('S9', m, n, a, b, q, r, nsim, method, metric, type = 'Correlation')

fig.fun('S10', m, n, a, b, q, r, nsim, method, metric, type = 'Non-parametric')
fig.fun('S11', m, n, a, b, q, r, nsim, method, metric, type = 'Non-parametric')
fig.fun('S12', m, n, a, b, q, r, nsim, method, metric, type = 'Non-parametric')
fig.fun('S13', m, n, a, b, q, r, nsim, method, metric, type = 'Non-parametric')


