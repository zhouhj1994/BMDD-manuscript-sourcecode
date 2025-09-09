library(foreach)

num.rej.fun <- function(pval.mat, curve.fdr.cutoffs) {
  qval.mat <- sapply(1 : ncol(pval.mat), function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  isna <- colSums(is.na(qval.mat)) == nrow(pval.mat)
  
  num.rej <- t(sapply(1 : length(curve.fdr.cutoffs), function(i)
    colSums(qval.mat <= curve.fdr.cutoffs[i], na.rm = TRUE)))
  num.rej[, isna] <- NA
  return(num.rej)
}

datanames <- c("IBD-1", "IBD-2", "IBD-3", "IBD-4")
n.dat <- length(datanames)
ind <- 1 : 6
n.met <- length(ind)
curve.fdr.cutoffs <- seq(0.01, 0.25, 0.01)
n.cut <- length(curve.fdr.cutoffs)
filenames <- list.files("output")
n.shuf <- length(filenames)

num.rej.list <- list()
exist.rej.list <- list()
for(i in 1 : n.dat) {
  num.rej <- NULL
  exist.rej <- NULL
  for(j in 1 : n.shuf) {
    pval.mat.list <- readRDS(paste0("output/", filenames[j]))
    tmp1 <- num.rej.fun(pval.mat.list[[i]][, ind], curve.fdr.cutoffs)
    tmp2 <- (tmp1 > 0) + 0
    num.rej <- rbind(num.rej, tmp1)
    exist.rej <- rbind(exist.rej, tmp2)
  }
  num.rej.list[[i]] <- num.rej
  exist.rej.list[[i]] <- exist.rej
}

output.list <- list()
for(i in 1 : n.dat) {
  output.raw <- exist.rej.list[[i]]
  output <- foreach (k = 1 : n.cut, .combine = 'rbind') %do% {
    res <- output.raw[seq(k, n.shuf * n.cut, n.cut), ]
    est <- colMeans(res, na.rm = TRUE)
    num.nan <- colSums(matrix(!is.na(res), nrow = n.shuf))
    est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
    
    res <- c(est, est.sd)
    res[which(is.nan(res))] <- NA
    res
  }
  output.list[[i]] <- output
}

saveRDS(output.list, 'fdr.list.rds')

