winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

preprocess.fun <- function(otu.tab, meta, prev.cut = 0.2, lib.cut = 1000, 
                           winsor.quan = 0.97) {
  keep.sam <- which(colSums(otu.tab) >= lib.cut)
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(meta[keep.sam, ])
  names(Z) <- names(meta)
  rownames(Z) <- rownames(meta)[keep.sam]
  keep.tax <- which(rowSums(Y > 0) / ncol(Y) >= prev.cut)
  Y <- Y[keep.tax, ]
  
  N <- colSums(Y)
  keep.sam1 <- which(colSums(Y) >= 1)
  Y1 <- Y[, keep.sam1]
  Z1 <- as.data.frame(Z[keep.sam1, ])
  names(Z1) <- names(Z)
  rownames(Z1) <- rownames(Z)[keep.sam1]
  
  Y1 <- winsor.fun(Y1, winsor.quan) 
  
  return(list(Y = Y1, Z = Z1, keep.sam = keep.sam[keep.sam1], keep.tax = keep.tax))
}

filenames <- list.files("IBD")
df <- as.data.frame(array(" ", c(4, 4)))
names(df) <- c("m", "n", "n (controls)", "n (cases)")
for(i in 1 : 4) {
  load(paste0("IBD/", filenames[i]))
  res <- preprocess.fun(otu_table, metadata)
  Y <- res$Y
  Z <- res$Z
  df[i, 1] <- nrow(Y)
  df[i, 2] <- ncol(Y)
  if(i != 3) {
    tab <- table(Z$grp)
    nam <- names(tab)
    df[i, 3] <- tab[which(nam == 'H')]
    df[i, 4] <- tab[-which(nam == 'H')]
  } else if(i == 3) {
    tab <- table(Z$grp)
    nam <- names(tab)
    df[i, 3] <- tab[which(nam == 'nonIBD')]
    df[i, 4] <- tab[-which(nam == 'nonIBD')]
  }
}

datanames <- c("IBD-1", "IBD-2", "IBD-3", "IBD-4")
rownames(df) <- datanames
library(xtable)
xtable(df)

