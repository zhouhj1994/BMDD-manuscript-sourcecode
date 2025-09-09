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

########################################################################
pval.mat.list <- readRDS("pval.mat.list.rds")
filenames <- list.files("IBD")

cutoff <- 0.1
ind <- c(1, 3, 5)
#ind <- c(2, 4, 6)
#ind <- 1 : 6
method <- c('LinDA', 'LinDA-BMDD', 'LinDA-SAVER')
#method <- c('ANCOMBC', 'ANCOMBC-BMDD', 'ANCOMBC-SAVER')
#method <- c('LinDA', 'ANCOMBC', 'LinDA-BMDD', 'ANCOMBC-BMDD', 'LinDA-SAVER', 'ANCOMBC-SAVER')

j <- 4

fun1 <- function(j) {
  pval.mat <- pval.mat.list[[j]][, ind]
  n.met <- ncol(pval.mat)
  qval.mat <- sapply(1 : n.met, function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  venn.list <- sapply(1 : n.met, function (i) 
    which(qval.mat[, i] <= cutoff))
  otu.all <- unique(names(unlist(venn.list)))
  upset.df <- data.frame(Name = otu.all)
  for(i in 1 : n.met) {
    upset.df[, i + 1] <- otu.all %in% names(venn.list[[i]]) + 0
  }
  names(upset.df)[2 : (n.met + 1)] <- method
  
  load(paste0("IBD/", filenames[j]))
  res <- preprocess.fun(otu_table, metadata)
  
  rej <- which(upset.df[, 2] == 1)
  length(rej)
  nam <- upset.df[rej, 1]
  rej1 <- tax_table[which(tax_table$otu %in% nam), ]
  table(rej1$phylum)
  table(rej1$class)
  table(rej1$family)
  
  rej <- which(upset.df[, 3] == 1)
  length(rej)
  nam <- upset.df[rej, 1]
  rej2 <- tax_table[which(tax_table$otu %in% nam), ]
  table(rej2$phylum)
  table(rej2$class)
  table(rej2$family)
  
  rej <- which(upset.df[, 4] == 1)
  length(rej)
  nam <- upset.df[rej, 1]
  rej3 <- tax_table[which(tax_table$otu %in% nam), ]
  table(rej3$phylum)
  table(rej3$class)
  table(rej3$family)
  
  rej <- which(upset.df[, 2] == 0 & upset.df[, 3] == 1 & upset.df[, 4] == 0)
  nam <- upset.df[rej, 1]
  tax_table[which(tax_table$otu %in% nam), ]
}










