otu.fun <- function(otu.tab, m, n, meta = NULL, otu.dist = NULL, otu.class = NULL) {
  keep.samp <- sort(order(colSums(otu.tab > 0), decreasing = TRUE)[1 : n])
  otu <- otu.tab[, keep.samp]
  keep.taxa <- sort(order(rowSums(otu > 0), decreasing = TRUE)[1 : m])
  otu <- otu[keep.taxa, ]
  if(!is.null(meta)) meta <- meta[keep.samp, ]
  if(!is.null(otu.dist)) otu.dist <- otu.dist[keep.taxa, keep.taxa]
  if(!is.null(otu.class)) otu.class <- otu.class[keep.taxa, ]
  return(list(otu = otu, meta = meta, otu.dist = otu.dist, otu.class = otu.class))
}

library(BMDD)

load('combo.RData')
meta <- cbind.data.frame(demo$sex, bmi.c)
otu.dist <- ape::cophenetic.phylo(tree)
otu.class <- otu.names.mat

m <- 102; n <- 80
data <- otu.fun(otu.tab, m, n, meta, otu.dist, otu.class)
U1 <- data$otu.class[, 1]
ind <- which(U1 %in% c('Other', 'Proteobacteria'))

otu <- data$otu[-ind, ]
Z1 <- as.vector(data$meta[, 1])
Z2 <- as.vector(scale(data$meta[, 2]))
U1 <- as.vector(data$otu.class[-ind, 1])
D <- data$otu.dist[-ind, -ind]
D.cen <- t(t(D) - colMeans(D))
pc1 <- eigen(t(D.cen) %*% D.cen)$vectors[, 1]
U2 <- as.vector(D.cen %*% pc1)

set.seed(666)
res.S1 <- bmdd(otu, type = 'count', trace = TRUE)
set.seed(666)
res.S2 <- bmdd(otu, type = 'count', trace = TRUE, 
               Z = data.frame(Z1), alp.eta = TRUE, pi.xi = TRUE, Z.standardizing = FALSE)
set.seed(666)
res.S3 <- bmdd(otu, type = 'count', trace = TRUE, 
               Z = data.frame(Z2), alp.eta = TRUE, pi.xi = TRUE, Z.standardizing = FALSE)
set.seed(666)
res.S4 <- bmdd(otu, type = 'count', trace = TRUE, 
               U = data.frame(U1), alp.kap = TRUE, pi.zeta = TRUE, U.standardizing = FALSE)
set.seed(666)
res.S5 <- bmdd(otu, type = 'count', trace = TRUE, 
               U = data.frame(U2), alp.kap = TRUE, pi.zeta = TRUE, U.standardizing = FALSE)

saveRDS(list(otu = otu, cova = list(Z1 = Z1, Z2 = Z2, U1 = U1, U2 = U2), 
             res = list(res.S1 = res.S1, res.S2 = res.S2, res.S3 = res.S3, 
                        res.S4 = res.S4, res.S5 = res.S5)), 'bmdd_fit_on_combo.rds')




