para.fun <- function(otu.tab, m) {
  otu.tab <- as.matrix(otu.tab)
  colnames(otu.tab) = NULL
  rownames(otu.tab) = NULL
  
  has.read <- rowSums(otu.tab > 0)
  ind.taxa <- sort(order(has.read, decreasing = TRUE)[1 : m])
  otu.tab.sel <- otu.tab[ind.taxa, ]
  
  lambda0 <- rep(NA, m)
  N <- colSums(otu.tab.sel)
  for(i in 1 : m) {
    fit <- glm(otu.tab.sel[i,] ~ N - 1, family = poisson(link = log))
    lambda0[i] <- fit$coefficients
  }
  out <- list(otu.tab.sel = otu.tab.sel, lambda0 = lambda0)
  return(out)
}

load("combo.RData")

m <- 500
res <- para.fun(otu.tab, m) 
Poisson.para <- list(lambda0 = res[[2]])
saveRDS(Poisson.para, "Poisson.para.rds")