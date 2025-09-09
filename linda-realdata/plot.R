library(ggplot2)
library(ggpubr)
library(dplyr)

pval.mat.list <- readRDS("pval.mat.list.rds")
fdr.list <- readRDS('fdr.list.rds')

num.rej.fun <- function(pval.mat, cutoffs) {
  qval.mat <- sapply(1 : ncol(pval.mat), function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  isna <- colSums(is.na(qval.mat)) == nrow(pval.mat)
  
  num.rej <- t(sapply(1 : length(cutoffs), function(i)
    colSums(qval.mat <= cutoffs[i], na.rm = TRUE)))
  num.rej[, isna] <- NA
  return(num.rej)
}

datasets <- c("IBD-1", "IBD-2", "IBD-3", "IBD-4")
n.dat <- length(datasets)

#ind <- c(1, 3, 5)
ind <- c(2, 4, 6)
n.met <- length(ind)
#methods <- factor(1 : n.met, labels = c('LinDA', 'LinDA-BMDD', 'LinDA-SAVER'))
methods <- factor(1 : n.met, labels = c('ANCOMBC', 'ANCOMBC-BMDD', 'ANCOMBC-SAVER'))

cutoffs <- seq(0.01, 0.25, 0.01)
n.cut <- length(cutoffs)

num.rej.list <- list()
for(i in 1 : n.dat) {
  num.rej.list[[i]] <- num.rej.fun(pval.mat.list[[i]][, ind], cutoffs)
}

cutoff <- rep(cutoffs, n.met * n.dat)
Method <- rep(rep(methods, each = n.cut), n.dat)
dataset <- rep(datasets, each = n.cut * n.met)
num.rej <- as.vector(do.call(cbind, num.rej.list))

fdr.mat <- do.call(cbind, fdr.list)
#ind <- seq(1,24,2)
ind <- seq(2,24,2)
fdr.val <- as.vector(fdr.mat[, c(c(1:6),c(1:6)+12,c(1:6)+24,c(1:6)+36)][,ind])
fdr.sd <- as.vector(fdr.mat[, c(c(1:6),c(1:6)+12,c(1:6)+24,c(1:6)+36)+6][,ind])

data.plot <- cbind.data.frame(dataset=dataset, Method=Method,
                              cutoff=cutoff,
                              value1=fdr.val,
                              value2=num.rej)

fun <- function(setup, ylab = c("", "")) {
  p1 <- ggplot(data = data.plot %>% filter(dataset == setup), 
               aes(x = cutoff, y = value1, group = Method)) +
    geom_line(aes(color = Method), size = 0.8) +
    #geom_point(aes(color = Method)) +
    labs(title = setup, x = "", y = ylab[1]) +
    theme_minimal(base_size = 15) +
    scale_colour_manual(values = c('lightgrey', "lightcoral", "lightblue")) +
    geom_abline(slope = 1, intercept = 0, color = "gray", 
                linetype = "dashed", size = 0.8)
  
  p2 <- ggplot(data = data.plot %>% filter(dataset == setup), 
               aes(x = cutoff, y = value2, group = Method)) +
    geom_line(aes(color = Method), size = 0.8) +
    #geom_point(aes(color = Method)) +
    labs(title = "", x = "", y = ylab[2]) +
    theme_minimal(base_size = 15) +
    scale_colour_manual(values = c('lightgrey', "lightcoral", "lightblue")) 
  
  return(list(p1, p2))
}

p1 <- fun('IBD-1', ylab = c('Empirical False Discovery Rate', 'Number of Discoveries'))
p2 <- fun('IBD-2')
p3 <- fun('IBD-3')
p4 <- fun('IBD-4')

big_plot <- ggarrange(plotlist = list(p1[[1]], p2[[1]], p3[[1]], p4[[1]],
                          p1[[2]], p2[[2]], p3[[2]], p4[[2]]), 
          nrow = 2, ncol = 4, common.legend = TRUE)

final_plot <- annotate_figure(
  big_plot,
  bottom = text_grob("Target FDR Level", size = 15)
)

pdf("curve_ancombc.pdf", width = 14, height = 8)
final_plot
dev.off()

#################################
library(VennDiagram)
pval.mat.list <- readRDS("pval.mat.list.rds")

fun <- function(pval.mat, setup) {
  n.met <- ncol(pval.mat)
  qval.mat <- sapply(1 : n.met, function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  venn.list <- sapply(1 : n.met, function (i) 
    names((which(qval.mat[, i] <= cutoff))))
  
  venn.plot <- venn.diagram(
    x = venn.list,
    category.names = method,
    filename = NULL, 
    output = TRUE,
    fill = c( "lightgrey", "lightcoral", "lightblue"),
    alpha = 0.5,
    cat.col = c("lightgrey", "lightcoral", "lightblue"),
    cat.cex = 1,
    cat.fontface = "bold",
    main = setup,
    main.cex = 1.5
  )
  pdf(paste0(setup, "_venn.pdf"), width = 8, height = 6)
  grid.draw(venn.plot)
  dev.off()
}

cutoff <- 0.1
ind <- c(1, 3, 5)
#ind <- c(2, 4, 6)
method <- c('LinDA', 'LinDA-BMDD', 'LinDA-SAVER')
#method <- c('ANCOMBC', 'ANCOMBC-BMDD', 'ANCOMBC-SAVER')

fun(pval.mat.list[[4]][, ind], 'IBD-4')



library(gridExtra)

fun2 <- function(pval.mat, setup) {
  n.met <- ncol(pval.mat)
  qval.mat <- sapply(1 : n.met, function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  venn.list <- sapply(1 : n.met, function (i) 
    names((which(qval.mat[, i] <= cutoff))))
  
  venn.plot <- venn.diagram(
    x = venn.list,
    category.names = method,
    filename = NULL, 
    output = TRUE,
    fill = c( "lightgrey", "lightcoral", "lightblue"),
    alpha = 0.5,
    cat.col = c("black", "red", "blue"),
    cat.cex = 1,
    cat.fontface = "bold",
    main = setup,
    main.cex = 1.5,
    cat.dist = c(0.03, -0.02, 0.03)
  )
  return(venn.plot)
}

p1 <- fun2(pval.mat.list[[1]][, ind], 'IBD-1')
p2 <- fun2(pval.mat.list[[2]][, ind], 'IBD-2')
p3 <- fun2(pval.mat.list[[3]][, ind], 'IBD-3')
p4 <- fun2(pval.mat.list[[4]][, ind], 'IBD-4')


pdf("venn.pdf", width = 8, height = 6)
grid.arrange(grobs = list(p1, p2, p3, p4), nrow = 2, ncol = 2,
             padding = unit(c(2, 3, 3, 3), "cm"))
dev.off()




