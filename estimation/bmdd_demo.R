library(BMDD)
library(DirichletMultinomial)

load('Daniel2018_AGP_UK.RData')
W0 <- feature.dat

m <- 200
n <- 1000

set.seed(666)
keep.samp <- sort(order(colSums(W0 > 0), decreasing = TRUE)[1:n])
W <- W0[,keep.samp]
keep.taxa <- sort(order(rowSums(W0 > 0), decreasing = TRUE)[1:m])
W <- W[keep.taxa,]

res.dmn <- dmn(t(W), k = 1)
res.bmdd <- bmdd(W, type = 'count', trace = TRUE)
saveRDS(list(W, res.bmdd, res.dmn),'res_bmdd_Daniel2018.rds')



################################################
library(ggplot2)
library(ggpubr)

W <- readRDS('res_bmdd_Daniel2018.rds')[[1]]
res.bmdd <- readRDS('res_bmdd_Daniel2018.rds')[[2]]
res.dmn <- readRDS('res_bmdd_Daniel2018.rds')[[3]]

theta <- as.vector(fitted(res.dmn))
theta.sum <- sum(theta)

beta <- res.bmdd$beta
Xhat <- t(t(beta) / colSums(beta))
X <- t(t(W) / colSums(W))

ind <- c(24, 48, 141, 173)

#i <- 1
#ind <- (1 : 4) + 4 * (i - 1)

plot.fun <- function(k) {
  j <- ind[k]
  x <- X[j, ]
  x <- x[x < quantile(x, 0.97)]
  
  x1 <- seq(min(x), max(x), len = 1000)
  y1 <- dbeta(x1, theta[j], theta.sum - theta[j])
  
  x2 <- Xhat[j, ]
  x2 <- x2[x2 < quantile(x2, 0.97)]
  
  adj <- c(0.5, 1, 0.5, 0.5)
  y.max <- c(30, 250, 20000, 15000)
  
  gp <- ggplot() +
    geom_histogram(aes(x = x, y = ..density..), #boundary = 0,
                   fill = 'steelblue', color = 'black') +
    xlab('') + ylab('') +
    theme_minimal(base_size = 25) +
    geom_line(aes(x = x1, y = y1, color = "Dirichlet distribution"), size = 1.2) +
    scale_y_continuous(limits = c(0, y.max[k])) + 
    stat_density(aes(x = x2, color = "BMDD"), size = 1.2, adjust = adj[k],
                 geom = "line", position = "identity") + 
    scale_color_manual(
      name = "",
      values = c("Dirichlet distribution" = "blue", "BMDD" = "red")
    ) +
    theme(legend.position = "top")
  
  if(k == 3) {
    gp <- gp +
      scale_x_continuous(
        breaks = c(0, 0.00025, 0.0005, 0.00075),  
        labels = c("0e+00", "2.5e−04", "5e−04", "7.5e−04")  
      ) 
  }
  
  return(gp)
}

p1 <- plot.fun(1)
p2 <- plot.fun(2)

p3 <- plot.fun(3)
p4 <- plot.fun(4)

big_plot <- ggarrange(plotlist = list(p1, p2, p3, p4), 
                      nrow = 2, ncol = 2, common.legend = TRUE)

final_plot <- annotate_figure(
  big_plot,
  bottom = text_grob("Proportion", size = 25),
  left = text_grob("Density", size = 25, rot = 90)
)

pdf('fit.pdf', width = 12, height = 8)
final_plot
dev.off()

