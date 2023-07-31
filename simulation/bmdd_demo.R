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
# saveRDS(list(W, res.bmdd, res.dmn),'res_bmdd_Daniel2018.rds')
# W <- readRDS('res_bmdd_Daniel2018.rds')[[1]]
# res.bmdd <- readRDS('res_bmdd_Daniel2018.rds')[[2]]
# res.dmn <- readRDS('res_bmdd_Daniel2018.rds')[[3]]

theta <- as.vector(fitted(res.dmn))
theta.sum <- sum(theta)

beta <- res.bmdd$beta
Xhat <- array(0, c(m, n))
for(j in 1 : 20) {
  Xhat.j <- matrix(rgamma(m * n, beta, 1), m)
  Xhat.j <- t(t(Xhat.j) / colSums(Xhat.j))
  Xhat <- Xhat + Xhat.j / 20
}
X <- t(t(W) / colSums(W))

i <- 10
ind <- ((i - 1) * 20 + 1) : (i * 20)
ind <- c(c(1, 3, 11, 15), c(4, 6, 11, 12) + 20, c(2, 5, 6, 8, 12) + 40,
         c(2, 7, 17, 19, 20) + 60, c(3, 4, 11, 12, 15) + 80,
         c(2, 5, 6, 11) + 100, c(1, 2, 7, 11, 15) + 120,
         c(1, 13, 15, 18) + 140, c(13, 16) + 160, c(2, 14) + 180)
par(mfrow = c(4, 5))

ind <- c(24, 105, 141, 173) #141, 15, 62, 194
pdf('fit2.pdf', width = 11, height = 8)
ind <- c(3, 48, 141, 173)
pdf('fit1.pdf', width = 11, height = 8)
ind <- c(24, 48, 141, 173)
pdf('fit.pdf', width = 11, height = 8)

par(mfrow = c(2, 2), cex = 1.2)
for(j in ind) {
  x <- X[j, ]
  x <- x[x < quantile(x, 0.97)]
  h <- hist(x, prob = TRUE, breaks = 20, main = '', xlab = '', ylab = '', 
            col = 'steelblue')
  xx <- seq(min(x), max(x), len = 1000)
  yy <- dbeta(xx, theta[j], theta.sum - theta[j])
  lines(xx, yy, lwd = 2, col = 'blue')
  x <- Xhat[j, ]
  x <- x[x < quantile(x, 0.97)]
  adj <- 0.5
  if(j == 48 | j == 3) adj <- 1
  lines(density(x, from = 0, adjust = adj), lwd = 2, col = 'red')
  if(j == ind[4])
    legend('topright', legend = c('BMDD', 'Dirichlet distribution'),
           col = c('red', 'blue'), lty = c(1, 1))
}
dev.off()






