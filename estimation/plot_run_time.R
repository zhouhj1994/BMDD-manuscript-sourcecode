library(ggplot2)
library(ggpubr)

run.time.plot <- function(s, m.vec, n.vec) {
  res <- read.table(paste0('output/bmdd_run_time/time', s, '.txt'))
  est <- colMeans(res)
  nsim <- nrow(res)
  est.sd <- sqrt(rowSums((t(res) - est) ^ 2) / ((nsim - 1) * nsim))
  
  lm <- length(m.vec)
  ln <- length(n.vec)
  n.vec <- factor(1 : ln, labels = n.vec)
  data <- cbind.data.frame(m = rep(m.vec, each = ln), n = rep(n.vec, lm),
                           value = est, sd = est.sd)
  ggplot(data, aes(x = m, y = est, group = n)) +
    geom_line(aes(color = n), linewidth = 1) +
    geom_point(aes(color = n, shape = n), size = 3) +
    theme_bw(base_size = 25) +
    xlab("m") +
    ylab("seconds per iteration") +
    theme(legend.key.width = unit(1, "cm"),
          plot.margin=unit(c(1, 1, 1, 1.5),"cm"))
}

s <- 1
m.vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000)
n.vec <- c(80, 160, 320)
plot1 <- run.time.plot(s, m.vec, n.vec)

pdf("bmdd_run_time_nocov.pdf",width = 11,height = 8)
print(plot1)
dev.off()

s <- 2
m.vec <- c(100, 200, 300, 400, 500)
n.vec <- c(80, 160, 320)
plot2 <- run.time.plot(s, m.vec, n.vec)

s <- 4
plot3 <- run.time.plot(s, m.vec, n.vec)

pdf("bmdd_run_time.pdf",width = 20,height = 6)
ggarrange(plotlist = list(plot1, plot2, plot3), ncol = 3, common.legend = TRUE, legend = 'right')
dev.off()
