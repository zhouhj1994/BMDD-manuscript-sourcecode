library(ggplot2)
library(ggpubr)

run.time.plot <- function(s, m.vec, n.vec) {
  res <- read.table(paste0('output/saver_run_time_2/time', s, '.txt'))[, -c(1 : 3)]
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
    ylab("seconds") +
    theme(legend.key.width = unit(1, "cm"),
          plot.margin=unit(c(1, 1, 1, 1.5),"cm"))
}

s <- 1
m.vec <- c(500, 1000, 1500, 2000, 2500, 3000)
n.vec <- c(80, 160, 320)
plot1 <- run.time.plot(s, m.vec, n.vec)

pdf("saver_run_time.pdf",width = 11,height = 8)
print(plot1)
dev.off()


