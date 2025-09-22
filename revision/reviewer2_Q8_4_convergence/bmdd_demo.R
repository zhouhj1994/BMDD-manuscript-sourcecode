source("bmdd4.R")
load('Daniel2018_AGP_UK.RData')
W0 <- feature.dat

m <- 200
n <- 1000

keep.samp <- sort(order(colSums(W0 > 0), decreasing = TRUE)[1:n])
W <- W0[,keep.samp]
keep.taxa <- sort(order(rowSums(W0 > 0), decreasing = TRUE)[1:m])
W <- W[keep.taxa,]

res.bmmd.Daniel2018 <- list()

set.seed(1)
res.bmmd.Daniel2018[[1]] <- bmdd4(W, type = 'count', trace = TRUE)

set.seed(11)
res.bmmd.Daniel2018[[2]] <- bmdd4(W, type = 'count', trace = TRUE)

saveRDS(res.bmmd.Daniel2018, 'res_bmmd_Daniel2018_list.rds')



res.bmdd <- res.bmmd.Daniel2018[[1]]
B <- length(res.bmdd)
beta.ave.mat.1 <- array(NA, c(B, m))
for(i in 1 : B) {
  beta <- res.bmdd[[i]]$beta
  beta <- t(t(beta) / colSums(beta))
  beta.ave.mat.1[i, ] <- rowMeans(beta)
}

res.bmdd <- res.bmmd.Daniel2018[[2]]
B <- length(res.bmdd)
beta.ave.mat.2 <- array(NA, c(B, m))
for(i in 1 : B) {
  beta <- res.bmdd[[i]]$beta
  beta <- t(t(beta) / colSums(beta))
  beta.ave.mat.2[i, ] <- rowMeans(beta)
}


library(ggplot2)
library(reshape2)


ind <- c(24, 48, 141, 173)

df1 <- as.data.frame(beta.ave.mat.1)[, ind]
df1$Iteration <- 1:nrow(beta.ave.mat.1)
df1$Source <- "Initialization 1"

df2 <- as.data.frame(beta.ave.mat.2)[, ind]
df2$Iteration <- 1:nrow(beta.ave.mat.2)
df2$Source <- "Initialization 2"

# Combine
df_all <- rbind(df1, df2)

# Melt to long format
df_long <- melt(df_all, id.vars = c("Iteration", "Source"))

# Plot
gp <- ggplot(df_long, aes(x = Iteration, y = value, color = Source)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal(base_size = 16) +   
  theme(
    axis.title = element_text(size = 18),     
    axis.text = element_text(size = 14),      
    strip.text = element_blank(),             
    legend.text = element_text(size = 14),
    legend.position = "top",
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(x = "Iteration", y = "Average posterior mean", color = NULL)

pdf("convergence.pdf", width = 11, height = 8)
gp
dev.off()



