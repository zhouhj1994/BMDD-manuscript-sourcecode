bmdd_fit_on_combo <- readRDS("bmdd_fit_on_combo.rds")
fit <- bmdd_fit_on_combo$res[[1]]
pi <- fit$pi
alp0 <- fit$alpha$alp0
alp1 <- fit$alpha$alp1

ind <- which(alp0 > alp1)
tmp <- alp0[ind]
alp0[ind] <- alp1[ind]
alp1[ind] <- tmp
pi[ind] <- 1 - pi[ind]

m <- 100
n <- 80


########################################################################
pdf("bmdd_fit_on_combo_plot.pdf", width = 11, height = 8)

df <- data.frame(alp0, alp1, pi)

library(ggplot2)
library(grid)

# Calculate statistics directly from df
means <- sapply(df, mean)
medians <- sapply(df, median)

# Base ggplot
p <- ggplot(df, aes(x = alp0, y = alp1, color = pi, size = pi)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_c() +
  labs(
    x = expression(alpha[0]),
    y = expression(alpha[1]),
    color = expression(pi),
    size = expression(pi)
  ) +
  theme_minimal() +
  theme(
    axis.title      = element_text(size = 18),  # axis label
    axis.text       = element_text(size = 14),  # tick labels
    legend.title    = element_text(size = 16),  # legend title
    legend.text     = element_text(size = 14),  # legend labels
    plot.margin     = margin(15, 15, 15, 15)   # margin around plot
  )

# Draw the plot
print(p)

# Add multi-line annotation using grid.text() with math symbols
grid.text(
  label = bquote(
    Mean(alpha[0]) == .(round(means["alp0"],2)) ~ "\n" ~
      Median(alpha[0]) == .(round(medians["alp0"],2)) ~ "\n" ~
      Mean(alpha[1]) == .(round(means["alp1"],2)) ~ "\n" ~
      Median(alpha[1]) == .(round(medians["alp1"],2)) ~ "\n" ~
      Mean(pi) == .(round(means["pi"],2)) ~ "\n" ~
      Median(pi) == .(round(medians["pi"],2))
  ),
  x = 0.95, y = 0.95, just = c("right","top"), gp = gpar(cex = 1.2)
)

dev.off()

############################################################

library(BMDD)

data.generation <- function(alp0, alp1, pi, m, n) {
  del <- matrix(rbinom(m * n, 1, pi), m)
  beta <- alp0 * (1 - del) + alp1 * del
  
  X <- matrix(rgamma(m * n, beta, 1), m)
  X <- t(t(X) / colSums(X))
  return(X)
}

set.seed(666)

X <- data.generation(alp0, alp1, pi, m, n)

res.BMDD.1 <- bmdd(X, type = 'proportion', trace = TRUE)

alp0.est <- res.BMDD.1$alpha$alp0
alp1.est <- res.BMDD.1$alpha$alp1
pi.est <- res.BMDD.1$pi

ind <- which(alp0.est > alp1.est)
tmp <- alp0.est[ind]
alp0.est[ind] <- alp1.est[ind]
alp1.est[ind] <- tmp
pi.est[ind] <- 1 - pi.est[ind]

###########################################################

set.seed(666)

X <- data.generation(alp0, alp1, pi=rep(0.5, m), m, n)

res.BMDD.2 <- bmdd(X, type = 'proportion', trace = TRUE)

alp0.est <- res.BMDD.2$alpha$alp0
alp1.est <- res.BMDD.2$alpha$alp1
pi.est <- res.BMDD.2$pi

ind <- which(alp0.est > alp1.est)
tmp <- alp0.est[ind]
alp0.est[ind] <- alp1.est[ind]
alp1.est[ind] <- tmp
pi.est[ind] <- 1 - pi.est[ind]

###########################################################

df <- data.frame(alp0, alp1, pi, alp0.est, alp1.est, pi.est)

library(ggplot2)
library(dplyr)
library(patchwork)

# Highlight points where |alpha1 - alpha0| < 0.5
df <- df %>%
  mutate(highlight = ifelse(abs(alp1 - alp0) < 0.5,
                            "|alpha0 - alpha1| < 0.5",
                            "|alpha0 - alpha1| > 0.5"))

# Define common colors
colors <- c("|alpha0 - alpha1| < 0.5" = "red",
            "|alpha0 - alpha1| > 0.5" = "grey")

# Text sizes
tick_size <- 16
axis_size <- 18
title_size <- 20
legend_size <- 16

# --- Plot 1 ---
p1 <- ggplot(df, aes(x = alp0, y = alp0.est, color = highlight)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = colors) +
  labs(title = expression(alpha[0]), x = "True value", y = "Estimated value", color = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = tick_size),
    axis.title = element_text(size = axis_size),
    plot.title = element_text(size = title_size, hjust = 0.5),
    legend.text = element_text(size = legend_size),
    legend.title = element_blank()
  )

# --- Plot 2 ---
p2 <- ggplot(df, aes(x = alp1, y = alp1.est, color = highlight)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = colors) +
  labs(title = expression(alpha[1]), x = "True value", y = "Estimated value", color = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = tick_size),
    axis.title = element_text(size = axis_size),
    plot.title = element_text(size = title_size, hjust = 0.5),
    legend.text = element_text(size = legend_size),
    legend.title = element_blank()
  )

# --- Plot 3 ---
p3 <- ggplot(df, aes(x = seq_along(pi.est), y = pi.est, color = highlight)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_color_manual(values = colors) +
  labs(title = expression(pi), x = "Index", y = "Estimated value", color = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = tick_size),
    axis.title = element_text(size = axis_size),
    plot.title = element_text(size = title_size, hjust = 0.5),
    legend.text = element_text(size = legend_size),
    legend.title = element_blank()
  )

# --- Save to PDF with shared legend ---
pdf("check_the_estimates.pdf", width = 22, height = 8)
(p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.title = element_blank())
dev.off()

