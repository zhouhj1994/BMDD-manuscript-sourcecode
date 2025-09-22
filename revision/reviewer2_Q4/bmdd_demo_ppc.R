library(ggplot2)
library(ggpubr)
library(dplyr)

cnt <- readRDS('res_bmdd_Daniel2018.rds')[[1]]
res.bmdd <- readRDS('res_bmdd_Daniel2018.rds')[[2]]

beta <- res.bmdd$beta

comp <- t(t(cnt) / colSums(cnt))

m <- nrow(cnt)
n <- ncol(cnt)
N <- colSums(cnt)

set.seed(123)
comp.post <- matrix(rgamma(m * n, beta, 1), m)
comp.post <- t(t(comp.post) / colSums(comp.post))
cont.post <- sapply(1 : n, function(i)rmultinom(1, N[i], comp.post[, i]))
comp.post <- t(t(cont.post) / colSums(cont.post))


ind <- c(24, 48, 141, 173)
y.max <- c(30, 250, 20000, 15000)


plot.fun <- function(k) {
  j <- ind[k]
  
  # Subset and truncate extreme values
  x1 <- comp[j, ]
  x1 <- x1[x1 < quantile(x1, 0.97)]
  x2 <- comp.post[j, ]
  x2 <- x2[x2 < quantile(x2, 0.97)]
  
  # Combine into a data frame
  df <- data.frame(
    value = c(x1, x2),
    group = c(rep("Observed", length(x1)), rep("Generated", length(x2)))
  )
  
  # x-axis range and span
  x_range <- range(df$value)
  x_span <- diff(x_range)
  
  # Compute summary stats
  stats <- df %>%
    group_by(group) %>%
    summarise(
      perc_zero = mean(value == 0) * 100,
      mean_val  = mean(value),
      sd_val    = sd(value),
      .groups = "drop"
    )
  
  # Minimal separation
  min_sep <- 0.18 * x_span
  
  # Always put group 1 left, group 2 right
  sep <- max(min_sep, diff(stats$mean_val) / 2)
  
  # Rightward shift to avoid left margin
  right_shift <- 0.24 * x_span
  
  stats <- stats %>%
    mutate(
      mean_str = formatC(mean_val, digits = 3, format = "g"),
      sd_str   = formatC(sd_val, digits = 3, format = "g"),
      label = paste0("0s = ", round(perc_zero, 1), "%\n",
                     "Mean = ", mean_str, "\n",
                     "SD = ", sd_str),
      # group 1 (Observed) on left, group 2 (Generated) on right
      x = c(stats$mean_val[stats$group == "Observed"] - sep,
            stats$mean_val[stats$group == "Generated"] + sep) + right_shift,
      y = y.max[k] * 0.9
    )
  
  # Plot
  gp <- ggplot(df, aes(x = value, fill = group)) +
    geom_histogram(aes(y = ..density..), 
                   alpha = 0.5, position = "identity", bins = 30, color = "black") +
    scale_y_continuous(limits = c(0, y.max[k])) +
    theme_minimal(base_size = 20) +
    labs(x = "", y = "") +
    theme(legend.title = element_blank()) +
    geom_text(
      data = stats,
      aes(x = x, y = y, label = label, color = group),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 1,
      size = 4,
      show.legend = FALSE
    )
  
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

pdf('ppc_hist.pdf', width = 12, height = 8)
final_plot
dev.off()

