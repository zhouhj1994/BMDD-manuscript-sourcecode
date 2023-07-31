library(ggplot2)
library(wesanderson) #wes_palette

#################################################
method.all <- c('BMDD', 'BMDD-a', 'BMDD-2', 'BMDD-2a',
                'mbDenoise', 'mbDenoise-2', 'mbImpute', 'mbImpute-2','pmr', 
                'SAVER', 'scImpute', 'ALRA', 
                'DMM', 'DMM-2', 'naive-1', 'naive-2', 'naive')
color.all <- c('#CC3333', '#FF6666', '#FF9999', '#FFCCCC',
               wes_palette('Darjeeling1', n = 5)[-1],wes_palette('Darjeeling2', n = 4),
               wes_palette('Chevalier1', n = 4), wes_palette('GrandBudapest2', n = 2)[2])
method.all <- method.all[1 : 8]
color.all <- color.all[1 : 8]

plot.fun <- function(val, sd, method, color) {
  val <- as.vector(t(val))
  sd <- as.vector(t(sd))
  
  model <- factor(1 : 5, labels = c('gamma', 'neg-binomial', 'Poisson', 
                                    'log-normal', 'log-multinormal'))
  metric <- factor(1 : 5, labels = c('SH', 'SP', 'HD', 'JS', 'BC'))
  s1 <- length(model); s2 <- length(method); s3 <- length(metric)
  method <- factor(1 : s2, labels = method)
  
  Model <- rep(model, each = s2 * s3)
  Method <- rep(rep(method, each = s3), s1)
  Metric <- rep(metric, s1 * s2)
  
  data <- cbind.data.frame(val, sd, Model, Method, Metric)
  
  plot1 <- ggplot(data, aes(x = Method, y = val, fill = Method)) + 
    geom_bar(stat = 'identity', position = "dodge") +
    #geom_errorbar(aes(ymin = val - 1.96 * sd, ymax = val + 1.96 * sd), width = .2) +
    scale_fill_manual(values = color) +
    labs(title = "", x = "", y = "") +
    #ggh4x::facet_grid2(Metric ~ Model, scales = "free_y", independent = "y") +
    ggh4x::facet_grid2(Metric ~ Model, scales = "free_y") +
    #facet_grid(Metric ~ Model) +
    theme_bw(base_size = 18) +
    theme(plot.margin = unit(c(1, 1, 1, 1.5), "cm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 
  return(plot1)
}

u <- c(1, 2, 4, 6, 8)
v <- c(1, 3, 5 : 8)
method <- method.all[v]
color <- color.all[v]
ind <- NULL
for(i in v) {
  ind <- c(ind, (i - 1) * 9 + u)
}

fun <- function(m, n, q, q1, p, p1, r, k, filename) {
  file <- paste0('result/output_cov_m_', m, '_n_', n, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data <- read.table(file)
  val <- data[seq(k, 10, 2), ind]
  sd <- data[seq(k, 10, 2), ind + ncol(data) / 2]
  val <- -log10(val)
  sd <- abs(-log10(exp(1)) / val) * sd
  
  pdf(filename, width = 16, height = 12)  
  print(plot.fun(val, sd, method, color))  
  dev.off()
}

fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0, k = 1, filename = 'S5_1_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0, k = 2, filename = 'S5_2_m_60_n_50.pdf')

#################################################
plot.fun.2 <- function(val, color) {
  val <- as.vector(t(val))
  
  model <- factor(1 : 5, labels = c('gamma', 'neg-binomial', 'Poisson', 
                                    'log-normal', 'log-multinormal'))
  metric <- factor(1 : 5, labels = c('SH', 'SP', 'HD', 'JS', 'BC'))
  setup <- factor(1 : 2, labels = c('S5.1', 'S5.2'))
  s1 <- length(model); s2 <- length(setup); s3 <- length(metric)
  
  Model <- rep(model, each = s2 * s3)
  Setup <- rep(rep(setup, each = s3), s1)
  Metric <- rep(metric, s1 * s2)
  
  data <- cbind.data.frame(val, Model, Setup, Metric)
  
  plot1 <- ggplot(data, aes(x = Setup, y = val, fill = Setup)) + 
    geom_bar(stat = 'identity', position = "dodge") +
    scale_fill_manual(values = color) + 
    labs(title = "", x = "", y = "") +
    ggh4x::facet_grid2(Metric ~ Model, scales = "free_y", independent = "y") +
    #ggh4x::facet_grid2(Metric ~ Model, scales = "free_y") +
    #facet_grid(Metric ~ Model) +
    theme_bw(base_size = 18) +
    theme(plot.margin = unit(c(1, 1, 1, 1.5), "cm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 
  return(plot1)
}

u <- c(1, 2, 4, 6, 8)

fun <- function(m, n, q, q1, p, p1, r) {
  file <- paste0('result/output_cov_m_', m, '_n_', n, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data <- read.table(file)
  val.bmdd <- data[, u]
  val.bmdd.2 <- data[, u + 18]
  r <- (val.bmdd - val.bmdd.2) / val.bmdd.2
  return(r)
}

r1 <- fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0)
r <- cbind(r1[seq(1, 10, 2), ], r1[seq(2, 10, 2), ])

pdf('BMDD_BMDD_cov.pdf', width = 16, height = 12)  
color <- wes_palette("Zissou1", 9, type = "continuous")
plot.fun.2(r, color[7:8])  
dev.off()




