library(ggplot2)
library(wesanderson) #wes_palette

#################################################
method.all <- c('BMDD', 'BMDD-2', 'BMDD-3', 'BMDD-4',
                'mbDenoise', 'mbDenoise-2', 'mbImpute', 'mbImpute-2','pmr', 
                'SAVER', 'scImpute', 'ALRA', 
                'DMM', 'DMM-2', 'naive-1', 'naive-2', 'naive')
color.all <- c('#CC3333', '#FF6666', '#FF9999', '#FFCCCC',
               wes_palette('Darjeeling1', n = 5)[-1],wes_palette('Darjeeling2', n = 4),
               wes_palette('Chevalier1', n = 4), wes_palette('GrandBudapest2', n = 2)[2])
method.all <- method.all[-c(3, 4, 6, 8)]
color.all <- color.all[-c(3, 4, 6, 8)]

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
v <- c(1, 3, 4, 5, 6, 7, 8, 9)
method <- method.all[v]
color <- color.all[v]
ind <- NULL
for(i in v) {
  ind <- c(ind, (i - 1) * 9 + u)
}

fun <- function(m, n, q, q1, p, p1, r, filename) {
  file <- paste0('result/output_m_', m, '_n_', n, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data <- read.table(file)
  val <- -log10(data[, ind])
  sd <- abs(-log10(exp(1)) / data[, ind]) * data[, ind + ncol(data) / 2]
  
  pdf(filename, width = 16, height = 12)  
  print(plot.fun(val, sd, method, color))  
  dev.off()
}

fun(m = 60, n = 50, q = 0, q1 = 0, p = 1, p1 = 0.5, r = 0, filename = 'S0_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0, q1 = 0, p = 0, p1 = 0, r = 0, filename = 'S1_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0, p1 = 0, r = 0, filename = 'S2_1_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0.4, q1 = 0.3, p = 0, p1 = 0, r = 0, filename = 'S2_2_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0.6, q1 = 0.3, p = 0, p1 = 0, r = 0, filename = 'S2_3_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0, filename = 'S3_m_60_n_50.pdf')
fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0.1, filename = 'S4_m_60_n_50.pdf')

#################################################
plot.fun.2 <- function(val, color) {
  val <- as.vector(t(val))
  
  model <- factor(1 : 5, labels = c('gamma', 'neg-binomial', 'Poisson', 
                                    'log-normal', 'log-multinormal'))
  metric <- factor(1 : 5, labels = c('SH', 'SP', 'HD', 'JS', 'BC'))
  setup <- factor(1 : 6, labels = c('S1', 'S2.1', 'S2.2', 'S2.3', 'S3', 'S4'))
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

###
fun <- function(m, n, q, q1, p, p1, r) {
  file <- paste0('result/output_m_', m, '_n_', n, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data <- read.table(file)
  val.bmdd <- data[, u]
  val.saver <- data[, u + 45]
  r <- (val.saver - val.bmdd) / val.bmdd
  return(r)
}
 
r1 <- fun(m = 60, n = 50, q = 0, q1 = 0, p = 0, p1 = 0, r = 0)
r2 <- fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0, p1 = 0, r = 0)
r3 <- fun(m = 60, n = 50, q = 0.4, q1 = 0.3, p = 0, p1 = 0, r = 0)
r4 <- fun(m = 60, n = 50, q = 0.6, q1 = 0.3, p = 0, p1 = 0, r = 0)
r5 <- fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0)
r6 <- fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0.1)
r <- cbind(r1, r2, r3, r4, r5, r6)

pdf('BMDD_SAVER_50.pdf', width = 16, height = 12)  
color <- wes_palette("Zissou1", 9, type = "continuous")
plot.fun.2(r, color[1:6])  
dev.off()

###
fun <- function(m, n, q, q1, p, p1, r) {
  file <- paste0('result/output_bmdd_saver_m_', m, '_n_', n, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data <- read.table(file)
  val.bmdd <- data[, u]
  val.saver <- data[, u + 18]
  r <- (val.saver - val.bmdd) / val.bmdd
  return(r)
}

r1 <- fun(m = 60, n = 500, q = 0, q1 = 0, p = 0, p1 = 0, r = 0)
r2 <- fun(m = 60, n = 500, q = 0.2, q1 = 0.3, p = 0, p1 = 0, r = 0)
r3 <- fun(m = 60, n = 500, q = 0.4, q1 = 0.3, p = 0, p1 = 0, r = 0)
r4 <- fun(m = 60, n = 500, q = 0.6, q1 = 0.3, p = 0, p1 = 0, r = 0)
r5 <- fun(m = 60, n = 500, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0)
r6 <- fun(m = 60, n = 500, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0.1)
r <- cbind(r1, r2, r3, r4, r5, r6)

pdf('BMDD_SAVER_500.pdf', width = 16, height = 12)  
color <- wes_palette("Zissou1", 9, type = "continuous")
plot.fun.2(r, color[1:6])  
dev.off()

###
fun <- function(m, n1, n2, q, q1, p, p1, r) {
  file1 <- paste0('result/output_m_', m, '_n_', n1, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data1 <- read.table(file1)
  val.bmdd1 <- data1[, u]
  
  file2 <- paste0('result/output_bmdd_saver_m_', m, '_n_', n2, '_q_', q, 
                  '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data2 <- read.table(file2)
  val.bmdd2 <- data2[, u]
  
  r <- (val.bmdd1 - val.bmdd2) / val.bmdd2
  return(r)
}

r1 <- fun(m = 60, n1 = 50, n2 = 500, q = 0, q1 = 0, p = 0, p1 = 0, r = 0)
r2 <- fun(m = 60, n1 = 50, n2 = 500, q = 0.2, q1 = 0.3, p = 0, p1 = 0, r = 0)
r3 <- fun(m = 60, n1 = 50, n2 = 500, q = 0.4, q1 = 0.3, p = 0, p1 = 0, r = 0)
r4 <- fun(m = 60, n1 = 50, n2 = 500, q = 0.6, q1 = 0.3, p = 0, p1 = 0, r = 0)
r5 <- fun(m = 60, n1 = 50, n2 = 500, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0)
r6 <- fun(m = 60, n1 = 50, n2 = 500, q = 0.2, q1 = 0.3, p = 0.3, p1 = 0.5, r = 0.1)
r <- cbind(r1, r2, r3, r4, r5, r6)

pdf('BMDD_BMDD.pdf', width = 16, height = 12)  
color <- wes_palette("Zissou1", 9, type = "continuous")
plot.fun.2(r, color[1:6])  
dev.off()

#################################################
plot.fun.3 <- function(val, color) {
  val <- as.vector(t(val))
  
  model <- factor(1 : 5, labels = c('gamma', 'neg-binomial', 'Poisson', 
                                    'log-normal', 'log-multinormal'))
  metric <- factor(1 : 5, labels = c('SH', 'SP', 'HD', 'JS', 'BC'))
  setup <- factor(1 : 3, labels = c('S2.1', 'S2.2', 'S2.3'))
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
file0 <- paste0('result/output_m_', 60, '_n_', 50, '_q_', 0, 
                '_q1_', 0, '_p_', 0, '_p1_', 0, '_r_', 0, ".txt")
data0 <- read.table(file0)

fun <- function(m, n, q, q1, p, p1, r) {
  file <- paste0('result/output_m_', m, '_n_', n, '_q_', q, 
                 '_q1_', q1, '_p_', p, '_p1_', p1, '_r_', r, ".txt")
  data <- read.table(file)
  
  val.bmdd <- data[, u]
  val.bmdd0 <- data0[, u]
  r <- (val.bmdd - val.bmdd0) / val.bmdd0
  return(r)
}

r2 <- fun(m = 60, n = 50, q = 0.2, q1 = 0.3, p = 0, p1 = 0, r = 0)
r3 <- fun(m = 60, n = 50, q = 0.4, q1 = 0.3, p = 0, p1 = 0, r = 0)
r4 <- fun(m = 60, n = 50, q = 0.6, q1 = 0.3, p = 0, p1 = 0, r = 0)
r <- cbind(r2, r3, r4)

pdf('BMDD_S2.pdf', width = 16, height = 12)  
color <- wes_palette("Zissou1", 9, type = "continuous")
plot.fun.3(r, color[2 : 4])  
dev.off()






