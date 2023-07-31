library(ggplot2)
library(wesanderson) #wes_palette

#################################################
plot.fun <- function(val, sd, method, color) {
  val <- as.vector(t(val))
  sd <- as.vector(t(sd))
  
  model <- factor(1 : 6, labels = c('S6 (stool)', 'S6 (vaginal)', 'S7 (stool)', 
                                    'S7 (vaginal)', 'S8 (stool)', 'S8 (vaginal)'))
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

file.fun <- function(data, m, n, q, q1, r) {
  paste0('result/output_nonpara_', data, 
         '_m_', m, '_n_', n, '_q_', q, '_q1_', q1, '_r_', r, ".txt")
}

out.stool.1 <- as.vector(unlist(read.table(file.fun(
  data = 'stool', m = 60, n = 50, q = 0, q1 = 0, r = 0))))
out.stool.2 <- as.vector(unlist(read.table(file.fun(
  data = 'stool', m = 60, n = 50, q = 0.2, q1 = 0.3, r = 0))))
out.stool.3 <- as.vector(unlist(read.table(file.fun(
  data = 'stool', m = 60, n = 50, q = 0.2, q1 = 0.3, r = 0.1))))

out.vaginal.1 <- as.vector(unlist(read.table(file.fun(
  data = 'vaginal', m = 60, n = 50, q = 0, q1 = 0, r = 0))))
out.vaginal.2 <- as.vector(unlist(read.table(file.fun(
  data = 'vaginal', m = 60, n = 50, q = 0.2, q1 = 0.3, r = 0))))
out.vaginal.3 <- as.vector(unlist(read.table(file.fun(
  data = 'vaginal', m = 60, n = 50, q = 0.2, q1 = 0.3, r = 0.1))))

method.all <- c('BMDD', 'BMDD-2', 'BMDD-3', 'BMDD-4',
                'mbDenoise', 'mbDenoise-2', 'mbImpute', 'mbImpute-2','pmr', 
                'SAVER', 'scImpute', 'ALRA', 
                'DMM', 'DMM-2', 'naive-1', 'naive-2', 'naive')
color.all <- c('#CC3333', '#FF6666', '#FF9999', '#FFCCCC',
               wes_palette('Darjeeling1', n = 5)[-1],wes_palette('Darjeeling2', n = 4),
               wes_palette('Chevalier1', n = 4), wes_palette('GrandBudapest2', n = 2)[2])
method.all <- method.all[-c(3, 4, 6, 8)]
color.all <- color.all[-c(3, 4, 6, 8)]

u <- c(1, 2, 4, 6, 8)
v <- c(1, 3, 4, 5, 6, 7, 8, 9)
method <- method.all[v]
color <- color.all[v]
ind <- NULL
for(i in v) {
  ind <- c(ind, (i - 1) * 9 + u)
}

out <- rbind(out.stool.1, out.vaginal.1, out.stool.2, out.vaginal.2, out.stool.3, out.vaginal.3)
val <- -log10(out[, ind])
sd <- abs(-log10(exp(1)) / out[, ind]) * out[, ind + ncol(out) / 2]

pdf('nonpara_stool_vaginal.pdf', width = 16, height = 12)  
plot.fun(val, sd, method, color)  
dev.off()




