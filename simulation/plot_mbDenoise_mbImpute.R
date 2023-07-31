library(ggplot2)
library(wesanderson) #wes_palette

#################################################
method.all <- c('BMDD-1', 'BMDD-2', 'BMDD', 'BMDD-4',
                'mbDenoise-1', 'mbDenoise', 'mbImpute-1', 'mbImpute','pmr', 
                'SAVER', 'scImpute', 'ALRA', 
                'DMM', 'DMM-2', 'naive-1', 'naive-2', 'naive')
color.all <- c('#CC3333', '#FF6666', '#FF9999', '#FFCCCC',
               wes_palette('Darjeeling1', n = 5)[-1],wes_palette('Darjeeling2', n = 4),
               wes_palette('Chevalier1', n = 4), wes_palette('GrandBudapest2', n = 2)[2])
color.all.copy <- color.all
color.all[3] <- color.all.copy[1]
color.all[6] <- color.all.copy[5]
color.all[8] <- color.all.copy[7]

plot.fun <- function(val, sd, method, color) {
  val <- as.vector(t(val))
  sd <- as.vector(t(sd))
  
  model <- factor(1 : 6, labels = c('mbDenoise 1', 'mbDenoise 2', 'mbDenoise 3', 
                                    'mbDenoise 6', 'mbImpute 1', 'mbImpute 2'))
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

file.fun.1 <- function(i){
  paste0('result/output_mbDenoise_simulation', i, '.txt')
}
file.fun.2 <- function(i){
  paste0('result/output_mbImpute_simulation', i, '.txt')
}

out.mbDenoise1 <- as.vector(unlist(read.table(file.fun.1(1))))
out.mbDenoise2 <- as.vector(unlist(read.table(file.fun.1(2))))
out.mbDenoise3 <- as.vector(unlist(read.table(file.fun.1(3))))
out.mbDenoise6 <- as.vector(unlist(read.table(file.fun.1(6))))
out.mbImpute1 <- as.vector(unlist(read.table(file.fun.2(1))))
out.mbImpute2 <- as.vector(unlist(read.table(file.fun.2(2))))

u <- c(1, 2, 4, 6, 8)
v <- c(3, 6, 8, 9 : 13)
method <- method.all[v]
color <- color.all[v]

ind <- NULL
for(i in v) {
  ind <- c(ind, (i - 1) * 9 + u)
}

out <- rbind(out.mbDenoise1, out.mbDenoise2, out.mbDenoise3,
             out.mbDenoise6, out.mbImpute1, out.mbImpute2)
val <- -log10(out[, ind])
sd <- abs(-log10(exp(1)) / out[, ind]) * out[, ind + ncol(out) / 2]

pdf('mbDenoise_mbImpute_simulations.pdf', width = 16, height = 12)  
plot.fun(val, sd, method, color)  
dev.off()

