winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

preprocess.fun <- function(otu.tab, prev.cut = 0.2, lib.cut = 1000, 
                           winsor.quan = 0.97) {
  keep.sam <- which(colSums(otu.tab) >= lib.cut)
  Y <- otu.tab[, keep.sam]
  
  keep.tax <- which(rowSums(Y > 0) / ncol(Y) >= prev.cut)
  Y <- Y[keep.tax, ]
  
  keep.sam1 <- which(colSums(Y) >= 1)
  Y1 <- Y[, keep.sam1]
  
  Y1 <- winsor.fun(Y1, winsor.quan) 
  
  return(list(Y = Y1, keep.sam = keep.sam[keep.sam1], keep.tax = keep.tax))
}

library(BMDD)

####################################################################

load('combo.RData')

otu.tab.combo <- otu.tab
otu.names.combo <- otu.names.mat[, 5]

dim(otu.tab.combo)
sum(otu.tab.combo == 0) / (dim(otu.tab.combo)[1] * dim(otu.tab.combo)[2])

genus.tab.combo <- rowsum(otu.tab.combo, group = otu.names.combo)
dim(genus.tab.combo)
sum(genus.tab.combo == 0) / (dim(genus.tab.combo)[1] * dim(genus.tab.combo)[2])


library(HMP16SData)
library(phyloseq)

phy <- V13() %>%
  as_phyloseq()

otu.tab.hmp <- otu_table(phy)
otu.names.hmp <- tax_table(phy)[, 6]

dim(otu.tab.hmp)
sum(otu.tab.hmp == 0) / (dim(otu.tab.hmp)[1] * dim(otu.tab.hmp)[2])

genus.tab.hmp <- rowsum(otu.tab.hmp, group = otu.names.hmp)
dim(genus.tab.hmp)
sum(genus.tab.hmp == 0) / (dim(genus.tab.hmp)[1] * dim(genus.tab.hmp)[2])

############################################

data.prep.combo <- preprocess.fun(genus.tab.combo, prev.cut = 0.05, lib.cut = 1000, winsor.quan = 0.97)
genus.tab.combo <- data.prep.combo$Y
dim(genus.tab.combo)
sum(genus.tab.combo == 0) / (dim(genus.tab.combo)[1] * dim(genus.tab.combo)[2])

fea.combo <- rownames(genus.tab.combo)


ind.stool <- which(sample_data(phy)$HMP_BODY_SUBSITE=="Stool")
genus.tab.hmp <- genus.tab.hmp[, ind.stool]
dim(genus.tab.hmp)
sum(genus.tab.hmp == 0) / (dim(genus.tab.hmp)[1] * dim(genus.tab.hmp)[2])

data.prep.hmp <- preprocess.fun(genus.tab.hmp, prev.cut = 0.05, lib.cut = 1000, winsor.quan = 0.97)
genus.tab.hmp <- data.prep.hmp$Y
dim(genus.tab.hmp)
sum(genus.tab.hmp == 0) / (dim(genus.tab.hmp)[1] * dim(genus.tab.hmp)[2])

fea.hmp <- rownames(genus.tab.hmp)


fea <- intersect(fea.combo, fea.hmp)
tab.combo <- genus.tab.combo[fea, ]
tab.hmp <- genus.tab.hmp[fea, ]

####################################################################

set.seed(666)
res.combo <- bmdd(tab.combo, type = 'count', trace = TRUE)

alp0.est <- res.combo$alpha$alp0
alp1.est <- res.combo$alpha$alp1
pi.est <- res.combo$pi

ind <- which(alp0.est > alp1.est)
tmp <- alp0.est[ind]
alp0.est[ind] <- alp1.est[ind]
alp1.est[ind] <- tmp
pi.est[ind] <- 1 - pi.est[ind]

alp0.est.combo <- alp0.est
alp1.est.combo <- alp1.est
pi.est.combo <- pi.est


set.seed(666)
res.hmp <- bmdd(tab.hmp, type = 'count', trace = TRUE)

alp0.est <- res.hmp$alpha$alp0
alp1.est <- res.hmp$alpha$alp1
pi.est <- res.hmp$pi

ind <- which(alp0.est > alp1.est)
tmp <- alp0.est[ind]
alp0.est[ind] <- alp1.est[ind]
alp1.est[ind] <- tmp
pi.est[ind] <- 1 - pi.est[ind]

alp0.est.hmp <- alp0.est
alp1.est.hmp <- alp1.est
pi.est.hmp <- pi.est


####################################################################

library(ggplot2)
library(ggrepel)

df <- data.frame(val1 = pi.est.combo, val2 = pi.est.hmp, name = fea)

p1 <- ggplot(df, aes(x = val1, y = val2, label = name)) +
  geom_point(color = "steelblue", alpha = 0.6) +  # point size by value
  geom_text_repel(
    size = 2.8,                # smaller text
    max.overlaps = Inf,        # force all labels
    force = 5,                 # stronger push
    max.iter = 10000,          # more iterations to resolve overlaps
    box.padding = 0.6,         # space around labels
    point.padding = 0.3,       # space around points
    segment.size = 0.2,        # thin leader lines
    segment.alpha = 0.6        # faint leader lines
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression(pi ~ " estimated from COMBO"),
    y = expression(pi ~ " estimated from HMP"),
    title = ""
  ) +
  theme(legend.position = "right")

pdf("pi_combo_hmp.pdf", width = 12, height = 9)
p1
dev.off()


####################################################################

library(ggpubr)

X.combo <- t(t(tab.combo) / colSums(tab.combo))
X.hmp <- t(t(tab.hmp) / colSums(tab.hmp))

plot.fun <- function(x, tit) {
  ggplot() +
    geom_histogram(aes(x = x, y = ..density..), #boundary = 0,
                   fill = 'steelblue', color = 'black') +
    xlab('') + ylab('') +
    ggtitle(tit) +
    theme_minimal(base_size = 25)
}


p1 <- plot.fun(X.combo[27, ], "Parabacteroides (COMBO)")
p2 <- plot.fun(X.hmp[27, ], "Parabacteroides (HMP)")

p3 <- plot.fun(X.combo[17, ], "Eubacterium (COMBO)")
p4 <- plot.fun(X.hmp[17, ], "Eubacterium (HMP)")

big_plot <- ggarrange(plotlist = list(p1, p2, p3, p4), 
                      nrow = 2, ncol = 2, common.legend = TRUE)

final_plot <- annotate_figure(
  big_plot,
  bottom = text_grob("Proportion", size = 25),
  left = text_grob("Density", size = 25, rot = 90)
)

pdf('hist_combo_hmp.pdf', width = 12, height = 8)
final_plot
dev.off()








