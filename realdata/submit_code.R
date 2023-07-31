library(MicrobiomeStat)
library(BMDD)
library(phyloseq)
otu_filter <- function(feature.dat, prev = 0.1, dep = 1000){ # otu * sam
  idx <- apply(feature.dat, 1, function(x) sum(x>0)>(ncol(feature.dat)*prev))
  idx2 <- colSums(feature.dat) > dep
  return(feature.dat[idx,idx2])
} 


## ==== real data analysis [# discoveries]====
dd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Data/noshuffle"
rd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Result/noshuffle"

## Dateset1 - CRC: YuJ_2015.Rdata
setwd(dd)
file = 'YuJ_2015.Rdata'
load(file)
otu.tab <- as.data.frame(as.matrix(otu_table(phy)))
meta.dat <- as.data.frame(as.matrix(sample_data(phy)))
meta.dat$grp <- as.factor(meta.dat$grp)
feature.dat <- otu_filter(otu.tab)
meta.dat <- meta.dat[colnames(feature.dat),]

bmdd.fit <- bmdd(W = feature.dat, type ='count')
prop.bmdd <- t(t(bmdd.fit$beta)/colSums(bmdd.fit$beta))
bmdd.obj  <- linda(feature.dat=prop.bmdd, meta.dat=meta.dat,
                   formula = '~grp', feature.dat.type = 'proportion')
bmdd.res <- bmdd.obj$output[[1]][,'padj',drop =F]

linda.obj  <- linda(feature.dat=feature.dat, meta.dat=meta.dat,
                    formula = '~grp', feature.dat.type = 'count')
linda.res <- linda.obj$output[[1]][,'padj',drop =F]
setwd(rd)
save(bmdd.fit, bmdd.res, linda.res, file= paste0(gsub('.Rdata','',file),'_res.RData'))




## Dateset2 - [1] USA vs Malawi: adults
setwd(dd)
file = 'USA_Malawi.Rdata'
load(file)
phy <- tax_glom(phy, taxrank = 'Species')
otu.tab <- as.data.frame(as.matrix(otu_table(phy)))
meta.dat <- as.data.frame(as.matrix(sample_data(phy)))
meta.dat$grp <- as.factor(meta.dat$grp)
meta.dat <- meta.dat[meta.dat$AGE_IN_YEARS!='None',]
meta.dat$AGE <- as.numeric(meta.dat$AGE)
meta.dat <- meta.dat[meta.dat$AGE >17,] # subset adults
otu.tab <- otu.tab[,rownames(meta.dat)]
feature.dat <- otu_filter(otu.tab)
meta.dat <- meta.dat[colnames(feature.dat),]

bmdd.fit  <- bmdd(W = feature.dat, type ='count')
prop.bmdd <- t(t(bmdd.fit$beta)/colSums(bmdd.fit$beta))
bmdd.obj  <- linda(feature.dat=prop.bmdd, meta.dat=meta.dat,
                   formula = '~grp', feature.dat.type = 'proportion')
bmdd.res <- bmdd.obj$output[[1]][,'padj',drop =F]

linda.obj  <- linda(feature.dat=feature.dat, meta.dat=meta.dat,
                    formula = '~grp', feature.dat.type = 'count')
linda.res <- linda.obj$output[[1]][,'padj',drop =F]
setwd(rd)
save(bmdd.fit, saver.fit, bmdd.res, linda.res, saver.res, file= paste0(gsub('.Rdata','',file),'_Species_adults_res.RData'))

## Dateset2 - [2] USA vs Malawi: infants
setwd(dd)
file = 'USA_Malawi.Rdata'
load(file)
phy <- tax_glom(phy, taxrank = 'Species')
otu.tab <- as.data.frame(as.matrix(otu_table(phy)))
meta.dat <- as.data.frame(as.matrix(sample_data(phy)))
meta.dat$grp <- as.factor(meta.dat$grp)
meta.dat <- meta.dat[meta.dat$AGE_IN_YEARS!='None',]
meta.dat$AGE <- as.numeric(meta.dat$AGE)
meta.dat <- meta.dat[meta.dat$AGE <=3,] # subset infants
otu.tab <- otu.tab[,rownames(meta.dat)]
feature.dat <- otu_filter(otu.tab)
meta.dat <- meta.dat[colnames(feature.dat),]

bmdd.fit  <- bmdd(W = feature.dat, type ='count')
prop.bmdd <- t(t(bmdd.fit$beta)/colSums(bmdd.fit$beta))
bmdd.obj  <- linda(feature.dat=prop.bmdd, meta.dat=meta.dat,
                   formula = '~grp', feature.dat.type = 'proportion')
bmdd.res <- bmdd.obj$output[[1]][,'padj',drop =F]

linda.obj  <- linda(feature.dat=feature.dat, meta.dat=meta.dat,
                    formula = '~grp', feature.dat.type = 'count')
linda.res <- linda.obj$output[[1]][,'padj',drop =F]
setwd(rd)
save(bmdd.fit, saver.fit, bmdd.res, linda.res, saver.res, file= paste0(gsub('.Rdata','',file),'_Species_infants_res.RData'))





## Dateset3 - CRC: t1d_mejialeon.Rdata
setwd(dd)
file = "t1d_mejialeon.Rdata"
load(file)
otu.tab <- otu_table
meta.dat <- metadata
meta.dat$grp <- as.factor(meta.dat$grp)
feature.dat <- otu_filter(otu.tab)
meta.dat <- meta.dat[colnames(feature.dat),]

bmdd.fit <- bmdd(W = feature.dat, type ='count')
prop.bmdd <- t(t(bmdd.fit$beta)/colSums(bmdd.fit$beta))
bmdd.obj  <- linda(feature.dat=prop.bmdd, meta.dat=meta.dat,
                   formula = '~grp', feature.dat.type = 'proportion')
bmdd.res <- bmdd.obj$output[[1]][,'padj',drop =F]

linda.obj  <- linda(feature.dat=feature.dat, meta.dat=meta.dat,
                    formula = '~grp', feature.dat.type = 'count')
linda.res <- linda.obj$output[[1]][,'padj',drop =F]
setwd(rd)
save(bmdd.fit, bmdd.res, linda.res, file= paste0(gsub('.Rdata','',file),'_res.RData'))




## ==== FDR calculation ====
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  file = args[1]
  iter = as.integer(args[2])
}

dd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Data/shuffle"
rd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Result/shuffle"

setwd(md)
load(file)
meta.dat <- as.data.frame(meta[,iter, drop =F])
colnames(meta.dat) <- 'grp'

bmdd.fit <- bmdd(W = feature.dat, type ='count')
prop.bmdd <- t(t(bmdd.fit$beta)/colSums(bmdd.fit$beta))
bmdd.obj  <- linda(feature.dat=prop.bmdd, meta.dat=meta.dat, 
                   formula = '~grp', feature.dat.type = 'proportion')
bmdd.res <- bmdd.obj$output[[1]][,'padj',drop =F]
linda.obj  <- linda(feature.dat=feature.dat, meta.dat=meta.dat, 
                    formula = '~grp', feature.dat.type = 'count')
linda.res <- linda.obj$output[[1]][,'padj',drop =F]
setwd(rd)
save(bmdd.fit,bmdd.res,linda.res, file= paste0(gsub('.Rdata','',file),'_',iter,'_res.RData'))






## === summary and plot====
library(dplyr)
library(tidyr)
library(tibble)
library(tidyverse)
library(ggvenn)
library(ggpubr)

wd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Data/noshuffle/"
rd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Result/noshuffle/"
sd <- "/Users/M216453/Documents/Mayo_project/2023_06_22_BMDD/submit/Result/shuffle/"


## ==== CRC ====
# venn
setwd(rd)
df.res <- NULL
plts <- list()
files <- list.files(pattern = 'USA_Malawi')
for(file in files){
  cat('[ ',file,' ]\n')
  setwd(rd)
  load(file)
  cutoff = 0.05
  plts[[file]] <- merge(bmdd.res, linda.res, by = 0, all = T) %>%
    column_to_rownames('Row.names') %>% dplyr::rename(bmdd=1, linda=2) %>%
    mutate(across(everything(), ~if_else(. <= cutoff, TRUE, FALSE))) %>% 
    ggplot() +
    geom_venn(aes(A = linda, B = bmdd),set_names = c('LinDA','BMDD+LinDA'), show_percentage = FALSE, fill_color = c('forestgreen','yellow')) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  df <- merge(bmdd.res, linda.res, by = 0, all = T) %>% 
    column_to_rownames('Row.names') %>% dplyr::rename(bmdd=1, linda=2)
  rownames(df) <- paste0('o',rownames(df))
  
  setwd(wd)
  load(gsub('_Species_adults_res.RData|_Species_infants_res.RData','.Rdata',file))
  tax <- as.data.frame(tax_table(phy))
  rownames(tax) <- paste0('o',rownames(tax))
  tax <- tax[rownames(df),, drop =F]
  df.res <- rbind(df.res, merge(tax[,sub('.*_','',gsub('_infants_res.RData|_adults_res.RData','',file)),drop =F] %>% dplyr::rename(taxa = 1), df, by = 0, all = T) %>% mutate(file = gsub('_res.RData','',file)))
}

ggarrange(plts$USA_Malawi_Species_adults_res.RData, plts$USA_Malawi_Species_infants_res.RData)

# FDR
setwd(sd)
files <- list.files(pattern = 'adults|infants')
bmdd.ct <- linda.ct <- c()
for(file in files){
  bmdd.res <- linda.res <-NULL
  load(file)
  bmdd.ct[file] <- sum(bmdd.res$padj <= cutoff)
  linda.ct[file] <- sum(linda.res$padj <= cutoff)
}
linda.df <- as.data.frame(linda.ct) %>% rownames_to_column('id') %>% mutate(data = gsub('_shuffle.*','',id), level = sub('.*_','',gsub('_infants.*|_adults.*','',id))) %>% mutate(class = sub('.*_','',data))
bmdd.df <- as.data.frame(bmdd.ct) %>% rownames_to_column('id') %>% mutate(data = gsub('_shuffle.*','',id), level = sub('.*_','',gsub('_infants.*|_adults.*','',id))) %>% mutate(class = sub('.*_','',data))
fdr.gut <- merge(aggregate(bmdd.ct~data + level + class, bmdd.df, function(x) mean(x>0)), aggregate(linda.ct~data, linda.df, function(x) mean(x>0)), by = 'data') 
fdr.gut[fdr.gut$level=='Species',]



### ===== YuJ_2015 =====
## venn
file = "YuJ_2015_res.RData"
setwd(rd)
load(file)
merge(linda.res, bmdd.res, by = 0, all=T) %>% 
  full_join(saver.res %>% rownames_to_column('Row.names')) %>%
  column_to_rownames('Row.names') %>% 
  dplyr::rename(linda=1, bmdd = 2) %>% 
  mutate(across(everything(), ~if_else(. <= 0.05, TRUE, FALSE))) %>% 
  ggplot() +
  geom_venn(aes(A = linda, B = bmdd),set_names = c('LinDA','BMDD+LinDA'), show_percentage = FALSE, fill_color = c('forestgreen','yellow')) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

## FDR
setwd(sd)
files <- list.files(pattern = 'YuJ_')
bmdd.ct <- linda.ct <- c()
for(file in files){
  load(file)
  bmdd.ct[file] <- sum(bmdd.res$padj <= cutoff)
  linda.ct[file] <- sum(linda.res$padj <= cutoff)
}
mean(linda.ct>0);mean(bmdd.ct>0)





####====== t1d_mejialeon ======
## venn
setwd(rd)
file = "t1d_mejialeon_res.RData"
load(file)
cutoff = 0.05
merge(linda.res, bmdd.res, by = 0, all=T) %>% 
  column_to_rownames('Row.names') %>% 
  dplyr::rename(linda=1, bmdd = 2) %>% 
  mutate(across(everything(), ~if_else(. <= cutoff, TRUE, FALSE))) %>% 
  ggplot() +
  geom_venn(aes(A = linda, B = bmdd),set_names = c('LinDA','BMDD+LinDA'), show_percentage = FALSE, fill_color = c('forestgreen','yellow')) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
df <- merge(bmdd.res, linda.res, by = 0, all = T) %>% 
  column_to_rownames('Row.names')%>% dplyr::rename(bmdd=1, linda=2)
df[df$bmdd <=cutoff&df$linda>cutoff,]
setwd(dd)
load('t1d_mejialeon.Rdata')
tax <- tax_table %>% column_to_rownames('otu')
tax <- tax[rownames(df),, drop =F]
df <- merge(tax, df, by = 0, all = T) %>% column_to_rownames('Row.names')
df1 <- df[df$bmdd <= cutoff & df$linda > cutoff,-1]
dim(df1[df1$genus %in% c('Prevotella'),])
dim(df1[df1$genus %in% c('Bacteroides'),])


## FDR
setwd(sd)
files <- list.files(pattern = 't1d_mejialeon')
linda.ct <- bmdd.ct <- NULL
for(file in files){
  linda.res <- bmdd.res <- NULL
  load(file)
  linda.ct <- c(linda.ct, sum(linda.res$padj <= cutoff))
  bmdd.ct <- c(bmdd.ct, sum(bmdd.res$padj <= cutoff))
}
mean(linda.ct>0);mean(bmdd.ct>0)



