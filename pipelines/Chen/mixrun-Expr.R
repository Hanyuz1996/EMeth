library(ggplot2)
library(gridExtra)
library(scales)
setwd('~/Hutch-Research/R_batch3')
load('truerho.Rdata')
load('samp_est.RData')
load('samp_gen.RData')
cellTypes <- c("Monocyte","Neutrophill","Tcell")

library(data.table)
tcel <- fread('tcel_gene_nor_combat_20151109.txt')
neut <- fread('neut_gene_nor_combat_20151109.txt')
mono <- fread('mono_gene_nor_combat_20151109.txt')
alldat <- list(mono,neut,tcel)

common <- intersect(colnames(tcel),intersect(colnames(neut),colnames(mono)))
samp_est <- unique(samp_est)
common_est <- intersect(samp_est,common)
common_gen <- intersect(samp_gen,common)

levels = list()
level1 = list('mono','tcel')
level2 = list('neut','tcel')
level3 = list('neut','mono')
level4 = list('neut',c('mono','tcel'))
level5 = list('tcel',c('neut','mono'))
level6 = list('mono',c('tcel','neut'))
levels = list(level1,level2,level3,level4,level5,level6)
labels = rep(c('mono','neut','tcel'),each = length(common_est))
source('selectGene.R')

mu.all <- matrix(NA,dim(mono)[1],3)
sd.all <- matrix(NA,dim(mono)[1],3)
unstable <- list()
for(i in 1:3){
  temp <- 2^(alldat[[i]][,-1])
  mu.all[,i] <- rowMeans(temp)
  sd.all[,i] <- apply(temp,1,sd)
  test <- mu.all[,i]/sd.all[,i]
  cat(quantile(sd.all[,i],0.95),quantile(mu.all[,i],0.95),quantile(test,0.95),'\n')
  level <- list(which((sd.all[,i] > quantile(sd.all[,i],0.95)) |
                      (mu.all[,i] > quantile(mu.all[,i],0.95)) |
                      (test < quantile(test,0.05))))
  unstable <- c(unstable,level)
}
unstable <- unlist(unique(unstable))
genesym <- mono[-unstable,1]
mono <- 2^(mono[-unstable,-1])
neut <- 2^(neut[-unstable,-1])
tcel <- 2^(tcel[-unstable,-1])
exprdat.est <- cbind(subset(mono,select = common_est),
                     subset(neut,select = common_est),
                     subset(tcel,select = common_est))
row.names(exprdat.est) <- unlist(genesym)

genes <- filterGene(exprdat.est, pv = 10^(-6),level = levels,label = labels)
genes <- unique(unlist(genes))
#genes <- mono[genes,1]
mono <- mono[genes,]
neut <- neut[genes,]
tcel <- tcel[genes,]

dim(mono)
table(is.na(mono[,]))

genesym <- unlist(genesym)[genes]
genesym <- genesym[-2]
mono <- mono[-2,]
neut <- neut[-2,]
tcel <- tcel[-2,]

mono <- cbind(genesym,mono)
neut <- cbind(genesym,neut)
tcel <- cbind(genesym,tcel)

truerho <- deconv_expr
rownames(truerho) <- samp_gen
int <- intersect(samp_gen,common_gen)

exprdat.est <- list(subset(mono,select = common_est),
                     subset(neut,select = common_est),
                     subset(tcel,select = common_est))
exprdat.gen <- list(subset(mono,select = int),
                    subset(neut,select = int),
                    subset(tcel,select = int))

mu =  matrix(NA,nrow(tcel),3)
mu_gen = matrix(NA,nrow(tcel),3)
s2  = matrix(NA,nrow(tcel),3)
s2_gen = matrix(NA,nrow(tcel),3)
colnames(mu) <- cellTypes
for(i in 1:3){
  mu[,i] <- rowMeans(exprdat.est[[i]][,-1])
  mu_gen[,i] <- rowMeans(exprdat.gen[[i]][,-1])
  s2[,i] <- apply(exprdat.est[[i]][,-1],1,sd)
  s2_gen[,i] <- apply(exprdat.gen[[i]][,-1],1,sd)
}

dir.create('~/Hutch-Research/figures/Chen')
setwd("~/Hutch-Research/figures/Chen")

utypes = intersect(cellTypes,colnames(deconv_expr))

rho_expr <- truerho[int,]
mix <- matrix(NA,nrow(neut),length(int))

pdf('testmean_vs_trainmean_expr_noisy.pdf',width = 20)
plist = list()
plist <- lapply(1:length(cellTypes),FUN=function(j){
  tempdata <- cbind(mu[,j],mu_gen[,j])
  colnames(tempdata) <- c('train_mean','test_mean')
  newplot <- ggplot(data = as.data.frame(tempdata),aes(x=train_mean,y=test_mean))+ geom_point()+geom_abline(intercept = 0,slope = 1) + 
             ggtitle(cellTypes[j]) + theme(text = element_text(size = 20))
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()

ind <- sample(length(int),10)
dir.create('check_individual')
for(i in 1:10){
  pdf(sprintf('check_individual/expr_noise_%s.pdf',int[ind[i]]),width = 15)
  plist = list()
  plist <- lapply(1:length(cellTypes),FUN=function(j){
    sampname <- int[ind[i]]
    tempdata <- cbind(mu[,j],subset(exprdat.gen[[j]],select = sampname))
    colnames(tempdata) <- c('train_mean',sampname)
    newplot <- ggplot(data = as.data.frame(tempdata),aes_string(x="train_mean",y=sampname))+ geom_point()+geom_abline(intercept = 0,slope = 1)+ 
               ggtitle(cellTypes[j]) + theme(text = element_text(size = 20))
  })
  grid.arrange(grobs = plist, ncol = 3)
  dev.off()
}


for(i in 1:length(int)){
  temp <- rho_expr[i,]
  samp <- int[i]
  mix[,i] <- temp[1]*exprdat.gen[[1]][[samp]]+ 
             temp[2]*exprdat.gen[[2]][[samp]]+
             temp[3]*exprdat.gen[[3]][[samp]]
}
mu <- cbind(mono[,1],as.data.frame(mu))
colnames(mu) <- c('Gene_Symbols',cellTypes)
fwrite(mu,'refChen.txt',sep='\t',col.names = TRUE)

mix <- cbind(as.matrix(mono[,1]),as.data.frame(mix))
colnames(mix) <- c('Gene_Symbols',int)
fwrite(mix,'mixChen.txt',sep = '\t', col.names = TRUE)

cbs <- fread('cibersortx_chen.txt')
sampname <- unlist(cbs[,1])
cbs <- cbs[,2:4]
cbsm <- as.matrix(cbs)
rownames(cbsm) <- sampname
rownames(deconv_expr) <- samp_gen
common = intersect(samp_gen,rownames(cbsm))
cbsm <- cbsm[common,]
deconv_expr <- deconv_expr[common,]

cbsm[,1] <- cbsm[,1]/1.4
cbsm[,2] <- cbsm[,2]/0.15
cbsm[,3] <- cbsm[,3]/0.4
cbsm <- cbsm/rowSums(cbsm)

save(common_est, file = 'common_est.RData')
save(common_gen,file = 'common_gen.RData')

