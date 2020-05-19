library(data.table)
#library(multcomp)

setwd('~/Hutch-Research/R_batch3')
source('step1-estMeanVar.R')
source('~/Hutch-Research/R_batch1/_lib.R')
subsize = dim(dat.est)[2]/3
labels = c(rep(cellTypes[1],subsize),rep(cellTypes[2],subsize),rep(cellTypes[3],subsize))
dat.est <- as.matrix(dat.est)
rownames(dat.est) <- cpgname

levels = list()
levels = list()
level1 = list('mono','tcel')
level2 = list('neut','tcel')
level3 = list('neut','mono')
level4 = list('neut',c('mono','tcel'))
level5 = list('tcel',c('neut','mono'))
level6 = list('mono',c('tcel','neut'))
levels = list(level1,level2,level3,level4,level5,level6)
labels = rep(c('mono','neut','tcel'),each = length(common_est))
genes <- filterGene(dat.est, pv = 10^(-6),level = levels,label = labels)
#load('Expr_Gene.RData')
genes <- unique(unlist(genes))
genes <- cpgname[genes]
#genes <- mono[genes,1]

save(genes, file = 'probes.RData')

