library(data.table)
#library(multcomp)

setwd('~/Hutch-Research/R_batch3')
source('MeanVar.R')
source('~/Hutch-Research/R_batch1/_lib.R')
subsize = 60
sam <- data.table(id = 1:dim(dat.est)[2],label = c(rep(cellTypes[1],subsize),rep(cellTypes[2],subsize),rep(cellTypes[3],subsize)))
dat.est <- as.matrix(dat.est)
rownames(dat.est) <- cpgname

uct <- unique(sam$label)
uct

levels = list()
for(k in 1:length(uct)){
   ct1 = uct[k]
   ct2 = paste("Non", ct1)

   level1 = list()

   level1[[ct1]] = ct1
   level1[[ct2]] = setdiff(uct,ct1)

   levels[[k]]   = level1
}
length(levels)

dir.create('~/Hutch-Research/figures/Chen_Filter')
path = '~/Hutch-Research/figures/Chen_Filter'
p2use = list()

for(kk in 1:length(levels)){
  cat(kk, date(), "\n")
  level1 = levels[[kk]]
  nms    = names(level1)
  nm1    = sprintf("%s_vs_%s", nms[1], nms[2])

  p2use[[nm1]] = probeSelect(level1, sam, dat.est, path=path, dataType='methylation')
}

lapply(p2use,dim)

save(p2use, file = 'Chen-p2use.RData')

quit(save = 'no')
