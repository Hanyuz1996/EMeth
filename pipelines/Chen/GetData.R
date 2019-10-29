library('data.table')
setwd('~/Hutch-Research/Data/Chen')

dat1 <- fread('mono.txt.gz')
dat2 <- fread('neut.txt.gz')
dat3 <- fread('tcel.txt.gz')
datM <- cbind(dat1,dat2[,-1],dat3[,-1])
#fwrite(datM,file = 'chen_methylation_beta.txt',sep = '\t')

label <- c(rep('Monocyte',dim(dat1[,-1])[2]),rep('Neutrophil',dim(dat2[,-1])[2]),
           rep('Tcell',dim(dat3[,-1])[2]))

samp <- data.table(sample = colnames(datM[,-1]),label = label)
#fwrite(samp,file = 'chen_sample_info.txt', sep = '\t')
