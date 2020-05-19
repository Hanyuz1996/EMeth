library('data.table')

dat1 <- fread('~/Hutch-Research/Data/Chen/mono.txt.gz')
dat2 <- fread('~/Hutch-Research/Data/Chen/neut.txt.gz')
dat3 <- fread('~/Hutch-Research/Data/Chen/tcel.txt.gz')


pert <- function(mat,type = 'methylation', C = perturbconst){
  id <- mat[,1]
  num <- mat[,-1]
  if(type == 'methylation'){
    num <- log((num)/(1-num))
  }
  s2.all <- apply(num,1,var,na.rm = TRUE)
  s2.mat <- matrix(rep(s2.all,each = ncol(num)),ncol = ncol(num),byrow = TRUE)
  temp <- matrix(rnorm(dim(num)[1]*dim(num)[2]),ncol = ncol(num),nrow = nrow(num))
  noise <- sqrt(s2.mat*C) * temp
  test_add <- num + noise
  if(type == 'methylation'){
    test_add = exp(test_add)/(1+exp(test_add))
  }
  res <- cbind(id,test_add)
  colnames(res) <- colnames(mat)
  rownames(res) <- rownames(mat)
  return(res)
}

dat1 <- pert(dat1)
dat2 <- pert(dat2)
dat3 <- pert(dat3)

datM <- cbind(dat1,dat2[,-1],dat3[,-1])
#fwrite(datM,file = 'chen_methylation_beta.txt',sep = '\t')

label <- c(rep('Monocyte',dim(dat1[,-1])[2]),rep('Neutrophil',dim(dat2[,-1])[2]),
           rep('Tcell',dim(dat3[,-1])[2]))

samp <- data.table(sample = colnames(datM[,-1]),label = label)
#fwrite(samp,file = 'chen_sample_info.txt', sep = '\t')

#pdf('../figures/Chen/VarM.pdf')
#plist = list()
#plist <- lapply(1:length(cellTypes),FUN=function(j){
#  newplot <- ggplot(data = as.data.frame(s2.all),aes_string(x=cellTypes[j]))+
#    geom_histogram(aes(y=..density..))+
#    geom_density(alpha=.2,fill='#FF6666')+
#    ggtitle(cellTypes[j])
#})
#grid.arrange(grobs = plist, ncol = 3)
#dev.off()
