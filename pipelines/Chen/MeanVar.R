library('data.table')
setwd('~/Hutch-Research/R_batch3')
source('GetData.R')

loadfromfile = TRUE

setwd('~/Hutch-Research/R_batch3')
dim(datM)
dim(samp)

datlist <- list(dat1[,-1],dat2[,-1],dat3[,-1])
cellTypes = unique(label)

#---------------------------------------------
# Roughtly Filter probe to keep those 
# with small variance
#---------------------------------------------

mu.all <- matrix(NA, nrow(datM),3)
sd.all <- matrix(NA, nrow(datM),3)
colnames(mu.all) <- cellTypes
colnames(sd.all) <- cellTypes
unstable <- list()

for(i in 1:3){
    mu.all[,i] <- rowMeans(datlist[[i]],na.rm = TRUE)
    sd.all[,i] <- apply(datlist[[i]],1,sd)
    unstable[[i]] <- c(unstable, which(sd.all[,i]>0.05))
}

unstable.all <- unique(unlist(unstable))
length(unstable.all)

datM <- cbind(dat1,dat2[,-1],dat3[,-1])
cpgname <- unlist(datM[-unstable.all,1])
datM <- datM[-unstable.all,-1]
rownames(datM) <- cpgname
mu.all <- mu.all[-unstable.all,]
sd.all <- sd.all[-unstable.all,]

if(use_new_probe){
  load('Chen-p2use.RData')
  nms = names(p2use)
  for(nm1 in nms){
    p1 = p2use[[nm1]]
    p2use[[nm1]] = p1[which(p1$id %in% cpgname),]
  }
  
  lapply(p2use, dim)
  
  genes = NULL
  for(lb1 in names(p2use)){
    genes = c(genes, p2use[[lb1]]$id)
  }
  
  length(genes)
  length(unique(genes))
  genes = unique(genes)
}else{
  load('_946probes.RData')
  genes = rownames(mu)
}
genes = intersect(genes, cpgname)

if(perturb = TRUE){
  mu.all = s2.all = matrix(NA,nrow = nrow(dat.est),ncol = 3)
  test <- subset(datM, rownames(datM) %in% genes)
  test <- log(test/(1-test))
  for(i in 1:3){
    s2.all[,i] <- apply(test,1,var, na.rm = TRUE)
  }
  s2.all.mixsample <- apply(test,1,var,na.rm=TRUE)
  test_add <- test + 
    diag(s2.all.mixsample) %*%  matrix(rnorm(nrow(test)*ncol(test)),nrow(test),ncol(test))  
  colnames(test_add) <- colnames(datM)
  rownames(test_add) <- genes
  
  pdf('../figures/Chen/VarM.pdf')
  plist = list()
  plist <- lapply(1:length(cellTypes),FUN=function(j){
    newplot <- ggplot(data = as.data.frame(s2.all),aes_string(x=cellTypes[j]))+
      geom_histogram(aes(y=..density..))+
      geom_density(alpha=.2,fill='#FF6666')+
      ggtitle(cell[j])
  })
  grid.arrange(grobs = plist, ncol = 3)
  dev.off()
  
}

#---------------------------------------------
# Separate dat and dat_gen
#---------------------------------------------

individual <- intersect(colnames(dat1),colnames(dat2))
individual <- intersect(individual,colnames(dat3))
samp <- sample(individual[-1],60)
gen <- setdiff(individual[-1], samp)
if(loadfromfile = TRUE){
  load('common_est.RData')
  load('common_gen.RData') 
  samp <- common_est
  gen <- common_gen
}

test_add <- exp(test_add)/(exp(test_add)+1)
colnames(test_add) <- colnames(datM)
rownames(test_add) <- genes
dat1 <- test_add[,1:196]
dat2 <- test_add[,197:393]
dat3 <- test_add[,394:525]

dat1.est <- subset(dat1,select = samp)
dat2.est <- subset(dat2,select = samp)
dat3.est <- subset(dat3,select = samp)

dat.est <- as.matrix(cbind(dat1.est,dat2.est,dat3.est))
dat.gen <- list(as.matrix(subset(dat1,select = gen)),
                as.matrix(subset(dat2,select = gen)),
                as.matrix(subset(dat3,select = gen)))
