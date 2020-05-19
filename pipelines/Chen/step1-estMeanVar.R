library('data.table')
setwd('~/Hutch-Research/R_batch3')
source('step0-GetData.R')

loadfromfile = TRUE
perturb = TRUE
use_new_probe = TRUE

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
dat1 <- dat1[-unstable.all,]
dat2 <- dat2[-unstable.all,]
dat3 <- dat3[-unstable.all,]


#---------------------------------------------
# Separate dat and dat_gen
#---------------------------------------------

individual <- intersect(colnames(dat1),colnames(dat2))
individual <- intersect(individual,colnames(dat3))
samp <- sample(individual[-1],60)
gen <- setdiff(individual[-1], samp)
if(loadfromfile){
  load('common_est.RData')
  load('common_gen.RData') 
  samp <- common_est
  gen <- common_gen
}

#colnames(test_add) <- colnames(datM)

dat1.est <- subset(dat1,select = samp)
dat2.est <- subset(dat2,select = samp)
dat3.est <- subset(dat3,select = samp)

dat.est <- as.matrix(cbind(dat1.est,dat2.est,dat3.est))
dat.gen <- list(as.matrix(subset(dat1,select = gen)),
                as.matrix(subset(dat2,select = gen)),
                as.matrix(subset(dat3,select = gen)))


