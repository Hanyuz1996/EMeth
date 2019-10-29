library(data.table)
library(stringr)
library(quadprog)
library(e1071)
library(ggplot2)
library(MASS)
library(gridExtra)
library(quantreg)
setwd('~/Hutch-Research/R_batch3')
source('~/Hutch-Research/R_batch1/_lib.R')
#source('GetData.R')
source('MeanVar.R')
#source('filterProbes.R')
load('truerho.RData')
truerho <- deconv_expr
perturb_gen <- TRUE

sampsize = length(common_gen)
use_new_probe = TRUE
penalty = 500
aber = FALSE

rownames(dat.est) = cpgname
for(i in 1:3){
rownames(dat.gen[[i]]) = cpgname
}
table(colnames(dat.gen[[1]]) == colnames(dat.gen[[2]]))
table(colnames(dat.gen[[1]]) == colnames(dat.gen[[3]]))

dat.est = dat.est[genes,]
dat.gen = list(dat.gen[[1]][genes,],dat.gen[[2]][genes,],dat.gen[[3]][genes,])
dim(dat.est)

mu = s2 = matrix(NA, nrow = nrow(dat.est), ncol = 3)
mu.gen = s2.gen = matrix(NA, nrow = nrow(dat.gen[[1]]), ncol = 3)
size.est = dim(dat.est)[2]/3
for(i in 1:3){
    mu[,i] <- rowMeans(dat.est[,((i-1)*size.est+1):(size.est*i)], na.rm = TRUE)
    s2[,i] <- apply(dat.est[,((i-1)*size.est+1):(size.est*i)],1,var, na.rm = TRUE)
    mu.gen[,i] <- rowMeans(dat.gen[[i]],na.rm=TRUE)
    s2.gen[,i] <- apply(dat.gen[[i]],1,var,na.rm = TRUE)
}

rownames(mu) = rownames(s2) = rownames(dat.est)
rownames(mu.gen) = rownames(s2.gen) = rownames(dat.est)

truerho <- deconv_expr
rownames(truerho) <- samp_gen
int <- intersect(samp_gen,common_gen)
rho_expr <- truerho[int,]
mixvalue <- matrix(NA,length(genes),length(int))

mix <- function(sampsize,lambda = 10){
  rho <- matrix(NA,sampsize,3)
  #mix <- matrix(NA,length(genes),sampsize)
  colnames(rho) <- cellTypes
  for(i in 1:length(int)){
  temp <- rho_expr[i,]
  temp <- temp/sum(temp)
  rho[i, ] <- temp
  mixvalue[, i] <- temp[1]* dat.gen[[1]][,i] + temp[2] * dat.gen[[2]][,i] + temp[3]*dat.gen[[3]][,i]
  stemp <- lapply(1:3, FUN = function(j){dat.gen[[j]][,i]*(1-dat.gen[[j]][,i])})
  stemp <- as.matrix(cbind(stemp[[1]],stemp[[2]],stemp[[3]]))
  print(dim(stemp))
  print(stemp[1:5,])
  stemp <- stemp %*% temp
  idic <- sample(nrow(stemp)*ncol(stemp),0.1*nrow(stemp)*ncol(stemp),replace=FALSE)
  stemp[idic] <- stemp[idic]*lambda
  mixvalue[,i] <- mixvalue[,i] + rnorm(nrow(mixvalue))*stemp
} 
  list(rho = rho, bulk = mixvalue)
}

samp <- mix(sampsize)
Y <- samp$bulk
true.rho <- samp$rho
colnames(true.rho) <- cellTypes
eta <- rep(0, ncol(Y))

eta[eta > 0.99] = 0.99


temp       = runif(sampsize * length(cellTypes)) * 0.2 -0.1
methods = c("WeightEM","OriEM","svr","ls")
rho     = array(data = NA, dim = c(ncol(Y),length(cellTypes),length(methods)),
                dimnames = list(1:ncol(Y),cellTypes,methods))
alpha = rep(1/length(cellTypes),length(cellTypes))
simsize = ncol(Y)

  temp       = runif(simsize * length(cellTypes)) * 0.2 -0.1
  rho_init   = matrix(temp,ncol = length(cellTypes))
  nu0_init   = runif(nrow(Y))
  sigma_c_init = 0.1
  lambda_init  = 2

 for(j in 1:ncol(Y)){
    if(j %% 50 == 0){ cat(j, date(), "\n") }
    y    = Y[,j]
    X    = as.data.frame(mu)
    Xmat = mu

    svrmodel        = svm(y ~ ., data = X, kernel = 'linear')
    temp            = (t(svrmodel$coefs) %*% svrmodel$SV)
    temp[temp < 0]  = 0
    cat('sum of coefficients',sum(temp))
    rho[j,,'svr']   = (1-eta[j])*temp/sum(temp)

    temp            = lm(y ~ .-1,data = X)$coefficients
    temp[temp < 0]  = 0
    rho[j,,'ls']    = (1-eta[j])*temp/sum(temp)
  }

  rho_init = rho[,,'ls']
  K = nrow(Y)
  Y_ab     = Y - mu %*% t(rho_init)
    if(aber){
    for(k in 1:K){
      nu0_init[k] = min(1,max(0,sum(eta * Y_ab[k,])/sum( eta^2)))
    }
    }

try({
  temp <- apply(s2,1,max)
  temp <- temp/median(temp)
  lb <- quantile(temp,0.15)
  ub <- quantile(temp,0.85)
  temp[temp<lb] <- lb
  temp[temp>ub] <- ub
  W <- matrix(rep(temp,ncol(Y)),ncol = ncol(Y),byrow = FALSE)
  hundrediter_weight = deconvEM(Y,eta,mu,aber = aber, V='w', weight = W,
                                0.5,rho_init,nu0_init,0.1, lambda_init = 2,
                                nu = penalty, maxiter = 100)
  rho[,,'WeightEM'] = hundrediter_weight$rho
  sigma_c_est_weight   = hundrediter_weight$sigma_c
  pi_est_weight      = hundrediter_weight$pi_a
  nu0_est_weight     = hundrediter_weight$nu0
  iter_weight        = hundrediter_weight$iter
  gamma_weight       = hundrediter_weight$gamma
})

try({
  hundrediter = deconvEM(Y,eta,mu,aber = aber, V='c', weight = matrix(1,5,5),
                            0.5,rho_init,nu0_init,sigma_c_init = 0.1, lambda_init = 2,
                            nu = penalty, maxiter = 100)
  rho[,,'OriEM'] = hundrediter$rho
  sigma_c_est = hundrediter$sigma_c
  pi_est      = hundrediter$pi_a
  nu0_est     = hundrediter$nu0
  iter        = hundrediter$iter
  gamma       = hundrediter$gamma
})

dir.create('~/Hutch-Research/figures/Chen')
setwd("~/Hutch-Research/figures/Chen")

utypes = intersect(cellTypes,colnames(true.rho))

cormat <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- utypes

err <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(err) <- methods
rownames(err) <- utypes
for(i in 1:length(utypes)){
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){
    cor(rho[,utypes[i],methods[j]],true.rho[,utypes[i]])
  })
  err[i,] <- sapply(1:length(methods), FUN = function(j){
    mean((rho[,utypes[i],methods[j]] - true.rho[,utypes[i]])^2)
  })
}

for(i in 1:length(utypes)){
  pdf(sprintf('%s_methylation_vs_true_mix.pdf',utypes[i]))
  plist = list()
  plist <- lapply(1:length(methods), FUN = function(j){
    tempdata = cbind(rho[,utypes[i],methods[j]],true.rho[,utypes[i]],eta)
    colnames(tempdata) <- c("methylation","truemix","eta")
    newplot <- ggplot(data = as.data.frame(tempdata), aes(x=methylation,y=truemix,color=eta))+ xlim(0,1) + ylim(0,1) +
      geom_point() + geom_abline(intercept = 0,slope = 1) + ggtitle(methods[j])
  })
  grid.arrange(grobs = plist,ncol=2)
  dev.off()
}
 pdf('correlation.pdf')
  plist = list()
  plist <- lapply(1:length(utypes),FUN = function(i){
  tempdata = data.frame(methods,correlation = cormat[utypes[i],] )
  corplot <- ggplot(tempdata,aes(methods,correlation))+geom_col()+ggtitle(utypes[i])
  })
  grid.arrange(grobs = plist, ncol = 2)
  dev.off()

  pdf('RootedMSE.pdf')
  plist = list()
  plist <- lapply(1:length(utypes),FUN = function(i){
    tempdata = data.frame(methods, rootedMSE = sqrt(err[utypes[i],]) )
    corplot <- ggplot(tempdata,aes(methods, rootedMSE))+geom_col()+ggtitle(utypes[i])
  })
  grid.arrange(grobs = plist, ncol = 2)
  dev.off()

if(1<0){
pdf('~/Hutch-Research//figures/Chen/Methy/testmean_vs_trainmean_methy_noisy.pdf',width = 20)
plist = list()
 plist <- lapply(1:length(cellTypes),FUN=function(j){
  tempdata <- cbind(mu[,j],mu.gen[,j])
  colnames(tempdata) <- c('train_mean','test_mean')
  newplot <- ggplot(data = as.data.frame(tempdata),aes(x=train_mean,y=test_mean))+ xlim(0,1)+ylim(0,1)+geom_point()+geom_abline(intercept = 0,slope = 1) + ggtitle(cellTypes[j]) 
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()

ind <- sample(ncol(dat.gen[[1]]),10)
dir.create('~/Hutch-Research//figures/Chen/Methy/check_individual')
for(i in 1:10){
  pdf(sprintf('~/Hutch-Research//figures/Chen/Methy/check_individual/methy_noisy_%s.pdf',colnames(dat.gen[[1]])[ind[i]]))
  plist = list()
 plist <- lapply(1:length(cellTypes),FUN=function(j){
  tempdata <- cbind(mu[,j],dat.gen[[j]][,ind[i]])
  sampname <- colnames(dat.gen[[j]])[ind[i]]
  colnames(tempdata) <- c('train_mean',sampname)
  newplot <- ggplot(data = as.data.frame(tempdata),aes_string(x="train_mean",y=sampname))+ xlim(0,1)+ylim(0,1)+geom_point()+geom_abline(intercept = 0,slope = 1) + ggtitle(cellTypes[j])
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()
}

cell <- c(cellTypes,'Average')
samp2mean <- matrix(NA, ncol(dat.gen[[1]]),4)
for(i in 1:ncol(dat.gen[[1]])){
   samp2mean[i,1] <- sqrt(mean((dat.gen[[1]][,i]-mu[,1])^2))
   samp2mean[i,2] <- sqrt(mean((dat.gen[[2]][,i]-mu[,2])^2))
   samp2mean[i,3] <- sqrt(mean((dat.gen[[3]][,i]-mu[,3])^2))
   samp2mean[i,4] <- sqrt(mean(samp2mean[i,1:3]^2))
}
samp2mean <- as.data.frame(samp2mean)
colnames(samp2mean) <- cell
pdf('~/Hutch-Research//figures/Chen/Methy/hist_MSE_noisy.pdf')
plist = list()
plist <- lapply(1:length(cell),FUN=function(j){
newplot <- ggplot(data = samp2mean,aes_string(x=cell[j]))+
             geom_histogram(aes(y=..density..))+
             geom_density(alpha=.2,fill='#FF6666')+
             ggtitle(cell[j])
})
grid.arrange(grobs = plist, ncol = 2)
dev.off()
}
  
print(cormat)
print(err)
quit(save = 'no')
