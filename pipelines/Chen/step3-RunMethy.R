library(data.table)
library(stringr)
library(quadprog)
library(e1071)
library(ggplot2)
library(MASS)
library(gridExtra)
library(quantreg)
setwd('~/Hutch-Research/R_batch3')
perturbconst = 5
source('~/Desktop/EMeth/source/_lib.R')
source('~/Desktop/EMeth/pipelines/Chen/step1-estMeanVar.R')
source('~/Desktop/EMeth/pipelines/Chen/step2-filterProbes.R')

sampsize = dim(dat.gen[[1]])[2]
use_new_probe = TRUE
penalty = 500
aber = FALSE
if(use_new_probe){
  load('probes.RData')
}

genes = intersect(genes, cpgname)
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
for(i in 1:3){
  mu[,i] <- rowMeans(dat.est[,((i-1)*subsize+1):(subsize*i)], na.rm = TRUE)
  s2[,i] <- apply(dat.est[,((i-1)*subsize+1):(subsize*i)],1,var, na.rm = TRUE)
  mu.gen[,i] <- rowMeans(dat.gen[[i]],na.rm=TRUE)
  s2.gen[,i] <- apply(dat.gen[[i]],1,var,na.rm = TRUE)
}
rownames(mu) = rownames(s2) = rownames(dat.est)
rownames(mu.gen) = rownames(s2.gen) = rownames(dat.est)

load('truerho.Rdata')
truerho <- deconv_expr
mix <- function(sampsize,lambda = 10){
  rho <- matrix(NA,sampsize,3)
  mix <- matrix(NA,length(genes),sampsize)
  colnames(rho) <- cellTypes
  for(i in 1:sampsize){
    temp <- truerho[i,]
    rho[i, ] <- temp
    mix[, i] <- temp[1]* dat.gen[[1]][,i] + temp[2] * dat.gen[[2]][,i] + temp[3]*dat.gen[[3]][,i]
    stemp <- lapply(1:3, FUN = function(j){dat.gen[[j]][,i]*(1-dat.gen[[j]][,i])})
    stemp <- as.matrix(cbind(stemp[[1]],stemp[[2]],stemp[[3]]))
    print(dim(stemp))
    print(stemp[1:5,])
    stemp <- stemp %*% temp
    idic <- sample(nrow(stemp)*ncol(stemp),0.1*nrow(stemp)*ncol(stemp),replace=FALSE)
    stemp[idic] <- stemp[idic]*lambda
    mix[,i] <- mix[,i] + rnorm(nrow(mix))*stemp
  } 
  list(rho = rho, bulk = mix)
}

samp <- mix(sampsize,lambda = 100)
Y <- samp$bulk
deconv_expr <- samp$rho
colnames(deconv_expr) <- cellTypes
eta <- rep(0, ncol(Y))

eta[eta > 0.99] = 0.99


temp       = runif(sampsize * length(cellTypes)) * 0.2 -0.1
methods = c("LaplaceEM","OriEM","svr","ls","rls","qp")
rho     = array(data = NA, dim = c(ncol(Y),length(cellTypes),length(methods)),
                dimnames = list(1:ncol(Y),cellTypes,methods))
alpha = rep(1/length(cellTypes),length(cellTypes))
simsize = ncol(Y)

temp       = runif(simsize * length(cellTypes)) * 0.2 -0.1
rho_init   = matrix(temp,ncol = length(cellTypes))
nu0_init   = runif(nrow(Y))
sigma_c_init = 0.1
lambda_init  = 2
pi_a   = rep(0.5,simsize)
C = c(0.1,1/sqrt(10),1,sqrt(10),10)
for(j in 1:ncol(Y)){
  if(j %% 50 == 0){ cat(j, date(), "\n") }
  y    = Y[,j]
  X    = as.data.frame(mu)
  Xmat = mu
  
  cv_svr = rep(0,5)
  svrmodel1       = svm(y~., data = X,kernel = 'linear', cost = 0.1, cross= 5)
  cv_svr[1]       = mean(svrmodel1$MSE)
  svrmodel2       = svm(y~., data = X,kernel = 'linear', cost = 1/sqrt(10), cross= 5)
  cv_svr[2]       = mean(svrmodel2$MSE)
  svrmodel3       = svm(y~., data = X,kernel = 'linear', cost = 1, cross= 5)
  cv_svr[3]       = mean(svrmodel2$MSE)
  svrmodel4       = svm(y~., data = X,kernel = 'linear', cost = sqrt(10), cross= 5)
  cv_svr[4]       = mean(svrmodel2$MSE)
  svrmodel5       = svm(y~., data = X,kernel = 'linear', cost = 10)
  cv_svr[5]       = mean(svrmodel5$MSE)
  best_svr        = which.min(cv_svr)
  svrmodel        = svm(y~., data = X, kernel = 'linear', cost = C[best_svr])
  temp            = (t(svrmodel$coefs) %*% svrmodel$SV)
  temp[temp < 0]  = 0
  rho[j,,'svr']   = (1-eta[j])*temp/sum(temp)
  
  temp            = lm(y ~ .-1,data = X)$coefficients
  temp[temp < 0]  = 0
  rho[j,,'ls']    = (1-eta[j])*temp/sum(temp)
  
  temp            = rlm(y ~ .-1,data = X)$coefficients
  temp[temp < 0]  = 0
  rho[j,,'rls']   = (1-eta[j])*temp/sum(temp)
  
  A = rbind(diag(rep(1,length(cellTypes))),rep(-1,length(cellTypes)))
  b = c(rep(0,length(cellTypes)),-1+eta[j])
  D = t(Xmat) %*% Xmat
  d = t(t(Xmat) %*% y)
  rho[j,,'qp']   = (solve.QP(D,d,t(A),b)$solution)
}

rho_init = rho[,,'ls']
K = nrow(Y)
Y_ab     = Y - mu %*% t(rho_init)
if(aber){
  for(k in 1:K){
    nu0_init[k] = min(1,max(0,sum(eta * Y_ab[k,])/sum( eta^2)))
  }
}
print('LaplaceEM')
hundrediter_laplace = cv.emeth(Y,eta,mu,aber = aber, V='c', init = 'default',
                               family = 'laplace', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
rho[,,'LaplaceEM'] = hundrediter_laplace[[1]]$rho

print('OriEM')
hundrediter = cv.emeth(Y,eta,mu,aber = aber, V='c', init = 'default',
                       family = 'normal', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
rho[,,'OriEM'] = hundrediter[[1]]$rho
path = sprintf('~/Hutch-Research/figures/Chen/%s',perturbconst)
dir.create(path)
setwd(path)

utypes = colnames(deconv_expr)


cormat <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- utypes

err <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(err) <- methods
rownames(err) <- utypes
for(i in 1:length(utypes)){
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){
    cor(rho[,utypes[i],methods[j]],deconv_expr[,utypes[i]])
  })
  err[i,] <- sapply(1:length(methods), FUN = function(j){
    mean((rho[,utypes[i],methods[j]] - deconv_expr[,utypes[i]])^2)
  })
}

for(i in 1:length(utypes)){
  pdf(sprintf('%s_methylation_vs_true_mix.pdf',utypes[i]),width = 15)
  plist = list()
  plist <- lapply(1:length(methods), FUN = function(j){
    tempdata = cbind(rho[,utypes[i],methods[j]],deconv_expr[,utypes[i]],eta)
    colnames(tempdata) <- c("methylation","truemix","eta")
    newplot <- ggplot(data = as.data.frame(tempdata), aes(x=methylation,y=truemix))+ xlim(0,1) + ylim(0,1) +
      geom_point() + geom_abline(intercept = 0,slope = 1) + ggtitle(methods[j])
  })
  grid.arrange(grobs = plist,ncol=3)
  dev.off()
}
pdf('correlation.pdf',width=15)
plist = list()
plist <- lapply(1:length(utypes),FUN = function(i){
  tempdata = data.frame(methods,correlation = cormat[utypes[i],] )
  corplot <- ggplot(tempdata,aes(methods,correlation))+geom_col()+ggtitle(utypes[i])+theme(axis.text.x = element_text(angle = 90, hjust = 1))
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()

pdf('RootedMSE.pdf',width=15)
plist = list()
plist <- lapply(1:length(utypes),FUN = function(i){
  tempdata = data.frame(methods, rootedMSE = sqrt(err[utypes[i],]) )
  corplot <- ggplot(tempdata,aes(methods, rootedMSE))+geom_col()+ggtitle(utypes[i])+theme(axis.text.x = element_text(angle = 90, hjust = 1))
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()

pdf('testmean_vs_trainmean.pdf',width = 20)
plist = list()
plist <- lapply(1:length(cellTypes),FUN=function(j){
  tempdata <- cbind(mu[,j],mu.gen[,j])
  colnames(tempdata) <- c('train_mean','test_mean')
  newplot <- ggplot(data = as.data.frame(tempdata),aes(x=train_mean,y=test_mean))+ xlim(0,1)+ylim(0,1)+geom_point()+geom_abline(intercept = 0,slope = 1) + ggtitle(cellTypes[j]) 
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()

ind <- sample(ncol(dat.gen[[1]]),10)
dir.create('check_individual')
for(i in 1:10){
  pdf(sprintf('check_individual/%s.pdf',colnames(dat.gen[[1]])[ind[i]]))
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
pdf('hist_MSE.pdf')
plist = list()
plist <- lapply(1:length(cell),FUN=function(j){
  newplot <- ggplot(data = samp2mean,aes_string(x=cell[j]))+
    geom_histogram(aes(y=..density..))+
    geom_density(alpha=.2,fill='#FF6666')+
    ggtitle(cell[j])
})
grid.arrange(grobs = plist, ncol = 3)
dev.off()
print(cormat)
print(err)
save(cormat,file = sprintf('cormat_%d.RData',perturbconst))
save(err,file = sprintf('err_%d,RData',perturbconst))
setwd('~/Hutch-Research/R_batch3')
samp_est <- colnames(dat.est)
save(samp_est,file = 'samp_est.RData')
samp_gen <- colnames(dat.gen[[1]])
save(samp_gen,file = 'samp_gen.RData')
save(deconv_expr,file = 'truerho.Rdata')
save(rho,file = 'methyrho.RData')
