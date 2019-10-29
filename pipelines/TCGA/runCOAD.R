library(data.table)
library(scales)
library(stringr)
library(quadprog)
library(e1071)
library(ggplot2)
library(MASS)
library(gridExtra)
library(quantreg)
setwd("~/Hutch-Research/R_batch2")
#source('Expr_Methy_Data_SKCM.R')

use_abs_eta = TRUE
purity_correction =  TRUE
discard_high_purity = FALSE
discard_mast_cells = FALSE 
testalgorithm = FALSE
printplot = TRUE
aber = TRUE

#-----------------------------------------------------------
# Compare absolute purity and estiamted purity
#-----------------------------------------------------------
if(1<0){
print(cor(abs_purity$absolute_extract_purity,est_purity$AscatPurity))
pdf("../../../../figures/abs_vs_est_purity.pdf", width = 13.5, height = 6)
eta_plot = data.frame(cbind(abs_purity$absolute_extract_purity,est_purity$AscatPurity))
colnames(eta_plot) = c('abs_purity','inferred_purity')
etascatter <- ggplot(data = eta_plot, aes(x=abs_purity, y=inferred_purity)) + xlim(0,1) + ylim(0,1) +
                     geom_point() + geom_abline(intercept = 0,slope = 1)
print(etascatter)
dev.off()
}

#-----------------------------------------------------------
# Compare results with different methods and with expression data
#-----------------------------------------------------------
#Y = 2^Y/(2^Y+1)

if(use_abs_eta){
  print("Absolute Purity")
  setwd("~/Hutch-Research/R_batch2")
  source("COAD.R")
  eta = eta_abs[which(colnames(ys) %in% rownames(deconv_expr))]
  Y   = ys[,intersect(colnames(ys),rownames(deconv_expr))]
  ref = mu[rownames(Y),]
  #printplot = FALSE
  if(purity_correction){
    temp <- rownames(deconv_expr)
    deconv_expr <- diag(1-eta) %*% deconv_expr
    rownames(deconv_expr) <- temp
  }
  ref[ref < 0.05] = 0.05
  ref[ref > 0.95] = 0.95
  mu = ref
  #penalty = dim(Y)[1]*(10^seq(-2,1,1)) 
  penalty = 1000
}

if(discard_high_purity){
  high = which(eta > 0.8)
  Y = Y[,-high]
  eta = eta[-high]
  deconv_expr = deconv_expr[-high,]
}

if(testalgorithm){
  print("test algorithm")
  Y = ref %*% t(deconv_expr[,1:7]) + 
       matrix(rnorm(nrow(Y)*ncol(Y),0,sd=0.05),nrow = nrow(Y), ncol = ncol(Y))
  eta = rep(0.001,ncol(Y))
}

if(discard_mast_cells){
  cat('dim of Y before discarding other cells', dim(Y))
  print(summary(other))
  high = which(other > quantile(other,0.3))
  cat('30% quantile of other cell types', quantile(other,0.3))
  cat('length of high: ', length(high))
  Y = Y[,-high]
  eta = eta[-high]
  deconv_expr = deconv_expr[-high,]
  print(dim(Y))
}

eta[eta > 0.99] = 0.99

print(dim(ref))
methods = c("MaxVarEM","LaplaceEM","BinomEM","OriEM","svr","ls","rls","qp")
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
    
    svrmodel        = svm(y ~ ., data = X, kernel = 'linear', cost = c(1,sqrt(10),10,10^{3/2},100),cross = 5)
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

try({
  hundrediter_laplace = deconvEM_laplace(Y,eta,mu,aber = aber, V='c', weight = matrix(1,5,5),
                                0.5,rho_init,nu0_init,sigma_c_init = 0.1, lambda_init = 2,
                                nu = penalty, maxiter = 100, verbose = TRUE)
  rho[,,'LaplaceEM'] = hundrediter_laplace$rho
  sigma_c_est_laplace   = hundrediter_laplace$sigma_c
  pi_est_laplace      = hundrediter_laplace$pi_a
  nu0_est_laplace     = hundrediter_laplace$nu0
  iter_laplace        = hundrediter_laplace$iter
  gamma_laplace       = hundrediter_laplace$gamma
})

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
  rho[,,'MaxVarEM'] = hundrediter_weight$rho
  sigma_c_est_weight   = hundrediter_weight$sigma_c
  pi_est_weight      = hundrediter_weight$pi_a
  nu0_est_weight     = hundrediter_weight$nu0
  iter_weight        = hundrediter_weight$iter
  gamma_weight       = hundrediter_weight$gamma
})

try({
  hundrediter_het = deconvEM(Y,eta,mu,aber = aber, V='b', weight = matrix(1,5,5),
                                0.5,rho_init,nu0_init,sigma_c_init = 0.1, lambda_init = 2,
                                nu = penalty, maxiter = 100)
  rho[,,'BinomEM'] = hundrediter_het$rho
  sigma_c_est_het       = hundrediter_het$sigma_c
  pi_est_het      = hundrediter_het$pi_a
  nu0_est_het     = hundrediter_het$nu0
  iter_het        = hundrediter_het$iter
  gamma_het       = hundrediter_het$gamma
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
setwd('~/Hutch-Research/R_batch2')
rho_COAD = rho
deconv_expr_COAD = deconv_expr
save(rho_COAD, file = 'rho_inf_COAD.RData')
save(deconv_expr_COAD, file = 'deconv_expr_COAD.RData')

dir.create('~/Hutch-Research/figures/MethyVSExpr/COAD')
setwd("~/Hutch-Research/figures/MethyVSExpr/COAD")

utypes = intersect(cellTypes,colnames(deconv_expr))

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
  pdf(sprintf('%s_express_methy_discard_high_purity.pdf',utypes[i]))
  plist = list()
  plist <- lapply(1:length(methods), FUN = function(j){
    tempdata = cbind(rho[,utypes[i],methods[j]],deconv_expr[,utypes[i]],eta)
    colnames(tempdata) <- c("methylation","expression","eta")
    newplot <- ggplot(data = as.data.frame(tempdata), aes(x=methylation,y=expression,color=eta))+ xlim(0,0.3) + ylim(0,0.3) +
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
  
print(cormat)
print(err)

colnames(cormat)[3] <- 'BinomEM'
colnames(err)[3] <- 'BinomEM'
methods <- c('BinomEM','OriEM','ls','svr')
cormat <- subset(cormat,select = methods)
err <- subset(err,select = methods)
methods <- factor(methods,levels = methods)

OneMinusCorr <- 1-matrix(cormat,ncol = 1, byrow = FALSE)
RMSE <- matrix(err,ncol = 1, byrow = FALSE )
CellType <- rep(cellTypes,length(methods))
Methods <- rep(methods,each = length(cellTypes))
res <- cbind.data.frame(OneMinusCorr,RMSE,CellType,Methods)
pdf('Comparison.pdf')
complot<- ggplot(res, aes(x=OneMinusCorr,y=RMSE, color =  Methods))+ggtitle('COAD')+geom_point()+scale_y_continuous(trans = log2_trans(),
    breaks = trans_breaks('log10',function(x) 10^x),
    labels = trans_format('log10',math_format(10^.x)))+
    geom_text(label = res[,3])+xlim(0.1,1.05)
print(complot)
dev.off()
cormat_COAD <- cormat
err_COAD <- err
save(cormat_COAD,file='cormat_COAD.RData')
save(err_COAD,file = 'err_COAD.RData')

quit(save = 'no')
