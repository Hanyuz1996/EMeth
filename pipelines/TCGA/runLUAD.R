library(data.table)
library(scales)
library(stringr)
library(quadprog)
library(e1071)
library(ggplot2)
library(MASS)
library(gridExtra)
library(quantreg)
setwd("~/test/pipelines/TCGA")
#source('ICeDT.R')
#source('Expr_Methy_Data_LUAD.R')

use_abs_eta = TRUE
purity_correction =  TRUE
discard_high_purity = FALSE
discard_mast_cells = FALSE 
testalgorithm = FALSE
printplot = FALSE
aber = TRUE

#-----------------------------------------------------------
# Compare absolute purity and estiamted purity
#-----------------------------------------------------------
if(printplot){
  print(cor(abs_purity$absolute_extract_purity,est_purity$AscatPurity))
  pdf("~/Hutch-Research/figures/abs_vs_est_purity.pdf", width = 13.5, height = 6)
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
  source("LUAD.R")
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
  penalty = dim(Y)[1]*(10^seq(-2,1,1)) 
  #penalty  
  pi_a_init = rep(0.5,ncol(Y))
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
methods = c("LaplaceEM","OriEM","svr","ls","rls","qp")
rho     = array(data = NA, dim = c(ncol(Y),length(cellTypes),length(methods)),
                dimnames = list(1:ncol(Y),cellTypes,methods))
alpha = rep(1/length(cellTypes),length(cellTypes))
simsize = ncol(Y)
C = c(0.1,1/sqrt(10),1,sqrt(10),10)

for(j in 1:ncol(Y)){
  if(j %% 10 == 0){ cat(j, date(), "\n") }
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
  temp        = (t(svrmodel$coefs) %*% svrmodel$SV)
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

print('LaplaceEM')
hundrediter_laplace = cv.emeth(Y,eta,mu,aber = aber, V='c', init = 'default',
                               family = 'laplace', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
rho[,,'LaplaceEM'] = hundrediter_laplace[[1]]$rho

print('OriEM')
hundrediter = cv.emeth(Y,eta,mu,aber = aber, V='c', init = 'default',
                       family = 'normal', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
rho[,,'OriEM'] = hundrediter[[1]]$rho

####################################################################
# Suppose to run ICeDt
#Yexpr  = fread('~/Hutch-Research/Data/Real/cleaned_expression_LUAD_data.txt')
#muexpr = fread('~/Hutch-Research/Data/Real/LM22.txt')
#samname <- colnames(Yexpr)
#for(i in 1:length(samname)){samname[i] <- str_replace(samname[i],"X","")}
#colnames(Yexpr) <- samname
#genesymbol <- unlist(Yexpr[,1])
#Yexpr = subset(Yexpr,select = colnames(Y))
#rownames(Yexpr) <- genesymbol
#gene = intersect(genesymbol,muexpr$'Gene symbol')
#Yexpr = as.matrix(exp(subset(Yexpr, rownames(Yexpr)  %in%  gene)))
#muexpr = as.matrix(subset(muexpr,unlist(muexpr$'Gene symbol') %in% gene)[,-1])
#rownames(Yexpr) = rownames(muexpr) = gene
#icedt_LUAD = ICeDT(Y=as.matrix(Yexpr),Z=muexpr, tumorPurity = eta, refVar = NULL)

setwd('~/')
rho_LUAD = rho
deconv_expr_LUAD = deconv_expr
save(rho_LUAD, file = 'rho_inf_LUAD.RData')
save(deconv_expr_LUAD, file = 'deconv_expr_LUAD.RData')
#save(icedt_LUAD, file = 'icedt_LUAD.RData')
dir.create('~/Hutch-Research/figures/MethyVSExpr/LUAD')
setwd("~/Hutch-Research/figures/MethyVSExpr/LUAD")

utypes = intersect(cellTypes,colnames(deconv_expr))

cormat <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- utypes

err <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(err) <- methods
rownames(err) <- utypes

rss <- matrix(NA,nrow = length(utypes),ncol = length(methods))
colnames(rss) <- methods
rownames(rss) <- utypes


for(i in 1:length(utypes)){
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){
    cor(rho[,utypes[i],methods[j]],deconv_expr[,utypes[i]])
  })
  err[i,] <- sapply(1:length(methods), FUN = function(j){
    sqrt(mean((rho[,utypes[i],methods[j]] - deconv_expr[,utypes[i]])^2))
  }) 
  rss[i,] <- sapply(1:length(methods), FUN = function(j){
    temp <- lm(deconv_expr[,utypes[i]]~rho[,utypes[i],methods[j]])
    return(sum(temp$residuals^2))
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
print(rss)

OneMinusCorr <- 1-matrix(cormat,ncol = 1, byrow = FALSE)
RMSE <- matrix(err,ncol = 1, byrow = FALSE )
CellType <- rep(cellTypes,length(methods))
Methods <- rep(methods,each = length(cellTypes))
res <- cbind.data.frame(OneMinusCorr,RMSE,CellType,Methods)
pdf('Comparison.pdf')
complot<- ggplot(res, aes(x=OneMinusCorr,y=RMSE, color =  Methods))+ggtitle('LUAD')+geom_point()+scale_y_continuous(trans = log2_trans(),
                                                                                                                    breaks = trans_breaks('log10',function(x) 10^x),
                                                                                                                    labels = trans_format('log10',math_format(10^.x)))+
  geom_text(label = res[,3])+xlim(0.1,1.05)
print(complot)
dev.off()
cormat_LUAD <- cormat
err_LUAD <- err
save(cormat_LUAD,file='cormat_LUAD.RData')
save(err_LUAD,file = 'err_LUAD.RData')

quit(save = 'no')

