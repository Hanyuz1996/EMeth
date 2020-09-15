library(nnls)
library(data.table)
library(MASS)
library(quadprog)
library(ggplot2)
library(e1071)

dir = getwd()
source("EstimatingMethylation.R")
set.seed(1750905)

#---------------------------------------------------------
# Generate data from estimated distribution
#---------------------------------------------------------
gen_methy_beta <- function(mu,alpha,sample.size, pi=0, cellnum = 100, noise){
    cellNum = cellnum
    rho = matrix(NA,nrow = sample.size, ncol = length(cellTypes))
    beta = matrix(NA, nrow= nrow(mu), ncol = length(cellTypes))
    Y = V = matrix(0, nrow = nrow(mu),ncol = sample.size)

    rho   = matrix(0.2*runif(sample.size * length(cellTypes))-0.1,ncol = length(cellTypes))
    rho   = rho + matrix(rep(alpha,sample.size),nrow = sample.size)
    rho   = rho/rowSums(rho)
    eta   = runif(sample.size)*0.1
    rho   = diag(1-eta) %*% rho
    colnames(rho) = cellTypes
    
    nu0   = runif(nrow(mu))
    nu0.m = matrix(rep(nu0,times = sample.size),ncol = sample.size,byrow = FALSE)
   
    Y = mu %*% t(rho) + nu0.m %*% diag(eta)
    V = (mu*(1-mu)) %*% t(rho)
    V = V + matrix(rep(nu0*(1-nu0),sample.size),ncol=sample.size,byrow=FALSE) %*% diag(eta)
    V = V / cellNum
    
    pi.m = rep(runif(sample.size,pi-0.05,pi+0.05),nrow(Y))
    idic         = runif(n = nrow(mu)*sample.size) < pi.m
    mixerr       = rnorm(n = nrow(mu)*sample.size) * as.vector(t(sqrt(V)))
    mixerr[idic] = sqrt(noise) * mixerr[idic]
    mixerr       = matrix(mixerr,nrow = nrow(mu), ncol = sample.size, byrow = TRUE)
    
    list(bulk_sample = Y + mixerr, mix = rho, V = V,eta = eta, nu0 = nu0, idic = idic)
}

runsim <- function(simsize,simnoise,simpi,reptime,aber = FALSE, penalty = penalty, cellnum,maxiter = 50 ){
  
  genes2use = rownames(mu_gen)
  mean2pert = sample(genes2use, floor(simpi*length(genes2use)))
  mugenpert = mu_gen
  
  for(cpg in mean2pert){
    cell_temp = sample(cellTypes,7)
    mugenpert[cpg,cell_temp] = runif(7)
  }
  
  simdata   = gen_methy_beta(mugenpert,alpha,sample.size = simsize, pi = simpi, cellnum = cellnum,noise = simnoise)
  print("Data generated")
  Y         = simdata$bulk_sample
  eta       = simdata$eta
  nu0       = simdata$nu0
  rho.true  = simdata$mix
  idic.true = simdata$idic
  V.true    = simdata$V
  nu0cor    = rep(0,length(cellTypes))
  for(i in 1:length(cellTypes)){
     nu0cor[i] = cor(nu0,mu[,i])
  }

  sigma_c   = 0.1
  sigma_a   = simnoise * sigma_c
  
  #-------------------------------------------------------------
  # Estimating rho_qi by several methods: wls; ls; rls; qp
  #-------------------------------------------------------------
  methods = c("LaplaceEM","OriEM","svr","ls","rls","qp")
  rho     = array(data = NA, dim = c(nrow(rho.true),ncol(rho.true),length(methods)+1),
                  dimnames = list(1:nrow(rho.true),cellTypes,c(methods,"true")))
  err     = array(data = NA, dim = c(nrow(rho.true),ncol(rho.true),length(methods)),
                  dimnames = list(1:nrow(rho.true),cellTypes,methods))
  
  rho[,,length(methods)+1]=rho.true
  
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
  
  #temp <- apply(s2,1,max)
  #temp <- temp/median(temp)
  #lb <- quantile(temp,0.15)
  #ub <- quantile(temp,0.85)
  #temp[temp<lb] <- lb
  #temp[temp>ub] <- ub
  #W <- matrix(rep(temp,ncol(Y)),ncol = ncol(Y),byrow = FALSE)

  print('LaplaceEM')
  hundrediter_laplace = cv.emeth(Y,eta,mu,aber = aber, V='c', init = 'default',
                                 family = 'laplace', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
  rho[,,'LaplaceEM'] = hundrediter_laplace[[1]]$rho
  
  print('OriEM')
  hundrediter = cv.emeth(Y,eta,mu,aber = aber, V='c', init = 'default',
                         family = 'normal', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
  rho[,,'OriEM'] = hundrediter[[1]]$rho
        
    for(k in 1:length(methods)){
      err[,,k] = (rho[,,k]-rho.true)^2
    }
    err.mean = sqrt(apply(err,c(2,3),mean))

    cormat <- matrix(NA,nrow = length(cellTypes),ncol = length(methods))
    colnames(cormat) <- methods
    rownames(cormat) <- cellTypes

    rss <- matrix(NA,nrow = length(cellTypes),ncol = length(methods))
    colnames(rss) <- methods
    rownames(rss) <- cellTypes

    print(table(is.na(rho.true)))
    print(table(is.na(rho)))

    for(i in 1:length(cellTypes)){{
      cormat[i,] <- sapply(1:length(methods), FUN = function(j){{
      cor(rho[,i,methods[j]],rho.true[,i])
    }})
      rss[i,] <- sapply(1:length(methods), FUN = function(j){{
      temp <- lm(rho.true[,i]~rho[,i,methods[j]])
      return(sum(temp$residuals^2))
    }})
}}
#  dir.create(sprintf('../figures/%s',Sys.Date()))
#  dir.create(sprintf('../figures/%s/%d',Sys.Date(),reptime))
   
#  for(i in 1:length(cellTypes)){{
#  pdf(sprintf('../figures/%s/%d/%s-%s-%s-estimate-vs-true.pdf',Sys.Date(),reptime,cellnum, simnoise, cellTypes[i]))
#  plist = list()
#  plist <- lapply(1:length(methods), FUN = function(j){{
#    tempdata = cbind(res.rho[,i,methods[j]],rho.true[,i])
#    colnames(tempdata) <- c("methylation","true")
#    newplot <- ggplot(data = as.data.frame(tempdata), aes(x=methylation,y=true))+ xlim(0,0.3) + ylim(0,0.3) +
#      geom_point() + geom_abline(intercept = 0,slope = 1) + ggtitle(methods[j])
#  }})
#  grid.arrange(grobs = plist,ncol=2)
#  dev.off()
#}}

#print(nu0cor)
  
    list(rho = rho, err = err.mean, cor = cormat,
         laplace = hundrediter_laplace[[1]],
         ori = hundrediter[[1]], 
         true = simdata,
         rss = rss,
         nu0cor = nu0cor)
}






