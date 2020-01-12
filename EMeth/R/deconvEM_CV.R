deconvEM_CV <-
function(Y,eta,mu,aber = TRUE, V = 'c', weight = matrix(1,5,5),
                        pi_a_init, rho_init, nu0_init = rep(0,50), sigma_c_init, lambda_init = 10,
                        nu = penlist, folds = 5, usesubset = TRUE, maxiter=50, verbose = FALSE){
  if(length(nu) == 1){
    result = deconvEM(Y,eta,mu,aber,V,weight = weight, pi_a_init,rho_init,nu0_init,
                      sigma_c_init,lambda_init,nu=nu,maxiter)
    return(list(result,nu))
  }
  
  rdnumber <- runif(nrow(Y))
  losslist <- matrix(NA, nrow = folds, ncol = length(nu))
  for(i in 1:folds){
    testind <- which(rdnumber < i/folds & rdnumber > (i-1)/folds )
    trainind <- setdiff(1:nrow(Y),testind)
    sampind <- 1:ncol(Y)
    if(usesubset){
      nsamp <- min(50, ncol(Y))
      sampind <- sample(ncol(Y),nsamp)
    }
    
    Y_train <- Y[trainind,sampind]
    Y_test  <- Y[testind,sampind] 
    pi_a_train <- pi_a_init[sampind]
    eta_train <- eta[sampind]
    mu_train<- mu[trainind,]
    weight_train <- weight[trainind,sampind]
    mu_test <- mu[testind,]
    
    rho_init_train <- rho_init[sampind,]
    nu0_init_train <- nu0_init[trainind]
    nu0test <- rep(0,length(testind))

    for(j in 1:length(nu)){
    #  print('In cross validation')
      temp   <- deconvEM(Y_train,eta_train,mu_train,aber,V, weight = weight_train
                         ,pi_a_train,rho_init_train,nu0_init_train,
                         sigma_c_init,lambda_init,nu=nu[j],maxiter = 10)
      rhotest<- temp$rho
      pred   <- mu_test %*% t(rhotest) 
      Y_ab     = Y_test - pred
      
      nu0test <- lapply(1:length(testind),FUN=function(l){
        min(1,max(0,sum(eta_train * Y_ab[l,])/sum(eta_train^2)))
      })
      nu0test <- unlist(nu0test)
      nu0test.m <- matrix(rep(nu0test,times = length(sampind)),ncol = length(sampind),byrow = FALSE)
      losslist[i,j] <-  mean((Y_ab - nu0test.m %*% diag(eta_train))^2)
      if(is.na(losslist[i,j])){
        cat("Is this NA? rhotest", table(is.na(rhotest)),'\n')
        cat("Is this NA? nu0test", table(is.na(nu0test)),'\n')
        cat("Sum of eta", sum(eta_train))
      }
    }
  }
  
  avaloss <- apply(losslist,2,mean,na.rm = TRUE)
  choosenu <- nu[which((avaloss) == min(avaloss))]
  
  print('finish cross validation')
  result = deconvEM(Y,eta,mu,aber,V,weight = weight,pi_a_init,rho_init,nu0_init,
                    sigma_c_init,lambda_init,nu=choosenu,maxiter, verbose = verbose)
    
  return(list(result,choosenu,losslist))
}
