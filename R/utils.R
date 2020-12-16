deconv.init <-
  function(Y,eta,mu,aber = TRUE){
    Q <- ncol(mu)
    simsize <- ncol(Y)
    temp       = runif(simsize * Q) * 0.2 -0.1
    rho_init   = matrix(0.5,ncol = Q, nrow = ncol(Y))
    nu0_init   = runif(nrow(Y))
    sigma_c_init = 0.1
    lambda_init  = 2
    pi_a_init   = rep(0.5,simsize)
    
    for(j in 1:ncol(Y)){
      if(j %% 50 == 0){ cat(j, date(), "\n") }
      y    = Y[,j]
      X    = as.data.frame(mu)
      Xmat = mu
      
      temp            = lm(y ~ .-1,data = X)$coefficients
      temp[temp < 0]  = 0
      rho_init[j,]    = (1-eta[j])*temp/sum(temp)
    }
    
    K = nrow(Y)
    Y_ab     = Y - mu %*% t(rho_init)
    if(aber){
      for(k in 1:K){
        nu0_init[k] = min(1,max(0,sum(eta * Y_ab[k,])/sum(eta^2)))
      }
    }
    
    init.paras <- list(weight = matrix(1,5,5),
                       pi_a_init = pi_a_init, rho_init = rho_init, nu0_init = nu0_init, 
                       sigma_c_init=0.1 , lambda_init = 10)
    return(init.paras)
  }

deconvEM <-
  function(Y,eta,mu,aber = TRUE, V = 'c', weight = matrix(1,5,5),
           pi_a_init, rho_init, nu0_init = rep(0,50), sigma_c_init=0.1 , lambda_init = 10,
           nu = 0, maxiter=50, verbose = FALSE ){
    
    Q         = ncol(mu)
    K         = nrow(mu)
    I         = ncol(Y)
    X         = as.data.frame(mu) 
    
    #Initialization
    pi_a      = pi_a_init
    
    if(! aber){
      nu0_init = rep(0,K)
      eta = rep(0,I)
    }
    nu0       = nu0_init
    rho       = rho_init
    sigma_c   = sigma_c_init
    
    if(V == 'b'){
      mu00 = mu*(1-mu)
      mu00[mu00 < 0.05] = 0.05 
      Vf = mu00 %*% t(rho_init)
      Vf = Vf + matrix(rep(nu0_init*(1-nu0_init),ncol(Y)),ncol=ncol(Y),byrow=FALSE) %*% diag(eta)
      Vc = sigma_c_init * sigma_c_init * Vf
      Va = lambda_init * Vc
    }else if(V == 'c'){
      Vc = matrix(sigma_c_init^2, K,I)
      Va = matrix(lambda_init * sigma_c_init^2, K,I)
    }else if(V == 'w'){
      Vf = weight 
      Vc = Vf
      Va = lambda_init * Vf
    }
    
    lambda    = lambda_init
    nu0       = nu0_init
    nu0.m     = matrix(rep(nu0,times = I),ncol = I,byrow = FALSE)
    Mu        = mu %*% t(rho) + nu0.m  %*% diag(eta)
    
    for(iter in 1:maxiter){
      
      if(iter %% 10 == 0 & verbose){
        cat("-------------------\n")
        cat(iter, date(), "\n") 
      }
      
      # Save old
      pi_a_old = pi_a
      rho_old  = rho
      Vc_old    = Vc
      Va_old    = Va
      lambda_old = lambda
      nu0_old  = nu0
      Va[Va < 0.001] = 0.001
      Vc[Vc < 0.001] = 0.001
      
      # Truncation to avoid extreme values
      e = Y - Mu
      e[which(e > 0.99)] = 0.99
      e[which(e < -0.99)] = -0.99
      
      # E step
      pdf_a    =  dnorm(e/sqrt(Va)) %*% diag(pi_a) / sqrt(Va)
      pdf_c    =  dnorm(e/sqrt(Vc)) %*% diag(1-pi_a) /sqrt(Vc)
      
      gamma    = pdf_a/(pdf_a+pdf_c)
      
      # M step
      V_renorm = Vc/norm(Vc)
      W        = ( gamma + lambda*(1 - gamma))/(lambda * V_renorm) 
      
      pi_a     = colMeans(gamma)
      Y_c      = Y - nu0.m %*% diag(eta)
      nulist   = rep(0,I)
      
      for(i in 1:I){
        C      = norm(t(mu) %*% diag(W[,i]) %*% mu)
        A      = rbind(rep(-1,Q),diag(rep(1,Q)))
        b      = c(-1+eta[i],rep(0,Q))*sqrt(C)
        D      = (t(mu) %*% diag(W[,i]) %*% mu + nu * diag(ncol(mu)) )/(C)
        d      = t(t(mu) %*% diag(W[,i]) %*% Y_c[,i])/sqrt(C)
        temp   = solve.QP(D,d,t(A),b,meq = 1)$solution/sqrt(C)
        rho[i,]= temp
      }
      
      Y_ab     = Y - mu %*% t(rho)
      if(aber){
        for(k in 1:K){
          nu0[k] = min(1,max(0,sum(W[k,]*eta * Y_ab[k,])/sum(W[k,] * eta^2)))
        }
      }  
      nu0.m    = matrix(rep(nu0,times = I),ncol = I,byrow = FALSE)
      Mu      = mu %*% t(rho) + nu0.m %*% diag(eta)
      
      if(V == 'b'){
        Vf = mu00 %*% t(rho)
        Vf = Vf + matrix(rep(nu0*(1-nu0),ncol(Y)),ncol=ncol(Y),byrow=FALSE) %*% diag(eta)
        Vc = sigma_c * sigma_c * Vf
        Va = lambda * Vc
      }else if(V == 'c'){
        Vf = matrix(1, K,I)
      }else if(V == 'w'){
        Vf = weight
      }
      
      sigma_a = (sum(gamma*((Y-Mu)^2)/Vf)/sum(gamma))
      sigma_c = (sum((1-gamma)*(Y-Mu)^2/Vf)/sum(1-gamma))
      
      lambda = sigma_a/sigma_c
      
      Va = Vf * sigma_a
      Vc = Vf * sigma_c
      
      # Check convergence
      if(max(abs(rho_old-rho))<1e-4){break}
    }
    
    list(rho = rho, sigma_c = sigma_c, lambda = lambda, nu0 = nu0, pi_a = pi_a, gamma = gamma, weights = W, iter = iter)
  }


deconvEM_laplace <-
  function(Y,eta,mu,aber = TRUE, V = 'c', weight = matrix(1,5,5),
           pi_a_init, rho_init, nu0_init = rep(0,50), sigma_c_init, lambda_init = 10,
           nu = 0, maxiter=50, verbose = FALSE ){
    
    Q         = ncol(mu)
    K         = nrow(mu)
    I         = ncol(Y)
    X         = as.data.frame(mu) 
    
    #Initialization
    pi_a      = pi_a_init
    
    if(! aber){
      nu0_init = rep(0,K)
      eta  = rep(0,I)
    }
    nu0       = nu0_init
    rho       = rho_init
    sigma_c   = sigma_c_init
    
    if(V == 'b'){
      mu00 = mu*(1-mu)
      mu00[mu00 < 0.05] = 0.05 
      Vf = mu00 %*% t(rho_init)
      Vf = Vf + matrix(rep(nu0_init*(1-nu0_init),ncol(Y)),ncol=ncol(Y),byrow=FALSE) %*% diag(eta)
      Vc = sigma_c_init * Vf
      Va = lambda_init * Vc
    }else if(V == 'c'){
      Vc = matrix(sigma_c_init, K,I)
      Va = matrix(lambda_init * sigma_c_init, K,I)
    }else if(V == 'w'){
      Vf = weight 
      Vc = Vf
      Va = lambda_init * Vf
    }
    
    lambda    = lambda_init
    nu0       = nu0_init
    nu0.m     = matrix(rep(nu0,times = I),ncol = I,byrow = FALSE)
    Mu        = mu %*% t(rho) + nu0.m  %*% diag(eta)
    
    for(iter in 1:maxiter){
      
      if(iter %% 10 == 0 & verbose){
        cat("-------------------\n")
        cat(iter, date(), "\n") 
      }
      
      # Save old
      pi_a_old = pi_a
      rho_old  = rho
      Vc_old    = Vc
      Va_old    = Va
      lambda_old = lambda
      nu0_old  = nu0
      Va[Va < 0.001] = 0.001
      Vc[Vc < 0.001] = 0.001
      
      # Truncation to avoid extreme values
      e = Y - Mu
      e[which(e > 0.99)] = 0.99
      e[which(e < -0.99)] = -0.99
      
      # E step
      pdf_a    = dlaplace(e/Va) %*% diag(pi_a) / Va
      pdf_c    = dlaplace(e/Vc) %*% diag(1-pi_a)/Vc
      
      gamma    = pdf_a/(pdf_a+pdf_c)
      #print(dim(gamma))
      
      # M step
      V_renorm = Vc/norm(Vc)
      W        = ( gamma + lambda*(1 - gamma))/(lambda * V_renorm) 
      #print(dim(W))
      
      pi_a     = colMeans(gamma)
      Y_c      = Y - nu0.m %*% diag(eta)
      nulist   = rep(0,I)
      
      for(i in 1:I){
        y <- Y_c[,i]
        X <- mu
        
        opt <- optim(rho[i,],fn=function(a){
          sum(abs(y-X %*% a)*W[,i])+nu*sum(a^2)},gr = 'Brent')
        temp = opt$par
        temp[temp < 0] = 0
        
        rho[i,]= (1-eta[i])*temp/sum(temp)
      }
      
      Y_ab     = Y - mu %*% t(rho)
      if(aber){
        for(k in 1:K){
          y <- Y_ab[k,]
          opt <- optim(nu0[k],fn=function(a){
            sum(abs(y-eta * a)*W[k,])},gr = 'Brent')
          temp = opt$par
          nu0[k] = max(0,min(1,temp))
        }
      }  
      nu0.m    = matrix(rep(nu0,times = I),ncol = I,byrow = FALSE)
      
      Mu      = mu %*% t(rho) + nu0.m %*% diag(eta)
      
      if(V == 'b'){
        Vf = mu00 %*% t(rho)
        Vf = Vf + matrix(rep(nu0*(1-nu0),ncol(Y)),ncol=ncol(Y),byrow=FALSE) %*% diag(eta)
        Vc = sigma_c * sigma_c * Vf
        Va = lambda * Vc
      }else if(V == 'c'){
        Vf = matrix(1, K,I)
      }else if(V == 'w'){
        Vf = weight %*% t(rho^2)
        Vf = Vf + matrix(rep(nu0*(1-nu0),ncol(Y)),ncol=ncol(Y),byrow=FALSE) %*% diag(eta)
      }
      
      sigma_a = (sum(gamma*(abs(Y-Mu))/Vf)/sum(gamma))
      sigma_c = (sum((1-gamma)*(abs(Y-Mu))/Vf)/sum(1-gamma))
      
      lambda = sigma_a/sigma_c
      
      Va = Vf * sigma_a
      Vc = Vf * sigma_c
      
      # Check convergence
      if(max(abs(rho_old-rho))<1e-4){break}
    }
    
    list(rho = rho, sigma_c = sigma_c, lambda = lambda,nu0 = nu0, pi_a = pi_a, gamma = gamma, weights = W, iter = iter)
  }

deconvEM_CV <-
  function(Y,eta,mu,aber = TRUE, V = 'c', weight = matrix(1,5,5),
           pi_a_init, rho_init, nu0_init = rep(0,50), sigma_c_init, lambda_init = 10,
           nu, folds = 5, usesubset = TRUE, maxiter=50, verbose = FALSE){
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
      mu_test <- mu[testind,]
      
      rho_init_train <- rho_init[sampind,]
      nu0_init_train <- nu0_init[trainind]
      nu0test <- rep(0,length(testind))
      
      for(j in 1:length(nu)){
        print('In cross validation')
        if(verbose){
          cat("penalty parameter", nu[j])
        }
        temp   <- deconvEM(Y_train,eta_train,mu_train,aber,V, 
                           ,pi_a_train,rho_init_train,nu0_init_train,
                           sigma_c_init,lambda_init,nu=nu[j],maxiter = 10,verbose = verbose)
        rhotest<- temp$rho
        pred   <- mu_test %*% t(rhotest) 
        Y_ab     = Y_test - pred
        
        nu0test <- lapply(1:length(testind),FUN=function(l){
          min(1,max(0,sum(eta_train * Y_ab[l,])/sum(eta_train^2)))
        })
        nu0test <- unlist(nu0test)
        nu0test.m <- matrix(rep(nu0test,times = length(sampind)),ncol = length(sampind),byrow = FALSE)
        losslist[i,j] <-  mean((Y_ab - nu0test.m %*% diag(eta_train))^2)
      }
    }
    
    avaloss <- apply(losslist,2,mean,na.rm = TRUE)
    choosenu <- nu[which((avaloss) == min(avaloss))]
    
    print('finish cross validation')
    result = deconvEM(Y,eta,mu,aber,V,weight = weight,pi_a_init,rho_init,nu0_init,
                      sigma_c_init,lambda_init,nu=choosenu,maxiter, verbose = verbose)
    
    return(list(result,choosenu,losslist))
  }
deconvEM_CV_laplace <-
  function(Y,eta,mu,aber = TRUE, V = 'c', weight = matrix(1,5,5),
           pi_a_init, rho_init, nu0_init = rep(0,50), sigma_c_init, lambda_init = 10,
           nu, folds = 5, usesubset = TRUE, maxiter=50, verbose = FALSE){
    if(length(nu) == 1){
      result = deconvEM_laplace(Y,eta,mu,aber,V,weight = weight, pi_a_init,rho_init,nu0_init,
                                sigma_c_init,lambda_init,nu=nu,maxiter,verbose)
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
      #weight_train <- weight[trainind,sampind]
      mu_test <- mu[testind,]
      
      rho_init_train <- rho_init[sampind,]
      nu0_init_train <- nu0_init[trainind]
      nu0test <- rep(0,length(testind))
      
      for(j in 1:length(nu)){
        print('In cross validation')
        if(verbose){
          cat("penalty parameter", nu[j])
        }
        temp   <- deconvEM_laplace(Y_train,eta_train,mu_train,aber,V,
                                   ,pi_a_train,rho_init_train,nu0_init_train,
                                   sigma_c_init,lambda_init,nu=nu[j],maxiter = 10,verbose)
        rhotest<- temp$rho
        pred   <- mu_test %*% t(rhotest) 
        Y_ab     = Y_test - pred
        
        nu0test <- lapply(1:length(testind),FUN=function(l){
          min(1,max(0,sum(eta_train * Y_ab[l,])/sum(eta_train^2)))
        })
        nu0test <- unlist(nu0test)
        nu0test.m <- matrix(rep(nu0test,times = length(sampind)),ncol = length(sampind),byrow = FALSE)
        losslist[i,j] <-  mean((Y_ab - nu0test.m %*% diag(eta_train))^2)
      }
    }
    
    avaloss <- apply(losslist,2,mean,na.rm = TRUE)
    choosenu <- nu[which((avaloss) == min(avaloss))]
    
    print('finish cross validation')
    result = deconvEM_laplace(Y,eta,mu,aber,V,weight = weight,pi_a_init,rho_init,nu0_init,
                              sigma_c_init,lambda_init,nu=choosenu,maxiter, verbose = verbose)
    
    return(list(result,choosenu,losslist))
  }

dlaplace <-
  function(y,m=0,s=1){
    return(exp(-abs(y-m)/s)/(2*s))
  }
