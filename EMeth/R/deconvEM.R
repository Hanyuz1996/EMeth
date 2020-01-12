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
    cat(pi_a, "\n")
    cat(nu0[1:5],"\n")
    cat(rho[1,],"\n")
    cat(lambda,"\n")
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
    
    if(verbose){
        cat(sprintf("%s NAs in W",length(which(is.na(W)))),"\n")
    }

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
    #print(nu0[1:5])
    nu0.m    = matrix(rep(nu0,times = I),ncol = I,byrow = FALSE)
    
    #cat('dim of mu', dim(mu), '\n')
    #cat('dim of rho', dim(rho),'\n')

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
  cat('sigma_c_square', sigma_c,'\n')
  
  list(rho = rho, sigma_c = sigma_c, lambda = lambda, nu0 = nu0, pi_a = pi_a, gamma = gamma, weights = W, iter = iter)
}
