deconv.init <-
function(Y,eta,mu,aber = TRUE){
  temp       = runif(simsize * length(cellTypes)) * 0.2 -0.1
  rho_init   = matrix(0.5,ncol = length(cellTypes), nrow = ncol(Y))
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
