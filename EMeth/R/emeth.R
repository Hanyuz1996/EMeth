emeth <-
function(Y,eta,mu,aber,V, init = 'default', family = 'laplace',
                  nu = 0, maxiter = 50, verbose = FALSE){
    if(init == 'default'){
      initlist = deconv.init(Y,eta,mu,aber)
    } 
    else{
      initlist = init
    }
    weight = initlist$weight
    pi_a_init = initlist$pi_a_init
    rho_init = initlist$rho_init
    nu0_init = initlist$nu0_init
    sigma_c_init = initlist$sigma_c_init
    lambda_init = initlist$lambda_init
    
    if(family == 'laplace'){
      result = deconvEM_laplace(Y,eta,mu,aber,V,weight = weight,pi_a_init,rho_init,nu0_init,
                        sigma_c_init,lambda_init,nu=nu, maxiter = maxiter, verbose = verbose)
      return(result)
    }
    else if(family == 'normal'){
      result = deconvEM(Y,eta,mu,aber,V,weight = weight,pi_a_init,rho_init,nu0_init,
                        sigma_c_init,lambda_init,nu=nu, maxiter = maxiter, verbose = verbose)
      return(result)
    }
    else{
      stop("Must specify the family from laplace or normal!")
    }
}
