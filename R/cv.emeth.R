cv.emeth <-
function(Y,eta,mu,aber,V, init = 'default', nu, family = 'laplace',
                     folds = 5, usesubset = TRUE, maxiter = 50, verbose = FALSE){
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
    result = deconvEM_CV_laplace(Y,eta,mu,aber = aber, V = 'c', weight,
                                 pi_a_init, rho_init, nu0_init, sigma_c_init, lambda_init,
                                 nu, folds, usesubset, maxiter, verbose = verbose)
  }
  else if(family == 'normal'){
    result = deconvEM_CV(Y,eta,mu,aber = aber, V = 'c', weight,
                         pi_a_init, rho_init, nu0_init, sigma_c_init, lambda_init,
                         nu, folds, usesubset, maxiter, verbose = verbose)
  }
  else{
    stop("Must specify the family from laplace or normal!")
  }
  
  if(median(result[[1]]$pi_a > 0.5)){
    warning("maximum proportion of aberrant probes is large, suggesting inappropriate references\n")
  }
  
  return(result)
}
