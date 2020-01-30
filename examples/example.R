# load example data
load('./data/etaexample.RData') #eta: tumor purity
load('./data/muexample.RData') #mu: reference data
load('./data/nu0example.RData') #nu0: true methylation of unknown cell type 
load('./data/Yexample.RData') #Y: sample methylation
load('./data/mixexample.RData') #rho.true: true proportion of each cell type
penalty= (dim(mu)[1])*(10^seq(-2,1,0.5)) 
source('../source/_lib.R')
cellTypes = colnames(mu)
  
print('LaplaceEM')
hundrediter_laplace = cv.emeth(Y,eta,mu,aber = TRUE, V='c', init = 'default',
                               family = 'laplace', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
rho.laplace = hundrediter_laplace[[1]]$rho

print('NormalEM')
hundrediter = cv.emeth(Y,eta,mu,aber = TRUE, V='c', init = 'default',
                       family = 'normal', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
rho.normal = hundrediter[[1]]$rho