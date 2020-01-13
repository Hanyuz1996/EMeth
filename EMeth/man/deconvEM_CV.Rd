\name{deconvEM_CV}
\alias{deconvEM_CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ~~ cross validation for EMeth algorithm ~~
}
\description{
  ~~ Automatic tuning for ridge penalty in EMeth algorithm. ~~
}
\usage{
deconvEM_CV(Y, eta, mu, aber = TRUE, V = "c", weight = matrix(1, 5, 5), pi_a_init, rho_init, nu0_init = rep(0, 50), sigma_c_init, lambda_init = 10, nu = penlist, folds = 5, usesubset = TRUE, maxiter = 50, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
 ~~ Bulk tissue sample for which cell type decomposition is to be estimated. Usually a matrix of size K*I, K is the number of probes and I is the number of samples. ~~
}
  \item{eta}{
    ~~ Tumor purity vector. A vector input of size I = number of samples. ~~
}
  \item{mu}{
    ~~ Reference data. A matrix of size K*Q, K is the number of probes and Q is the number of cell types. ~~
}
  \item{aber}{
     ~~Logic variable. Indicate if there are tumor cells or not.~~
}
  \item{V}{
     ~~ Variance structure, accept "c" = Constant variance, "b" = binomial variance structure, "w" = user defined weight. ~~
}
  \item{weight}{
    ~~ user defined weight of variance. ~~
}
  \item{pi_a_init}{
     ~~ A vector of size I. The probability of being abberant.~~
}
  \item{rho_init}{
     ~~ Initialization of a I*Q matrix, I = number of samples, Q = number of cell types.~~
}
  \item{nu0_init}{
     ~~ Initializati of tumor cell types.~~
}
  \item{sigma_c_init}{
     ~~ Initialization of overdispersion parameter for consistent probes. ~~
}
  \item{lambda_init}{
    ~~ Initialization of ratio of overdispersion parameter for aberrant probes versus consistent probes.~~
}
  \item{nu}{
     ~~ Tuning parameter for ridge penalty. ~~
}
  \item{folds}{
     ~~ Specify how many folds to use in cross validation. ~~
}
  \item{usesubset}{
    ~~ Logic variable. If TRUE, use a subset of the probes to perform the cross validation. ~~
}
  \item{maxiter}{
    ~~ Maximum time of iterations. ~~
}
  \item{verbose}{
     ~~Print out if there is anything printed during running.~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, eta, mu, aber = TRUE, V = "c", weight = matrix(1, 
    5, 5), pi_a_init, rho_init, nu0_init = rep(0, 50), sigma_c_init, 
    lambda_init = 10, nu = penlist, folds = 5, usesubset = TRUE, 
    maxiter = 50, verbose = FALSE) 
{
    if (length(nu) == 1) {
        result = deconvEM(Y, eta, mu, aber, V, weight = weight, 
            pi_a_init, rho_init, nu0_init, sigma_c_init, lambda_init, 
            nu = nu, maxiter)
        return(list(result, nu))
    }
    rdnumber <- runif(nrow(Y))
    losslist <- matrix(NA, nrow = folds, ncol = length(nu))
    for (i in 1:folds) {
        testind <- which(rdnumber < i/folds & rdnumber > (i - 
            1)/folds)
        trainind <- setdiff(1:nrow(Y), testind)
        sampind <- 1:ncol(Y)
        if (usesubset) {
            nsamp <- min(50, ncol(Y))
            sampind <- sample(ncol(Y), nsamp)
        }
        Y_train <- Y[trainind, sampind]
        Y_test <- Y[testind, sampind]
        pi_a_train <- pi_a_init[sampind]
        eta_train <- eta[sampind]
        mu_train <- mu[trainind, ]
        weight_train <- weight[trainind, sampind]
        mu_test <- mu[testind, ]
        rho_init_train <- rho_init[sampind, ]
        nu0_init_train <- nu0_init[trainind]
        nu0test <- rep(0, length(testind))
        for (j in 1:length(nu)) {
            temp <- deconvEM(Y_train, eta_train, mu_train, aber, 
                V, weight = weight_train, pi_a_train, rho_init_train, 
                nu0_init_train, sigma_c_init, lambda_init, nu = nu[j], 
                maxiter = 10)
            rhotest <- temp$rho
            pred <- mu_test \%*\% t(rhotest)
            Y_ab = Y_test - pred
            nu0test <- lapply(1:length(testind), FUN = function(l) {
                min(1, max(0, sum(eta_train * Y_ab[l, ])/sum(eta_train^2)))
            })
            nu0test <- unlist(nu0test)
            nu0test.m <- matrix(rep(nu0test, times = length(sampind)), 
                ncol = length(sampind), byrow = FALSE)
            losslist[i, j] <- mean((Y_ab - nu0test.m \%*\% diag(eta_train))^2)
            if (is.na(losslist[i, j])) {
                cat("Is this NA? rhotest", table(is.na(rhotest)), 
                  "\n")
                cat("Is this NA? nu0test", table(is.na(nu0test)), 
                  "\n")
                cat("Sum of eta", sum(eta_train))
            }
        }
    }
    avaloss <- apply(losslist, 2, mean, na.rm = TRUE)
    choosenu <- nu[which((avaloss) == min(avaloss))]
    print("finish cross validation")
    result = deconvEM(Y, eta, mu, aber, V, weight = weight, pi_a_init, 
        rho_init, nu0_init, sigma_c_init, lambda_init, nu = choosenu, 
        maxiter, verbose = verbose)
    return(list(result, choosenu, losslist))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
