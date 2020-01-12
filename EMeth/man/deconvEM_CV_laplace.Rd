\name{deconvEM_CV_laplace}
\alias{deconvEM_CV_laplace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deconvEM_CV_laplace(Y, eta, mu, aber = TRUE, V = "c", weight = matrix(1, 5, 5), pi_a_init, rho_init, nu0_init = rep(0, 50), sigma_c_init, lambda_init = 10, nu = penlist, folds = 5, usesubset = TRUE, maxiter = 50, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{aber}{
%%     ~~Describe \code{aber} here~~
}
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{pi_a_init}{
%%     ~~Describe \code{pi_a_init} here~~
}
  \item{rho_init}{
%%     ~~Describe \code{rho_init} here~~
}
  \item{nu0_init}{
%%     ~~Describe \code{nu0_init} here~~
}
  \item{sigma_c_init}{
%%     ~~Describe \code{sigma_c_init} here~~
}
  \item{lambda_init}{
%%     ~~Describe \code{lambda_init} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{folds}{
%%     ~~Describe \code{folds} here~~
}
  \item{usesubset}{
%%     ~~Describe \code{usesubset} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
        result = deconvEM_laplace(Y, eta, mu, aber, V, weight = weight, 
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
            temp <- deconvEM_laplace(Y_train, eta_train, mu_train, 
                aber, V, weight = weight_train, pi_a_train, rho_init_train, 
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
    result = deconvEM_laplace(Y, eta, mu, aber, V, weight = weight, 
        pi_a_init, rho_init, nu0_init, sigma_c_init, lambda_init, 
        nu = choosenu, maxiter, verbose = verbose)
    return(list(result, choosenu, losslist))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
