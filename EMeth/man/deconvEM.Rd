\name{deconvEM}
\alias{deconvEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   ~~deconvEM ~~
}
\description{
  ~~ Implementation of EMeth algorithm of normal distribution version.  ~~
}
\usage{
deconvEM(Y, eta, mu, aber = TRUE, V = "c", weight = matrix(1, 5, 5), pi_a_init, rho_init, nu0_init = rep(0, 50), sigma_c_init = 0.1, lambda_init = 10, nu = 0, maxiter = 50, verbose = FALSE)
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
    5, 5), pi_a_init, rho_init, nu0_init = rep(0, 50), sigma_c_init = 0.1, 
    lambda_init = 10, nu = 0, maxiter = 50, verbose = FALSE) 
{
    Q = ncol(mu)
    K = nrow(mu)
    I = ncol(Y)
    X = as.data.frame(mu)
    pi_a = pi_a_init
    if (!aber) {
        nu0_init = rep(0, K)
        eta = rep(0, I)
    }
    nu0 = nu0_init
    rho = rho_init
    sigma_c = sigma_c_init
    if (V == "b") {
        mu00 = mu * (1 - mu)
        mu00[mu00 < 0.05] = 0.05
        Vf = mu00 \%*\% t(rho_init)
        Vf = Vf + matrix(rep(nu0_init * (1 - nu0_init), ncol(Y)), 
            ncol = ncol(Y), byrow = FALSE) \%*\% diag(eta)
        Vc = sigma_c_init * sigma_c_init * Vf
        Va = lambda_init * Vc
    }
    else if (V == "c") {
        Vc = matrix(sigma_c_init^2, K, I)
        Va = matrix(lambda_init * sigma_c_init^2, K, I)
    }
    else if (V == "w") {
        Vf = weight
        Vc = Vf
        Va = lambda_init * Vf
    }
    lambda = lambda_init
    nu0 = nu0_init
    nu0.m = matrix(rep(nu0, times = I), ncol = I, byrow = FALSE)
    Mu = mu \%*\% t(rho) + nu0.m \%*\% diag(eta)
    for (iter in 1:maxiter) {
        if (iter\%\%10 == 0 & verbose) {
            cat("-------------------\n")
            cat(iter, date(), "\n")
            cat(pi_a, "\n")
            cat(nu0[1:5], "\n")
            cat(rho[1, ], "\n")
            cat(lambda, "\n")
        }
        pi_a_old = pi_a
        rho_old = rho
        Vc_old = Vc
        Va_old = Va
        lambda_old = lambda
        nu0_old = nu0
        Va[Va < 0.001] = 0.001
        Vc[Vc < 0.001] = 0.001
        e = Y - Mu
        e[which(e > 0.99)] = 0.99
        e[which(e < -0.99)] = -0.99
        pdf_a = dnorm(e/sqrt(Va)) \%*\% diag(pi_a)/sqrt(Va)
        pdf_c = dnorm(e/sqrt(Vc)) \%*\% diag(1 - pi_a)/sqrt(Vc)
        gamma = pdf_a/(pdf_a + pdf_c)
        V_renorm = Vc/norm(Vc)
        W = (gamma + lambda * (1 - gamma))/(lambda * V_renorm)
        if (verbose) {
            cat(sprintf("\%s NAs in W", length(which(is.na(W)))), 
                "\n")
        }
        pi_a = colMeans(gamma)
        Y_c = Y - nu0.m \%*\% diag(eta)
        nulist = rep(0, I)
        for (i in 1:I) {
            C = norm(t(mu) \%*\% diag(W[, i]) \%*\% mu)
            A = rbind(rep(-1, Q), diag(rep(1, Q)))
            b = c(-1 + eta[i], rep(0, Q)) * sqrt(C)
            D = (t(mu) \%*\% diag(W[, i]) \%*\% mu + nu * diag(ncol(mu)))/(C)
            d = t(t(mu) \%*\% diag(W[, i]) \%*\% Y_c[, i])/sqrt(C)
            temp = solve.QP(D, d, t(A), b, meq = 1)$solution/sqrt(C)
            rho[i, ] = temp
        }
        Y_ab = Y - mu \%*\% t(rho)
        if (aber) {
            for (k in 1:K) {
                nu0[k] = min(1, max(0, sum(W[k, ] * eta * Y_ab[k, 
                  ])/sum(W[k, ] * eta^2)))
            }
        }
        nu0.m = matrix(rep(nu0, times = I), ncol = I, byrow = FALSE)
        Mu = mu \%*\% t(rho) + nu0.m \%*\% diag(eta)
        if (V == "b") {
            Vf = mu00 \%*\% t(rho)
            Vf = Vf + matrix(rep(nu0 * (1 - nu0), ncol(Y)), ncol = ncol(Y), 
                byrow = FALSE) \%*\% diag(eta)
            Vc = sigma_c * sigma_c * Vf
            Va = lambda * Vc
        }
        else if (V == "c") {
            Vf = matrix(1, K, I)
        }
        else if (V == "w") {
            Vf = weight
        }
        sigma_a = (sum(gamma * ((Y - Mu)^2)/Vf)/sum(gamma))
        sigma_c = (sum((1 - gamma) * (Y - Mu)^2/Vf)/sum(1 - gamma))
        lambda = sigma_a/sigma_c
        Va = Vf * sigma_a
        Vc = Vf * sigma_c
        if (max(abs(rho_old - rho)) < 1e-04) {
            break
        }
    }
    cat("sigma_c_square", sigma_c, "\n")
    list(rho = rho, sigma_c = sigma_c, lambda = lambda, nu0 = nu0, 
        pi_a = pi_a, gamma = gamma, weights = W, iter = iter)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
