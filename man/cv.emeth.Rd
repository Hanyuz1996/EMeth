\name{cv.emeth}
\alias{cv.emeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cv.emeth: cross validation for emeth.
}
\description{
An implementation of cross validation for the ridge penalty parameter in EMeth.
}
\usage{
cv.emeth(Y, eta, mu, aber, V, init = "default", nu, family = "laplace", folds = 5, usesubset = TRUE, maxiter = 50, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Y}{
     Matrix of size K*I. Methylation of bulk sample for which the cell type decomposition is to be estimated. Y is usually a K*I matrix where K is the number of probes used and I is the number os samples. 
}
  \item{eta}{
      Vector of size I. Tumor purity of each sample. 
}
  \item{mu}{
      Matrix of size K*Q. Reference matrix, provided by literature. A sample reference data in the example directory is provided, that is obtained from six previous literatures. 
}
  \item{aber}{
     Logic variable: if there is unknown aberrant cell type. 
}
  \item{V}{
     string: default to be 'c' which stands for constant weight for all probes. It might be  'b' for binomial variance structure or 'w' for specific weight structure of variance.
}
  \item{init}{
     If init is a string 'default', we will adopt the default random initialization of all parameters for the algorithm. Otherwise one can provide a list of initialization of parameters.
}
  \item{family}{
     string: accept 'normal' or 'laplace' to specify what likelihood will be used in the algorithm.
}
  \item{nu}{
    a nonnegative number that stands for the penalty. 
}
  \item{folds}{
     Specify the number of folds for the cross validation.
}
  \item{usesubset}{
     Logic variable, if it is true, a random sampled subset of all probes are used to perform cross-validation.
}
    \item{maxiter}{
    max time of iteration of the EM algorithm, default to be 50.
}
  \item{verbose}{
    logic variable. If TRUE then will print additional information in iteration of EMeth.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{result }{The result of EMeth algorithm after fitting the best selected parameter.}
  \item{choosenu }{The chosen parameter of nu.}
  \item{losslist}{A matrix saving the loss for each fold and each choice of nu.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hanyu Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\examples{
## See examples folder.
}