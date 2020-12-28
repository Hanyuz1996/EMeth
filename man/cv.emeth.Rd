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
     Matrix of size K*I. Methylation of bulk sample for which the cell type decomposition is to be estimated. Y is usually a K*I matrix where K is the number of probes used and I is the number of samples. 
}
  \item{eta}{
      Vector of size I. Tumor purity of each sample. 
}
  \item{mu}{
      Matrix of size K*Q. Reference matrix, provided by literature. A sample reference data of immune cells is provided in the example directory. 
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
    nonnegative numbers that stand for the penalties, and the opitmal penalty value will be chosen by cross-validation. 
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

\value{

  \item{result }{The result of EMeth algorithm using the penalty value selected by cross-validation. It is a list and documentation of its entries can be found in the help file for function \code{emeth}}
  \item{choosenu }{The chosen value of nu, the penalty.}
  \item{losslist}{A matrix saving the loss for each fold and each choice of nu.}
}
\author{
Hanyu Zhang
}

\keyword{ methods }
\examples{
## See examples folder.
}