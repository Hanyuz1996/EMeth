\name{emeth}
\alias{emeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
emeth: cell type decomposition from DNA methylation data based on EM-type algorithm and penalized likelihood maximization. 
}
\description{
This function is an implementation of EMeth algorithm. It provides estimation of cell type proportion for each sample in bulk tissue samples Y. It needs a reference matrix mu to conduct the estimation. EMeth has two families, normal and laplace, depending on what likelihood is used. It also has a penalization term (ridge penalty).
}
\usage{
emeth(Y, eta, mu, aber, V, init = "default", family = "laplace", nu = 0, maxiter = 50, verbose = FALSE)
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
  \item{maxiter}{
    max time of iteration of the EM algorithm, default to be 50.
}
  \item{verbose}{
    logic variable. If TRUE then will print additional information in iteration of EMeth.
}
}
\value{
  \item{rho }{estimated proportion of cell types in each sample.}
  \item{sigma_c }{estimated variance.}
  \item{lambda}{estimated ratio between aberrant variance and consistent variance.}
  \item{nu0}{estimated methylation level of unknown cell type.}
  \item{pi_a}{estimates of the proportion of aberrant probes in each sample.}
  \item{gamma}{estimated probability of each probe in each sample to be aberrant.}
  \item{weights}{the working weights used in the iteration of EMeth.}
  \item{iter}{Number of iterations.}
}
\author{
  Hanyu Zhang 
}
\examples{
## See examples folder
}
\keyword{ methods }
