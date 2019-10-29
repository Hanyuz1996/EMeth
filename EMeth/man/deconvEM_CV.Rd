\name{deconvEM_CV}
\alias{deconvEM_CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{deconvEM_CV
%%  ~~function to do ... ~~
}
\description{
Performs a cross-validation on the ridge penalty hyperparameter for EMeth deconvolution algorithm. ~~
}
\usage{
deconvEM_CV(Y, eta, mu, aber = TRUE, V = "c", weight = matrix(1, 5, 5), pi_a_init, rho_init, nu0_init = rep(0, 50), sigma_c_init, lambda_init = 10, nu = penlist, folds = 5, usesubset = TRUE, maxiter = 50, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
    mixture sample: K*I matrix with K probes and I samples.
}
  \item{eta}{
    purity of each tumor sample: a vector with lenght I.
}
  \item{mu}{
    reference dataset: K*Q matrix, where Q is the number of cell types.
}
  \item{aber}{
    logic variable: TRUE for there are aberrant cells and FALSE for no aberrant cells.
}
  \item{V}{
    Variance structure: can accept 'b' for BinomEM, 'c' for OriEM. A more convenient way to set specific weight structure.
}
  \item{weight}{
    Weight structure: accept specific weight matrix for MaxVarEM. Users can also design specific weights.
}
  \item{pi_a_init}{
   initial value of mixture proportion pi_a.
}
  \item{rho_init}{
    initial value of cell type proportions. I*Q matrix.
}
  \item{nu0_init}{
    initial value for aberrant cell type methylation.
}
  \item{sigma_c_init}{
    initial value for variance of consistent cell types.
}
  \item{lambda_init}{
    initial value for ratio sigma_a/sigma_c.
}
  \item{nu}{
     penalty parameter: can accept one single value or a vector, recommend with same order of number of probes.
}
  \item{folds}{
    number folds to perform cross-validation, recommend: folds = 5.
}
  \item{usesubset}{
    logic variable: use all probes or a random subset of probes to do cross validation.
}
  \item{maxiter}{
     maximum number of iteration times.
}
  \item{verbose}{
     logic varible: will the algorithm print outputs during iterations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  This function returns a list. 
  Result: See the documentation for deconvEM function to get a detailed explanation of the "result" value.
  choosenu: returns the value of nu choosed.
  losslist: returns the weighted L2 loss for each nu for each sample.
}
\references{
 TBDone
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
\examples{
See github pages for a concrete example of usage.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ deconv}% use one of  RShowDoc("KEYWORDS")
\keyword{ EMeth }% __ONLY ONE__ keyword per line
